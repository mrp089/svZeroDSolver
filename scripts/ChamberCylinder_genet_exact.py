"""Build the EXACT Genet 2023 valve + 2-stage Windkessel circulation
(eq 35 piecewise valves, eq 36b-c two-stage C-R-C-R Windkessel) around the
ChamberCylinder, using only existing svZeroDSolver blocks:
  PressureBC(Pat) -> upstream(tiny) -[mitral PiecewiseValve]-> ventricle(ChamberCylinder)
    -[aortic PiecewiseValve]-> arterial(BloodVesselRC: Car,Rp) -> distal(BloodVesselRC: Cd,Rd) -> PressureBC(Pvs)

Genet Table 1 valve/circulation parameters (SI):
  K_at^-1 = 1.1111e5, K_ar^-1 = 7.6923e4, K_iso^-1 = 2.0e9  (Pa.s.m^-3)
  Cvalve = 9.0e-9 (m^3/Pa); Rp = 1.35e7, Rd = 1.0949e8 (Pa.s.m^-3)
  Car = 2.115e-10, Cd = 2.0158e-8 (m^3/Pa)
"""
import json

# Genet Table 1
K_at_inv = 1.1111e5
K_ar_inv = 7.6923e4
K_iso_inv = 2.0e9
Cvalve = 9.0e-9
Rp = 1.35e7
Rd = 1.0949e8
Car = 2.115e-10
Cd = 2.0158e-8

# Genet Table 1 geometry + material (baseline) — EXACT values, no fitting.
GEOM = dict(Ri=0.019, Re=0.033, length=0.0571, alpha_endo=60.0, alpha_epi=-60.0)
MAT = dict(C1=7.0, C2=0.0, C3=700.0, C4=2.0, C5=50.0, C6=4.0, gamma=70.0,
           k_s=1e8, k_0=260e3, mu=70.0)
# Genet enforces incompressibility with a mixed u/p Lagrange-multiplier pressure
# (mixed=1, the faithful default here). kappa is only the fallback displacement
# penalty (mixed=0); a low-order penalty locks volumetrically at coarse meshes,
# so it is NOT used for the reproduction. kappa is retained for the fallback.
KAPPA = 1e7
SIGMA0 = 65e3  # Genet Table 1 "Maximum active stress" sigma_0

# Atrial and venous pressures are prescribed inputs Genet does not tabulate; the
# MEDISIM PhysioBlocks reference (physioblocks/physioblocks, references/
# full_configurations/spherical_heart_sim.jsonc) gives them. P_at is the
# `atrial.blood_pressure` waveform: a `rescale_two_phases_function` between
# min=450 Pa and max=900 Pa, with the ATRIAL KICK (450->900) in LATE diastole and
# P_at HELD at 900 through end-diastole into early systole (then dropping back to
# 450). The 450 floor during mid-diastole gives the slow E-wave; the late kick
# gives the A-wave -> the two-phase filling of Fig. 5. (An earlier hand-made kick
# that instead dropped to 450 at end-diastole did NOT hold EDV and mismatched.)
# Systemic venous P_vs = 1600 Pa sets the afterload floor -> exact peak pressure.
PVS_PHYSIOBLOCKS = 1600.0
PAT_MIN, PAT_MAX = 450.0, 900.0
# PhysioBlocks reference_function (cycle fraction, value) + phase per interval
# (0 = diastole, 1 = systole) and diastole_scaling_factor alpha.
PAT_REF = [(0.0, PAT_MAX), (0.02, PAT_MAX), (0.07, PAT_MIN),
           (0.84, PAT_MIN), (0.9, PAT_MAX)]
PAT_PHASES = [1, 1, 0, 1]
PAT_ALPHA = 0.8


def atrial_pressure(period=0.8):
    """PhysioBlocks `atrial.blood_pressure` (rescale_two_phases_function),
    evaluated for a cardiac `period`. Stretches the phase-0 (diastole) vs
    phase-1 (systole) intervals of the reference waveform so it fills `period`
    while preserving the kick shape. Returns (times[s], pressures[Pa])."""
    ref, phases, alpha = PAT_REF, PAT_PHASES, PAT_ALPHA
    rp = ref[-1][0] - ref[0][0]
    beta = period / rp
    d0 = sum(ref[i][0] - ref[i - 1][0]
             for i in range(1, len(ref)) if phases[i - 1] == 0)
    d1 = rp - d0
    s0 = 1.0 + alpha * (beta - 1.0) * rp / d0
    s1 = 1.0 + (1.0 - alpha) * (beta - 1.0) * rp / d1
    t = [ref[0][0]]
    for i in range(1, len(ref)):
        w = (ref[i][0] - ref[i - 1][0]) * (s0 if phases[i - 1] == 0 else s1)
        t.append(t[-1] + w)
    return t, [v[1] for v in ref]


# Back-compat aliases (the atrial-kick waveform for the default 0.8 s period).
PAT_KICK_T, PAT_KICK_P = atrial_pressure(0.8)


# Activation: the PhysioBlocks BCS nu(t) trapezoid (activation_mode=1, set in the
# vv dict below). `tsys` anchors its systole onset (the nu=0 upstroke), retimed to
# Fig 5 against the VOLUME curve: tsys=0.12 s aligns the ejection downstroke
# (130->80 mL crossings within ~5 ms of Fig 5) and the torsion rise. (The pressure
# upstroke then sits ~40 ms late -- the model has essentially no isovolumic phase,
# a separate valve/afterload-dynamics point.) `tdias` is unused in mode 1 (the
# trapezoid carries its own relaxation ramp). NOTE (documented discrepancy): the
# PhysioBlocks spherical-sim activation relaxes with nu=-20 s^-1, slower than
# whatever Genet used for the cylindrical Fig 5 -> diastolic filling is delayed
# (the residual full-cycle volume mismatch); kept for fidelity to the reference.
def build(P_at=900.0, P_vs=PVS_PHYSIOBLOCKS, sigma_max=SIGMA0, bcs_alpha=12.0, ne=12,
          tsys=0.12, tdias=0.40, steepness=0.02, integrator="stiff",
          rho_infty=0.5, ncycle=8, aortic_Rmax=None, active_model=1,
          n0_flat=True, atrial_kick=True, mixed=True):
    # Integrator note: Genet's temporal scheme is the non-dissipative midpoint
    # (rho_infty=1) made *stable* by energy-preserving algorithmic stresses + the
    # Chapelle sqrt(k_c) internal-variable update. The plain midpoint alone
    # (rho_infty=1, without those) is UNDER-damped for the stiff k_s=1e8 series
    # spring and rings (spurious HF pressure oscillations). So the default here is
    # the L-stable "stiff" integrator (rho=0), which damps that mode cleanly; a
    # faithful energy-preserving scheme needs the bespoke integrator (not built).
    aortic_Rmax = aortic_Rmax if aortic_Rmax is not None else K_iso_inv
    vv = dict(GEOM); vv.update(MAT); vv["kappa"] = KAPPA
    # BCS activation nu(t): PhysioBlocks trapezoid (activation_mode=1) between
    # alpha_min=-20 (diastole) and alpha_max=+35 (systole) -- the MEDISIM
    # `active_law.activation.{min,max}`. Its own timing supersedes tsys/tdias
    # (kept only for the mode-0 tanh fallback). C_valve is NOT put on the
    # ventricle (that over-buffers Pv); it is the aortic-root compliance in the
    # circulation below (Genet Eq 36 / PhysioBlocks capacitance_valve), so
    # c_valve=0 here.
    vv.update(dict(sigma_max=sigma_max, alpha_max=35.0, alpha_min=-20.0,
                   activation_mode=1.0, tsys=tsys, tdias=tdias, steepness=steepness,
                   num_elements=ne, active_model=active_model, bcs_alpha=bcs_alpha,
                   c_valve=0.0, mixed=1.0 if mixed else 0.0))
    # n0(e_c): Frank-Starling reduction factor. Genet/Chapelle leave it a general
    # 0<=n0<=1 factor with no formula -> the non-fitted default is n0=1 (full
    # recruitment), realized by making the Gaussian effectively flat.
    if n0_flat:
        vv["n0_center"] = 0.0
        vv["n0_width"] = 1.0e3
    cfg = {
        "boundary_conditions": [
            {"bc_name": "ATRIUM", "bc_type": "PRESSURE",
             "bc_values": ({"P": list(PAT_KICK_P), "t": list(PAT_KICK_T)}
                           if atrial_kick else {"P": [P_at, P_at], "t": [0.0, 0.8]})},
            {"bc_name": "VENOUS", "bc_type": "PRESSURE",
             "bc_values": {"P": [P_vs, P_vs], "t": [0.0, 0.8]}},
        ],
        "simulation_parameters": {
            "number_of_cardiac_cycles": ncycle,
            "number_of_time_pts_per_cardiac_cycle": 400,
            "cardiac_period": 0.8, "steady_initial": False,
            "output_variable_based": True, "absolute_tolerance": 1e-9,
            "maximum_nonlinear_iterations": 50, "output_all_cycles": False,
            "integrator": integrator, "rho_infty": rho_infty,
        },
        "initial_condition": {"pressure_all": P_at, "volume:ventricle": 6.5e-5},
        "vessels": [
            {"boundary_conditions": {"inlet": "ATRIUM"}, "vessel_id": 0,
             "vessel_length": 1.0, "vessel_name": "upstream",
             "zero_d_element_type": "BloodVessel",
             "zero_d_element_values": {"R_poiseuille": 1000.0, "C": 1e-12}},
            {"boundary_conditions": {}, "vessel_id": 1, "vessel_length": 1.0,
             "vessel_name": "ventricle", "zero_d_element_type": "ChamberCylinder",
             "zero_d_element_values": vv},
            # arterial stage: aortic-root node carries C_valve + C_ar in parallel
            # (PhysioBlocks capacitance_valve + circulation_aorta_proximal, both
            # at aorta_proximal), then R_p to the distal stage (eq 36a-b).
            {"boundary_conditions": {}, "vessel_id": 2, "vessel_length": 1.0,
             "vessel_name": "arterial", "zero_d_element_type": "BloodVesselRC",
             "zero_d_element_values": {"Rpd": Rp, "Cp": Car + Cvalve}},
            # distal stage: Cd, Rd -> venous (eq 36c)
            {"boundary_conditions": {"outlet": "VENOUS"}, "vessel_id": 3,
             "vessel_length": 1.0, "vessel_name": "distal",
             "zero_d_element_type": "BloodVesselRC",
             "zero_d_element_values": {"Rpd": Rd, "Cp": Cd}},
        ],
        "valves": [
            {"type": "PiecewiseValve", "name": "mitral",
             "params": {"Rmax": K_iso_inv, "Rmin": K_at_inv,
                        "upstream_block": "upstream", "downstream_block": "ventricle"}},
            {"type": "PiecewiseValve", "name": "aortic",
             "params": {"Rmax": aortic_Rmax, "Rmin": K_ar_inv,
                        "upstream_block": "ventricle", "downstream_block": "arterial"}},
        ],
        "junctions": [
            {"junction_name": "J_ar", "junction_type": "NORMAL_JUNCTION",
             "inlet_blocks": ["arterial"], "outlet_blocks": ["distal"]},
        ],
    }
    return cfg
