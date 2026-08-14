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

# Atrial and venous pressures are prescribed inputs Genet does not tabulate.
# P_at: with the mixed-u/p wall (locking-free) a CONSTANT P_at = 0.9 kPa (the
# end-diastolic pressure read from Fig. 5) fills the ventricle to EDV~135 mL on
# its own -- so no atrial "kick" is needed (`atrial_kick=False` default). The
# earlier PhysioBlocks kick waveform was a crutch for the penalty wall's
# under-fill; on the mixed wall it produced a spurious sharp late filling step
# that does not match Fig. 5's smooth E-wave/diastasis/A-wave. (The kick arrays
# are kept for reference / experimentation only.) Systemic venous P_vs = 1600 Pa
# (PhysioBlocks) sets the afterload floor -> exact peak pressure.
PVS_PHYSIOBLOCKS = 1600.0
PAT_KICK_T = [0.0, 0.50, 0.66, 0.72, 0.80]
PAT_KICK_P = [450.0, 450.0, 900.0, 900.0, 450.0]


# Activation timing read from Fig. 5 (Genet does not tabulate it): systole onset
# t_sys=0.11 s (pressure begins rising) and relaxation onset t_dias=0.40 s
# (pressure leaves the plateau) -> ~290 ms systole, matching the figure. (Earlier
# ad-hoc 0.13/0.45 held the ventricle contracted ~45 ms too long.)
def build(P_at=900.0, P_vs=0.0, sigma_max=SIGMA0, bcs_alpha=12.0, ne=12,
          tsys=0.11, tdias=0.40, steepness=0.02, integrator="stiff",
          rho_infty=0.5, ncycle=8, aortic_Rmax=None, active_model=1,
          n0_flat=True, atrial_kick=False, mixed=True):
    # Integrator note: Genet's temporal scheme is the non-dissipative midpoint
    # (rho_infty=1) made *stable* by energy-preserving algorithmic stresses + the
    # Chapelle sqrt(k_c) internal-variable update. The plain midpoint alone
    # (rho_infty=1, without those) is UNDER-damped for the stiff k_s=1e8 series
    # spring and rings (spurious HF pressure oscillations). So the default here is
    # the L-stable "stiff" integrator (rho=0), which damps that mode cleanly; a
    # faithful energy-preserving scheme needs the bespoke integrator (not built).
    aortic_Rmax = aortic_Rmax if aortic_Rmax is not None else K_iso_inv
    vv = dict(GEOM); vv.update(MAT); vv["kappa"] = KAPPA
    vv.update(dict(sigma_max=sigma_max, alpha_max=30.0, alpha_min=-30.0,
                   tsys=tsys, tdias=tdias, steepness=steepness,
                   num_elements=ne, active_model=active_model, bcs_alpha=bcs_alpha,
                   c_valve=Cvalve, mixed=1.0 if mixed else 0.0))
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
            # arterial stage: Car, Rp  (eq 36b)
            {"boundary_conditions": {}, "vessel_id": 2, "vessel_length": 1.0,
             "vessel_name": "arterial", "zero_d_element_type": "BloodVesselRC",
             "zero_d_element_values": {"Rpd": Rp, "Cp": Car}},
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
