"""Build the EXACT Genet 2023 valve + 2-stage Windkessel circulation around the
ChamberCylinder, using only existing svZeroDSolver blocks:
  PressureBC(Pat) -> upstream(tiny) -[mitral PiecewiseValve]-> ventricle(ChamberCylinder)
    -[aortic PiecewiseValve]-> arterial(BloodVesselRC: Car,Rp) -> distal(BloodVesselRC: Cd,Rd) -> PressureBC(Pvs)

This is the Caruel et al. (2013, Biomech Model Mechanobiol 13:897-914) reduced
circulation that Genet inherits. Caruel Eq 12 is the piecewise valve conductance
law (filling K_at / isovolumic K_p / ejection K_ar), realized here by the
svZeroDSolver PiecewiseValve (a diode resistor, Regazzoni 2022 Eq 16/22, with
R=1/K). Caruel's "two-stage Windkessel" is the C-R-C-R ladder
  C_p Pdot_ar + (P_ar - P_d)/R_p = Q            (proximal node)
  C_d Pdot_d  + (P_d - P_ar)/R_p = (P_vs - P_d)/R_d   (distal node)
which BloodVesselRC (C at inlet, R to outlet) reproduces exactly with the two
stages in series: proximal Cp=C_ar, Rpd=R_p; distal Cp=C_d, Rpd=R_d.

Genet Table 1 valve/circulation parameters (SI):
  K_at^-1 = 1.1111e5, K_ar^-1 = 7.6923e4, K_iso^-1 = 2.0e9  (Pa.s.m^-3)
    -> valve conductances K_at=9e-6, K_ar=1.3e-5, K_p=5e-10 == Caruel Table 1 EXACT
  Rp = 1.35e7, Rd = 1.0949e8 (Pa.s.m^-3); Car = 2.115e-10, Cd = 2.0158e-8 (m^3/Pa)
  Cvalve = 9.0e-9 (m^3/Pa) -- a Genet add-on ABSENT from Caruel; see below.

DISCREPANCY FOUND & CORRECTED (circulation): the proximal compliance is C_ar
ALONE. An earlier build put C_ar + Cvalve at the aortic-root node, but Cvalve
(9e-9) is 43x C_ar, which blows the proximal time constant R_p*C_p from ~3e-3 s
up to 0.12 s. Caruel calibrates R_p*C_p ~ 1e-2 s ("R_p conditions the peak
ventricular pressure"); the oversized C_p made the ejection pressure rise far too
slowly (Fig 5 peak-pressure timing 331 ms vs data ~227 ms; pressure RMS 0.52 kPa).
Using C_p = C_ar recovers the fast isovolumic->ejection pressure rise (peak ~239 ms;
pressure RMS 0.19, volume RMS 3.3 mL). Caruel has NO valve capacitance and our
PiecewiseValve is a pure resistor, so Cvalve has no node in this topology; its
intended role is an OPEN QUESTION for Genet (kept as a constant for provenance).
"""
import json

# Genet Table 1
K_at_inv = 1.1111e5
K_ar_inv = 7.6923e4
K_iso_inv = 2.0e9
Cvalve = 9.0e-9  # Genet add-on, NOT used (absent from Caruel; see module docstring)
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
# full_configurations/spherical_heart_sim.jsonc) gives their SHAPE and levels:
# P_at is the `atrial.blood_pressure` waveform, a two-level trace between min=450
# and max=900 Pa. The 450 Pa floor through mid-diastole sets the E-wave (which is
# actually driven by ventricular relaxation dropping P_v below the floor); an
# atrial KICK (450->900 Pa) in late diastole drives the A-wave, and P_at is held at
# 900 through end-diastole into early systole (the mitral valve shuts once P_v
# rises, so the fall back to 450 in early systole does not affect flow).
# Systemic venous P_vs sets the mean-arterial (afterload) floor. Caruel Table 1
# uses P_vs = 1000 Pa; we use that (it is the source value Genet inherits). An
# earlier build used PhysioBlocks' 1600 Pa to pin the peak pressure, but with the
# corrected proximal compliance the peak is set by R_p (Caruel), so 1600 is no
# longer needed and 1000 both is faithful and lowers ESV toward Fig 5 (78->76 mL).
P_VS = 1000.0
PAT_MIN, PAT_MAX = 450.0, 900.0
# The kick ONSET is the one free timing of P_at (Genet tabulates neither P_at nor
# its timing). It is calibrated to the Fig 5 volume curve: kick_onset=0.64 s (rise
# over 0.05 s) lands the mitral A-wave inflow in the Fig 5 700-750 ms window so the
# ventricle refills to EDV=135 mL exactly at end-diastole. (The PhysioBlocks
# spherical-sim reference placed the kick at ~0.75 of the cycle, which -- with the
# mitral/compliance filling lag -- pushed the A-wave past end-diastole into the
# next cycle; that is a reference-timing artifact, not Genet's Fig 5.)
PAT_KICK_ONSET = 0.64
PAT_KICK_RISE = 0.05


def atrial_pressure(period=0.8, kick_onset=PAT_KICK_ONSET, kick_rise=PAT_KICK_RISE):
    """Atrial pressure P_at(t) on the real cardiac timeline (period seconds),
    calibrated to Fig 5. PAT_MIN floor through mid-diastole; kick PAT_MIN->PAT_MAX
    over [kick_onset, kick_onset+kick_rise] driving the A-wave; held at PAT_MAX
    through end-diastole into early systole, relaxing back to PAT_MIN after systole
    onset (mitral shut). Returns (times[s], pressures[Pa]); periodic in `period`."""
    f = period / 0.8
    t = [0.0, 0.048 * f, 0.088 * f, kick_onset, kick_onset + kick_rise, period]
    P = [PAT_MAX, PAT_MAX, PAT_MIN, PAT_MIN, PAT_MAX, PAT_MAX]
    return t, P


# The atrial-kick waveform for the default 0.8 s period.
PAT_KICK_T, PAT_KICK_P = atrial_pressure(0.8)


# Activation: the PhysioBlocks BCS nu(t) trapezoid (activation_mode=1, set in the
# vv dict below), retimed to Fig 5 with TWO anchors (see get_activation): the
# nu=0 upstroke is placed at `tsys` and the nu=0 downstroke at `tdias`, with the
# PhysioBlocks ramp/plateau proportions preserved between them. Both are read
# directly off the Fig 5 PRESSURE curve (the direct readout of active tension):
#   tsys  = 0.11 s -> contraction onset (P starts rising, Fig 5 ~110 ms)
#   tdias = 0.40 s -> relaxation onset  (P starts falling, Fig 5 ~400 ms; ESV@409)
# This fixes the earlier "stays contracted too long": the PhysioBlocks shape has a
# FIXED 360 ms systole, so a single-shift retiming put the down-crossing at ~480 ms
# (~80 ms late). nu(t) is the one ingredient Genet does not tabulate, so these two
# timings are the only calibrated quantities (open question for the author).
# NOTE: with the corrected proximal compliance (C_p=C_ar) the model now has a real
# isovolumic contraction (P_v rises to the aortic diastolic ~9 kPa with both valves
# shut before the aortic valve opens ~170 ms). Small residuals remain: a slight ESV
# under-ejection (~2-4 mL) and a diastolic E-wave/untwist that recover a bit fast.
def build(P_at=900.0, P_vs=P_VS, sigma_max=SIGMA0, bcs_alpha=12.0, ne=12,
          tsys=0.11, tdias=0.40, alpha_max=35.0, alpha_min=-20.0, alpha_r=0.0,
          steepness=0.02, integrator="stiff", rho_infty=0.5, ncycle=8,
          aortic_Rmax=None, active_model=1, atrial_kick=True, mixed=True,
          bcs_relax=True):
    # Integrator note: Genet's temporal scheme is the non-dissipative midpoint
    # (rho_infty=1) made *stable* by energy-preserving algorithmic stresses + the
    # Chapelle sqrt(k_c) internal-variable update. The plain midpoint alone
    # (rho_infty=1, without those) is UNDER-damped for the stiff k_s=1e8 series
    # spring and rings (spurious HF pressure oscillations). So the default here is
    # the L-stable "stiff" integrator (rho=0), which damps that mode cleanly; a
    # faithful energy-preserving scheme needs the bespoke integrator (not built).
    aortic_Rmax = aortic_Rmax if aortic_Rmax is not None else K_iso_inv
    vv = dict(GEOM); vv.update(MAT); vv["kappa"] = KAPPA
    # BCS activation nu(t): PhysioBlocks trapezoid (activation_mode=1). Its two
    # rate plateaus (alpha_max in systole, alpha_min in diastole -- the MEDISIM
    # `active_law.activation.{min,max}`, default +35/-20) and the two-anchor
    # retiming (tsys/tdias) are builder parameters. alpha_r is the load-dependent
    # relaxation time constant (Caruel 2013 w/m0, bcs_relax=1); alpha_r=0 uses the
    # instantaneous limit w=m0(e_c). c_valve=0 on the ChamberCylinder: Cvalve is a
    # Genet add-on that has no place in the Caruel valve+Windkessel topology (see
    # docstring).
    vv.update(dict(sigma_max=sigma_max, alpha_max=alpha_max, alpha_min=alpha_min,
                   activation_mode=1.0, tsys=tsys, tdias=tdias, steepness=steepness,
                   num_elements=ne, active_model=active_model, bcs_alpha=bcs_alpha,
                   c_valve=0.0, mixed=1.0 if mixed else 0.0,
                   bcs_relax=1.0 if bcs_relax else 0.0, alpha_r=alpha_r))
    # n0(e_c): the Frank-Starling recruitment factor is the fixed PhysioBlocks
    # piecewise-linear curve baked into the kernel (frank_starling()); it is NOT a
    # builder input (the kernel's n0_center/n0_width inputs are vestigial/ignored).
    # Over the operating range e_c in [-0.06, 0.29] it rides the ascending
    # limb/plateau (n0 ~ 0.4 -> 1.0), so n0 ~ 1 at peak contraction.
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
            # proximal stage (Caruel first Windkessel eq): C_p=C_ar at the aortic
            # root, then R_p to the distal stage. C_p is C_ar ALONE -- Cvalve is
            # NOT added here (see module docstring: it destroys the proximal time
            # constant and is absent from Caruel).
            {"boundary_conditions": {}, "vessel_id": 2, "vessel_length": 1.0,
             "vessel_name": "arterial", "zero_d_element_type": "BloodVesselRC",
             "zero_d_element_values": {"Rpd": Rp, "Cp": Car}},
            # distal stage (Caruel second Windkessel eq): Cd, Rd -> venous
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
