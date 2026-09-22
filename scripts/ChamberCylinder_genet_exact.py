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
import math

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
# Genet enforces incompressibility with a mixed u/p Lagrange-multiplier pressure,
# the (only) formulation of the ChamberCylinder block.
SIGMA0 = 65e3  # Genet Table 1 "Maximum active stress" sigma_0

# ---------------------------------------------------------------------------
# ECG-parametric input functions -- the AUTHORS' actual code (Genet collaborators,
# cardiac_input_functions.py: `nagumo_piecewise_linear` + `atrial_pressure_
# piecewise_linear`). All timings DERIVE from ECG intervals and scale with the
# heart-beat duration T (activation via the Bazett QT = QTc*sqrt(T)); the values
# at T=0.8 s reproduce the authors' Fig-5 waveforms. This supersedes the earlier
# hand-timed estimates (tsys/tdias/kick were fitted to the digitized curves; they
# are now the authors' derived quantities).
# ---------------------------------------------------------------------------
MMHG = 134.0  # authors' Pa-per-mmHg conversion (cardiac_input_functions.py)

# Activation nu(t) [ActiveLaw/Activation]: delay + QRS duration + Bazett-corrected
# QT + final repolarization ramp; alpha_max/min are the +-35 plateaus.
ACT_DELAY = 0.10
ACT_QRSD = 0.080
ACT_QTC = 0.330
ACT_RAMP = 0.010
ACT_AMAX, ACT_AMIN = 35.0, -35.0


def activation_times(period=0.8):
    """Derived activation anchors (authors' ActivationParameters): returns
    (tsys, tdias) = the nu=0 upstroke (depolarization_time = delay+ramp) and the
    nu=0 downstroke (depolarization_end). The kernel builds the internal knees
    (t_max, t_plateau, t_finish) from tsys, tdias, ACT_QRSD and ACT_RAMP."""
    QT = ACT_QTC * math.sqrt(period)      # Bazett corrected QT
    ST = QT - ACT_QRSD
    tsys = ACT_RAMP + ACT_DELAY           # depolarization_time (nu=0 up)
    tdias = ACT_QRSD + tsys + ST          # depolarization_end  (nu=0 down)
    return tsys, tdias


# Atrial pressure P_at(t) [Atria/PressureLaw]: two-level trace, peak 7 mmHg, floor
# 7-2.5 = 4.5 mmHg; the KICK onset/end are set by the PQ interval and AV-node delay.
PAT_MAX = 7.0 * MMHG            # 938 Pa
PAT_MIN = PAT_MAX - 2.5 * MMHG  # 603 Pa
ATR_PQ = 0.14                   # PQ interval
ATR_AVNODE = 0.080              # AV-node delay
# Systemic venous pressure P_vs: NOT tabulated in genet23 (it only appears in the
# Eq 36 distal node). FITTED here to the Fig 5 end-systolic volume: P_vs=350 Pa
# (~2.6 mmHg) gives ESV=74.0 mL (data 74). This is the ONE remaining fitted
# quantity -- Caruel Table 1 uses 1000 Pa (~7.5 mmHg), which overshoots ESV to 76.4.
# P_vs sets only the afterload/ESV floor (~1.8 mL per 500 Pa); EDV, peak P and
# torsion are insensitive to it.
P_VS = 350.0  # FIT to Fig 5 ESV (genet23 does not tabulate P_vs; Caruel = 1000)


def atrial_pressure(period=0.8):
    """Authors' `atrial_pressure_piecewise_linear`: PAT_MAX held from end-diastole
    into early systole, linear drop to PAT_MIN over [0.05, 0.10] s, floor through
    mid-diastole, then the atrial KICK PAT_MIN->PAT_MAX over
    [T-PQ+ramp, T-AVnode+ramp]. Returns (times[s], pressures[Pa])."""
    kick_on = period - ATR_PQ + ACT_RAMP       # 0.670 at T=0.8
    kick_end = period - ATR_AVNODE + ACT_RAMP  # 0.730 at T=0.8
    t = [0.0, 0.05, 0.10, kick_on, kick_end, period]
    P = [PAT_MAX, PAT_MAX, PAT_MIN, PAT_MIN, PAT_MAX, PAT_MAX]
    return t, P


# Derived anchors + atrial waveform for the default 0.8 s period.
TSYS_08, TDIAS_08 = activation_times(0.8)
PAT_KICK_T, PAT_KICK_P = atrial_pressure(0.8)


# Activation: the authors' ECG-derived piecewise-linear nu(t) (`nagumo`,
# activation_mode=1). The kernel builds the exact waveform from the anchors tsys
# (nu=0 upstroke) and tdias (nu=0 downstroke) plus act_qrs (QRS width) and act_ramp
# (final repolarization ramp); tsys/tdias default to activation_times(0.8) =
# (0.110, 0.405) s and scale with the heart rate via Bazett. alpha_max/min = +-35.
# NOTE: with the corrected proximal compliance (C_p=C_ar) the model has a real
# isovolumic contraction (P_v rises to the aortic diastolic ~9 kPa with both valves
# shut before the aortic valve opens ~170 ms).
def build(P_at=PAT_MAX, P_vs=P_VS, sigma_max=SIGMA0, bcs_alpha=12.0, ne=12,
          tsys=TSYS_08, tdias=TDIAS_08, alpha_max=ACT_AMAX, alpha_min=ACT_AMIN,
          act_qrs=ACT_QRSD, act_ramp=ACT_RAMP,
          integrator="stiff", rho_infty=0.5, ncycle=8,
          aortic_Rmax=None, atrial_kick=True, density=1000.0):
    # Integrator note: Genet's temporal scheme is the non-dissipative midpoint
    # (rho_infty=1) made *stable* by energy-preserving algorithmic stresses + the
    # Chapelle sqrt(k_c) internal-variable update. The plain midpoint alone
    # (rho_infty=1, without those) is UNDER-damped for the stiff k_s=1e8 series
    # spring and rings (spurious HF pressure oscillations). So the default here is
    # the L-stable "stiff" integrator (rho=0), which damps that mode cleanly; a
    # faithful energy-preserving scheme needs the bespoke integrator (not built).
    # Aortic valve closed resistance: genet23 Eq 36 (arterial line) has ONLY
    # K_ar<P_v-P_ar>_+, i.e. NO closed-leak on the aortic side (conductance 0 when
    # P_v<P_ar). So the aortic Rmax is effectively infinite (a hard one-way valve),
    # NOT K_iso^-1. (The K_iso isovolumic leak in Eq 36 is on the ATRIAL side and is
    # carried by the mitral valve's Rmax=K_iso^-1.) Effect on Fig 5 is negligible.
    aortic_Rmax = aortic_Rmax if aortic_Rmax is not None else 1.0e12
    vv = dict(GEOM); vv.update(MAT)
    # BCS activation nu(t): authors' ECG-derived nagumo (the kernel's only
    # activation law), built from tsys/tdias (nu=0 up/down) + act_qrs/act_ramp;
    # +-35 plateaus.
    # Relaxation: the authors' alpha_min=-35 already IS the diastolic relaxation
    # rate, so the kernel uses a FIXED relaxation (w=1), not the Caruel
    # length-dependent m0 amplification (which would double-count it and fire the
    # E-wave ~30 ms early).
    # Cavity/valve compliance Cvalve is a Genet add-on that has no place in the
    # Caruel valve+Windkessel topology (see docstring), so it is not modeled.
    # Full dynamics (genet23 Eqs 8,18,45): the model is fully DYNAMICAL -- the
    # inertia force rho_0*u_ddot (consistent mass matrix + velocity companion DOFs)
    # with density rho_0 = 1 kg/L = 1000 (Table 1). For cardiac parameters the
    # inertia is small (~1e-4 of the internal/pressure forces), so it barely shifts
    # P/V and only modestly the twist. (genet23's own non-dissipative midpoint
    # scheme is singular on this DAE, so the L-stable "stiff" integrator is used;
    # the O(zeta_dot^2) convective term D2u is omitted, as in the kernel note.)
    vv.update(dict(sigma_max=sigma_max, alpha_max=alpha_max, alpha_min=alpha_min,
                   tsys=tsys, tdias=tdias, act_qrs=act_qrs, act_ramp=act_ramp,
                   num_elements=ne, bcs_alpha=bcs_alpha, density=density))
    # n0(e_c): the Frank-Starling recruitment factor is the fixed PhysioBlocks
    # piecewise-linear curve baked into the kernel (frank_starling()); it is not a
    # builder or block input. Over the operating range e_c in [-0.06, 0.29] it
    # rides the ascending limb/plateau (n0 ~ 0.4 -> 1.0), so n0 ~ 1 at peak.
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
            "number_of_time_pts_per_cardiac_cycle": 800,  # dt=1 ms (genet23 sec 3.1)
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
