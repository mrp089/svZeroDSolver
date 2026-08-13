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
# kappa is NOT a Genet parameter (Genet uses a Lagrange-multiplier pressure; we
# use a penalty). Large enough that J~1. Documented as a numerical choice.
KAPPA = 1e7
SIGMA0 = 65e3  # Genet Table 1 "Maximum active stress" sigma_0

# Atrial and venous pressures are prescribed inputs Genet does not tabulate; the
# MEDISIM PhysioBlocks reference (references/full_configurations/
# spherical_heart_sim.jsonc) gives the concrete values. P_at is a waveform
# (baseline 450 Pa, pre-systolic atrial kick to 900 Pa peaking just before
# systole onset t_sys~0.13 s), and the systemic venous pressure P_vs = 1600 Pa.
PVS_PHYSIOBLOCKS = 1600.0
PAT_KICK_T = [0.0, 0.03, 0.06, 0.13, 0.17, 0.80]
PAT_KICK_P = [450.0, 450.0, 900.0, 900.0, 450.0, 450.0]


def build(P_at=900.0, P_vs=0.0, sigma_max=SIGMA0, bcs_alpha=12.0, ne=12,
          tsys=0.13, tdias=0.45, steepness=0.02, integrator="stiff",
          ncycle=8, aortic_Rmax=None, active_model=1, n0_flat=True,
          atrial_kick=False):
    aortic_Rmax = aortic_Rmax if aortic_Rmax is not None else K_iso_inv
    vv = dict(GEOM); vv.update(MAT); vv["kappa"] = KAPPA
    vv.update(dict(sigma_max=sigma_max, alpha_max=30.0, alpha_min=-30.0,
                   tsys=tsys, tdias=tdias, steepness=steepness,
                   num_elements=ne, active_model=active_model, bcs_alpha=bcs_alpha,
                   c_valve=Cvalve))
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
            "integrator": integrator,
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
