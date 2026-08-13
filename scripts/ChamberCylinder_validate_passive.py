# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
# University of California, and others. SPDX-License-Identifier: BSD-3-Clause
"""Passive-response validation for the ChamberCylinder block.

Compares the block's purely passive pressure-volume curve against an
independent closed-form-style incompressible thick-walled cylinder inflation
that uses the *same* Genet (2023) isotropic strain-energy law (Eq. 25).

KEY RESULT (volumetric locking): the displacement penalty (`mixed=0`) is
*stiffer* than the incompressible analytical at coarse meshes and converges to
it only at ne~32 - classic volumetric locking (linear elements + fully
integrated near-incompressible bulk term). The **mixed u/p** formulation
(`mixed=1`, Genet's Lagrange-multiplier pressure) matches the analytical
mesh-independently from ne=8. Run this script to see both columns.

IMPORTANT: the cavity pressure is applied by loading the ventricle inlet
*directly* (via a junction), not through a smooth `ValveTanh`. A tanh valve
sits half-open at zero pressure drop (R ~ (Rmin+Rmax)/2), which throttles a
constant-pressure inflation test - an artifact of the test harness.

Run:  PYTHONPATH=build/python python scripts/ChamberCylinder_validate_passive.py
"""
import numpy as np
import pysvzerod

Ri, Re, L0 = 0.019, 0.033, 0.0571            # Genet Table 1 geometry (m)
C1, C2, C3, C4, C5, C6 = 7.0, 0.0, 700.0, 2.0, 50.0, 4.0  # Genet Table 1 passive

trapz = getattr(np, "trapezoid", None) or (
    lambda y, x: np.sum(0.5 * (y[1:] + y[:-1]) * (x[1:] - x[:-1])))


def analytical_iso(a, n=600):
    """Cavity pressure & volume for an incompressible cylinder inflated so the
    inner radius is `a`, isotropic part of Genet W^e, closed-end axial balance."""
    def cav(a, lz):
        b = np.sqrt(a**2 + (Re**2 - Ri**2) / lz)
        r = np.linspace(a, b, n)
        R = np.sqrt(np.clip((r**2 - a**2) * lz + Ri**2, 1e-12, None))
        lt = r / R
        lr = 1.0 / (lt * lz)
        I1 = lr**2 + lt**2 + lz**2
        dW = C1 + 2 * C3 * C4 * (I1 - 3) * np.exp(np.minimum(C4 * (I1 - 3)**2, 300))
        st = 2 * dW * lt**2          # lam_theta * dW/dlam_theta
        sr = 2 * dW * lr**2          # lam_r     * dW/dlam_r
        return trapz((st - sr) / r, r), b
    def axial(a, lz):
        b = np.sqrt(a**2 + (Re**2 - Ri**2) / lz)
        r = np.linspace(a, b, n)
        R = np.sqrt(np.clip((r**2 - a**2) * lz + Ri**2, 1e-12, None))
        lt = r / R
        lr = 1.0 / (lt * lz)
        I1 = lr**2 + lt**2 + lz**2
        dW = C1 + 2 * C3 * C4 * (I1 - 3) * np.exp(np.minimum(C4 * (I1 - 3)**2, 300))
        sz = 2 * dW * lz**2 - 2 * dW * lr**2
        P, _ = cav(a, lz)
        return 2 * np.pi * trapz(sz * r, r) - P * np.pi * a**2
    lo, hi, flo = 0.5, 1.8, axial(a, 0.5)
    for _ in range(60):
        m = 0.5 * (lo + hi)
        fm = axial(a, m)
        if flo * fm <= 0:
            hi = m
        else:
            lo, flo = m, fm
    lz = 0.5 * (lo + hi)
    P, _ = cav(a, lz)
    return P, np.pi * a**2 * (L0 * lz) * 1e6


def model_V(P_pa, ne=16, C5v=C5, mixed=0):
    """Purely passive block volume at prescribed cavity pressure P_pa.

    The pressure is applied directly to the ventricle inlet through a junction
    (no valve), so the cavity sees the full pressure and equilibrates.
    `mixed=1` selects the mixed u/p (locking-free) incompressibility."""
    vv = {"Ri": Ri, "Re": Re, "length": L0, "alpha_endo": 60.0, "alpha_epi": -60.0,
          "C1": C1, "C2": C2, "C3": C3, "C4": C4, "C5": C5v, "C6": C6,
          "kappa": 1e7, "gamma": 70.0, "sigma_max": 0.0, "alpha_max": 0.0,
          "alpha_min": 0.0, "tsys": 0.13, "tdias": 0.45, "steepness": 0.02,
          "num_elements": ne, "active_model": 0, "c_valve": 0.0,
          "mixed": float(mixed)}
    cfg = {
        "boundary_conditions": [
            {"bc_name": "ATRIUM", "bc_type": "PRESSURE",
             "bc_values": {"P": [P_pa, P_pa], "t": [0, 0.8]}},
            {"bc_name": "OUT", "bc_type": "FLOW",
             "bc_values": {"Q": [0.0, 0.0], "t": [0, 0.8]}}],
        "simulation_parameters": {"number_of_cardiac_cycles": 6,
                                  "number_of_time_pts_per_cardiac_cycle": 150,
                                  "cardiac_period": 0.8, "integrator": "genalpha",
                                  "absolute_tolerance": 1e-9,
                                  "output_variable_based": True},
        "initial_condition": {"pressure_all": 0.0, "volume:ventricle": 6.5e-5},
        "vessels": [
            {"boundary_conditions": {"inlet": "ATRIUM"}, "vessel_id": 0,
             "vessel_length": 1.0, "vessel_name": "upstream",
             "zero_d_element_type": "BloodVessel",
             "zero_d_element_values": {"R_poiseuille": 1e-3, "C": 1e-12}},
            {"boundary_conditions": {"outlet": "OUT"}, "vessel_id": 1,
             "vessel_length": 1.0, "vessel_name": "ventricle",
             "zero_d_element_type": "ChamberCylinder", "zero_d_element_values": vv}],
        "valves": [],
        "junctions": [{"junction_name": "J", "junction_type": "NORMAL_JUNCTION",
                       "inlet_blocks": ["upstream"], "outlet_blocks": ["ventricle"]}],
    }
    try:
        out = pysvzerod.simulate(cfg)
    except Exception:
        return float("nan")
    d = out[out["name"] == "volume:ventricle"].sort_values("time")["y"].to_numpy()
    return d[-1] * 1e6


if __name__ == "__main__":
    print("Passive P-V: analytical incompressible cylinder (Genet W^e, isotropic")
    print("part) vs the ChamberCylinder block, isotropic (C5=0) so it is directly")
    print("comparable to the analytical. Penalty locks; mixed u/p matches.\n")
    print(f"{'a[mm]':>6} {'P[kPa]':>8} {'V_analyt':>9} | penalty ne=8/16/32 | mixed ne=8/16/32")
    for a_mm in [24, 27]:
        P, Va = analytical_iso(a_mm / 1000.0)
        pen = [model_V(P, ne=n, C5v=0.0, mixed=0) for n in (8, 16, 32)]
        mix = [model_V(P, ne=n, C5v=0.0, mixed=1) for n in (8, 16, 32)]
        sp = " ".join(f"{v:5.1f}" for v in pen)
        sm = " ".join(f"{v:5.1f}" for v in mix)
        print(f"{a_mm:6.1f} {P/1000:8.3f} {Va:9.1f} | {sp} | {sm}")
    print("\nPenalty needs ne~32 to reach the analytical (volumetric locking);")
    print("mixed u/p matches it from ne=8 (mesh-independent, locking-free).")
