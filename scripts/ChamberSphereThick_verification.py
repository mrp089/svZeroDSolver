# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
# University of California, and others. SPDX-License-Identifier: BSD-3-Clause
"""Verification of the ChamberSphereThick block.

Run from the repository root with the built Python module on the path::

    PYTHONPATH=<build>/python python scripts/ChamberSphereThick_verification.py

Four checks:

1. Static passive inflation against an independent reference. With
   ``density = 0``, ``gamma = 0`` and ``sigma_max = 0`` the wall equation of
   the block is exactly algebraic, so at every time step the cavity pressure
   must equal the static inflation pressure of the current endocardial
   displacement. The reference takes the *energy* route (the derivative of
   the Mooney-Rivlin strain energy along the stretch, integrated over the
   wall) rather than the block's tensorial stress route, so it checks the
   material evaluation, the kinematics and the quadrature together. It is
   run at wall-thickness ratios up to Re/Ri = 3.

2. Thin-shell limit. As d0/r0 -> 0 the block must reproduce ChamberSphere,
   including the dynamics, the viscosity and the active contraction.

3. Quadrature convergence in ``num_elements`` on a thick wall.

4. Transmural heterogeneity, the feature ChamberSphere cannot represent.
"""

import copy
import json
import os
import sys

import numpy as np
import svzerod

THIS_DIR = os.path.abspath(os.path.dirname(__file__))
BASE = json.load(
    open(os.path.join(THIS_DIR, "..", "tests", "cases", "chamber_sphere.json"))
)
W1, W2 = 10e3, 40.0  # the mooney_rivlin material of that case
GAUSS_X, GAUSS_W = np.polynomial.legendre.leggauss(40)


def series(res, name):
    """Time history of one output variable."""
    s = res[res["name"] == name]
    return s["time"].to_numpy(), s["y"].to_numpy()


def run(chamber, **sim):
    cfg = copy.deepcopy(BASE)
    cfg["chambers"][0] = chamber
    cfg["simulation_parameters"].update(sim)
    return svzerod.simulate(cfg)


def thick(Ri, Re, ne=10, gamma=10.0, density=1e3, sigma_max=185e3):
    ch = copy.deepcopy(BASE["chambers"][0])
    ch["type"] = "ChamberSphereThick"
    ch["values"] = {
        "Ri": Ri,
        "Re": Re,
        "gamma": gamma,
        "num_elements": ne,
        "density": density,
    }
    ch["active_stress"]["sigma_max"] = sigma_max
    return ch


def thin(r0, d0, eta=10.0, rho=1e3, sigma_max=185e3):
    ch = copy.deepcopy(BASE["chambers"][0])
    ch["values"] = {"rho": rho, "thick0": d0, "radius0": r0, "eta": eta}
    ch["active_stress"]["sigma_max"] = sigma_max
    return ch


def dW_dlambda(lam):
    """d/dlambda of the Mooney-Rivlin energy on C = diag(l^-4, l^2, l^2).

    With I1 = l^-4 + 2 l^2 and I2 = 2 l^-2 + l^4 (J = 1, so the isochoric
    invariants are the invariants).
    """
    return W1 * (-4.0 / lam**5 + 4.0 * lam) + W2 * (-4.0 / lam**3 + 4.0 * lam**3)


def p_static(rho_i, Ri, Re):
    """Exact passive inflation pressure of the thick sphere.

    Stationarity of the total potential over the one-parameter incompressible
    family r^3 = R^3 + r_i^3 - R_i^3 gives
    P = int dW/dlambda dlambda/drho R^2 dR / r_i^2.
    """
    ri = Ri + rho_i
    R = 0.5 * (Ri + Re) + 0.5 * (Re - Ri) * GAUSS_X
    w = GAUSS_W * 0.5 * (Re - Ri)
    r = np.cbrt(R**3 + ri**3 - Ri**3)
    return np.sum(w * dW_dlambda(r / R) * (ri**2 / (R * r**2)) * R**2) / ri**2


def check_static():
    print("1. STATIC PASSIVE INFLATION (density = 0, gamma = 0, sigma_max = 0)")
    print("   The wall equation is then exactly algebraic: P = P_ref(rho).")
    for Ri, Re in [(0.040, 0.045), (0.025, 0.045), (0.020, 0.060)]:
        res = run(
            thick(Ri, Re, ne=10, gamma=0.0, density=0.0, sigma_max=0.0),
            number_of_time_pts_per_cardiac_cycle=200,
        )
        _, rho_i = series(res, "radius:ventricle")
        _, pres = series(res, "pressure:inlet_valve:ventricle")
        ref = np.array([p_static(x, Ri, Re) for x in rho_i])
        err = np.abs(pres - ref)
        print(
            f"   Re/Ri = {Re / Ri:.2f}   |P - P_ref|: median {np.median(err):.2e} Pa,"
            f" max {err.max():.2e} Pa   (|P_ref| up to {np.abs(ref).max():.0f} Pa)"
        )


def check_thin_limit():
    print("\n2. THIN-SHELL LIMIT (vs ChamberSphere, dynamics + contraction)")
    r0 = 0.05
    prev = None
    for d0 in [0.004, 0.002, 0.001, 0.0005]:
        a = run(thin(r0, d0))
        b = run(thick(r0 - d0 / 2, r0 + d0 / 2, ne=8))
        _, ya = series(a, "radius:ventricle")
        _, yb = series(b, "radius:ventricle")
        _, Pa = series(a, "pressure:inlet_valve:ventricle")
        _, Pb = series(b, "pressure:inlet_valve:ventricle")
        # cavity volume change, each from its own cavity radius
        Va = 4 / 3 * np.pi * ((r0 + ya) ** 3 - r0**3)
        Vb = 4 / 3 * np.pi * ((r0 - d0 / 2 + yb) ** 3 - (r0 - d0 / 2) ** 3)
        eV = np.abs(Va - Vb).max() / np.abs(Va).max()
        eP = np.abs(Pa - Pb).max() / np.abs(Pa).max()
        rate = (
            ""
            if prev is None
            else f"   rates dV {np.log2(prev[0] / eV):.2f}, P {np.log2(prev[1] / eP):.2f}"
        )
        print(f"   d0/r0 = {d0 / r0:6.3f}   rel diff dV {eV:.3e}, P {eP:.3e}{rate}")
        prev = (eV, eP)


def check_quadrature():
    print("\n3. QUADRATURE CONVERGENCE (Ri = 0.025, Re = 0.045, full dynamics)")
    Ri, Re = 0.025, 0.045
    ref = run(thick(Ri, Re, ne=64))
    _, Vref = series(ref, "volume:ventricle")
    _, Pref = series(ref, "pressure:inlet_valve:ventricle")
    for ne in [1, 2, 4, 8, 16]:
        r = run(thick(Ri, Re, ne=ne))
        _, V = series(r, "volume:ventricle")
        _, P = series(r, "pressure:inlet_valve:ventricle")
        print(
            f"   num_elements = {ne:3d}   rel diff to 64:  V {np.abs(V - Vref).max() / np.abs(Vref).max():.3e}"
            f"   P {np.abs(P - Pref).max() / np.abs(Pref).max():.3e}"
        )


def check_transmural():
    print("\n4. TRANSMURAL HETEROGENEITY (Ri = 0.025, Re = 0.045)")
    Ri, Re = 0.025, 0.045
    r = run(thick(Ri, Re, ne=4))
    _, rho_i = series(r, "radius:ventricle")
    for label, k in [("most dilated", int(np.argmax(rho_i)))]:
        ri = Ri + rho_i[k]
        lam = [np.cbrt(R**3 + ri**3 - Ri**3) / R for R in (Ri, 0.5 * (Ri + Re), Re)]
        print(
            f"   {label}: lambda endo {lam[0]:.4f}, mid {lam[1]:.4f}, epi {lam[2]:.4f}"
            f"   (ChamberSphere has one value for the whole wall)"
        )


if __name__ == "__main__":
    sys.exit(
        [check_static(), check_thin_limit(), check_quadrature(), check_transmural()]
        and None
    )
