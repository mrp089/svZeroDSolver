# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
# University of California, and others. SPDX-License-Identifier: BSD-3-Clause
"""Regenerate the Genet (2023) Fig. 5-9 comparison figures for ChamberCylinder.

Overlays the current model (exact Table-1 parameters, PhysioBlocks n0(e_c),
exact PiecewiseValve + 2-stage Windkessel, force-velocity, C_valve=0) on the
digitized Genet Fig. 5 (pressure/volume/PV-loop/twist) and the sensitivity
Figs. 6/8/9 (twist vs aspect ratio, peak P & twist vs fiber angle, peak P vs
wall volume). Writes fig5.png and sensitivity.png to the working directory.

Run:  PYTHONPATH=build/python python scripts/ChamberCylinder_figures.py
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pysvzerod

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ChamberCylinder_genet_exact as G  # noqa: E402

kPa = 1e-3


def run(geom=None, ne=12, ncycle=8):
    # Mixed-u/p wall (default) + constant P_at=0.9 kPa (no kick needed) + venous
    # P_vs=1.6 kPa afterload floor; activation timing read from Fig. 5.
    cfg = G.build(bcs_alpha=12.0, ne=ne, ncycle=ncycle, P_vs=G.PVS_PHYSIOBLOCKS)
    vv = cfg["vessels"][1]["zero_d_element_values"]
    vv["c_valve"] = 0.0  # exact C_valve over-buffers the peak (see discrepancies doc)
    if geom:
        vv.update(geom)
    cfg["simulation_parameters"]["absolute_tolerance"] = 1e-6
    try:
        out = pysvzerod.simulate(cfg)
    except Exception as e:
        print("  FAIL", str(e)[:60])
        return None

    def s(n):
        d = out[out["name"] == n].sort_values("time")
        return d["time"].to_numpy(), d["y"].to_numpy()

    t, V = s("volume:ventricle")
    _, Pv = s("pressure:ventricle:aortic")
    _, b = s("beta:ventricle")
    L = vv["length"]
    m = t >= t.max() - 0.8 - 1e-9
    return dict(t=(t[m] - t[m].min()) * 1000, V=V[m] * 1e6, P=Pv[m] * kPa,
                tw=b[m] * L * 180 / np.pi, Ppk=Pv[m].max() * kPa,
                EF=(1 - V[m].min() / V[m].max()) * 100,
                twpk=(b[m] * L * 180 / np.pi).max())


# Digitized Genet Fig. 5 (1 ms resolution; see scripts/genet23_fig5_digitized.csv).
def load_fig5():
    """Load the digitized Genet Fig. 5 traces (skips the '#' provenance header)."""
    import io
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "genet23_fig5_digitized.csv")
    with open(path) as f:
        body = "".join(ln for ln in f if not ln.lstrip().startswith("#"))
    return np.genfromtxt(io.StringIO(body), delimiter=",", names=True)


_FIG5 = load_fig5()
pP_t, pP = _FIG5["time_ms"], _FIG5["pressure_kPa"]
pV_t, pV = _FIG5["time_ms"], _FIG5["volume_mL"]
pT_t, pT = _FIG5["time_ms"], _FIG5["torsion_deg"]


def fig5():
    r = run()
    print(f"Fig 5: EDV={r['V'].max():.0f} ESV={r['V'].min():.0f} EF={r['EF']:.0f}% "
          f"Ppk={r['Ppk']:.1f} twist={r['twpk']:.0f}")
    fig, ax = plt.subplots(1, 4, figsize=(20, 4.4))
    ax[0].plot(pP_t, pP, 'k-', lw=2, label="Genet Fig 5 (digitized)")
    ax[0].plot(r['t'], r['P'], 'C3-', lw=1.6, label="model")
    ax[0].set(title="pressure", xlabel="t [ms]", ylabel="P [kPa]"); ax[0].legend()
    ax[1].plot(pV_t, pV, 'k-', lw=2); ax[1].plot(r['t'], r['V'], 'C3-', lw=1.6)
    ax[1].set(title="volume", xlabel="t [ms]", ylabel="V [mL]")
    ax[2].plot(pV, pP, 'k-', lw=2)
    ax[2].plot(r['V'], r['P'], 'C3-', lw=1.6)
    ax[2].set(title="P-V loop", xlabel="V [mL]", ylabel="P [kPa]")
    ax[3].plot(pT_t, pT, 'k-', lw=2); ax[3].plot(r['t'], r['tw'], 'C3-', lw=1.6)
    ax[3].set(title="torsion", xlabel="t [ms]", ylabel="torsion [deg]")
    fig.suptitle("Genet Fig 5 (black, digitized) vs ChamberCylinder mixed u/p (red): "
                 "exact params, constant P_at, exact valves + 2-stage Windkessel, "
                 "force-velocity")
    plt.tight_layout(rect=[0, 0, 1, 0.95]); plt.savefig("fig5.png", dpi=110); plt.close()
    print("saved fig5.png")


def sensitivity():
    Vc = np.pi * 0.019**2 * 0.0571
    Vw = np.pi * 0.0571 * (0.033**2 - 0.019**2)

    def geom_AR(AR):
        Ri = (Vc / np.pi / (2 * AR))**(1 / 3)
        L = 2 * AR * Ri
        return dict(Ri=Ri, Re=np.sqrt(Ri**2 + Vw / np.pi / L), length=L)

    def geom_wall(W):
        return dict(Ri=0.019, Re=np.sqrt(0.019**2 + W / np.pi / 0.0571), length=0.0571)

    ang, ars, wal = [30, 45, 60, 75, 90], [1.0, 1.5, 2.0], [117, 130, 143]
    pFib_tw = [16, 22.5, 20, 13, 8]; pFib_P = [11.5, 12.6, 12.7, 12.5, 11.5]
    pAsp_tw = [13.5, 20, 26.5]; pWal_P = [12.2, 12.8, 13.3]
    fib = [run(dict(alpha_endo=a, alpha_epi=-a), ne=8, ncycle=6) for a in ang]
    asp = [run(geom_AR(a), ne=8, ncycle=6) for a in ars]
    wl = [run(geom_wall(w * 1e-6), ne=8, ncycle=6) for w in wal]
    gv = lambda L, k: [x[k] if x else np.nan for x in L]  # noqa: E731

    fig, ax = plt.subplots(1, 4, figsize=(20, 4.4))
    ax[0].plot(ang, pFib_P, 'ks-', lw=2, label="Genet"); ax[0].plot(ang, gv(fib, 'Ppk'), 'C3o-', label="model")
    ax[0].set(title="Fig 9: peak P vs fiber angle", xlabel="helix +-alpha [deg]", ylabel="peak P [kPa]"); ax[0].legend()
    ax[1].plot(ang, pFib_tw, 'ks-', lw=2); ax[1].plot(ang, gv(fib, 'twpk'), 'C3o-')
    ax[1].set(title="Fig 9: peak twist vs fiber angle", xlabel="helix +-alpha [deg]", ylabel="twist [deg]")
    ax[2].plot(ars, pAsp_tw, 'ks-', lw=2); ax[2].plot(ars, gv(asp, 'twpk'), 'C3o-')
    ax[2].set(title="Fig 6: peak twist vs aspect ratio", xlabel="aspect ratio", ylabel="twist [deg]")
    ax[3].plot(wal, pWal_P, 'ks-', lw=2, label="Genet"); ax[3].plot(wal, gv(wl, 'Ppk'), 'C3o-', label="model")
    ax[3].set(title="Fig 8: peak P vs wall volume", xlabel="wall volume [mL]", ylabel="peak P [kPa]"); ax[3].legend()
    fig.suptitle("Genet sensitivity (black) vs ChamberCylinder (red): PhysioBlocks n0, force-velocity")
    plt.tight_layout(rect=[0, 0, 1, 0.95]); plt.savefig("sensitivity.png", dpi=110); plt.close()
    print("saved sensitivity.png")


if __name__ == "__main__":
    fig5()
    sensitivity()
