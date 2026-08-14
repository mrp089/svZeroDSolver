# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
# University of California, and others. SPDX-License-Identifier: BSD-3-Clause
"""Spatial (mesh) and temporal (time-step) convergence of the ChamberCylinder
Genet reproduction. Targets: peak systolic pressure and minimum volume (ESV).

Runs the full mixed-u/p beat + circulation (ChamberCylinder_genet_exact.build)
sweeping the number of through-wall elements `ne` (npts fixed) and the number of
time points per cycle `npts` (ne fixed), each to the Windkessel limit cycle. The
mixed-u/p wall is ~2nd-order and near mesh-independent; the L-stable stiff
integrator is temporally converged by dt = 2 ms (halving moves peak P ~0.05 Pa,
ESV ~0.01 mL). Coarser dt (>=4 ms) and finer meshes (ne>=16) hit stiff-solver
robustness limits, so the usable window is ne<=12, dt<=2 ms, whose coarsest
corner (ne=3, dt=2 ms) is within ~0.5% of the extrapolated limit.

Run:  PYTHONPATH=build/python python scripts/ChamberCylinder_convergence.py
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
G_PK, G_MV = 12.7989, 74.0037  # Genet Fig 5 (digitized) targets
NCYCLE = 15                    # reach the Windkessel limit cycle


def metrics(ne, npts):
    cfg = G.build(bcs_alpha=12.0, ne=ne, ncycle=NCYCLE, P_vs=G.PVS_PHYSIOBLOCKS)
    cfg["simulation_parameters"]["number_of_time_pts_per_cardiac_cycle"] = npts
    cfg["simulation_parameters"]["absolute_tolerance"] = 1e-9
    out = pysvzerod.simulate(cfg)
    P = out[out["name"] == "pressure:ventricle:aortic"]["y"].to_numpy() * kPa
    V = out[out["name"] == "volume:ventricle"]["y"].to_numpy() * 1e6
    return P.max(), V.min()


def sweep(label, cases, run):
    xs, pk, mv = [], [], []
    print(f"{label}:")
    for c in cases:
        try:
            p, m = run(c)
            xs.append(c); pk.append(p); mv.append(m)
            print(f"  {c:>7}: peakP={p:.4f} kPa  minV(ESV)={m:.4f} mL")
        except Exception as e:  # stiff-solver robustness edge
            print(f"  {c:>7}: FAIL ({str(e)[:40]})")
    return np.array(xs, float), np.array(pk), np.array(mv)


if __name__ == "__main__":
    ne, ne_pk, ne_mv = sweep("SPATIAL (npts=400, dt=2ms)",
                             [3, 4, 6, 8, 12, 16, 24],
                             lambda ne: metrics(ne, 400))
    npts, np_pk, np_mv = sweep("TEMPORAL (ne=12)", [200, 400, 800, 1600],
                               lambda n: metrics(12, n))
    dt = 800.0 / npts

    fig, ax = plt.subplots(2, 2, figsize=(12, 8))
    ax[0, 0].plot(1.0 / ne, ne_pk, 'C0o-'); ax[0, 0].axhline(G_PK, ls=':', c='k')
    ax[0, 0].set(title="spatial: peak pressure", xlabel="h = 1/ne", ylabel="peak P [kPa]")
    ax[0, 1].plot(1.0 / ne, ne_mv, 'C1o-'); ax[0, 1].axhline(G_MV, ls=':', c='k')
    ax[0, 1].set(title="spatial: min volume (ESV)", xlabel="h = 1/ne", ylabel="min V [mL]")
    ax[1, 0].plot(dt, np_pk, 'C0s-'); ax[1, 0].axhline(G_PK, ls=':', c='k')
    ax[1, 0].set(title="temporal: peak pressure", xlabel="dt [ms]", ylabel="peak P [kPa]")
    ax[1, 1].plot(dt, np_mv, 'C1s-'); ax[1, 1].axhline(G_MV, ls=':', c='k')
    ax[1, 1].set(title="temporal: min volume (ESV)", xlabel="dt [ms]", ylabel="min V [mL]")
    for a in ax.flat:
        a.grid(alpha=0.3)
    fig.suptitle("ChamberCylinder convergence (mixed u/p): dotted = Genet Fig 5")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig("convergence.png", dpi=120)
    print("saved convergence.png")
