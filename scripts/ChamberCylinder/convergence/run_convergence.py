# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
# University of California, and others. SPDX-License-Identifier: BSD-3-Clause
"""Convergence, cost and robustness of ChamberCylinder in the genet23 setup.

Runs the genet23 baseline (ChamberCylinder_genet_exact.build, fifth beat) over
the numerical parameters a user chooses: the number of elements through the
wall (ne), the time step (dt), the integrator and the Newton tolerance. Errors
are those of the fifth beat against a finer run of the same setup, evaluated at
the run's own output times.

Usage, from the repository root, with the Python module of a build on the
PYTHONPATH (the published numbers use the block branch, which has exact
tangents; this branch's complex-step build gives the same solutions, slower):

  python scripts/ChamberCylinder/convergence/run_convergence.py run [--workers 4]
  python scripts/ChamberCylinder/convergence/run_convergence.py cost
  python scripts/ChamberCylinder/convergence/run_convergence.py report

"run" caches every simulation in scripts/ChamberCylinder/output/convergence/
(not tracked), "cost" times one beat of each configuration serially, and
"report" writes results.json and the figures in this folder.
"""
import argparse
import json
import multiprocessing as mp
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
KIT = os.path.dirname(HERE)
sys.path.insert(0, KIT)
import ChamberCylinder_genet_exact as G  # noqa: E402

OUT = os.path.join(KIT, "output", "convergence")
PERIOD = 0.8  # [s]
MMHG = 133.322387415
RI, RE, LEN = G.GEOM["Ri"], G.GEOM["Re"], G.GEOM["length"]

# Integrators: "stiff" is generalized-alpha with rho_infty = 0 plus a line
# search and sub-stepping; "genalpha" is the plain scheme.
INTEGRATORS = {"stiff": ("stiff", 0.0), "ga05": ("genalpha", 0.5),
               "ga0": ("genalpha", 0.0)}
LABELS = {"stiff": "stiff", "ga05": r"gen.-$\alpha$, $\rho_\infty=0.5$",
          "ga0": r"gen.-$\alpha$, $\rho_\infty=0$"}
ATOL = 1e-7  # tolerance of the studies (1e-9 fails for ne >= 16, see README)
NE_GRID = [1, 2, 3, 4, 6, 8, 12, 16]
DT_GRID = [4.0, 2.0, 1.0, 0.5, 0.25]
NE_SPACE = NE_GRID + [24, 32, 48, 64]
DT_TIME = DT_GRID + [0.125, 0.0625]
ATOLS = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10]
# Waveforms of the last beat whose errors are reported
WAVES = {"V": "cavity volume", "P": "cavity pressure", "tw": "twist",
         "h": "wall thickness", "Q": "aortic flow"}
STATE = ("V", "P", "tw", "h")  # waveforms of the summary errors "max" and "rms"


def spec(integ, ne, dt, atol=ATOL, ncycle=G.NCYCLE, all_cycles=False):
    """One simulation: integrator key, elements, time step [ms], tolerance."""
    return dict(integ=integ, ne=ne, dt=dt, atol=atol, ncycle=ncycle,
                all_cycles=all_cycles)


def key(s):
    k = f"{s['integ']}_ne{s['ne']}_dt{s['dt']:g}_atol{s['atol']:.0e}_nc{s['ncycle']}"
    return k + ("_all" if s["all_cycles"] else "")


def studies():
    """The accuracy studies (run in parallel)."""
    S = {}
    # full grid of the cheap configurations, all integrators
    S["grid"] = [spec(i, ne, dt) for i in INTEGRATORS for ne in NE_GRID for dt in DT_GRID]
    # reference of the grid: fine in space and time
    S["reference"] = [spec(i, 32, 0.125) for i in ("stiff", "ga05")]
    # spatial series at a fixed step, against ne = 64
    S["space"] = [spec("stiff", ne, 0.5) for ne in NE_SPACE]
    # temporal series at fixed meshes, against dt = 0.0625 ms
    S["time"] = [spec(i, ne, dt) for i, ne in (("stiff", 8), ("ga05", 8), ("stiff", 16))
                 for dt in DT_TIME]
    # Newton tolerance, against the tightest tolerance that runs
    S["tol"] = [spec(i, ne, 1.0, atol=a) for i in ("stiff", "ga05") for ne in (8, 32)
                for a in ATOLS]
    # beats to a periodic state
    S["periodic"] = [spec("stiff", 8, 1.0, ncycle=30, all_cycles=True)]
    return S


def configure(s):
    """The svZeroDSolver configuration of a run."""
    integrator, rho = INTEGRATORS[s["integ"]]
    cfg = G.build(ne=s["ne"], ncycle=s["ncycle"], integrator=integrator, rho_infty=rho)
    sp = cfg["simulation_parameters"]
    sp["number_of_time_pts_per_cardiac_cycle"] = int(round(PERIOD * 1e3 / s["dt"])) + 1
    sp["absolute_tolerance"] = s["atol"]
    # output every 1 ms (every step for coarser steps), so all runs share times
    sp["output_interval"] = max(1, int(round(1.0 / s["dt"])))
    sp["output_all_cycles"] = s["all_cycles"]
    return cfg


def simulate(s):
    """Run one configuration; return its waveforms (or the error) and timing."""
    import svzerod
    cfg = configure(s)
    t0 = time.perf_counter()
    try:
        out = svzerod.simulate(cfg)
    except Exception as e:  # noqa: BLE001 - a failure is a result here
        return dict(s, ok=False, error=str(e), elapsed=time.perf_counter() - t0)
    elapsed = time.perf_counter() - t0

    def get(name):
        d = out[out["name"] == name].sort_values("time")
        return d["time"].to_numpy(), d["y"].to_numpy()

    t, V = get("volume:ventricle")
    res = dict(t=t, V=V * 1e6, P=get("pressure:ventricle:aortic")[1] / MMHG,
               Q=get("flow:ventricle:aortic")[1] * 1e6,
               tw=get("beta:ventricle")[1] * LEN * 180 / np.pi,
               h=((RE + get(f"rho_{s['ne']}:ventricle")[1])
                  - (RI + get("rho_0:ventricle")[1])) * 1e3)
    return dict(s, ok=True, elapsed=elapsed, **res)


def path(s):
    return os.path.join(OUT, key(s) + ".npz")


def save(r):
    os.makedirs(OUT, exist_ok=True)
    meta = {k: r[k] for k in ("integ", "ne", "dt", "atol", "ncycle", "all_cycles",
                              "ok", "elapsed")}
    meta["error"] = r.get("error", "")
    arrays = {k: r[k] for k in ("t",) + tuple(WAVES) if k in r}
    np.savez_compressed(path(r), meta=json.dumps(meta), **arrays)


def load(s):
    """Cached result of a run, or None."""
    if not os.path.exists(path(s)):
        return None
    d = np.load(path(s))
    r = json.loads(str(d["meta"]))
    r.update({k: d[k] for k in d.files if k != "meta"})
    return r


def cost_estimate(s):
    return s["ncycle"] * PERIOD / s["dt"] * s["ne"] ** 2


def cmd_run(workers):
    todo = {key(s): s for ss in studies().values() for s in ss if load(s) is None}
    order = sorted(todo.values(), key=cost_estimate, reverse=True)
    print(f"{len(order)} runs, {workers} workers", flush=True)
    with mp.get_context("spawn").Pool(workers) as pool:
        for r in pool.imap_unordered(simulate, order):
            save(r)
            msg = "ok" if r["ok"] else r["error"][:70]
            print(f"{key(r):40s} {r['elapsed']:8.1f} s  {msg}", flush=True)


# --- timing -------------------------------------------------------------------
def cost_specs():
    """One beat per configuration: the grid for all integrators, the large
    meshes at 1 ms."""
    S = [spec(i, ne, dt, ncycle=1) for i in INTEGRATORS for ne in NE_GRID for dt in DT_GRID]
    S += [spec(i, ne, 1.0, ncycle=1) for i in INTEGRATORS for ne in NE_SPACE[len(NE_GRID):]]
    return S


def cmd_cost():
    """Time one beat of each configuration, serially (best of 3 below 2 s)."""
    import svzerod
    fn = os.path.join(OUT, "cost.json")
    done = json.load(open(fn)) if os.path.exists(fn) else {}
    for s in cost_specs():
        k = key(s)
        if k in done:
            continue
        cfg = configure(s)
        times, err = [], ""
        for _ in range(3):
            t0 = time.perf_counter()
            try:
                svzerod.simulate(cfg)
            except Exception as e:  # noqa: BLE001
                err = str(e)
            times.append(time.perf_counter() - t0)
            if err or times[-1] > 2.0:
                break
        done[k] = dict(s, seconds=min(times), ok=not err, error=err)
        print(f"{k:40s} {min(times):8.3f} s {err[:60]}", flush=True)
        with open(fn, "w") as f:
            json.dump(done, f, indent=1)


# --- analysis -----------------------------------------------------------------
def last_beat(r):
    """Waveforms of the last beat, times in integer microseconds from its start."""
    t = r["t"]
    m = t >= t.max() - PERIOD - 1e-9
    us = np.round((t[m] - t[m].min()) * 1e6).astype(np.int64)
    return us, {k: r[k][m] for k in WAVES}


def scalars(w):
    V, P, tw = w["V"], w["P"], w["tw"]
    edv, esv = V.max(), V.min()
    # stroke work: area of the PV loop
    sw = abs(np.sum(0.5 * (P[1:] + P[:-1]) * np.diff(V)))
    return dict(EDV=edv, ESV=esv, EF=(edv - esv) / edv, Pmax=P.max(), twmax=tw.max(),
                SW=sw, hmax=w["h"].max())


def errors(r, ref):
    """Errors of a run against a reference at the run's output times: waveforms
    max |x - x_ref| / range(x_ref) (key k) and RMS(x - x_ref) / range(x_ref)
    (key k_rms), scalars relative. "max" and "rms" are the largest of the
    volume, pressure, twist and wall thickness errors."""
    t, w = last_beat(r)
    tr, wr = last_beat(ref)
    common, i, j = np.intersect1d(t, tr, return_indices=True)
    assert len(common) == len(t), "reference lacks the run's output times"
    e = {}
    for k in WAVES:
        x, xr = w[k][i], wr[k][j]
        e[k] = float(np.max(np.abs(x - xr)) / np.ptp(xr))
        e[k + "_rms"] = float(np.sqrt(np.mean((x - xr) ** 2)) / np.ptp(xr))
    s = scalars({k: w[k][i] for k in WAVES})
    sr = scalars({k: wr[k][j] for k in WAVES})
    for k in s:
        e[k] = float(abs(s[k] - sr[k]) / abs(sr[k]))
    e["max"] = max(e[k] for k in STATE)
    e["rms"] = max(e[k + "_rms"] for k in STATE)
    return e


def periodicity(r):
    """Per beat: EDV, ESV, peak pressure, peak twist and the largest change of the
    waveforms from the previous beat (relative to their range)."""
    t = r["t"]
    beats = []
    prev = None
    for n in range(int(round(t.max() / PERIOD))):
        m = (t >= n * PERIOD - 1e-9) & (t <= (n + 1) * PERIOD + 1e-9)
        w = {k: r[k][m] for k in WAVES}
        row = dict(beat=n + 1, **{k: float(v) for k, v in scalars(w).items()})
        if prev is not None and len(prev["V"]) == len(w["V"]):
            row["change"] = max(float(np.max(np.abs(w[k] - prev[k])) / np.ptp(w[k]))
                                for k in STATE)
        beats.append(row)
        prev = w
    return beats


def cmd_report():
    import convergence_report  # noqa: F401 - figures and tables
    convergence_report.main(sys.modules[__name__])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("command", choices=("run", "cost", "report"))
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    if a.command == "run":
        cmd_run(a.workers)
    elif a.command == "cost":
        cmd_cost()
    else:
        sys.path.insert(0, HERE)
        cmd_report()


if __name__ == "__main__":
    main()
