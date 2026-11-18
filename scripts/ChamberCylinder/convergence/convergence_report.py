# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
# University of California, and others. SPDX-License-Identifier: BSD-3-Clause
"""Tables and figures of the ChamberCylinder convergence study.

Called by `run_convergence.py report`. Writes results.json and
figures/{convergence,cost,waveforms,valve_opening,periodicity}.{pdf,svg} (SVG
for display on GitHub) and prints the tables of the README.
"""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
# Okabe-Ito, fixed order; every series also has its own marker
COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00"]
MARKERS = "os^Dv"
TARGETS = (1e-2, 1e-3)  # accuracy targets of the recommendation
METRICS = {"max": "max error / range", "rms": "RMS error / range"}
KW = dict(lw=1.5, ms=5.5, mfc="white", mew=1.3)


def fit_order(x, e):
    """Least-squares slope of log e over log x."""
    x, e = np.asarray(x, float), np.asarray(e, float)
    m = e > 0
    return float(np.polyfit(np.log(x[m]), np.log(e[m]), 1)[0]) if m.sum() >= 2 else None


def collect(RC):
    S = RC.studies()
    res = {}

    def get(s):
        r = RC.load(s)
        if r is None:
            raise SystemExit(f"missing run {RC.key(s)}: run `run_convergence.py run`")
        return r

    # grid against the fine reference (stiff; generalized-alpha agrees to 8e-4)
    ref = get(RC.spec("stiff", 32, 0.125))
    ref_ga = get(RC.spec("ga05", 32, 0.125))
    res["reference"] = dict(ga05_vs_stiff=RC.errors(ref_ga, ref),
                            scalars=RC.scalars(RC.last_beat(ref)[1]))
    grid = []
    for s in S["grid"]:
        r = get(s)
        row = dict(integ=s["integ"], ne=s["ne"], dt=s["dt"], ok=bool(r["ok"]),
                   error_msg=r.get("error", ""))
        if r["ok"]:
            row["errors"] = RC.errors(r, ref)
        grid.append(row)
    res["grid"] = grid

    # spatial series against ne = 64
    sref = get(RC.spec("stiff", 64, 0.5))
    res["space"] = [dict(ne=s["ne"], errors=RC.errors(get(s), sref)) for s in S["space"][:-1]]

    # temporal series against dt = 0.0625 ms
    res["time"] = {}
    for integ, ne in (("stiff", 8), ("ga05", 8), ("stiff", 16)):
        tref = get(RC.spec(integ, ne, RC.DT_TIME[-1]))
        res["time"][f"{integ}_ne{ne}"] = [
            dict(dt=dt, errors=RC.errors(get(RC.spec(integ, ne, dt)), tref))
            for dt in RC.DT_TIME[:-1]]

    # tolerance against the tightest tolerance that runs
    tol = {}
    for integ in ("stiff", "ga05"):
        for ne in (8, 32):
            runs = [(a, get(RC.spec(integ, ne, 1.0, atol=a))) for a in RC.ATOLS]
            a_ref, r_ref = min(((a, r) for a, r in runs if r["ok"]), key=lambda ar: ar[0])
            tol[f"{integ}_ne{ne}"] = dict(reference_atol=a_ref, rows=[
                dict(atol=a, ok=bool(r["ok"]),
                     errors=RC.errors(r, r_ref) if r["ok"] and a != a_ref else None,
                     error_msg=r.get("error", ""))
                for a, r in runs])
    res["tolerance"] = tol

    res["periodic"] = RC.periodicity(get(S["periodic"][0]))

    # serial timing of one beat
    fn = os.path.join(RC.OUT, "cost.json")
    cost = json.load(open(fn)) if os.path.exists(fn) else {}
    res["cost"] = [dict(integ=c["integ"], ne=c["ne"], dt=c["dt"], seconds=c["seconds"],
                        ok=c["ok"]) for c in cost.values()]
    return res


def cost_of(res, integ, ne, dt):
    for c in res["cost"]:
        if (c["integ"], c["ne"], c["dt"]) == (integ, ne, dt) and c["ok"]:
            return c["seconds"]
    return None


def recommend(res):
    """Cheapest grid configuration per integrator, metric and accuracy target."""
    out = []
    for integ in ("stiff", "ga0", "ga05"):
        for metric in METRICS:
            for tgt in TARGETS:
                cands = [(cost_of(res, integ, g["ne"], g["dt"]), g["ne"], g["dt"],
                          g["errors"][metric]) for g in res["grid"]
                         if g["integ"] == integ and g["ok"] and g["errors"][metric] <= tgt]
                cands = [c for c in cands if c[0] is not None]
                if cands:
                    c, ne, dt, e = min(cands)
                    out.append(dict(integ=integ, metric=metric, target=tgt, ne=ne, dt=dt,
                                    error=e, seconds_per_beat=c))
    return out


# --- figures ------------------------------------------------------------------
def style(a, xlabel, ylabel, title):
    a.set(xlabel=xlabel, ylabel=ylabel)
    a.set_title(title, fontsize=10.5)
    a.grid(True, which="major", color="0.9", lw=0.6)
    a.minorticks_off()


def slope(a, x0, x1, e0, p, label, below=True):
    """Guide line e0 (x/x0)^p from x0 to x1, labeled at its middle."""
    x = np.array([x0, x1])
    a.loglog(x, e0 * (x / x0) ** p, ":", color="0.45", lw=1.2)
    xm = np.sqrt(x0 * x1)
    a.annotate(label, (xm, e0 * (xm / x0) ** p), color="0.35", fontsize=9, ha="center",
               textcoords="offset points", xytext=(0, -13 if below else 5))


def pareto(points):
    """Points (cost, error, label) not beaten by a cheaper, more accurate one."""
    front, best = [], np.inf
    for c, e, lab in sorted(points):
        if e < best:
            front.append((c, e, lab))
            best = e
    return front


def fig_convergence(plt, RC, res):
    fig, axs = plt.subplots(2, 3, figsize=(15.5, 8.4), sharey="row")
    for row, (metric, ylabel) in enumerate(METRICS.items()):
        suffix = "" if metric == "max" else "_rms"
        # space
        a = axs[row, 0]
        ne = [r["ne"] for r in res["space"]]
        for i, (k, lab) in enumerate(RC.WAVES.items()):
            a.loglog(ne, [r["errors"][k + suffix] for r in res["space"]], "-",
                     color=COLORS[i], marker=MARKERS[i], label=lab, **KW)
        slope(a, 4, 32, 1e-3 if metric == "max" else 5e-4, -2, r"$n_e^{-2}$")
        a.set_xticks(ne, [str(n) for n in ne], fontsize=8)
        style(a, "elements through the wall $n_e$", ylabel,
              rf"({'ad'[row]}) Space: $\Delta t$ = 0.5 ms, vs. $n_e$ = 64")
        # time
        for col, (key_, title) in enumerate((
                ("stiff_ne8", "Time: stiff"),
                ("ga05_ne8", r"Time: gen.-$\alpha$, $\rho_\infty$ = 0.5")), start=1):
            a = axs[row, col]
            rows = res["time"][key_]
            dt = [r["dt"] for r in rows]
            for i, (k, lab) in enumerate(RC.WAVES.items()):
                a.loglog(dt, [r["errors"][k + suffix] for r in rows], "-",
                         color=COLORS[i], marker=MARKERS[i], label=lab, **KW)
            slope(a, 2.0, 0.25, 4e-4 if metric == "max" else 6e-5, 2, r"$\Delta t^2$")
            a.set_xticks(dt, [f"{d:g}" for d in dt], fontsize=8)
            style(a, r"time step $\Delta t$ [ms]", ylabel,
                  f"({'abcdef'[3 * row + col]}) {title}, $n_e$ = 8, vs. 0.0625 ms")
    axs[0, 0].legend(fontsize=8, frameon=False, loc="lower left")
    fig.tight_layout()
    return fig


def fig_cost(plt, RC, res):
    fig, axs = plt.subplots(1, 3, figsize=(16.0, 4.5))
    integs = ("stiff", "ga05", "ga0")
    a = axs[0]
    for i, integ in enumerate(integs):
        rows = sorted((c["ne"], c["seconds"]) for c in res["cost"]
                      if c["integ"] == integ and c["dt"] == 1.0 and c["ok"])
        a.loglog(*zip(*rows), "-", color=COLORS[i], marker=MARKERS[i],
                 label=RC.LABELS[integ], **KW)
    c2 = cost_of(res, "ga0", 2, 1.0)
    slope(a, 2, 64, 0.5 * c2, 1, r"$\propto n_e$")
    a.set_xticks(RC.NE_SPACE, [str(n) for n in RC.NE_SPACE], fontsize=8)
    style(a, "elements through the wall $n_e$", "run time per beat [s]",
          r"(a) Cost: $\Delta t$ = 1 ms (800 steps per beat)")
    a.legend(fontsize=8, frameon=False, loc="upper left")
    for a, (metric, ylabel) in zip(axs[1:], METRICS.items()):
        for i, integ in enumerate(integs):
            pts = [(cost_of(res, integ, g["ne"], g["dt"]), g["errors"][metric],
                    f"{g['ne']}/{g['dt']:g}")
                   for g in res["grid"] if g["integ"] == integ and g["ok"]]
            pts = [p for p in pts if p[0] is not None]
            c, e, _ = zip(*pts)
            a.loglog(c, e, MARKERS[i], color=COLORS[i], ms=3.5, alpha=0.35, ls="none")
            front = pareto(pts)
            fc, fe, fl = zip(*front)
            a.loglog(fc, fe, "-", color=COLORS[i], marker=MARKERS[i],
                     label=RC.LABELS[integ], **KW)
            if integ == "ga0":
                for x, y, lab in front:
                    a.annotate(lab, (x, y), textcoords="offset points", xytext=(4, -10),
                               fontsize=7, color="0.25")
        for tgt in TARGETS:
            a.axhline(tgt, color="0.6", lw=0.8, ls="--")
        style(a, "run time per beat [s]", ylabel.replace("error", "error (V, P, twist, h)"),
              f"({'bc'[list(METRICS).index(metric)]}) Best accuracy for the cost "
              r"($n_e$/$\Delta t$ [ms])")
        a.legend(fontsize=8, frameon=False, loc="upper right")
    fig.tight_layout()
    return fig


def fig_waveforms(plt, RC, res):
    fig, axs = plt.subplots(1, 3, figsize=(15.0, 4.3))
    ref = RC.load(RC.spec("stiff", 32, 0.125))
    tr, wr = RC.last_beat(ref)
    axs[0].plot(wr["V"], wr["P"], "k-", lw=2.4, label="reference")
    axs[1].plot(tr / 1e3, wr["tw"], "k-", lw=2.4, label="reference")
    axs[2].plot(tr / 1e3, wr["h"], "k-", lw=2.4, label="reference")
    for i, ne in enumerate((1, 2, 4, 8)):
        t, w = RC.last_beat(RC.load(RC.spec("stiff", ne, 1.0)))
        kw = dict(color=COLORS[i], lw=1.3, ls=(0, (4, 2)), label=f"$n_e$ = {ne}")
        axs[0].plot(w["V"], w["P"], **kw)
        axs[1].plot(t / 1e3, w["tw"], **kw)
        axs[2].plot(t / 1e3, w["h"], **kw)
    style(axs[0], "cavity volume [mL]", "cavity pressure [mmHg]", "(a) PV loop")
    style(axs[1], "time in beat [ms]", r"twist [$\degree$]", "(b) Twist")
    style(axs[2], "time in beat [ms]", "wall thickness [mm]", "(c) Wall thickness")
    axs[0].legend(fontsize=8, frameon=False, loc="center")
    fig.suptitle(r"Stiff integrator, $\Delta t$ = 1 ms; reference $n_e$ = 32, "
                 r"$\Delta t$ = 0.125 ms (fifth beat)", fontsize=10)
    fig.tight_layout()
    return fig


def fig_valve_opening(plt, RC, res):
    fig, axs = plt.subplots(1, 2, figsize=(11.0, 4.3), sharey=True)
    for a, integ in zip(axs, ("stiff", "ga05")):
        for i, dt in enumerate((2.0, 1.0, 0.5)):
            t, w = RC.last_beat(RC.load(RC.spec(integ, 8, dt)))
            m = (t >= 155e3) & (t <= 185e3)
            a.plot(t[m] / 1e3, w["Q"][m], "-", color=COLORS[i], marker=MARKERS[i], ms=4,
                   mfc="white", mew=1.1, lw=1.2, label=f"{dt:g} ms")
        t, w = RC.last_beat(RC.load(RC.spec(integ, 8, 0.0625)))
        m = (t >= 155e3) & (t <= 185e3)
        a.plot(t[m] / 1e3, w["Q"][m], "k-", lw=1.8, label="0.0625 ms")
        style(a, "time in beat [ms]", "aortic flow [mL/s]",
              f"({'ab'[integ != 'stiff']}) {RC.LABELS[integ]}, $n_e$ = 8")
        a.legend(fontsize=8, frameon=False, loc="upper left")
    fig.suptitle("Aortic valve opening: the abrupt valve switch excites a fast wall "
                 "transient", fontsize=10)
    fig.tight_layout()
    return fig


def fig_periodicity(plt, RC, res):
    fig, axs = plt.subplots(1, 2, figsize=(10.0, 4.0))
    rows = res["periodic"]
    axs[0].semilogy([r["beat"] for r in rows[1:]], [r["change"] for r in rows[1:]], "-",
                    color=COLORS[0], marker=MARKERS[0], **KW)
    for tgt in TARGETS:
        axs[0].axhline(tgt, color="0.6", lw=0.8, ls="--")
    style(axs[0], "beat", "max change from previous beat / range (V, P, twist, h)",
          r"(a) Approach to the periodic state: $n_e$ = 8, $\Delta t$ = 1 ms")
    # the end-diastolic volume is set by the prescribed atrial pressure and is
    # constant from the first beat
    esv = [r["ESV"] for r in rows]
    axs[1].plot([r["beat"] for r in rows], esv, "-", color=COLORS[1], marker=MARKERS[1], **KW)
    axs[1].axhline(esv[-1], color="0.6", lw=0.8, ls="--")
    axs[1].annotate(f"beat 5 (as in genet23): {esv[4]:.2f} mL", (5, esv[4]),
                    textcoords="offset points", xytext=(10, -14), fontsize=8, color="0.25")
    axs[1].annotate(f"periodic: {esv[-1]:.2f} mL", (30, esv[-1]), ha="right",
                    textcoords="offset points", xytext=(0, -14), fontsize=8, color="0.25")
    style(axs[1], "beat", "end-systolic volume [mL]", "(b) End-systolic volume")
    fig.tight_layout()
    return fig


def save(fig, name):
    os.makedirs(FIG, exist_ok=True)
    for ext in ("pdf", "svg"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), bbox_inches="tight")


# --- tables -------------------------------------------------------------------
def fmt(e):
    return f"{e:.1e}" if e is not None else "fail"


def print_tables(RC, res):
    for metric in METRICS:
        print(f"\n## Grid: {metric} error (V, P, twist, h) vs ne 32, dt 0.125 ms")
        for integ in ("stiff", "ga05", "ga0"):
            print(f"\n[{integ}]\n| ne | " + " | ".join(f"{d:g} ms" for d in RC.DT_GRID) + " |")
            for ne in RC.NE_GRID:
                cells = []
                for dt in RC.DT_GRID:
                    g = next(g for g in res["grid"]
                             if (g["integ"], g["ne"], g["dt"]) == (integ, ne, dt))
                    cells.append(fmt(g["errors"][metric]) if g["ok"] else "fail")
                print(f"| {ne} | " + " | ".join(cells) + " |")
    cols = ["V", "P", "tw", "h", "Q", "V_rms", "P_rms", "tw_rms", "h_rms", "Q_rms",
            "EDV", "ESV", "EF", "Pmax", "twmax", "SW"]
    print("\n## Space (dt 0.5 ms vs ne 64)\n| ne | " + " | ".join(cols) + " |")
    for r in res["space"]:
        print(f"| {r['ne']} | " + " | ".join(fmt(r["errors"][k]) for k in cols) + " |")
    for k, rows in res["time"].items():
        print(f"\n## Time {k} (vs 0.0625 ms)\n| dt | " + " | ".join(cols) + " |")
        for r in rows:
            print(f"| {r['dt']:g} | " + " | ".join(fmt(r["errors"][c]) for c in cols) + " |")
    print("\n## Orders")
    for k, v in res["orders"].items():
        print(k, {q: (round(o, 2) if o is not None else None) for q, o in v.items()})
    print("\n## Tolerance (max, rms vs tightest)")
    for k, v in res["tolerance"].items():
        print(k, "ref", v["reference_atol"], [
            (r["atol"], (fmt(r["errors"]["max"]), fmt(r["errors"]["rms"])) if r["errors"]
             else ("ref" if r["ok"] else "fail")) for r in v["rows"]])
    print("\n## Periodicity: beat EDV ESV Pmax twmax change")
    for r in res["periodic"]:
        print(r["beat"], round(r["EDV"], 3), round(r["ESV"], 3), round(r["Pmax"], 3),
              round(r["twmax"], 3), fmt(r.get("change")))
    print("\n## Recommendation")
    for v in res["recommend"]:
        print(v)
    print("\n## Reference: ga05 vs stiff at ne 32, 0.125 ms:",
          {k: fmt(v) for k, v in res["reference"]["ga05_vs_stiff"].items()})
    print("\n## Cost (s per beat)")
    for c in sorted(res["cost"], key=lambda c: (c["integ"], c["ne"], -c["dt"])):
        print(c["integ"], c["ne"], c["dt"], round(c["seconds"], 3), c["ok"])


def main(RC):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9.5, "svg.fonttype": "none"})
    res = collect(RC)
    res["orders"] = {}
    rows = [r for r in res["space"] if 4 <= r["ne"] <= 32]
    res["orders"]["space_ne4-32"] = {
        k: fit_order([r["ne"] for r in rows], [r["errors"][k] for r in rows])
        for k in list(RC.WAVES) + [k + "_rms" for k in RC.WAVES]}
    for key_, rows in res["time"].items():
        rows = [r for r in rows if r["dt"] <= 1.0]
        res["orders"][f"time_{key_}_dt<=1"] = {
            k: fit_order([r["dt"] for r in rows], [r["errors"][k] for r in rows])
            for k in list(RC.WAVES) + [k + "_rms" for k in RC.WAVES]}
    res["recommend"] = recommend(res)
    with open(os.path.join(HERE, "results.json"), "w") as f:
        json.dump(res, f, indent=1, default=float)
    figs = [("convergence", fig_convergence), ("waveforms", fig_waveforms),
            ("valve_opening", fig_valve_opening), ("periodicity", fig_periodicity)]
    if res["cost"]:
        figs.append(("cost", fig_cost))
    for name, fn in figs:
        fig = fn(plt, RC, res)
        save(fig, name)
        plt.close(fig)
    print_tables(RC, res)
