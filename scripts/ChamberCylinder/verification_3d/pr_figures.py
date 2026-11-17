"""Key figures of the 0D vs 3D verification, for the pull request.

Writes figures/convergence.{pdf,svg} and figures/pv_loop.{pdf,svg} from the
results of run_study.py and run_active_study.py (CASE/output/summary.json,
CASE/runs/*/coupled.json). SVG is for display on GitHub.

Usage: python pr_figures.py
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
OUT = os.path.join(HERE, "figures")
MMHG = 133.322387415
# Okabe-Ito, colorblind-safe, fixed order
COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"]


def levels(case):
    """3D levels of a study, coarse to fine: nr, h [mm], errors."""
    with open(os.path.join(HERE, case, "output", "summary.json")) as f:
        s = json.load(f)
    return [dict(r, nr=int(k[2:]), h_mm=r["h"] * 1e3) for k, r in
            sorted(s["three_d"].items(), key=lambda kv: -kv[1]["h"])], s


def convergence(plt):
    fig, axs = plt.subplots(1, 3, figsize=(13.5, 4.4))
    style = dict(lw=1.6, ms=6, mfc="white", mew=1.4)

    def curve(a, case, key, label, color, marker, ls="-"):
        lev = [r for r in levels(case)[0] if r["errors"].get(key)]
        a.loglog([r["h_mm"] for r in lev], [r["errors"][key] for r in lev], ls,
                 marker=marker, color=color, label=label, **style)

    def slope(a, h0, e0, p):
        h = np.array([h0, h0 / 4])
        a.loglog(h, e0 * (h / h0) ** p, ":", color="0.45", lw=1.2)
        a.annotate(f"$h^{p}$", (h[1], e0 * (h[1] / h0) ** p), color="0.35",
                   textcoords="offset points", xytext=(-4, -12), fontsize=9)

    # (a) passive inflation, circumferential fibers: cavity volume
    a = axs[0]
    walls = [("passive_isotropic", "neo-Hookean"),
             ("passive_anisotropic", "Holzapfel-Ogden, one fiber family"),
             ("passive_anisotropic_full", "Holzapfel-Ogden, all terms"),
             ("passive_anisotropic_full_closed", "all terms, free length")]
    for i, (case, label) in enumerate(walls):
        curve(a, case, "V", label, COLORS[i], "osdv"[i], "--" if "closed" in case else "-")
    slope(a, 7.0, 2.5e-5, 4)
    a.set(title="(a) Passive inflation: cavity volume")

    # (b) helical fibers, 2 L0 tube, mid-length section
    a = axs[1]
    for i, (key, label) in enumerate((("V", "cavity volume"), ("rho", "radial displacement"),
                                      ("beta", "twist rate"))):
        curve(a, "passive_helical_L2", key, label, COLORS[i], "os^"[i])
    slope(a, 14.0, 3e-3, 3)
    a.set(title=r"(b) Helical fibers $\pm60^\circ$: mid-length section")

    # (c) active contraction, 3D cylinder coupled to the 0D circulation
    a = axs[2]
    for case, ls, steps in (("active_prescribed", "-", 400), ("active_prescribed_dt2", "--", 800)):
        for i, (key, label) in enumerate((("SV", "stroke volume"), ("V", "cavity volume"),
                                          ("P", "cavity pressure"))):
            curve(a, case, key, f"{label}, {steps} steps" if steps == 400 else None,
                  COLORS[i], "o^s"[i], ls)
    a.plot([], [], "k--", label="800 steps per cycle")
    slope(a, 14.0, 2e-3, 3)
    a.set(title="(c) Coupled PV loop, first cycle")

    for a, loc in zip(axs, ("lower right", "lower right", "upper left")):
        a.set(xlabel="mesh size $h = (R_e - R_i)/n_r$ [mm]", ylabel="max relative error, 3D vs 0D")
        a.set_xticks([0.875, 1.75, 3.5, 7.0, 14.0], ["0.875\n(nr 16)", "1.75\n(8)", "3.5\n(4)",
                                                       "7\n(2)", "14\n(1)"])
        a.minorticks_off()
        a.set_xlim(0.7, 18.0)
        a.grid(True, which="major", color="0.9", lw=0.6)
        a.legend(fontsize=8, loc=loc, frameon=False)
        a.title.set_fontsize(10.5)
    fig.tight_layout()
    return fig


def pv_loop(plt):
    import run_active_study as RA
    case = "active_prescribed"
    root = os.path.join(HERE, case)
    lev, summary = levels(case)
    k = summary["cycle"]
    with open(os.path.join(root, "chamber_cylinder_0d.json")) as f:
        period = float(json.load(f)["simulation_parameters"]["cardiac_period"])
    ref, _ = RA.reference(os.path.join(root, "chamber_cylinder_0d.json"),
                          os.path.join(root, "output"))
    fig, axs = plt.subplots(1, 2, figsize=(10.0, 4.2))

    def cycle(t):
        m = RA.cycle_mask(t, period, k)
        return m & (t >= RA.STARTUP) if k == 1 else m

    m = cycle(ref["t"])
    axs[0].plot(ref["V"][m] * 1e6, ref["P"][m] / MMHG, "k-", lw=2.6, label="0D block")
    axs[1].axhline(0.0, color="k", lw=1.0)
    for i, r in enumerate(lev):
        with open(os.path.join(root, "runs", f"struct-p2_nr{r['nr']}", "coupled.json")) as f:
            s = json.load(f)
        t = np.asarray(s["t"])
        c = cycle(t)
        V = np.asarray(s["V"])[c]
        P = np.interp(t, s["t_0d"], s["P"])[c]
        P0 = np.interp(t[c], ref["t"], ref["P"])
        kw = dict(color=COLORS[i], lw=1.4, ls=(0, (4, 2)) if i < len(lev) - 1 else "-")
        axs[0].plot(V * 1e6, P / MMHG, label=f"3D, nr {r['nr']}", **kw)
        axs[1].plot(t[c] - (k - 1) * period, (P - P0) / MMHG, label=f"3D, nr {r['nr']}", **kw)
    axs[0].set(xlabel="cavity volume [mL]", ylabel="cavity pressure [mmHg]",
               title="(a) Pressure-volume loop, 0D block vs. coupled 3D")
    axs[1].set(xlabel="time [s]", ylabel="pressure difference, 3D − 0D [mmHg]",
               title="(b) Cavity pressure difference")
    for a in axs:
        a.grid(True, color="0.9", lw=0.6)
        a.legend(fontsize=8, frameon=False)
        a.title.set_fontsize(10.5)
    fig.tight_layout()
    return fig


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"svg.fonttype": "none", "font.size": 9})
    os.makedirs(OUT, exist_ok=True)
    for name, make in (("convergence", convergence), ("pv_loop", pv_loop)):
        fig = make(plt)
        for ext in ("pdf", "svg"):
            fig.savefig(os.path.join(OUT, f"{name}.{ext}"))
        plt.close(fig)
        print("wrote", os.path.join(OUT, name + ".{pdf,svg}"))


if __name__ == "__main__":
    main()
