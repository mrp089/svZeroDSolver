"""Run a ChamberSphere 0D simulation and render a two-panel animation.

Left panel: pressure-volume loop of the sphere chamber with a moving marker
at the current time. Right panel: the "beating" sphere drawn as a circle
whose radius and wall thickness evolve with the simulation.

- Inner radius is recovered from the simulated sphere volume.
- Wall thickness is recovered from incompressibility and the undeformed
  (radius0, thick0) values declared in the input file.

Usage:
    python chamber_sphere_animate.py INPUT.json [--vessel NAME] [--out FILE]
                                     [--color C] [--fps N]
"""

import argparse
import json
import os
import sys

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
import pysvzerod
from matplotlib.patches import Circle

PA_PER_MMHG = 133.322387415
M3_PER_ML = 1.0e-6


def find_chamber_sphere(config, requested_name=None):
    spheres = [
        v for v in config.get("vessels", [])
        if v.get("zero_d_element_type") == "ChamberSphere"
    ]
    if not spheres:
        raise ValueError("No ChamberSphere vessel found in input file.")
    if requested_name is not None:
        for v in spheres:
            if v["vessel_name"] == requested_name:
                return v
        raise ValueError(f"ChamberSphere '{requested_name}' not found.")
    if len(spheres) > 1:
        names = ", ".join(v["vessel_name"] for v in spheres)
        raise ValueError(
            f"Multiple ChamberSphere vessels present ({names}); "
            "pick one with --vessel."
        )
    return spheres[0]


def extract_series(result, name):
    d = result[result["name"] == name].sort_values("time")
    return d["time"].to_numpy(), d["y"].to_numpy()


def chamber_pressure_name(config, vessel_name):
    """Return a variable name giving the chamber's internal pressure.

    The ChamberSphere residual enforces Pin = Pout, so either node pressure
    adjacent to the chamber carries the chamber pressure. We look for a
    downstream neighbor (valve or vessel) first, falling back to upstream.
    """
    candidates = []
    for valve in config.get("valves", []):
        p = valve.get("params", {})
        if p.get("upstream_block") == vessel_name:
            candidates.append(f"pressure:{vessel_name}:{valve['name']}")
        if p.get("downstream_block") == vessel_name:
            candidates.append(f"pressure:{valve['name']}:{vessel_name}")
    for j in config.get("junctions", []):
        inlets = j.get("inlet_vessels", []) or j.get("inlet_blocks", [])
        outlets = j.get("outlet_vessels", []) or j.get("outlet_blocks", [])
        if vessel_name in inlets:
            candidates.append(f"pressure:{vessel_name}:{j['junction_name']}")
        if vessel_name in outlets:
            candidates.append(f"pressure:{j['junction_name']}:{vessel_name}")
    return candidates


def thickness_from_volume(volume_m3, radius0, thick0):
    """Current wall thickness under the thin-wall incompressibility limit.

    Thin-wall wall volume V_wall = 4 pi R^2 t = const = 4 pi R0^2 t0,
    hence t = t0 (R0 / R)^2. R is recovered from the sphere volume as
    R = (3V / 4 pi)^(1/3).
    """
    radius = np.cbrt(3.0 * volume_m3 / (4.0 * np.pi))
    thickness = thick0 * (radius0 / radius) ** 2
    return radius, thickness


def build_animation(times, pressures_mmHg, volumes_mL, radius, thickness,
                    color, close_loop=True,
                    pv_trace_volumes=None, pv_trace_pressures=None):
    fig, (ax_pv, ax_sphere) = plt.subplots(1, 2, figsize=(12, 5))

    # ---- PV loop ----
    # Draw the trace from the full-resolution arrays if provided so the loop
    # doesn't look polygonal when the animation frames are strided.
    trace_v = pv_trace_volumes if pv_trace_volumes is not None else volumes_mL
    trace_p = pv_trace_pressures if pv_trace_pressures is not None else pressures_mmHg
    if close_loop:
        trace_v = np.append(trace_v, trace_v[0])
        trace_p = np.append(trace_p, trace_p[0])
    ax_pv.plot(trace_v, trace_p, "-", color=color, lw=1.2)
    (marker,) = ax_pv.plot([], [], "o", color="black", ms=8)
    ax_pv.set_xlabel("Volume [mL]")
    ax_pv.set_ylabel("Pressure [mmHg]")
    ax_pv.set_title("Pressure-Volume Loop")
    pad_v = 0.05 * (volumes_mL.max() - volumes_mL.min() + 1e-12)
    pad_p = 0.05 * (pressures_mmHg.max() - pressures_mmHg.min() + 1e-12)
    ax_pv.set_xlim(volumes_mL.min() - pad_v, volumes_mL.max() + pad_v)
    ax_pv.set_ylim(pressures_mmHg.min() - pad_p, pressures_mmHg.max() + pad_p)
    ax_pv.grid(True, alpha=0.3)

    # ---- Sphere (cm) ----
    # Thin-wall convention: `radius` is the mid-wall radius; draw the wall
    # straddling it symmetrically, inner = R - t/2, outer = R + t/2.
    radius_cm = radius * 100.0
    thickness_cm = thickness * 100.0
    inner_cm = radius_cm - 0.5 * thickness_cm
    outer_cm = radius_cm + 0.5 * thickness_cm
    lim = 1.15 * outer_cm.max()
    ax_sphere.set_xlim(-lim, lim)
    ax_sphere.set_ylim(-lim, lim)
    ax_sphere.set_aspect("equal")
    ax_sphere.set_title("Sphere")
    ax_sphere.set_xlabel("x [cm]")
    ax_sphere.set_ylabel("y [cm]")
    ax_sphere.grid(True, alpha=0.3)
    ax_sphere.set_axisbelow(True)

    inner = Circle((0, 0), inner_cm[0], facecolor="white",
                   edgecolor="black", lw=1.0, zorder=2)
    outer = Circle((0, 0), outer_cm[0], facecolor=color,
                   edgecolor="black", lw=1.0, zorder=1)
    ax_sphere.add_patch(outer)
    ax_sphere.add_patch(inner)

    time_text = ax_sphere.text(
        0.02, 0.97, "", transform=ax_sphere.transAxes,
        va="top", ha="left", fontsize=10,
    )

    def update(i):
        marker.set_data([volumes_mL[i]], [pressures_mmHg[i]])
        inner.set_radius(inner_cm[i])
        outer.set_radius(outer_cm[i])
        time_text.set_text(
            f"t = {times[i]:.3f} s\n"
            f"R = {radius_cm[i]:.2f} cm\n"
            f"h = {thickness_cm[i]:.2f} cm"
        )
        return marker, inner, outer, time_text

    fig.tight_layout()
    return fig, update


def save_animation(fig, update, n_frames, out_path, fps):
    anim = animation.FuncAnimation(
        fig, update, frames=n_frames, interval=1000 / fps, blit=False,
    )
    ext = os.path.splitext(out_path)[1].lower()
    if ext in (".mp4", ".mov", ".m4v"):
        writer = animation.FFMpegWriter(
            fps=fps,
            codec="libx264",
            extra_args=[
                "-preset", "slow",
                "-crf", "16",
                "-pix_fmt", "yuv420p",
            ],
        )
    elif ext == ".gif":
        writer = animation.PillowWriter(fps=fps)
    else:
        raise ValueError(f"Unsupported output extension: {ext}")
    anim.save(out_path, writer=writer, dpi=220)
    plt.close(fig)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("input", help="Path to 0D solver input JSON.")
    p.add_argument("--vessel", default=None,
                   help="Name of ChamberSphere vessel (auto if single).")
    p.add_argument("--out", default=None,
                   help="Output video path (default: <input>_sphere.mp4).")
    p.add_argument("--color", default="crimson",
                   help="Wall fill color (matplotlib color spec).")
    p.add_argument("--speed", type=float, default=1.0,
                   help="Playback speed multiplier (1.0 = real time).")
    p.add_argument("--fps", type=float, default=30.0,
                   help="Target playback frame rate (default: 30).")
    p.add_argument("--stride", type=int, default=None,
                   help="Use every Nth time sample (default: auto to hit "
                        "--fps at real-time duration).")
    p.add_argument("--all-cycles", action="store_true",
                   help="Animate every cardiac cycle instead of only the "
                        "last one (forces output_all_cycles=True).")
    args = p.parse_args(argv)

    with open(args.input) as f:
        config = json.load(f)

    vessel = find_chamber_sphere(config, args.vessel)
    vname = vessel["vessel_name"]
    radius0 = vessel["zero_d_element_values"]["radius0"]
    thick0 = vessel["zero_d_element_values"]["thick0"]

    config.setdefault("simulation_parameters", {})
    config["simulation_parameters"]["output_variable_based"] = True
    config["simulation_parameters"]["output_all_cycles"] = bool(args.all_cycles)
    result = pysvzerod.simulate(config)

    t, dvolume = extract_series(result, f"volume:{vname}")
    if t.size == 0:
        raise RuntimeError(
            f"No output for 'volume:{vname}'. Check vessel name."
        )
    # The 'volume' state is the change from the undeformed sphere volume;
    # add it back so we plot the physical total volume.
    volume = dvolume + (4.0 / 3.0) * np.pi * radius0 ** 3

    # Pick a chamber pressure signal; Pin = Pout by construction.
    pressure = None
    for pname in chamber_pressure_name(config, vname):
        tp, pressure = extract_series(result, pname)
        if pressure.size:
            break
    if pressure is None or pressure.size == 0:
        raise RuntimeError(f"Could not locate chamber pressure for {vname}.")

    radius, thickness = thickness_from_volume(volume, radius0, thick0)

    # Keep the full-resolution PV trace; frame stride only affects animation.
    full_volumes_mL = volume / M3_PER_ML
    full_pressures_mmHg = pressure / PA_PER_MMHG

    duration = float(t[-1] - t[0])
    if duration <= 0:
        raise RuntimeError("Simulation time vector has zero duration.")

    if args.stride is None:
        # Downsample to hit ~fps frames per real-time second.
        stride = max(1, int(round(args.speed * len(t) / (duration * args.fps))))
    else:
        stride = max(1, args.stride)

    sl = slice(None, None, stride)
    t = t[sl]
    radius = radius[sl]
    thickness = thickness[sl]
    volumes_mL = full_volumes_mL[sl]
    pressures_mmHg = full_pressures_mmHg[sl]

    duration = float(t[-1] - t[0])
    fps = args.speed * len(t) / duration

    out_path = args.out or os.path.splitext(args.input)[0] + "_sphere.mp4"
    fig, update = build_animation(
        t, pressures_mmHg, volumes_mL, radius, thickness, args.color,
        close_loop=not args.all_cycles,
        pv_trace_volumes=full_volumes_mL,
        pv_trace_pressures=full_pressures_mmHg,
    )
    save_animation(fig, update, len(t), out_path, fps)
    print(
        f"Wrote {out_path} ({len(t)} frames, {fps:.1f} fps, "
        f"{len(t) / fps:.2f} s real time)"
    )


if __name__ == "__main__":
    sys.exit(main())
