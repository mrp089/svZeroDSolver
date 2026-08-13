# ChamberCylinder vs Genet et al. (2023): parameter provenance & discrepancy hunt

Goal: reproduce Genet et al. 2023 (Fig. 5) with **no fitted parameters** — every
value is either taken exactly from the paper's Table 1, from a cited reference,
or explicitly flagged as *not specified by the sources*. Any remaining mismatch
is then a traceable modeling difference, not a tuning artifact.

References: **[G]** Genet, Diaz, Chapelle, Moireau, *Reduced LV dynamics modeling
based on a cylindrical assumption*, Int. J. Numer. Meth. Biomed. Engng. 2023
e3711. **[6]** Chapelle, Le Tallec, Moireau, Sorine, *Energy-preserving muscle
tissue model*, IJMCE 2012. **[32]** Sainte-Marie, Chapelle, Cimrman, Sorine,
*Modeling and estimation of the cardiac electromechanical activity*, Comput.
Struct. 2006. **[33]** Kimmig, Chapelle, Moireau, *Thermodynamic properties of
muscle contraction models*, AMSES 2019.

## 1. Parameter provenance

### A. Exact from Genet Table 1 (used verbatim)

| Group | Parameters (SI) |
|---|---|
| Geometry | `Ri=19 mm, Re=33 mm, L=57.1 mm`, ref. chamber vol 65 mL, wall vol 130 mL |
| Microstructure | helix angle `α_endo=+60°, α_epi=−60°` (linear through wall) |
| Passive `W^e` (Eq. 25) | `C1=7, C2=0, C3=700, C5=50 Pa`; `C4=2, C6=4` (dimensionless); `ρ0=1 kg/L`; viscosity `γ=70` |
| Active (Eq. 33) | `k_s=1e8 Pa`, **`σ0=65e3 Pa`**, `μ=70 Pa·s`, `α=12.0`, `k0=260e3 Pa` |
| Valves (Eq. 35) | `K_at⁻¹=1.1111e5`, `K_ar⁻¹=7.6923e4`, `K_iso⁻¹=2.0e9` Pa·s·m⁻³ |
| Circulation (Eq. 36) | `C_valve=9.0e-9`, `C_ar=2.115e-10`, `C_d=2.0158e-8` m³/Pa; `R_p=1.35e7`, `R_d=1.0949e8` Pa·s·m⁻³ |

> **Correction applied this pass:** earlier runs used `σ0 = 45–55 kPa`; the exact
> Table 1 value is **65 kPa**. All results below use 65 kPa.

### B. Cited to a reference but given **no value/formula** in Genet

| Parameter | What Genet says | Reference content | Choice made (flagged, not fitted) |
|---|---|---|---|
| `n0(e_c)` Frank–Starling | "a function accounting for the Frank–Starling mechanism", cites [6] | [6] Remark 4: `n0` is a reduction factor `0≤n0≤1`, no closed form. **Ref [13] Caruel et al. 2014 Fig 2(a)** plots the actual `n0(e_c)`: piecewise-linear, 0 at `e_c≈−0.4`, rising to a plateau of 1 over `e_c∈[1,1.3]`, back to 0 at `e_c≈2.5` | **piecewise-linear digitized from Caruel Fig 2(a)** (breakpoints `(−0.4,0),(0.3,0.4),(1,1),(1.3,1),(2.5,0)`), in `frank_starling()`. See §3 OPEN #1: this resolves the over-ejection but exposes an `e_c`-scale question. |
| `ν(t)` activation | "prescribed as detailed in [6,33]", triggered when `[Ca²⁺]>c_th` | [6],[32]: `ν=+k_ATP` while depolarised, `−k_RS` while repolarised (`|ν|₊=k_ATP·1_{Ca>C}`). Rate values not tabulated for this LV. | tanh systole/diastole window; onset `t_sys≈0.13 s`, relaxation `t_dias≈0.45 s`, period 0.8 s **read from [G] Fig 5's own timeline**; rate `|ν|≈30/s`. Rate magnitude is immaterial to the peak (cancels at steady state). |
| `P_at` atrial pressure | Table 1 lists it but the cell is only a citation to [6] | [32]: prescribed low pressure + pre-systolic **atrial kick**; no number for this LV | `P_at ≈ 0.9 kPa`, read from [G] Fig 5 end-diastolic pressure (the paper's stated static end-diastolic loading). Atrial kick omitted. |
| `P_vs` venous pressure | not in Table 1; appears in Eq. 36c | [32] uses a venous `P_ve` in its Windkessel; no value for [G] | `P_vs = 0` (documented). |

### C. Numerical choices (not physical parameters of [G])

| Choice | [G] | Here | Note |
|---|---|---|---|
| Incompressibility | Lagrange-multiplier (exact `J=1`), mixed `P_k/P_{k-1}` | penalty `κ(J−1)J C⁻¹`, `κ=1e7` | penalty; `κ` not a [G] parameter |
| Spatial discretization | single high-order `P_k` element | linear `P1`, `ne≈12` + 3-pt Gauss | converged (see §3) |
| Time integration | energy-preserving midpoint + Chapelle internal-var scheme | generalized-α or `ConsistentStiffIntegrator` | |
| Nonlinear tangent | analytic | **complex-step** (machine accuracy) | |
| Residual tolerance | — | `abs_tol=1e-6` | `1e-9` is `~1e-13` *relative* for `~1e4 Pa`; unreachable |

## 2. Result with all-exact parameters (zero fitting)

Exact valves + 2-stage Windkessel, `σ0=65 kPa`, force-velocity `α=12`, `ne=12`,
`P_at=0.9 kPa`, `P_vs=0`. Two `n0` variants shown: `n0=1` (initial placeholder)
and the **Caruel Fig 2(a) `n0(e_c)`** (current).

| Quantity | Genet Fig 5 | `n0=1` | **Caruel `n0`** | Status |
|---|---|---|---|---|
| Diastolic pressure | ~0.9 kPa | 0.81 | 0.85 kPa | **match** |
| EDV | 137 mL | 124 | 123 mL | good (~10% low) |
| ESV | 74 mL | 24 | **62 mL** | Caruel `n0` fixes it |
| EF | 46% | 81% | **50%** | Caruel `n0` fixes it |
| Peak twist | ~20° | 64° | **31°** | Caruel `n0` fixes it |
| Peak systolic pressure | 12.8 kPa | 9.5 | **5.0 kPa** | now too *low* — see OPEN #1 |

## 3. Discrepancy log (systematic)

### RESOLVED — passive material is correct (not a discrepancy)
An earlier measurement suggested the passive wall was ~4× too stiff (67 mL at
0.9 kPa vs 137). **This was a test artifact, not the model.** A `ValveTanh`
inflation test sits half-open at zero pressure drop (`R≈(R_min+R_max)/2≈1e9`),
throttling a constant-pressure fill so the wall never equilibrates. With the
cavity pressure applied **directly** (no valve), the block reaches ~128–137 mL
at 0.9 kPa, matching both Genet Fig 5 and an independent analytical
incompressible thick-wall cylinder inflation with the same `W^e`
(`scripts/ChamberCylinder_validate_passive.py`). Ruled out along the way:
parameters (exact), the `W^e` formula (matches Eq. 25 term-by-term), the fiber
term (`C5=0` unchanged), active coupling (`active_model` 0 vs 1 identical),
mesh (`ne` 6→24), and penalty `κ` (1e5→1e8). In the full simulation, diastole
has a real pressure gradient so EDV=124 mL ≈ paper's 137.

### OPEN #1 — largely RESOLVED by Frank–Starling; residual is an `e_c`-scale question
The over-ejection was **Frank–Starling**. Ruled out first: valve law (hard vs
smooth identical), `C_valve` (ESV≈24–37 for 0→9e-9), force-velocity
(ESV≈24–37 for `α` 0→12). Implementing the actual `n0(e_c)` from Caruel 2014
Fig 2(a) (§1.B) moves ESV **24→62 mL** (paper 74), EF **81→50%** (paper 46),
twist **64→31°** (paper 20) — the ventricle now stops near reference volume as
the paper does.

**New, sharper open item:** peak pressure drops to **5.0 kPa** (paper 12.8),
now `C_valve`-independent. Cause: the model's contractile strain `e_c` operates
at **−0.11…0.22 (mean 0.095)** — the *rising limb* of Caruel's curve, where
`n0≈0.2–0.35`, so effective tension is only ~30% of `σ0`. The paper's 12.8 kPa
peak implies operation near the **plateau** (`n0≈1`, `e_c≈1`). So either the
model's `e_c` is defined on a different scale than Caruel's abscissa, or the
curve's abscissa is a normalized sarcomere length rather than `e_c=√I4−1`.
**This is now the single most valuable question for M. Genet** (see requests
list): what is the abscissa of the `n0` curve, and at what `e_c` should a
well-filled beat sit?

### OPEN #2 — peak twist ~3× (64° vs 20°)
Global twist `β` is an unconstrained DOF (no basal/apical rotation constraint —
a limitation [G] notes for its own model). Operating-point sensitive.

### OPEN #3 — atrial kick / P_at, P_vs
`P_at` is a prescribed waveform with a pre-systolic atrial kick ([32]); we use a
constant `P_at=0.9 kPa`. Not expected to affect ESV (see OPEN #1) but affects the
diastolic filling detail and EDV.

## 4. Reproduce

- Passive validation: `PYTHONPATH=build/python python scripts/ChamberCylinder_validate_passive.py`
- Exact-parameter circulation builder: `scripts/ChamberCylinder_genet_exact.py` —
  all `build()` defaults are the Table-1 / flagged values above.
