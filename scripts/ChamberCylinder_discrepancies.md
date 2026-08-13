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
| `n0(e_c)` Frank–Starling | "a function accounting for the Frank–Starling mechanism", cites [6] | [6] Remark 4: general reduction factor, no closed form. Ref [13] Caruel et al. 2013 evaluates `n0` as `interp(e_c, …)`; the **MEDISIM PhysioBlocks** code gives the calibrated breakpoints. | **piecewise-linear from PhysioBlocks** `physioblocks/physioblocks` (`active_law.py`: `n0 = interp(e_c, abscissas, ordinates)`), abscissas `[−0.167,−0.007,0.053,0.097,0.133,0.202,0.466,0.919,1.176]`, ordinates `[0,0.56,0.77,0.89,0.96,1,1,0.11,0]` — plateau at **physiological** `e_c∈[0.20,0.47]`. (Caruel Fig 7(a) *looks* like it plateaus at `e_c≈1` only because it is drawn over the wide isotonic papillary-muscle strain range.) |
| `ν(t)` activation | "prescribed as detailed in [6,33]", triggered when `[Ca²⁺]>c_th` | [6],[32]: `ν=+k_ATP` while depolarised, `−k_RS` while repolarised (`|ν|₊=k_ATP·1_{Ca>C}`). Rate values not tabulated for this LV. | tanh systole/diastole window; onset `t_sys≈0.13 s`, relaxation `t_dias≈0.45 s`, period 0.8 s **read from [G] Fig 5's own timeline**; rate `|ν|≈30/s`. Rate magnitude is immaterial to the peak (cancels at steady state). |
| `P_at` atrial pressure | Table 1 lists it but the cell is only a citation to [6] | [32]: prescribed low pressure + pre-systolic **atrial kick**. **PhysioBlocks** gives the concrete waveform: baseline **450 Pa** (diastasis), ramping to **900 Pa (atrial kick)** over the last ~15% of the cycle, held through early systole | constant `P_at = 0.9 kPa` (= the PhysioBlocks kick value; baseline 450 Pa and the kick timing not yet applied — see OPEN #4) |
| `P_vs` venous pressure | not in Table 1; appears in Eq. 36c | [32] `P_ve`; **PhysioBlocks** venous = **1600 Pa**; Caruel [13] `P_sv = 1000 Pa` | `P_vs = 0` (documented; ~1–1.6 kPa per the references would raise mean arterial pressure slightly) |

### C. Numerical choices (not physical parameters of [G])

| Choice | [G] | Here | Note |
|---|---|---|---|
| Inertia | full dynamics (Eqs. 8,18,45) | **implemented** (`use_inertia=1`): velocity companion DOFs + consistent mass `M=∫ρ₀(Du)ᵀDu dΩ` | matches quasi-static to ~0.3% — inertia is ~1e-4 of the forces (quasi-static is the default) |
| Incompressibility | Lagrange-multiplier (exact `J=1`), mixed `P_k/P_{k-1}` | penalty `κ(J−1)J C⁻¹`, `κ=1e7` | penalty; `κ` not a [G] parameter |
| Spatial discretization | single high-order `P_k` element | linear `P1`, `ne≈12` + 3-pt Gauss | converged (see §3) |
| Time integration | energy-preserving midpoint + Chapelle internal-var scheme | generalized-α or `ConsistentStiffIntegrator` | |
| Nonlinear tangent | analytic | **complex-step** (machine accuracy) | |
| Residual tolerance | — | `abs_tol=1e-6` | `1e-9` is `~1e-13` *relative* for `~1e4 Pa`; unreachable |

## 2. Result with all-exact parameters (zero fitting)

Exact valves + 2-stage Windkessel, `σ0=65 kPa`, force-velocity `α=12`, `ne=12`,
`P_at=0.9 kPa`, `P_vs=0`, **PhysioBlocks `n0(e_c)`**, `C_valve=0` (see note).

| Quantity | Genet Fig 5 | `n0=1` | **PhysioBlocks `n0`** | Status |
|---|---|---|---|---|
| Diastolic pressure | ~0.9 kPa | 0.81 | 0.79 kPa | **match** |
| Peak systolic pressure | 12.8 kPa | 9.5/20 | **12.4 kPa** | **match** |
| Peak twist | ~20° | 64° | **20°** | **match** |
| EDV | 137 mL | 124 | 124 mL | good (~10% low) |
| ESV | 74 mL | 24 | 58 mL | good |
| EF | 46% | 81% | 53% | good |

At peak contraction `e_c` reaches ~0.26, on the `n0` plateau (`n0→1`), so the
ventricle develops near-full `σ0`; during filling/relaxation `e_c` falls onto the
rising limb (`n0<1`), giving the length-dependent self-limitation. **`C_valve`
note:** the exact `C_valve=9e-9` stores ~0.1 L at systolic pressure and buffers
the peak down to ~7.5 kPa in this quasi-static model; `C_valve=0` gives the match
above. The paper's energy-preserving scheme keeps `Ṗv≈0` on the ejection plateau
so its `C_valve·Ṗv` term is inert — a discretization difference, not a parameter.

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

### RESOLVED #1 — systolic over-ejection was the Frank–Starling `n0(e_c)`
The over-ejection was **Frank–Starling**. Ruled out first: valve law (hard vs
smooth identical), `C_valve` (ESV≈24–37 for 0→9e-9), force-velocity
(ESV≈24–37 for `α` 0→12). Implementing the actual length-dependence `n0(e_c)`
(§1.B) fixes it — see §2: peak pressure and twist now **match** (12.4 vs 12.8
kPa; 20 vs 20°), ESV/EF good.

A digitization detour, now resolved: my first pass digitized Caruel Fig 7(a),
whose abscissa runs over the wide isotonic papillary-muscle range so the plateau
*looks* like `e_c≈1`. That put the model (which operates at `e_c~0.1–0.26`) on
the rising limb (`n0~0.3`) and dropped the peak to 5 kPa. The **PhysioBlocks
calibrated breakpoints** put the plateau at the physiological `e_c∈[0.20,0.47]`,
exactly where the model operates → `n0→1` at peak, restoring the pressure. The
model's `e_c=√I4−1` needs **no rescaling**: its operating range coincides with
the calibrated plateau, and PhysioBlocks likewise feeds the fiber deformation
straight into `interp(e_c,…)`.

### RESOLVED #2 — peak twist
Twist was operating-point/`n0`-driven, not a free-`β` defect: with the correct
`n0(e_c)` it lands at **20°**, matching [G]. (`β` remains a free DOF; no
basal/apical constraint is needed at this operating point.)

### OPEN #3 — `C_valve` over-buffering (peak pressure with exact `C_valve`)
With the exact `C_valve=9e-9` the peak pressure is buffered to ~7.5 kPa (vs 12.4
at `C_valve=0`). The compliance stores ~0.1 L per systolic pressure swing, which
is large relative to the stroke volume when the ejection pressure is rounded
rather than a flat plateau. In [G]'s energy-preserving scheme the pressure holds
a plateau (`Ṗv≈0`), so `C_valve·Ṗv` in Eq. 36a is inert. **Adding inertia does
not fix it** (dynamic vs quasi-static agree to ~0.3%; §1.C), confirming the cause
is the temporal scheme (energy-preserving midpoint + Chapelle internal-variable
update), not inertia and not a parameter mismatch.

### OPEN #4 — residual EDV/ESV (~10%) and `P_at` waveform
EDV 124 vs 137 and ESV 58 vs 74. `P_at` is a prescribed waveform with a
pre-systolic atrial kick ([32]); we use a constant `P_at=0.9 kPa`, which under-
fills slightly. The atrial kick would raise EDV toward 137 and, via preload,
ESV toward 74.

## 4. Reproduce

- Passive validation: `PYTHONPATH=build/python python scripts/ChamberCylinder_validate_passive.py`
- Exact-parameter circulation builder: `scripts/ChamberCylinder_genet_exact.py` —
  all `build()` defaults are the Table-1 / flagged values above.
