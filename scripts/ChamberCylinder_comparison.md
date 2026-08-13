# ChamberCylinder vs. Genet et al. (2023): model comparison

Genet, Diaz, Chapelle, Moireau, *"Reduced left ventricular dynamics modeling
based on a cylindrical assumption"*, Int. J. Numer. Meth. Biomed. Engng. (2023)
e3711. Circulation/active-law details from Ref. 32 (Sainte-Marie et al., Comput.
Struct. 2006) and Ref. 6 (Chapelle et al., 2012).

This documents what the `ChamberCylinder` block reproduces faithfully, what is
simplified, and what a full 1:1 reproduction would additionally require. The
guiding constraint was to reuse svZeroDSolver's DAE machinery ("the governing
equations from the paper, but not the numerical scheme that is complicated").

## Physics / continuum model

| Ingredient | Genet et al. 2023 | ChamberCylinder | Status |
|---|---|---|---|
| Cylindrical kinematics `ζ=[ρ,β,φ,ε,η]` (Eq. 1, 13–16) | full: 3 radial fields + 2 scalars | full, identical | **faithful** |
| Passive deviatoric `W^e` (Eq. 25) | transversely isotropic, `C1..C6`, isochoric `Ī1,Ī2,Ī4` | identical | **faithful** |
| Myofiber orientation `e_F(R)` (Eq. 22) | helix angle linear endo→epi | identical | **faithful** |
| Viscous stress `Σ^v=γĖ` (Eq. 28) | included | included | **faithful** |
| Active fiber stress `Σ^a=σ_1D e_F⊗e_F` (Eq. 29–33) | Bestel–Clément–Sorine `e_c,τ_c,k_c` ODEs, force-velocity, Frank–Starling | same (`active_model=1`); force-velocity needs the stiff integrator | **faithful** (see notes) |
| Pressure loading, inner surface + lids (Eq. 37–43, 63) | full | identical | **faithful** |
| Rigid-body pins `φ(Rᵢ)=η(Rᵢ)=0` (Sec. 2.4) | included | included | **faithful** |
| Prestress | neglected | neglected | **match** |
| Incompressibility | Lagrange-multiplier pressure field `p`, `∫(J−1)p̂ dΩ=0`, mixed `Pₖ/Pₖ₋₁` | near-incompressible **penalty** `Σ^b=κ(J−1)J C⁻¹` | **simplified**: `J` within ~0.4% of 1; no pressure field, no inf-sup element |
| Inertia `P_a` (Eq. 18, mass matrix A8) | included | **neglected** (quasi-static) | **simplified**: dynamics come from viscosity + active ODEs + circulation |

### Active-law notes
- Equations match Chapelle 2012 Eq. (9)/(10) exactly (`τ_c`,`k_c` ODEs, `σ_1D=(τ_c+μė_c)/(1+e_fib)`).
- **Rigid series-spring form**: `σ_1D` is written via `τ_c+μė_c` (= fiber tension) rather than `k_s(e_fib−e_c)`, keeping the paper's stiff `k_s=1e8` out of the mechanical residual (it stays in the `e_c` ODE). Numerically equivalent for such large `k_s`.
- **Frank–Starling `n_0(e_c)`** is left deliberately general in Chapelle 2012 (Remark 4 — a reduction factor, function of the history of `e_c`, no fixed formula). Implemented here as a tunable Gaussian force-length curve, default ≈1.
- **Activation `ν(t)`** uses the same tanh systole/diastole window as the simple model (the paper prescribes `ν` from a Ca²⁺ threshold — a comparable prescription, not identical).
- A simpler `active_model=0` (single activated fiber stress, ChamberSphere-style) is also provided; it is *not* the paper's law.

## Circulation / boundary conditions

| Ingredient | Genet et al. 2023 | ChamberCylinder test setup | Status |
|---|---|---|---|
| Valves (Eq. 35) | piecewise-linear diode `K_at,K_iso,K_ar` | `ValveTanh` (`Rmin=1/K`, `Rmax=1/K_iso`) | **approximated** (smooth vs piecewise) |
| Arterial circulation (Eq. 36) | 2-stage Windkessel `C_valve, C_ar, R_p, C_d, R_d` | `WindkesselBC` RCR (`Rp=R_p, C=C_d, Rd=R_d`) + block `c_valve` | **approximated**: `C_ar` (2.1e‑10) dropped as negligible vs `C_d` (2e‑8) → 2-stage ≈ RCR |
| Atrial pressure `P_at` | prescribed low pressure with pre-systolic atrial kick (Ref. 32) | prescribed **constant**, value **tuned** to reach EDV | **simplified**: no atrial kick; `P_at` not tabulated in the paper so tuned |
| Venous `P_vs` | prescribed | `Pd≈0` | **simplified** (not given in paper) |

## Numerics

| Ingredient | Genet et al. 2023 | ChamberCylinder | Status |
|---|---|---|---|
| Spatial discretization | single high-order element (`Pₖ`, k≈3), high quadrature (5–11), `Pₖ₋₁` pressure | multiple **linear** `P₁` elements (default 10) + 3-pt Gauss, penalty (no pressure DOF) | **different**: both converge; linear needs ~10–15 elements |
| Radial integration measure | `2πL ∫ … R dR` (physical) | same (`R dR`) | **faithful** (paper's App. A prints `dR`; `R` restored for volume/incompressibility consistency) |
| Time integration | **energy-preserving** midpoint (equilibrium) + backward (incompressibility) + Chapelle internal-variable scheme (Eq. 59–60) | generalized-α (default) **or** `ConsistentStiffIntegrator` (L-stable gen-α + backtracking-line-search Newton + adaptive sub-stepping) | **simplified**: not energy-preserving; the stiff integrator adds *robustness*, not the paper's exact scheme |
| Nonlinear tangent | consistent (analytic) | **finite-difference** local kernels | **simplified**: residual exact so converged solution is exact, but FD tangent limits fine-mesh + force-velocity convergence |
| Force-velocity `α|ė_c|` term | handled by the bespoke energy-consistent scheme | needs `integrator:"stiff"`; a small `|ė_c|` smoothing (ε≈1e-2) regularizes the kink | **simplified**: converges at ne≤6–8; fails at ne=12 (FD-tangent accuracy) |

## Results vs. Fig. 5–9 (with genet's valve + Windkessel circulation, force-velocity on)

| Quantity | Genet Fig. 5 | ChamberCylinder | Match |
|---|---|---|---|
| EDV / ESV / EF | 137 mL / 74 mL / 46% | ~139 / 79 / 43% | **good** |
| P–V loop (shape, volume range) | rectangular, 74–137 mL | matches | **good** |
| Diastolic pressure | ~0.9 kPa | ~1 kPa | **good** |
| Peak systolic pressure | 12.8 kPa | ~15.5 kPa | ~20% high (overshoot) |
| Peak twist | ~20° | ~45° (with FV; 63° without) | ~2× (free twist) |
| Fig. 6 twist ∝ aspect ratio | yes; P,V unaffected | yes; P,V unaffected | **trend match** |
| Fig. 8 peak P & EF ↑ with wall vol | yes; twist ~flat | yes | **trend match** |
| Fig. 9 twist peaks ±45°, →0 at ±90° | yes | yes | **trend match** |
| Fig. 9 pressure optimum ~±60° | yes | broad/washed-out at this operating point | partial |

### Root causes of the residual quantitative gaps
1. **Peak-pressure overshoot (~20%)** — the ventricle ejects too fast, so the characteristic-resistance drop `R_p·Q` inflates `P_v`. The force-velocity term reduces this; the remainder is the FD tangent + absence of the energy-consistent scheme.
2. **Twist magnitude (~2×)** — the global twist `β` is unconstrained (no basal/apical rotation constraint — a limitation the paper explicitly notes for its own model). It is also operating-point sensitive (higher at smaller volume). Force-velocity lowers it (63°→45°).
3. **Passive filling stiffer** — reaching EDV≈137 mL needs a higher `P_at` (~3–4 kPa) than genet (~0.9 kPa); a passive operating-point difference amplified by the exponential fiber term.

## What a full 1:1 reproduction would additionally require
1. Mixed displacement–pressure (Lagrange-multiplier) incompressibility with an inf-sup-stable element pair, replacing the penalty.
2. Inertia (consistent mass matrix; velocity-companion DOFs for the second-order system).
3. The paper's **energy-preserving** temporal discretization (midpoint + backward) with the Chapelle `√k_c`, `τ_c/√k_c` internal-variable update — which is what carries the stiff force-velocity term to fine meshes.
4. **Analytic consistent tangents** for the stress/active kernels (replacing finite differences) — needed for robust fine-mesh + force-velocity convergence.
5. The exact piecewise valve law (Eq. 35) and full 2-stage Windkessel (Eq. 36), with genet's `P_at` (including the pre-systolic atrial kick) and `P_vs`.
6. High-order (`Pₖ`) spatial elements to match the paper's single-element discretization.
7. Optionally, a basal/apical twist constraint (or accept the free-`β` over-prediction, as the paper does).
