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
| Valves (Eq. 35) | piecewise-linear diode `K_at,K_iso,K_ar` | **`PiecewiseValve`** (exact diode, `R=1/K`): mitral `Rmin=1/K_at, Rmax=1/K_iso`; aortic `Rmin=1/K_ar, Rmax=1/K_iso`) | **faithful** (exact piecewise law; hard switch converges as well as the smooth one) |
| Arterial circulation (Eq. 36) | 2-stage Windkessel `C_valve, C_ar, R_p, C_d, R_d` | **exact 2-stage** = two `BloodVesselRC` in series (proximal `Cp=C_ar, Rpd=R_p`; distal `Cp=C_d, Rpd=R_d`) → venous `P_vs`; block `c_valve=C_valve` | **faithful** (`C_ar` now included; `C_valve` in the cavity mass balance) |
| Atrial pressure `P_at` | prescribed low pressure with pre-systolic atrial kick (Ref. 32) | prescribed **constant**, value **tuned** to reach EDV | **simplified**: no atrial kick; `P_at` not tabulated in the paper (Table 1 cites Ref. 6) so tuned |
| Venous `P_vs` | prescribed | `Pd≈0` | **simplified** (not given in paper) |

## Numerics

| Ingredient | Genet et al. 2023 | ChamberCylinder | Status |
|---|---|---|---|
| Spatial discretization | single high-order element (`Pₖ`, k≈3), high quadrature (5–11), `Pₖ₋₁` pressure | multiple **linear** `P₁` elements (default 10) + 3-pt Gauss, penalty (no pressure DOF) | **different**: both converge; linear needs ~10–15 elements |
| Radial integration measure | `2πL ∫ … R dR` (physical) | same (`R dR`) | **faithful** (paper's App. A prints `dR`; `R` restored for volume/incompressibility consistency) |
| Time integration | **energy-preserving** midpoint (equilibrium) + backward (incompressibility) + Chapelle internal-variable scheme (Eq. 59–60) | generalized-α (default) **or** `ConsistentStiffIntegrator` (L-stable gen-α + backtracking-line-search Newton + adaptive sub-stepping) | **simplified**: not energy-preserving; the stiff integrator adds *robustness*, not the paper's exact scheme |
| Nonlinear tangent | consistent (analytic) | **complex-step** local kernels (kernel templated on scalar type; imaginary perturbation `h≈1e-30`, `∂C=Im(f)/h`) | **faithful**: machine-accuracy `dC_dy`/`dC_dydot` with no finite-difference cancellation |
| Force-velocity `α|ė_c|` term | handled by the bespoke energy-consistent scheme | `integrator:"stiff"`; a small `|ė_c|` smoothing (ε≈1e-2) regularizes the kink | **works at fine mesh**: with complex-step tangents and an appropriately scaled residual tolerance (`abs_tol≈1e-6` — `1e-9` is `~1e-13` *relative* for `~1e4 Pa` and unreachable), α=12 converges through ne≥12 |

## Results vs. Fig. 5–9 (exact `PiecewiseValve` + exact 2-stage Windkessel, force-velocity `α=12` on)

| Quantity | Genet Fig. 5 | ChamberCylinder | Match |
|---|---|---|---|
| EDV | 137 mL | ~147 mL (`P_at` tuned) | **good** |
| Peak systolic pressure | 12.8 kPa | ~13–15 kPa | **good** (exact valve holds pressure; the smooth `ValveTanh` leaked to ~4–8 kPa) |
| Diastolic pressure | ~0.9 kPa | ~2 kPa (`P_at` tuned high — stiff filling) | fair |
| ESV / EF | 74 mL / 46% | ~31 mL / 79% with faithful `C_valve=9e-9`; **~69 mL / 53% with `C_valve=0`** | **`C_valve`-dependent** (see note) |
| P–V loop shape | rectangular (flat ejection plateau) | rounded (ejection pressure declines) | **partial** |
| Peak twist | ~20° | ~67° | ~3× (free twist) |
| Fig. 6 twist ∝ aspect ratio | yes; P,V unaffected | yes; P,V unaffected | **trend match** |
| Fig. 8 peak P & EF ↑ with wall vol | yes; twist ~flat | yes | **trend match** |
| Fig. 9 twist peaks ±45°, →0 at ±90° | yes | yes | **trend match** |
| Fig. 9 pressure optimum ~±60° | yes | broad/washed-out at this operating point | partial |

**On the exact valve:** replacing the smooth `ValveTanh` with the exact piecewise `PiecewiseValve` (Eq. 35) markedly *improves* the peak-pressure and filling match — the smooth valve leaked during the isovolumic phases (peak pressure only ~4–8 kPa and poor filling at the same parameters), while the hard diode holds pressure to ~13 kPa. The exact valve converges as robustly as the smooth one.

**On ESV / `C_valve`:** with the paper's faithful cavity compliance `C_valve = 9e-9`, the model over-ejects (ESV ~31 vs 74). Reason: the paper's ejection is a flat pressure *plateau* (`Ṗv ≈ 0`, so the `C_valve·Ṗv` term in Eq. 36a is inert), whereas our P–V loop is *rounded* (the compliant afterload sags during ejection), so `C_valve·Ṗv` actively over-drains the cavity. Setting `C_valve = 0` removes this coupling and recovers ESV ~69 / EF ~53%. The rounding itself is the afterload/ESPVR operating-point gap, not a valve issue.

### Root causes of the residual quantitative gaps
1. **Rounded (non-rectangular) P–V loop → ESV too low.** The ejection pressure declines instead of holding a plateau, because the compliant Windkessel afterload sags as blood is ejected and the ventricular ESPVR lets it keep ejecting. The rounding then couples with the faithful cavity compliance `C_valve` (the `C_valve·Ṗv` term) to over-drain the cavity (ESV 31 vs 74). The paper's flat plateau makes that term inert. This is an afterload/ESPVR operating-point difference, **not** a valve or tangent issue.
2. **Twist magnitude (~3×).** The global twist `β` is unconstrained (no basal/apical rotation constraint — a limitation the paper explicitly notes for its own model) and is operating-point sensitive.
3. **Passive filling stiffer.** Reaching EDV≈137 mL needs a higher `P_at` (~2 kPa) than genet (~0.9 kPa); a passive operating-point difference amplified by the exponential fiber term.

## What a full 1:1 reproduction would additionally require
1. Mixed displacement–pressure (Lagrange-multiplier) incompressibility with an inf-sup-stable element pair, replacing the penalty.
2. Inertia (consistent mass matrix; velocity-companion DOFs for the second-order system).
3. The paper's **energy-preserving** temporal discretization (midpoint + backward) with the Chapelle `√k_c`, `τ_c/√k_c` internal-variable update. With complex-step tangents and a properly scaled residual tolerance the force-velocity term already converges at fine meshes (below), so this scheme's remaining role is the *rectangular* pressure plateau (energy consistency), which would in turn fix the `C_valve`/ESV interaction.
4. ~~Analytic consistent tangents~~ — **done** via complex-step differentiation (machine accuracy; the kernel is templated on the scalar type and evaluated with an imaginary perturbation). Fine-mesh + force-velocity (`α=12`, ne≥12) now converges, once the residual tolerance is scaled to the SI magnitudes (`abs_tol≈1e-6`, not `1e-9`).
5. ~~Exact piecewise valve law (Eq. 35) and 2-stage Windkessel (Eq. 36)~~ — **done** (`PiecewiseValve` diodes + two `BloodVesselRC` stages). Still missing: genet's pre-systolic **atrial kick** in `P_at` and an independently-set `P_vs`.
6. High-order (`Pₖ`) spatial elements to match the paper's single-element discretization.
7. Optionally, a basal/apical twist constraint (or accept the free-`β` over-prediction, as the paper does).
