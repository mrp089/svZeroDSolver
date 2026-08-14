# ChamberCylinder — model ingredient audit (provenance & discrepancies)

Every ingredient of the ChamberCylinder + circulation model, its **source**, and
**where it could differ from Genet et al. (2023)**. Companion to the discrepancy
*hunt* log in `ChamberCylinder_discrepancies.md`; this file is the *provenance*
map. Cross-checked against the actual kernel (`src/model/ChamberCylinder.cpp/.h`),
the builder (`ChamberCylinder_genet_exact.py`), and the MEDISIM PhysioBlocks
config (`references/full_configurations/spherical_heart_sim.jsonc`).

Provenance tags: **[T1]** Genet Table 1 (exact); **[Eq]** Genet formulation
equation; **[PB]** MEDISIM PhysioBlocks (Genet's group's code — secondary source
for inputs Genet does not tabulate); **[F5]** read from the digitized Fig 5;
**[?]** my choice / not sourced (discrepancy risk).

## Ranked discrepancies (what is left)

1. **Time integration [Eq, not matched].** Genet: energy-preserving midpoint +
   Chapelle √k_c internal-variable scheme. Here: L-stable generalized-α
   (`rho_infty=0`). This adds numerical damping → the model's diastole (refill +
   untwist) is faster/more abrupt than Genet's. Dominant remaining difference.
2. **BCS activation ν(t) [PB, now implemented].** The exact PhysioBlocks
   trapezoid (ν ∈ [−20, +35], `activation_mode=1`) is used, its onset retimed to
   the Fig 5 **volume curve** (`tsys=0.12`; the ejection downstroke overlaps Fig 5
   to ~5 ms, and torsion aligns — the pressure upstroke then sits ~40 ms late, the
   model having essentially no isovolumic phase). Documented discrepancy: the
   **−20 s⁻¹ diastolic rate relaxes slower** than whatever Genet used for the
   cylindrical Fig 5 → filling is delayed → the diastolic **volume** fit is worse
   (RMS_V ≈ 20 vs ≈ 7 mL for an ad-hoc ±30 tanh), and ESV over-ejects (68 vs 74).
   The PhysioBlocks *spherical-sim* activation ≠ Genet's Fig-5 activation (which
   the paper does not tabulate). Kept for fidelity to the cited reference.
3. **Prescribed atrial P_at [PB].** Exact PhysioBlocks waveform is implemented,
   but its 450 Pa mid-diastole floor equilibrates the wall ~6 mL below Genet's
   diastasis (static(450 Pa)≈116 vs ~122 mL).
4. **k_0 (crossbridge rate) [T1 vs PB].** I use 260e3 (Table 1); PhysioBlocks
   uses 273e3. ~5%, minor.
5. **Frank–Starling n0(e_c) [PB].** Genet leaves n0 a general factor; the 9-point
   curve is taken from PhysioBlocks (matches exactly). An *inferred* input.

**RESOLVED — C_valve.** Exact 9e-9 [T1] now used as the **aortic-root compliance**
(Genet Eq 36 / PhysioBlocks `capacitance_valve` at `aorta_proximal`, in parallel
with C_ar). The earlier "over-buffering" (~7.5 vs 12.8 kPa) was a **placement bug**
— C_valve had been applied to the ventricular mass balance (−C_valve·dP_v/dt);
moved into the circulation it gives peak 12.9 kPa with a flat ejection plateau.

Everything else below is either exact [T1] / [Eq] or a converged numerical choice.

---

## A. Geometry & fiber architecture
| item | value | source | possible discrepancy |
|---|---|---|---|
| inner radius Ri | 0.019 m | **[T1]** | none |
| outer radius Re | 0.033 m | **[T1]** | none |
| length L | 0.0571 m | **[T1]** | none |
| helix angle | +60°(endo) → −60°(epi), **linear** in R | **[T1]** angle; **[Eq]** profile | Genet's exact transmural profile; linear is the standard assumption and reproduces the Fig 9 fiber-angle optimum |

## B. Kinematics (reduced cylindrical model)
- DOFs ζ = [ρ, β, φ, ε, η]; deformation gradient F (Eq 13), C = FᵀF, Green–Lagrange
  E; fiber invariant Ī4 = C̄:(e_f⊗e_f). **[Eq]** (Genet Eqs 13–21, 31).
- **Verified**: the implemented I4 equals Genet Eq 31's (1+e_fib)² term-by-term.
- Discrepancy: the through-wall basis (linear ρ, φ, η per element) is my FE choice,
  not Genet's spectral/basis choice — but the targets are mesh-converged (2nd order),
  so this is not a source of the residual.

## C. Passive material Wᵉ
| coef | value | source |
|---|---|---|
| C1, C2 | 7, 0 Pa | **[T1]** |
| C3, C4 | 700 Pa, 2 | **[T1]** |
| C5, C6 | 50 Pa, 4 | **[T1]** |
- Form: Wᵉ = C1(Ī1−3) + C2(Ī2−3) + C3(exp(C4(Ī1−3)²)−1) + C5(exp(C6(Ī4−1)²)−1),
  reduced invariants Īk = J^(−2/3)·Ik. **[Eq]** (Genet Eq 25).
- **Verified** against an independent incompressible thick-wall cylinder inflation
  with the same Wᵉ (`ChamberCylinder_validate_passive.py`): agree to <0.5 mL at
  converged mesh. No discrepancy in form or parameters.

## D. Incompressibility
- **Mixed u/p** (element-constant hydrostatic pressure p_e, weak ∫(J−1)=0), Σᵇ =
  p_e J C⁻¹. **[Eq]** (Genet enforces incompressibility with a Lagrange-multiplier
  pressure field). Penalty fallback (κ=1e7) locks at low order — **not** used.
- Discrepancy: my P0–P1 pairing vs Genet's exact mixed element; spatially converged
  (order ≈2), so negligible.

## E. Active model (Bestel–Clément–Sorine) — main parameter-provenance area
| item | value | source | note |
|---|---|---|---|
| structure τc/kc ODEs, σ_1D | Eqs 30–33 | **[Eq]** | σ_1D=(τc+μėc)/√I4 **verified** vs Kimmig 2019 Eq 31 |
| max active stress σ0 | 65e3 Pa | **[T1]** | PB(spherical) uses 90e3 — different heart; T1 is correct here |
| series stiffness k_s | 1e8 Pa | **[T1] = [PB]** | match |
| series damping μ | 70 Pa·s | **[T1] = [PB]** (`damping_coef`) | match |
| crossbridge rate k_0 | 260e3 | **[T1]**; PB 273e3 | ~5% discrepancy (source difference) |
| force–velocity α | 12 | **[T1] = [PB]** (`destruction_rate`) | match |
| Frank–Starling n0(e_c) | 9-pt piecewise-linear | **[PB]** | matches PB exactly; Genet leaves n0 general (inferred input) |
| **activation ν(t)** | **+35 / −20 trapezoid** (`activation_mode=1`), onset retimed to the volume curve (`tsys=0.12`) | **[PB]** | exact PhysioBlocks waveform; its −20 relaxation is slow → worse diastolic-volume fit (ranked #2) |

## F. Viscosity (parallel)
- Σᵛ = γ Ė, γ = 70 Pa·s. **[T1] = [PB]** (`damping_parallel`). **[Eq]** Genet Eq 28.

## G. Cavity pressure loading & volume
- Follower pressure Pv on the inner surface (ρ-equation) and lids (ε-equation),
  V = π(Ri+ρi)²(1+ε)L. **[Eq]** (Genet Eq 63). No parameter; the follower
  linearization is analytic. Discrepancy: none identified.

## H. Boundary conditions
- Rigid pins φ(Ri)=η(Ri)=0 (base/apex tethering). **[Eq]**. Discrepancy: the exact
  tethering set — at this operating point β is free and twist matches (20°), so no
  extra basal constraint is needed.

## I. Inertia
- density ρ0 = 1000 kg/m³ **[T1]**, but the beat runs quasi-dynamically
  (`use_inertia=0`). Consistent-mass inertia was tested: ~0.3% effect. Discrepancy:
  inertia omitted — negligible.

## J. Valves (Eq 35)
| item | value | source |
|---|---|---|
| mitral K_at⁻¹ | 1.1111e5 | **[T1]** |
| aortic K_ar⁻¹ | 7.6923e4 | **[T1]** |
| closed K_iso⁻¹ | 2.0e9 | **[T1]** |
- Implemented as `PiecewiseValve` (hard diode). Discrepancy: Genet's Eq 35 is a
  piecewise-*linear* resistance; the hard switch is the limiting case and matched
  the beat.

## K. Windkessel (Eq 36)
| item | value | source | note |
|---|---|---|---|
| C_ar | 2.115e-10 m³/Pa | **[T1]** | |
| C_d | 2.0158e-8 | **[T1]** | at `aorta_distal` |
| R_p | 1.35e7 Pa·s/m³ | **[T1]** | |
| R_d | 1.0949e8 | **[T1]** | |
| C_valve | 9.0e-9 | **[T1]** | **aortic-root compliance** at `aorta_proximal`, in parallel with C_ar (RESOLVED — was mis-placed on the ventricle) |

## L. Prescribed pressures
| item | value | source | discrepancy |
|---|---|---|---|
| atrial P_at | rescale_two_phases, 450↔900 Pa, late kick | **[PB]** | mid-diastole floor → diastasis ~6 mL low |
| venous P_vs | 1600 Pa | **[PB]** | Genet not tabulated (Caruel 1000); sets afterload floor → exact peak |

## M. Activation timing
- t_sys = 0.11 s, t_dias = 0.40 s (~290 ms systole). **[F5]** (read from Fig 5:
  pressure-rise and plateau-end). PB reference waveform uses different fractions.
  Discrepancy: read from the figure, not tabulated by Genet.

## N. Time integration
- L-stable generalized-α (`ConsistentStiffIntegrator`, ρ∞=0). **[?]** vs Genet's
  energy-preserving midpoint + Chapelle √k_c. **The dominant remaining difference**
  (diastolic dynamics). See discrepancies OPEN #6.

## O. Numerics (converged; see ChamberCylinder_convergence.py)
- Spatial: 2nd-order in ne; converged. Temporal: converged at dt=2 ms (halving
  moves peak P 0.05 Pa, ESV 0.011 mL). Coarsest **stable** = ne=3, dt=2 ms →
  within ~0.5% of the extrapolated limit; used for further studies.
- Robustness edges: mixed u/p + stiff integrator fail to solve at ne≥16 and dt≥4 ms.
- Tolerance 1e-9; ncycle≥12 to reach the limit cycle (Windkessel transient).
