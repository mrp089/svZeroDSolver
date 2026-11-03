# ChamberCylinder: genet23 reproduction and verification

Scripts and data behind the `ChamberCylinder` block, the thick-walled
cylindrical left ventricle of Genet et al. (2023), "Reduced left ventricular
dynamics modeling based on a cylindrical assumption", Int. J. Numer. Meth.
Biomed. Engng., e3711 (genet23).

| File | Content |
| --- | --- |
| `ChamberCylinder_genet_exact.py` | svZeroDSolver configuration of the genet23 baseline: Table 1 parameters, valve law (Eq 35) and two-stage Windkessel (Eq 36), activation and atrial-pressure inputs |
| `ChamberCylinder_figures.py` | runs the cases of genet23 Figs 5-9 and overlays them on the digitized traces |
| `genet23_fig{5..9}_digitized.csv` | numeric traces digitized from the published Figs 5-9 |
| `ChamberCylinder_verification.py` | verification against independent solutions of the same continuum model (symmetry, static inflation, free vibration, an independent closed-loop re-implementation); no data from genet23 |
| `ChamberCylinder_provenance.tex` | provenance of every parameter and input, the formulation audit, and the analysis of the remaining torsion discrepancy |
| `ChamberCylinder.yaml` | the continuous governing equations the block discretizes (documentation, not a Jacobian-generator input) |

Only digitized numbers are stored here; the paper's figures, text and the
authors' code are not.

## Running

Build the Python module (see the main README), then from the repository root:

```bash
PYTHONPATH=build/python python scripts/ChamberCylinder/ChamberCylinder_figures.py
PYTHONPATH=build/python python scripts/ChamberCylinder/ChamberCylinder_verification.py T1 T2 T3 T4 plot
```

Both write to `scripts/ChamberCylinder/output/` (not tracked). The full
verification takes a few minutes; its fast subset runs in the test suite
(`tests/test_chamber_cylinder_verification.py`). Build the document from this
folder with `pdflatex ChamberCylinder_provenance.tex` (twice).

The torsion diagnostics in the document (formulation switches B1-B4, the
literature screen and the damping fingerprint) used temporary kernel switches
that are not part of the block. They are archived in commit `bbd53847` on
branch `claude/chamber-cylinder-diagnostics-archive`; figures that need them
appear as placeholders when the document is built from this folder.
