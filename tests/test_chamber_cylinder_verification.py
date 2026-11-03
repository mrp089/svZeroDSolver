"""Verification of the ChamberCylinder block against independent solutions of the
same continuum model (no reference data from genet23's figures). Fast subset of
scripts/ChamberCylinder_verification.py, which also holds the derivations."""
import math
import os
import sys

import numpy as np
import pytest

pytest.importorskip("scipy")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import ChamberCylinder_verification as CV  # noqa: E402

GEOM = {k: CV.G.GEOM[k] for k in ("Ri", "Re", "length", "alpha_endo", "alpha_epi")}
DEG = 180.0 / math.pi


def test_transmural_shear_vanishes_and_mirror_symmetry():
    """A rotation by pi about a radial axis maps body, fibres and loads onto
    themselves and (phi, eta) onto (-phi, -eta), so phi' = eta' = 0; mirroring
    the fibre angle flips the twist and nothing else."""
    r = CV.t1_symmetry(ne=3, ncycle=2)
    assert r["shear_phi_max_deg"] < 1e-9
    assert r["shear_eta_max_mm"] < 1e-9
    assert max(r["mirror_dV_rel"], r["mirror_dP_rel"], r["mirror_deps_rel"]) < 1e-8
    assert r["mirror_beta_sum_rel"] < 1e-8


@pytest.mark.parametrize("C5", [50.0, 0.0])
def test_static_inflation_extension_torsion(C5):
    """Passive static state vs. the semi-analytical solution (exact
    incompressibility, virtual work in (r_i, lambda, beta))."""
    mat = dict(CV.PASSIVE, C5=C5)
    ref = CV.static_reference(1000.0, mat, GEOM)
    blk = CV.static_block(1000.0, mat, GEOM, ne=12, Thold=10.0, dt=1e-2)
    assert ref["ok"]
    assert abs(blk["V"] / ref["V"] - 1) < 1e-3
    assert abs(blk["eps"] - (ref["q"][1] - 1)) < 5e-4
    assert abs(blk["beta"] - ref["q"][2]) * GEOM["length"] * DEG < 0.02


def test_torsional_free_vibration():
    """Small-amplitude torsion of the free tube vs. the closed-form frequency and
    decay (twist inertia with momentum-free rigid rotation, shear stiffness,
    tissue viscosity and the fibre dashpot of the unactivated active branch)."""
    wd, s, _ = CV.twist_reference(1e4, 50.0, 70.0, 1000.0)
    t, b = CV.vibration_block("twist", dt=1e-4, T=0.15, gamma=50.0, mu=70.0)
    fit = CV.damped_fit(t, b, wd, s)
    assert abs(fit["omega"] / wd - 1) < 1e-3
    assert abs(fit["sigma"] / s - 1) < 5e-3


def test_unit_scaling_invariance():
    """The same torsion test in (mm, kPa, s) units gives the same history."""
    _, b1 = CV.vibration_block("twist", dt=1e-4, T=0.05)
    _, b2 = CV.vibration_block("twist", scale=(1e3, 1e-3), dt=1e-4, T=0.05)
    assert np.abs(b1 - b2).max() < 1e-6 * np.abs(b1).max()
