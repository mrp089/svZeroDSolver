import json
import os

import numpy as np
import pytest

from .utils import execute_pysvzerod, RTOL_PRES

this_file_dir = os.path.abspath(os.path.dirname(__file__))


def test_steady_flow_calibration():
    testfile = os.path.join(this_file_dir, "cases", "steadyFlow_calibration.json")

    result, _ = execute_pysvzerod(testfile, "calibrator")

    calibrated_parameters = result["vessels"][0]["zero_d_element_values"]

    assert np.isclose(
        np.mean(calibrated_parameters["R_poiseuille"]), 100, rtol=RTOL_PRES
    )
    assert np.isclose(np.mean(calibrated_parameters["C"]), 0.0001, rtol=RTOL_PRES)
    assert np.isclose(np.mean(calibrated_parameters["L"]), 1.0, rtol=RTOL_PRES)
    assert np.isclose(
        np.mean(calibrated_parameters["stenosis_coefficient"]), 0.0, rtol=RTOL_PRES
    )


def test_chamber_sphere_calibration():
    """Calibrate all twelve ChamberSphere parameters via the normal calibrator.

    ChamberSphere is calibrated through the standard point-wise calibrator using
    ``ChamberSphere::update_gradient``. The fixture
    ``chamber_sphere_calibration.json`` provides a full-state observation set
    (``y``/``dy`` for the ports and the internal variables radius, velo, stress,
    tau, volume) plus a ``t`` vector, constructed so every residual vanishes at
    ``_true_values``. The ``t`` vector lets the optimizer set ``model->time`` per
    observation, which makes the active-stress equation -- and hence the
    activation/timing parameters (sigma_max, alpha_max, alpha_min, tsys, tdias,
    steepness) -- identifiable alongside the six time-independent parameters. All
    twelve are recovered from a 20%-perturbed start.
    """
    testfile = os.path.join(
        this_file_dir, "cases", "chamber_sphere_calibration.json"
    )

    result, config = execute_pysvzerod(testfile, "calibrator")

    calibrated = result["vessels"][0]["zero_d_element_values"]
    for name in config["_calibrate_subset"]:
        true_value = config["_true_values"][name]
        assert np.isclose(calibrated[name], true_value, rtol=1e-6), (
            f"{name}: got {calibrated[name]}, expected {true_value}"
        )


@pytest.mark.parametrize(
    "test_case",
    [
        "0080_0001_calibrate_from_0d",
        "0104_0001_calibrate_from_0d",
        "0140_2001_calibrate_from_0d",
        # Calibrates only R_poiseuille via per-block ``calibrate`` fields,
        # starting from the reference values for C, L and
        # stenosis_coefficient and zeroed R_poiseuille. The calibrator should
        # recover R_poiseuille while leaving the other parameters untouched.
        "0104_0001_calibrate_R_only",
    ],
)
def test_calibration_vmr(test_case):
    """Calibrate a model from the vascular model repository and check that
    every parameter matches the corresponding reference."""
    test = os.path.join(
        this_file_dir, "cases", "vmr", "input", f"{test_case}.json"
    )
    model_id = test_case[:9]
    reference_file = os.path.join(
        this_file_dir, "cases", "vmr", "reference", f"{model_id}_optimal_from_0d.json"
    )
    with open(reference_file) as ff:
        reference = json.load(ff)

    result, _ = execute_pysvzerod(test, "calibrator")

    for i, vessel in enumerate(reference["vessels"]):
        for key, value in vessel["zero_d_element_values"].items():
            assert np.isclose(
                result["vessels"][i]["zero_d_element_values"][key],
                value,
                rtol=RTOL_PRES,
            )

    for i, junction in enumerate(reference["junctions"]):
        if "junction_values" in junction:
            for key, value in junction["junction_values"].items():
                assert np.allclose(
                    result["junctions"][i]["junction_values"][key],
                    value,
                    rtol=RTOL_PRES,
                )
