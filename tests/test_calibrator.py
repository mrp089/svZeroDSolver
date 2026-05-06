import json
import os
import pytest

import numpy as np

from .utils import execute_pysvzerod

this_file_dir = os.path.abspath(os.path.dirname(__file__))

# The calibrator should recover the input parameters used to generate the y
# observations. With y from the forward solver and dy reconstructed via the
# gen-alpha update relation, that recovery is bound only by the forward
# solver's nonlinear tolerance — the steady case drifts by about 1e-9 on R,
# the VMR cases match to machine precision.
RTOL = 1e-10
ATOL = 1e-9


def test_steady_flow_calibration():
    testfile = os.path.join(this_file_dir, "cases", "steadyFlow_calibration.json")

    result, _ = execute_pysvzerod(testfile, "calibrator")

    calibrated_parameters = result["vessels"][0]["zero_d_element_values"]

    assert np.isclose(
        np.mean(calibrated_parameters["R_poiseuille"]), 100, rtol=RTOL, atol=ATOL
    )
    assert np.isclose(
        np.mean(calibrated_parameters["C"]), 0.0001, rtol=RTOL, atol=ATOL
    )
    assert np.isclose(np.mean(calibrated_parameters["L"]), 1.0, rtol=RTOL, atol=ATOL)
    assert np.isclose(
        np.mean(calibrated_parameters["stenosis_coefficient"]),
        0.0,
        rtol=RTOL,
        atol=ATOL,
    )


@pytest.mark.parametrize("model_id", ["0080_0001", "0104_0001", "0140_2001"])
def test_calibration_vmr(model_id):
    """Test actual models from the vascular model repository."""
    with open(
        os.path.join(
            this_file_dir, "cases", "vmr", "input", f"{model_id}_calibrate_from_0d.json"
        )
    ) as ff:
        reference = json.load(ff)

    test = os.path.join(
        this_file_dir, "cases", "vmr", "input", f"{model_id}_calibrate_from_0d.json"
    )

    result, _ = execute_pysvzerod(test, "calibrator")

    for i, vessel in enumerate(reference["vessels"]):
        for key, value in vessel["zero_d_element_values"].items():
            assert np.isclose(
                result["vessels"][i]["zero_d_element_values"][key],
                value,
                rtol=RTOL,
                atol=ATOL,
            )

    for i, junction in enumerate(reference["junctions"]):
        if "junction_values" in junction:
            for key, value in junction["junction_values"].items():
                assert np.allclose(
                    result["junctions"][i]["junction_values"][key],
                    value,
                    rtol=RTOL,
                    atol=ATOL,
                )
