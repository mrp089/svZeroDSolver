import json
import os
import pytest

import numpy as np

from .utils import execute_pysvzerod, RTOL_PRES

this_file_dir = os.path.abspath(os.path.dirname(__file__))


def test_steady_flow_calibration():
    testfile = os.path.join(this_file_dir, "cases", "steadyFlow_calibration.json")

    result, _ = execute_pysvzerod(testfile, "calibrator")

    # New format uses dict with vessel name as key, and renamed parameters
    calibrated_parameters = result["vessels"]["branch0_seg0"]["values"]

    assert np.isclose(calibrated_parameters["resistance"], 100, rtol=RTOL_PRES)
    assert np.isclose(calibrated_parameters["capacitance"], 0.0001, rtol=RTOL_PRES)
    assert np.isclose(calibrated_parameters["inductance"], 1.0, rtol=RTOL_PRES)
    assert np.isclose(calibrated_parameters["stenosis"], 0.0, rtol=RTOL_PRES)


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

    # New format uses dict with vessel/junction name as key
    for vessel_name, vessel_config in reference["vessels"].items():
        for key, value in vessel_config["values"].items():
            np.isclose(
                result["vessels"][vessel_name]["values"][key],
                value,
                rtol=RTOL_PRES,
            )

    if "junctions" in reference:
        for junction_name, junction_config in reference["junctions"].items():
            if "values" in junction_config:
                for key, value in junction_config["values"].items():
                    np.allclose(
                        result["junctions"][junction_name]["values"][key],
                        value,
                        rtol=RTOL_PRES,
                    )
