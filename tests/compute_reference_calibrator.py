import os
import copy
import json
import pysvzerod
import argparse
import numpy as np
from scipy.interpolate import CubicSpline
import pdb

# Number of observations to extract in result refinement
NUM_OBS = 100

cases_solver = "cases_solver"
cases_calibrator = "cases_calibrator"
calibration_params = {
    "tolerance_gradient": 1e-5,
    "tolerance_increment": 1e-10,
    "maximum_iterations": 100,
    "calibrate_stenosis_coefficient": True,
    "set_capacitance_to_zero": False,
}


def refine(x: np.ndarray, y: np.ndarray, num: int) -> np.ndarray:
    """Refine a curve using cubic spline interpolation with derivative.

    Args:
        x: X-coordinates
        y: Y-coordinates
        num: New number of points of the refined data.


    Returns:
        new_y: New y-coordinates
        new_dy: New dy-coordinates
    """
    y = y.copy()
    y[-1] = y[0]
    x_new = np.linspace(x[0], x[-1], num)
    spline = CubicSpline(x, y, bc_type="periodic")
    new_y = spline(x_new)
    new_dy = spline.derivative()(x_new)
    return new_y.tolist(), new_dy.tolist()


def collect(out, res, name_vessel, name_other, is_out):
    for f in ["flow", "pressure"]:
        if is_out:
            name = ":".join([f, name_vessel, name_other])
            inout = "_out"
        else:
            name = ":".join([f, name_other, name_vessel])
            inout = "_in"
        i_vesel = res["name"] == name_vessel
        data = np.array(res[f + inout][i_vesel])
        time = np.array(res["time"][i_vesel])
        out["y"][name], out["dy"][name] = refine(time, data, NUM_OBS)


def collect_results_from_0d(inp, res):
    vessels = np.array(inp["vessels"])
    out = {"y": {}, "dy": {}}
    # collect connections between branches and junctions
    if "junctions" in inp:
        for jc in inp["junctions"]:
            for io in ["inlet", "outlet"]:
                for vs in vessels[jc[io + "_vessels"]]:
                    collect(
                        out, res, vs["vessel_name"], jc["junction_name"], io == "inlet"
                    )
    # collect connections between branches and boundary conditions
    for vs in inp["vessels"]:
        if "boundary_conditions" in vs:
            for bc_type, bc_name in vs["boundary_conditions"].items():
                collect(out, res, vs["vessel_name"], bc_name, bc_type == "outlet")

    # export results, time, and flow on same resampled time discretization
    flow = out["y"]["flow:INFLOW:branch0_seg0"]
    return out, np.linspace(0, np.max(res["time"]), len(flow)).tolist(), flow


def compare(geo, times):
    print("\n", geo, "\n\n")

    # run the estimation
    inp0, cali = estimate(geo, times)

    # compare 0D element values
    for i in range(len(inp0["vessels"])):
        for ele in cali["vessels"][i]["zero_d_element_values"].keys():
            sol = cali["vessels"][i]["zero_d_element_values"][ele]
            if ele in inp0["vessels"][i]["zero_d_element_values"]:
                ref = inp0["vessels"][i]["zero_d_element_values"][ele]
            else:
                ref = 0.0

            # calculate error
            if ref == 0.0:
                err = 1.0
            else:
                err = np.abs(sol / ref - 1.0)

            # print results
            out = str(i) + "\t" + ele[0]
            for j in [ref, sol]:
                out += "\t\t{:.1e}".format(j)
            out += "\t\t{:+d}".format(int(np.log(err)))
            # print(out)


def compute_ref_sol(test_name):
    # read input file
    with open(os.path.join(cases_solver, test_name)) as f:
        inp = json.load(f)

    # run forward simulation
    try:
        res = pysvzerod.simulate(inp)
    except RuntimeError as e:
        print("Simulation failed: ", e)
        return

    # collect results in format required for calibrator
    out, bc_time, bc_flow = collect_results_from_0d(inp, res)

    # replace inflow to match calibration data
    for bc in inp["boundary_conditions"]:
        if bc["bc_name"] == "INFLOW":
            bc["bc_values"]["t"] = bc_time
            bc["bc_values"]["Q"] = bc_flow

    # set all elements to zero
    for vessel in inp["vessels"]:
        for ele in vessel["zero_d_element_values"].keys():
            vessel["zero_d_element_values"][ele] = 0.0

    # add calibration parameters
    inp["calibration_parameters"] = calibration_params

    # only calibrate to last cycle
    inp["simulation_parameters"]["output_all_cycles"] = False

    # write output json
    fname_in = os.path.join(cases_calibrator, "input", test_name)
    with open(fname_in, "w") as f:
        json.dump(inp | out, f, indent=2)

    # run the calibration
    with open(fname_in) as ff:
        config = json.load(ff)
    try:
        cali = pysvzerod.calibrate(config)
    except RuntimeError as e:
        print("Calibration failed: ", e)
        return

    # write calibrated json
    fname_out = os.path.join(cases_calibrator, "results", test_name)
    with open(fname_out, "w") as f:
        json.dump(cali, f, indent=2)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute reference solution for a test case"
    )
    parser.add_argument("testname", help="name of the test case json file")
    args = parser.parse_args()
    compute_ref_sol(args.testname)
