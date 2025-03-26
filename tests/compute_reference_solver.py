import os
import json
import argparse
import pysvzerod

from .utils import cases_solver

def compute_ref_sol(testname):
    """
    compute reference solution for a test case

    :param testname: name of the test case json to compute reference solution for
    """
    # compute result
    test_filename = os.path.join(cases_solver, testname)
    result = pysvzerod.simulate(json.load(open(test_filename)))

    # save result
    result_filename = os.path.join(cases_solver, "results", "result_" + testname)

    # save to json
    with open(result_filename, "w") as f:
        f.write(result.to_json())

    # print for confirmation
    print(
        f"Reference solution for test case {testname} computed and saved to {result_filename}. Please verify that the results are as expected."
    )

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compute reference solution for a test case')
    parser.add_argument('testname', help='name of the test case json file')
    args = parser.parse_args()
    compute_ref_sol(args.testname)