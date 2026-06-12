@page jacobian Jacobian Generator for svZeroDSolver

This tool generates C++ code for block implementations in the svZeroDSolver framework using symbolic mathematics.

## Overview

The script `script/jacobian.py` reads block definitions from YAML files and generates C++ code for implementing the mathematical models in the solver. It uses symbolic differentiation (via SymPy) to automatically derive the necessary Jacobian matrices and system contributions.

## Usage

The script depends on SymPy, which is part of the project's `dev` dependency
group. After running `uv sync` (see the [Developer Guide](@ref developer_guide)),
run it inside the managed environment:

```bash
uv run python scripts/jacobian.py <yaml_file>
```

## YAML File Format

The YAML files define the mathematical model of a block with the following structure:

### Required Sections

- `variables`: List of variable names in the model
- `derivatives`: List of derivative names (must match variables with `_dt` suffix)
- `constants`: List of parameter names
- `residuals`: List of residual equations that define the system

### Optional Sections

- `time_dependent`: List of parameters that depend on time (e.g., activation functions)
- `helper_functions`: Python code defining helper functions used in the residuals
- `time_symbol`: Name of an in-cycle time symbol (e.g. `t`) used by the
  residuals. It is neither a state variable nor a calibration parameter (it gets
  no Jacobian column), but it is available to the residual/helper expressions.
  The generated `update_gradient` reads it from `model->time`. The helper
  `warp_signed` (a periodic phase wrap with unit derivative) is available for
  activation functions and is emitted as a C++ lambda.

### Example

```yaml
variables:
  - Pin
  - Qin
  - Pout
  - Qout

derivatives:
  - dPin_dt
  - dQin_dt
  - dPout_dt
  - dQout_dt

constants:
  - R
  - C
  - L
  - S

residuals:
  - Pin - Pout - (R + S * abs(Qin)) * Qin - L * dQout_dt
  - Qin - Qout - C * dPin_dt + C * (R + 2 * S * abs(Qin)) * dQin_dt
```

## Output

The script generates four C++ function implementations:

1. `update_constant` - Sets up constant matrix coefficients for the system
2. `update_time` - Updates time-dependent parameters
3. `update_solution` - Computes solution-dependent terms and Jacobians
4. `update_gradient` - Residual and its parameter-Jacobian (`dr/dalpha`) used by
   the calibrator. Parameters are read from the `alpha` vector and columns are
   addressed via `global_param_ids[ParamId::<name>]`, so the `constants` order
   must match the block's `ParamId` enum.

## Workflow

1. Create a YAML file defining your block's mathematical model
2. Run `jacobian.py` on this file to generate C++ code
3. Copy the generated code to your block implementation file
4. Complete the implementation with necessary boilerplate code

## Tips

- Ensure the number of variables equals the number of derivatives
- The number of residuals should be equal to the number of variables minus 2
- Use helper functions for complex expressions to improve readability
- Define time-dependent constants separately

## Examples

See the provided YAML examples:
- `ChamberSphere.yaml` - Spherical heart chamber model
- `BloodVessel.yaml` - Blood vessel model with optional stenosis
- `ClosedLoopCoronaryBC.yaml` - Coronary boundary condition model