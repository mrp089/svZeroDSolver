# Test Cases

This directory contains test cases for validating svZeroDSolver functionality. Each test case is a complete, working example that can also serve as a starting point for your own simulations.

For comprehensive documentation on creating input files, see the [Input File Guide](../../docs/pages/input_file_guide.md) in the documentation.

## File Naming Convention

Test files follow a naming pattern:

```
<flowType>_<vesselProperties>_<outletBC>.json
```

- **flowType**: `steadyFlow` (constant) or `pulsatileFlow` (time-varying)
- **vesselProperties**: R, RC, RL, RLC, stenosis, bifurcation, confluence
- **outletBC**: R (resistance), RCR, coronary, steadyPressure

## Test Categories

### Basic Single-Vessel Tests

| File | Description |
|------|-------------|
| `steadyFlow_R_R.json` | Single resistive vessel with steady inflow and resistance outlet |
| `steadyFlow_R_RCR.json` | Single resistive vessel with RCR (Windkessel) outlet |
| `steadyFlow_R_steadyPressure.json` | Single resistive vessel with fixed outlet pressure |
| `steadyFlow_R_coronary.json` | Single resistive vessel with coronary outlet |
| `steadyFlow_RC_R.json` | Vessel with resistance + capacitance |
| `steadyFlow_RL_R.json` | Vessel with resistance + inductance |
| `steadyFlow_RLC_R.json` | Vessel with resistance + inductance + capacitance |
| `steadyFlow_stenosis_R.json` | Vessel with stenosis (nonlinear pressure drop) |

### Pulsatile Flow Tests

| File | Description |
|------|-------------|
| `pulsatileFlow_R_RCR.json` | Pulsatile sinusoidal inflow with RCR outlet |
| `pulsatileFlow_R_coronary.json` | Pulsatile inflow with coronary outlet |
| `pulsatileFlow_CStenosis_steadyPressure.json` | Vessel with capacitance and stenosis |

### Network Topology Tests

| File | Description |
|------|-------------|
| `steadyFlow_bifurcationR_R1.json` | Symmetric bifurcation (equal downstream resistances) |
| `steadyFlow_bifurcationR_R2.json` | Asymmetric bifurcation (different resistances) |
| `steadyFlow_confluenceR_R.json` | Two vessels merging into one |
| `steadyFlow_blood_vessel_junction.json` | Bifurcation with BloodVesselJunction |

### Closed-Loop Heart Models

| File | Description |
|------|-------------|
| `closedLoopHeart_singleVessel.json` | Heart model with single aorta vessel |
| `closedLoopHeart_withCoronaries.json` | Heart model with coronary circulation |
| `coupledBlock_closedLoopHeart_*.json` | Heart models with external coupling |

### Valve and Chamber Models

| File | Description |
|------|-------------|
| `valve_tanh.json` | Two vessels connected by a tanh valve |
| `chamber_elastance_inductor.json` | Time-varying elastance chamber model |
| `chamber_sphere.json` | Thick-walled sphere chamber model |

### Output Format Tests

Files with these suffixes test different output options:

| Suffix | Simulation Parameter |
|--------|---------------------|
| `_variable` | `output_variable_based: true` |
| `_derivative` | `output_derivative: true` |
| `_mean` | `output_mean_only: true` |
| `_cycle_error` | `use_cycle_to_cycle_error: true` |

## Running Tests

```bash
# Run all solver tests
python -m pytest tests/test_solver.py -v

# Run a specific test
python -m pytest "tests/test_solver.py::test_solver[steadyFlow_R_R.json]" -v
```

## Reference Results

Reference solutions are stored in the `results/` subdirectory as JSON files (`result_<testname>.json`). These are used to validate that code changes don't affect simulation results.
