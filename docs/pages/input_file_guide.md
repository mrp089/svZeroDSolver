@page input_file_guide Input File Guide

[TOC]

This guide explains how to create JSON configuration files for svZeroDSolver. For working examples, see the test cases in `tests/cases/`.

# Overview

svZeroDSolver configuration files are JSON documents that define:
- Simulation parameters (time stepping, tolerances)
- Boundary conditions (inlets, outlets)
- Vessels (blood vessel segments)
- Junctions (connections between vessels)
- Optional: valves, chambers, external coupling blocks

# JSON Configuration Structure

A typical configuration file has the following structure:

```json
{
  "simulation_parameters": { ... },
  "boundary_conditions": [ ... ],
  "vessels": [ ... ],
  "junctions": [ ... ]
}
```

# Simulation Parameters

The `simulation_parameters` section controls time stepping and solver behavior:

```json
"simulation_parameters": {
  "number_of_cardiac_cycles": 10,
  "number_of_time_pts_per_cardiac_cycle": 100,
  "absolute_tolerance": 1e-8
}
```

## Required Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `number_of_cardiac_cycles` | int | Number of cardiac cycles to simulate |
| `number_of_time_pts_per_cardiac_cycle` | int | Time points per cycle (determines time step size) |

## Optional Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `absolute_tolerance` | float | 1e-8 | Absolute tolerance for the nonlinear solver |
| `cardiac_period` | float | auto | Cardiac cycle period (auto-detected from BC waveforms if not specified) |
| `steady_initial` | bool | true | Use steady-state initial conditions |
| `output_all_cycles` | bool | false | Output results for all cycles (not just last) |
| `output_mean_only` | bool | false | Output only cycle-averaged mean values |
| `output_derivative` | bool | false | Include time derivatives in output |
| `output_variable_based` | bool | false | Organize output by variable name |
| `use_cycle_to_cycle_error` | bool | false | Stop when cycle-to-cycle error is below threshold |
| `sim_cycle_to_cycle_percent_error` | float | 1.0 | Target cycle-to-cycle error percentage |

# Boundary Conditions

Boundary conditions define inlets and outlets. Each BC is an object in the `boundary_conditions` array:

```json
"boundary_conditions": [
  {
    "bc_name": "INFLOW",
    "bc_type": "FLOW",
    "bc_values": { ... }
  }
]
```

## Boundary Condition Types

### FLOW - Prescribed Flow Rate

Specifies flow rate as a function of time:

```json
{
  "bc_name": "INFLOW",
  "bc_type": "FLOW",
  "bc_values": {
    "Q": [0.0, 5.0, 5.0, 0.0],
    "t": [0.0, 0.1, 0.9, 1.0]
  }
}
```

For steady flow, use two identical values:
```json
"bc_values": {
  "Q": [5.0, 5.0],
  "t": [0.0, 1.0]
}
```

### PRESSURE - Prescribed Pressure

Specifies pressure as a function of time:

```json
{
  "bc_name": "OUTLET",
  "bc_type": "PRESSURE",
  "bc_values": {
    "P": [100.0, 100.0],
    "t": [0.0, 1.0]
  }
}
```

### RESISTANCE - Resistance Outlet

Models downstream vasculature as a simple resistance:

\f[
P = Q \cdot R + P_d
\f]

```json
{
  "bc_name": "OUTLET",
  "bc_type": "RESISTANCE",
  "bc_values": {
    "R": 1000.0,
    "Pd": 0.0
  }
}
```

| Parameter | Description |
|-----------|-------------|
| `R` | Resistance |
| `Pd` | Distal pressure |

### RCR - Windkessel Model

Three-element Windkessel model for downstream vasculature:

\f[
\begin{circuitikz}
\draw (0,0) to [R, l=$R_p$, *-] (2,0) to [R, l=$R_d$, -*] (4,0);
\draw (2,0) to [C, l=$C$] (2,-1.5) node[ground]{};
\end{circuitikz}
\f]

```json
{
  "bc_name": "OUTLET",
  "bc_type": "RCR",
  "bc_values": {
    "Rp": 100.0,
    "C": 0.001,
    "Rd": 1000.0,
    "Pd": 0.0
  }
}
```

| Parameter | Description |
|-----------|-------------|
| `Rp` | Proximal resistance |
| `C` | Capacitance |
| `Rd` | Distal resistance |
| `Pd` | Distal pressure |

### CORONARY - Coronary Outlet

Models coronary arteries with intramyocardial pressure effects:

```json
{
  "bc_name": "CORONARY",
  "bc_type": "CORONARY",
  "bc_values": {
    "Ra1": 100.0,
    "Ra2": 100.0,
    "Rv1": 100.0,
    "Ca": 0.0001,
    "Cc": 0.0001,
    "Pim": [0.0, 10.0, 0.0],
    "Pim_t": [0.0, 0.5, 1.0],
    "P_v": 0.0
  }
}
```

# Vessels

Blood vessels are defined in the `vessels` array. Each vessel has properties that determine its hemodynamic behavior.

```json
"vessels": [
  {
    "vessel_id": 0,
    "vessel_name": "aorta",
    "vessel_length": 10.0,
    "zero_d_element_type": "BloodVessel",
    "zero_d_element_values": {
      "R_poiseuille": 100.0
    },
    "boundary_conditions": {
      "inlet": "INFLOW",
      "outlet": "OUTLET"
    }
  }
]
```

## Vessel Properties

| Property | Type | Description |
|----------|------|-------------|
| `vessel_id` | int | Unique identifier for the vessel |
| `vessel_name` | string | Name of the vessel (used in output) |
| `vessel_length` | float | Length of the vessel segment |
| `zero_d_element_type` | string | Type of 0D element (usually "BloodVessel") |
| `zero_d_element_values` | object | Hemodynamic parameters |
| `boundary_conditions` | object | Connected BCs (for vessels at model boundaries) |

## Zero-D Element Values

The `zero_d_element_values` object specifies the vessel's hemodynamic properties:

| Parameter | Symbol | Description |
|-----------|--------|-------------|
| `R_poiseuille` | \f$R\f$ | Viscous resistance |
| `C` | \f$C\f$ | Capacitance (vessel wall compliance) |
| `L` | \f$L\f$ | Inductance (blood inertia) |
| `stenosis_coefficient` | \f$K_s\f$ | Stenosis coefficient for nonlinear pressure drop |

### Vessel Types by Parameters

- **R** (Resistance only): Simplest model, viscous resistance only
- **RC**: Resistance + Capacitance (wall compliance)
- **RL**: Resistance + Inductance (blood inertia)
- **RLC**: Full model with resistance, inductance, and capacitance

### Stenosis Model

When `stenosis_coefficient` is specified, the pressure drop includes a nonlinear term:

\f[
\Delta P = R \cdot Q + K_s \cdot Q \cdot |Q|
\f]

# Junctions

Junctions connect multiple vessels. They are defined in the `junctions` array:

```json
"junctions": [
  {
    "junction_name": "J0",
    "junction_type": "NORMAL_JUNCTION",
    "inlet_vessels": [0],
    "outlet_vessels": [1, 2]
  }
]
```

## Junction Types

### NORMAL_JUNCTION

Basic junction that enforces mass conservation and pressure continuity:

```json
{
  "junction_name": "J0",
  "junction_type": "NORMAL_JUNCTION",
  "inlet_vessels": [0],
  "outlet_vessels": [1, 2]
}
```

### BloodVesselJunction

Junction with additional resistance and inductance in each outlet branch:

```json
{
  "junction_name": "J0",
  "junction_type": "BloodVesselJunction",
  "inlet_vessels": [0],
  "outlet_vessels": [1, 2],
  "junction_values": {
    "R_poiseuille": [100.0, 200.0],
    "L": [0.0, 0.0],
    "stenosis_coefficient": [0.0, 0.0]
  }
}
```

## Network Topologies

### Bifurcation (1 inlet, multiple outlets)

One vessel splits into two or more vessels:

```json
"junctions": [
  {
    "junction_name": "bifurcation",
    "junction_type": "NORMAL_JUNCTION",
    "inlet_vessels": [0],
    "outlet_vessels": [1, 2]
  }
]
```

### Confluence (multiple inlets, 1 outlet)

Two or more vessels merge into one:

```json
"junctions": [
  {
    "junction_name": "confluence",
    "junction_type": "NORMAL_JUNCTION",
    "inlet_vessels": [0, 1],
    "outlet_vessels": [2]
  }
]
```

# Valves

Valves control flow direction between blocks. They are defined in the `valves` array:

```json
"valves": [
  {
    "name": "aortic_valve",
    "type": "ValveTanh",
    "params": {
      "Rmax": 100000.0,
      "Rmin": 100.0,
      "Steepness": 100.0,
      "upstream_block": "ventricle",
      "downstream_block": "aorta"
    }
  }
]
```

## ValveTanh

A valve with resistance that varies smoothly based on pressure gradient using a tanh function:

| Parameter | Description |
|-----------|-------------|
| `Rmax` | Maximum resistance (valve closed) |
| `Rmin` | Minimum resistance (valve open) |
| `Steepness` | Controls transition sharpness |
| `upstream_block` | Name of upstream block |
| `downstream_block` | Name of downstream block |

# Chambers

Cardiac chamber models are defined in the `chambers` array:

## ChamberElastanceInductor

Time-varying elastance model for ventricular contraction:

```json
"chambers": [
  {
    "name": "LV",
    "type": "ChamberElastanceInductor",
    "values": {
      "Emax": 1.0,
      "Emin": 0.1,
      "Vrd": 10.0,
      "Vrs": 5.0,
      "t_active": 0.3,
      "t_twitch": 0.5,
      "Impedance": 0.001
    }
  }
]
```

| Parameter | Description |
|-----------|-------------|
| `Emax` | Maximum elastance (end-systolic) |
| `Emin` | Minimum elastance (diastolic) |
| `Vrd` | Rest volume (diastole) |
| `Vrs` | Rest volume (systole) |
| `t_active` | Duration of active contraction |
| `t_twitch` | Total twitch duration |
| `Impedance` | Chamber impedance |

## ChamberSphere

Thick-walled sphere model with active and passive stress components.

# Complete Example

Here is a complete example of a simple configuration file:

```json
{
  "simulation_parameters": {
    "number_of_cardiac_cycles": 10,
    "number_of_time_pts_per_cardiac_cycle": 100,
    "absolute_tolerance": 1e-8
  },
  "boundary_conditions": [
    {
      "bc_name": "INFLOW",
      "bc_type": "FLOW",
      "bc_values": {
        "Q": [5.0, 5.0],
        "t": [0.0, 1.0]
      }
    },
    {
      "bc_name": "OUTLET",
      "bc_type": "RCR",
      "bc_values": {
        "Rp": 100.0,
        "C": 0.001,
        "Rd": 1000.0,
        "Pd": 0.0
      }
    }
  ],
  "vessels": [
    {
      "vessel_id": 0,
      "vessel_name": "aorta",
      "vessel_length": 10.0,
      "zero_d_element_type": "BloodVessel",
      "zero_d_element_values": {
        "R_poiseuille": 100.0,
        "C": 0.0001
      },
      "boundary_conditions": {
        "inlet": "INFLOW",
        "outlet": "OUTLET"
      }
    }
  ],
  "junctions": []
}
```

# Running a Simulation

## Python API

```python
import json
import pysvzerod

# Load configuration
with open("config.json") as f:
    config = json.load(f)

# Run simulation
result = pysvzerod.simulate(config)

# Result is a dictionary with flow and pressure time series
print(result.keys())
```

## Command Line

```bash
svzerodsolver config.json output.csv
```

# Additional Resources

- Test cases in `tests/cases/` provide working examples for various configurations
- [Developer Guide](@ref developer_guide) for contributing to svZeroDSolver
- [SimVascular documentation](https://simvascular.github.io/documentation/rom_simulation.html#0d-solver) for more details
