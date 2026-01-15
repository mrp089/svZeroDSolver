#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
# University of California, and others. SPDX-License-Identifier: BSD-3-Clause
"""Convert svZeroDSolver input files from old format to new simplified format.

Old format uses arrays for boundary_conditions, vessels, junctions, etc.
New format uses dictionaries with names as keys and standardized attribute naming.

Changes:
- boundary_conditions: array -> dict, bc_name -> key, bc_type -> type, bc_values -> values
- vessels: array -> dict, vessel_name -> key, zero_d_element_type -> type,
           zero_d_element_values -> values, remove vessel_id and vessel_length
- junctions: array -> dict, junction_name -> key, junction_type -> type,
             junction_values -> values, inlet_vessels/outlet_vessels -> inlet/outlet
- closed_loop_blocks: array -> dict with generated names, closed_loop_type -> type,
                      parameters -> values
- valves: array -> dict, params -> values
- chambers: already using correct format, just convert to dict
- external_solver_coupling_blocks: already using correct format, just convert to dict
- Add top-level connections array for BC-to-vessel connections
"""

import argparse
import json
import sys
from pathlib import Path


def convert_boundary_conditions(old_bcs: list) -> dict:
    """Convert boundary_conditions from array to dict format."""
    new_bcs = {}
    for bc in old_bcs:
        name = bc["bc_name"]
        new_bcs[name] = {
            "type": bc["bc_type"],
            "values": bc["bc_values"]
        }
    return new_bcs


def convert_vessels(old_vessels: list, connections: list) -> dict:
    """Convert vessels from array to dict format.

    Also extracts boundary condition connections to the connections list.
    Returns vessel_id_map for junction conversion.
    """
    new_vessels = {}
    vessel_id_map = {}

    for vessel in old_vessels:
        name = vessel["vessel_name"]
        vessel_id = vessel.get("vessel_id")

        if vessel_id is not None:
            vessel_id_map[vessel_id] = name

        new_vessels[name] = {
            "type": vessel["zero_d_element_type"],
            "values": vessel["zero_d_element_values"]
        }

        # Extract boundary condition connections
        if "boundary_conditions" in vessel:
            bcs = vessel["boundary_conditions"]
            if "inlet" in bcs:
                connections.append([bcs["inlet"], name])
            if "outlet" in bcs:
                connections.append([name, bcs["outlet"]])

    return new_vessels, vessel_id_map


def convert_junctions(old_junctions: list, vessel_id_map: dict) -> dict:
    """Convert junctions from array to dict format.

    Converts inlet_vessels/outlet_vessels (using IDs) to inlet/outlet (using names).
    Junction connections are encoded in the junction's inlet/outlet arrays,
    not in the global connections list.
    """
    new_junctions = {}

    for junction in old_junctions:
        name = junction["junction_name"]

        new_junction = {
            "type": junction["junction_type"]
        }

        # Convert junction_values if present
        if "junction_values" in junction:
            new_junction["values"] = junction["junction_values"]

        # Convert vessel IDs to names, or use block names directly
        if "inlet_vessels" in junction:
            inlet_names = [vessel_id_map[vid] for vid in junction["inlet_vessels"]]
            new_junction["inlet"] = inlet_names
        elif "inlet_blocks" in junction:
            new_junction["inlet"] = junction["inlet_blocks"]

        if "outlet_vessels" in junction:
            outlet_names = [vessel_id_map[vid] for vid in junction["outlet_vessels"]]
            new_junction["outlet"] = outlet_names
        elif "outlet_blocks" in junction:
            new_junction["outlet"] = junction["outlet_blocks"]

        new_junctions[name] = new_junction

    return new_junctions


def convert_closed_loop_blocks(old_blocks: list) -> dict:
    """Convert closed_loop_blocks from array to dict format.

    The outlet_blocks field is kept within the block definition
    and processed by the C++ code, not added to the global connections array.
    """
    new_blocks = {}

    for i, block in enumerate(old_blocks):
        # Generate a name for the closed loop block
        block_type = block["closed_loop_type"]
        if block_type == "ClosedLoopHeartAndPulmonary":
            name = "CLH"  # Standard name used in code
        else:
            name = f"closed_loop_{i}"

        new_block = {
            "type": block_type,
            "values": block.get("parameters", {})
        }

        # Copy cardiac_cycle_period if present
        if "cardiac_cycle_period" in block:
            new_block["cardiac_cycle_period"] = block["cardiac_cycle_period"]

        # Keep outlet_blocks as outlet (renamed)
        if "outlet_blocks" in block:
            new_block["outlet"] = block["outlet_blocks"]

        new_blocks[name] = new_block

    return new_blocks


def convert_valves(old_valves: list) -> dict:
    """Convert valves from array to dict format.

    Valve connections (upstream_block, downstream_block) are kept in the values
    and processed by the C++ code, not added to the global connections array.
    """
    new_valves = {}

    for valve in old_valves:
        name = valve["name"]
        # Rename 'params' to 'values', keeping upstream_block/downstream_block
        new_valves[name] = {
            "type": valve["type"],
            "values": valve.get("params", {})
        }

    return new_valves


def convert_chambers(old_chambers: list) -> dict:
    """Convert chambers from array to dict format."""
    new_chambers = {}

    for chamber in old_chambers:
        name = chamber["name"]
        new_chambers[name] = {
            "type": chamber["type"],
            "values": chamber.get("values", {})
        }

    return new_chambers


def convert_external_coupling(old_blocks: list) -> dict:
    """Convert external_solver_coupling_blocks from array to dict format."""
    new_blocks = {}

    for block in old_blocks:
        name = block["name"]
        new_blocks[name] = {
            "type": block["type"],
            "location": block["location"],
            "connected_block": block["connected_block"],
            "values": block.get("values", {})
        }
        # Copy periodic if present
        if "periodic" in block:
            new_blocks[name]["periodic"] = block["periodic"]

    return new_blocks


def convert_input_file(old_config: dict) -> dict:
    """Convert an svZeroDSolver input file from old to new format."""
    new_config = {}
    connections = []
    vessel_id_map = {}

    # Copy description and simulation_parameters as-is
    if "description" in old_config:
        new_config["description"] = old_config["description"]

    if "simulation_parameters" in old_config:
        new_config["simulation_parameters"] = old_config["simulation_parameters"]

    # Convert boundary_conditions
    if "boundary_conditions" in old_config:
        new_config["boundary_conditions"] = convert_boundary_conditions(
            old_config["boundary_conditions"]
        )

    # Convert vessels (and extract connections, build vessel_id_map)
    if "vessels" in old_config:
        new_config["vessels"], vessel_id_map = convert_vessels(
            old_config["vessels"], connections
        )

    # Convert junctions
    if "junctions" in old_config and old_config["junctions"]:
        new_config["junctions"] = convert_junctions(
            old_config["junctions"], vessel_id_map
        )

    # Convert closed_loop_blocks
    if "closed_loop_blocks" in old_config:
        new_config["closed_loop_blocks"] = convert_closed_loop_blocks(
            old_config["closed_loop_blocks"]
        )

    # Convert valves
    if "valves" in old_config:
        new_config["valves"] = convert_valves(old_config["valves"])

    # Convert chambers
    if "chambers" in old_config:
        new_config["chambers"] = convert_chambers(old_config["chambers"])

    # Convert external_solver_coupling_blocks
    if "external_solver_coupling_blocks" in old_config:
        new_config["external_solver_coupling_blocks"] = convert_external_coupling(
            old_config["external_solver_coupling_blocks"]
        )

    # Add connections
    if connections:
        new_config["connections"] = connections

    # Copy initial conditions as-is
    if "initial_condition" in old_config:
        new_config["initial_condition"] = old_config["initial_condition"]

    if "initial_condition_d" in old_config:
        new_config["initial_condition_d"] = old_config["initial_condition_d"]

    return new_config


def main():
    parser = argparse.ArgumentParser(
        description="Convert svZeroDSolver input files from old to new format"
    )
    parser.add_argument("input_file", help="Input JSON file (old format)")
    parser.add_argument(
        "-o", "--output",
        help="Output JSON file (new format). If not specified, prints to stdout"
    )
    parser.add_argument(
        "--inplace", "-i",
        action="store_true",
        help="Modify the input file in place"
    )

    args = parser.parse_args()

    # Read input file
    input_path = Path(args.input_file)
    with open(input_path) as f:
        old_config = json.load(f)

    # Convert
    new_config = convert_input_file(old_config)

    # Output
    output_json = json.dumps(new_config, indent=2)

    if args.inplace:
        with open(input_path, "w") as f:
            f.write(output_json)
            f.write("\n")
        print(f"Converted {input_path} in place", file=sys.stderr)
    elif args.output:
        output_path = Path(args.output)
        with open(output_path, "w") as f:
            f.write(output_json)
            f.write("\n")
        print(f"Wrote {output_path}", file=sys.stderr)
    else:
        print(output_json)


if __name__ == "__main__":
    main()
