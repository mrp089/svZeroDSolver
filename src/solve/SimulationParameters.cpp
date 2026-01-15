// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "SimulationParameters.h"

bool get_param_scalar(const svzero_json& data, const std::string& name,
                      const InputParameter& param, double& val) {
  if (data.contains(name)) {
    val = data[name];
  } else {
    if (param.is_optional) {
      val = param.default_val;
    } else {
      return true;
    }
  }
  return false;
}

bool get_param_vector(const svzero_json& data, const std::string& name,
                      const InputParameter& param, std::vector<double>& val) {
  if (data.contains(name)) {
    val = data[name].get<std::vector<double>>();
  } else {
    if (param.is_optional) {
      val = {param.default_val};
    } else {
      return true;
    }
  }
  return false;
}

bool has_parameter(
    const std::vector<std::pair<std::string, InputParameter>>& params,
    const std::string& name) {
  for (const auto& pair : params) {
    if (pair.first == name) {
      return true;
    }
  }
  return false;
}

int generate_block(Model& model, const svzero_json& block_params_json,
                   const std::string& block_type, const std::string_view& name,
                   bool internal, bool periodic) {
  // Generate block from factory
  auto block = model.create_block(block_type);

  // Read block input parameters
  std::vector<int> block_param_ids;
  int new_id;
  int err;

  // Check that all parameters defined for the current block are valid
  for (auto& el : block_params_json.items()) {
    // Ignore comments (starting with _)
    if (el.key()[0] == '_') {
      continue;
    }

    // Check if json input is a valid parameter for the current block
    if (!has_parameter(block->input_params, el.key())) {
      throw std::runtime_error("Unknown parameter " + el.key() +
                               " defined in " + block_type + " block " +
                               static_cast<std::string>(name));
    }
  }

  // The rest of this function reads the parameters for each block, adds them to
  // the model, and stores the corresponding param IDs in each block

  // Handle input parameters given as a list differently
  if (block->input_params_list) {
    for (const auto& block_param : block->input_params) {
      // todo: check error here
      for (double value : block_params_json[block_param.first]) {
        block_param_ids.push_back(model.add_parameter(value));
      }
    }
  } else {
    for (const auto& block_param : block->input_params) {
      // Time parameter is read at the same time as time-dependent value
      if (block_param.first.compare("time") == 0) {
        continue;
      }

      // Skip reading parameters that are not a number
      if (!block_param.second.is_number) {
        continue;
      }

      // Get vector parameter
      if (block_param.second.is_array) {
        // Get parameter vector
        std::vector<double> val;
        err = get_param_vector(block_params_json, block_param.first,
                               block_param.second, val);
        if (err) {
          throw std::runtime_error("Array parameter " + block_param.first +
                                   " is mandatory in " + block_type +
                                   " block " + static_cast<std::string>(name));
        }

        // Get time vector
        InputParameter t_param{false, true};
        std::vector<double> time;
        err = get_param_vector(block_params_json, "time", t_param, time);
        if (err) {
          throw std::runtime_error("Array parameter time is mandatory in " +
                                   block_type + " block " +
                                   static_cast<std::string>(name));
        }

        // Add parameters to model
        new_id = model.add_parameter(time, val, periodic);
      }

      // Get scalar parameter
      else {
        double val;
        err = get_param_scalar(block_params_json, block_param.first,
                               block_param.second, val);
        if (err) {
          throw std::runtime_error("Scalar parameter " + block_param.first +
                                   " is mandatory in " + block_type +
                                   " block " + static_cast<std::string>(name));
        }

        // Add parameter to model
        new_id = model.add_parameter(val);
      }
      // Store parameter IDs
      block_param_ids.push_back(new_id);
    }
  }

  // Add block to model (with parameter IDs)
  return model.add_block(block, name, block_param_ids, internal);
}

void validate_input(const svzero_json& config) {
  if (!config.contains("simulation_parameters")) {
    throw std::runtime_error("Define simulation_parameters");
  }
  if (!config.contains("boundary_conditions")) {
    throw std::runtime_error("Define at least one boundary condition");
  }
}

SimulationParameters load_simulation_params(const svzero_json& config) {
  DEBUG_MSG("Loading simulation parameters");
  SimulationParameters sim_params;
  const auto& sim_config = config["simulation_parameters"];
  sim_params.sim_coupled = sim_config.value("coupled_simulation", false);

  if (!sim_params.sim_coupled) {
    sim_params.sim_num_cycles = sim_config["number_of_cardiac_cycles"];
    sim_params.sim_pts_per_cycle =
        sim_config["number_of_time_pts_per_cardiac_cycle"];
    sim_params.sim_num_time_steps =
        (sim_params.sim_pts_per_cycle - 1) * sim_params.sim_num_cycles + 1;
    sim_params.use_cycle_to_cycle_error =
        sim_config.value("use_cycle_to_cycle_error", false);
    if (sim_params.use_cycle_to_cycle_error) {
      assert(sim_params.sim_num_cycles >=
             2);  // need at least two cycles to compute cycle-to-cycle error
      sim_params.sim_cycle_to_cycle_error =
          sim_config.value("sim_cycle_to_cycle_percent_error", 1.0) / 100;
    }
    sim_params.sim_external_step_size = 0.0;
  } else {
    sim_params.sim_num_cycles = 1;
    sim_params.sim_num_time_steps = sim_config["number_of_time_pts"];
    sim_params.sim_pts_per_cycle = sim_params.sim_num_time_steps;
    sim_params.sim_external_step_size =
        sim_config.value("external_step_size", 0.1);
  }
  sim_params.sim_abs_tol = sim_config.value("absolute_tolerance", 1e-8);
  sim_params.sim_nliter = sim_config.value("maximum_nonlinear_iterations", 30);
  sim_params.sim_steady_initial = sim_config.value("steady_initial", true);
  sim_params.sim_rho_infty = sim_config.value("rho_infty", 0.5);
  sim_params.output_variable_based =
      sim_config.value("output_variable_based", false);
  sim_params.output_interval = sim_config.value("output_interval", 1);
  sim_params.output_mean_only = sim_config.value("output_mean_only", false);
  sim_params.output_derivative = sim_config.value("output_derivative", false);
  sim_params.output_all_cycles = sim_config.value("output_all_cycles", false);
  sim_params.sim_cardiac_period = sim_config.value("cardiac_period", -1.0);
  DEBUG_MSG("Finished loading simulation parameters");
  return sim_params;
}

void load_simulation_model(const svzero_json& config, Model& model) {
  DEBUG_MSG("Loading model");
  // Create list to store block connections while generating blocks
  std::vector<std::tuple<std::string, std::string>> connections;

  // Move through different component names
  std::string component;

  // Create vessels
  DEBUG_MSG("Loading vessels");
  component = "vessels";
  std::set<std::string> vessel_names;
  if (config.contains(component)) {
    create_vessels(model, config, component, vessel_names);
  }

  // Create map for boundary conditions to boundary condition type
  component = "boundary_conditions";
  std::map<std::string, std::string> bc_type_map;
  for (auto& [bc_name, bc_config] : config[component].items()) {
    std::string bc_type = bc_config["type"];
    bc_type_map.insert({bc_name, bc_type});
  }

  // Create external coupling blocks
  component = "external_solver_coupling_blocks";
  if (config.contains(component)) {
    create_external_coupling(model, connections, config, component,
                             vessel_names, bc_type_map);
  }

  // Create boundary conditions
  component = "boundary_conditions";
  std::vector<std::string> closed_loop_bcs;
  create_boundary_conditions(model, config, component, bc_type_map,
                             closed_loop_bcs);

  // Create explicit junctions (e.g., BloodVesselJunction with parameters)
  component = "junctions";
  if (config.contains(component)) {
    create_junctions(model, connections, config, component);
  }

  // Collect top-level connections and analyze for auto-junction creation
  std::vector<std::tuple<std::string, std::string>> raw_connections;
  if (config.contains("connections")) {
    for (const auto& conn : config["connections"]) {
      std::string upstream = conn[0];
      std::string downstream = conn[1];
      raw_connections.push_back({upstream, downstream});

      // Update vessel types based on BC connections
      bool upstream_is_bc = bc_type_map.count(upstream) > 0;
      bool downstream_is_bc = bc_type_map.count(downstream) > 0;
      bool upstream_is_vessel = vessel_names.count(upstream) > 0;
      bool downstream_is_vessel = vessel_names.count(downstream) > 0;

      // If BC connects to vessel inlet, mark vessel as inlet type
      if (upstream_is_bc && downstream_is_vessel) {
        Block* vessel = model.get_block(downstream);
        if (vessel->vessel_type == VesselType::outlet) {
          vessel->update_vessel_type(VesselType::both);
        } else {
          vessel->update_vessel_type(VesselType::inlet);
        }
      }

      // If vessel connects to BC outlet, mark vessel as outlet type
      if (upstream_is_vessel && downstream_is_bc) {
        Block* vessel = model.get_block(upstream);
        if (vessel->vessel_type == VesselType::inlet) {
          vessel->update_vessel_type(VesselType::both);
        } else {
          vessel->update_vessel_type(VesselType::outlet);
        }
      }
    }
  }

  // Auto-create NORMAL_JUNCTIONs where blocks have multiple inlets or outlets
  // Build maps: block -> list of downstream blocks, block -> list of upstream
  // blocks
  std::map<std::string, std::vector<std::string>> downstream_map;
  std::map<std::string, std::vector<std::string>> upstream_map;
  for (const auto& conn : raw_connections) {
    downstream_map[std::get<0>(conn)].push_back(std::get<1>(conn));
    upstream_map[std::get<1>(conn)].push_back(std::get<0>(conn));
  }

  // Track which junctions we create to avoid duplicates
  std::set<std::string> auto_junctions;

  // Create junctions for blocks with multiple outlets (bifurcations)
  for (const auto& [block_name, downstreams] : downstream_map) {
    if (downstreams.size() > 1) {
      // Check if block already exists (could be an explicit junction)
      if (model.has_block(block_name)) {
        // Block exists, check if it's already a junction
        Block* block = model.get_block(block_name);
        if (block->block_class == BlockClass::junction) {
          // Already a junction, just add connections
          for (const auto& ds : downstreams) {
            connections.push_back({block_name, ds});
          }
          continue;
        }
      }
      // Need to create an auto-junction after this block
      std::string junction_name = "J_" + block_name + "_outlet";
      if (auto_junctions.find(junction_name) == auto_junctions.end()) {
        generate_block(model, {}, "NORMAL_JUNCTION", junction_name);
        auto_junctions.insert(junction_name);
        DEBUG_MSG("Auto-created junction " << junction_name);
      }
      // Connect: block -> junction -> each downstream
      connections.push_back({block_name, junction_name});
      for (const auto& ds : downstreams) {
        connections.push_back({junction_name, ds});
      }
    } else if (downstreams.size() == 1) {
      // Single downstream, direct connection
      connections.push_back({block_name, downstreams[0]});
    }
  }

  // Create junctions for blocks with multiple inlets (confluences)
  for (const auto& [block_name, upstreams] : upstream_map) {
    if (upstreams.size() > 1) {
      // Check if connections already go through an auto-junction
      bool already_handled = false;
      for (const auto& us : upstreams) {
        std::string potential_junction = "J_" + us + "_outlet";
        if (auto_junctions.find(potential_junction) != auto_junctions.end()) {
          // This upstream already has an outlet junction, check if it connects
          // to block_name
          already_handled = true;
          break;
        }
      }
      if (already_handled) {
        continue;  // Connections already added via outlet junction
      }

      // Check if block already exists and is a junction
      if (model.has_block(block_name)) {
        Block* block = model.get_block(block_name);
        if (block->block_class == BlockClass::junction) {
          continue;  // Already handled by explicit junction
        }
      }

      // Need to create an auto-junction before this block
      std::string junction_name = "J_" + block_name + "_inlet";
      if (auto_junctions.find(junction_name) == auto_junctions.end()) {
        generate_block(model, {}, "NORMAL_JUNCTION", junction_name);
        auto_junctions.insert(junction_name);
        DEBUG_MSG("Auto-created junction " << junction_name);
      }
      // Connect: each upstream -> junction -> block
      for (const auto& us : upstreams) {
        // Remove any existing direct connection that was added
        auto it = std::find(connections.begin(), connections.end(),
                            std::make_tuple(us, block_name));
        if (it != connections.end()) {
          connections.erase(it);
        }
        connections.push_back({us, junction_name});
      }
      connections.push_back({junction_name, block_name});
    }
    // Single upstream connections are already handled in the downstream loop
  }

  // Create closed-loop blocks
  component = "closed_loop_blocks";
  if (config.contains(component)) {
    create_closed_loop(model, connections, config, component, closed_loop_bcs);
  }

  // Create valves
  component = "valves";
  if (config.contains(component)) {
    create_valves(model, connections, config, component);
  }

  // Create chambers
  component = "chambers";
  if (config.contains(component)) {
    create_chambers(model, config, component);
  }

  // Create Connections
  for (auto& connection : connections) {
    auto ele1 = model.get_block(std::get<0>(connection));
    auto ele2 = model.get_block(std::get<1>(connection));
    model.add_node({ele1}, {ele2}, ele1->get_name() + ":" + ele2->get_name());
  }

  // Finalize model
  model.finalize();
}

void create_vessels(Model& model, const svzero_json& config,
                    const std::string& component,
                    std::set<std::string>& vessel_names) {
  // Loop all vessels (dictionary format: name -> {type, values})
  for (auto& [vessel_name, vessel_config] : config[component].items()) {
    const auto& vessel_values = vessel_config["values"];
    const std::string vessel_type = vessel_config["type"];
    vessel_names.insert(vessel_name);

    generate_block(model, vessel_values, vessel_type, vessel_name);
    DEBUG_MSG("Created vessel " << vessel_name);
  }
}

void create_boundary_conditions(Model& model, const svzero_json& config,
                                const std::string& component,
                                std::map<std::string, std::string>& bc_type_map,
                                std::vector<std::string>& closed_loop_bcs) {
  // Loop all boundary conditions (dictionary format: name -> {type, values})
  for (auto& [bc_name, bc_config] : config[component].items()) {
    std::string bc_type = bc_config["type"];
    const auto& bc_values = bc_config["values"];

    int block_id = generate_block(model, bc_values, bc_type, bc_name);

    // Keep track of closed loop blocks
    Block* block = model.get_block(block_id);

    if (block->block_type == BlockType::windkessel_bc) {
      model.update_has_windkessel_bc(true);
      double Rd = bc_values["resistance_distal"];
      double C = bc_values["capacitance"];
      double time_constant = Rd * C;
      model.update_largest_windkessel_time_constant(std::max(
          model.get_largest_windkessel_time_constant(), time_constant));
    }

    if (block->block_type == BlockType::closed_loop_rcr_bc) {
      if (bc_values["closed_loop_outlet"] == true) {
        closed_loop_bcs.push_back(bc_name);
      }
    } else if (block->block_class == BlockClass::closed_loop) {
      closed_loop_bcs.push_back(bc_name);
    }
    DEBUG_MSG("Created boundary condition " << bc_name);
  }
}

void create_external_coupling(
    Model& model,
    std::vector<std::tuple<std::string, std::string>>& connections,
    const svzero_json& config, const std::string& component,
    std::set<std::string>& vessel_names,
    std::map<std::string, std::string>& bc_type_map) {
  // Loop all external coupling blocks (dictionary format: name -> {...})
  for (auto& [coupling_name, coupling_config] : config[component].items()) {
    std::string coupling_type = coupling_config["type"];
    std::string coupling_loc = coupling_config["location"];
    bool periodic = coupling_config.value("periodic", true);
    const auto& coupling_values = coupling_config["values"];
    const bool internal = false;

    generate_block(model, coupling_values, coupling_type, coupling_name,
                   internal, periodic);

    // Determine the type of connected block
    std::string connected_block = coupling_config["connected_block"];
    std::string connected_type;
    int found_block = 0;
    if (connected_block == "ClosedLoopHeartAndPulmonary") {
      connected_type = "ClosedLoopHeartAndPulmonary";
      found_block = 1;
    } else {
      try {
        connected_type = bc_type_map.at(connected_block);
        found_block = 1;
      } catch (...) {
      }
      if (found_block == 0) {
        // Search for connected_block in the set of vessel names
        if (vessel_names.count(connected_block) > 0) {
          connected_type = "BloodVessel";
          found_block = 1;
        }
      }
      if (found_block == 0) {
        std::cout << "Error! Could not find connected type for block: "
                  << connected_block << std::endl;
        throw std::runtime_error("Terminating.");
      }
    }  // connected_block != "ClosedLoopHeartAndPulmonary"
    // Create connections
    if (coupling_loc == "inlet") {
      std::vector<std::string> possible_types = {"RESISTANCE",
                                                 "RCR",
                                                 "ClosedLoopRCR",
                                                 "SimplifiedRCR",
                                                 "CORONARY",
                                                 "ClosedLoopCoronaryLeft",
                                                 "ClosedLoopCoronaryRight",
                                                 "BloodVessel"};
      if (std::find(std::begin(possible_types), std::end(possible_types),
                    connected_type) == std::end(possible_types)) {
        throw std::runtime_error(
            "Error: The specified connection type for inlet"
            "external_coupling_block is invalid.");
      }
      connections.push_back({coupling_name, connected_block});
    } else if (coupling_loc == "outlet") {
      std::vector<std::string> possible_types = {
          "ClosedLoopRCR", "ClosedLoopHeartAndPulmonary", "BloodVessel"};
      if (std::find(std::begin(possible_types), std::end(possible_types),
                    connected_type) == std::end(possible_types)) {
        throw std::runtime_error(
            "Error: The specified connection type for outlet "
            "external_coupling_block is invalid.");
      }
      // Add connection only for closedLoopRCR and BloodVessel. Connection to
      // ClosedLoopHeartAndPulmonary will be handled in
      // ClosedLoopHeartAndPulmonary creation.
      if ((connected_type == "ClosedLoopRCR") ||
          (connected_type == "BloodVessel")) {
        connections.push_back({connected_block, coupling_name});
      }  // connected_type == "ClosedLoopRCR"
    }  // coupling_loc
    DEBUG_MSG("Created external coupling block " << coupling_name);
  }
}

void create_junctions(
    Model& model,
    std::vector<std::tuple<std::string, std::string>>& connections,
    const svzero_json& config, const std::string& component) {
  // Loop all junctions (dictionary format: name -> {type, values, inlet,
  // outlet})
  for (auto& [junction_name, junction_config] : config[component].items()) {
    std::string j_type = junction_config["type"];

    if (!junction_config.contains("values")) {
      generate_block(model, {}, j_type, junction_name);
    } else {
      generate_block(model, junction_config["values"], j_type, junction_name);
    }

    // Check for connections to inlets and outlets (using block names)
    if (junction_config.contains("inlet")) {
      for (const auto& block_name : junction_config["inlet"]) {
        connections.push_back({block_name, junction_name});
      }
    }
    if (junction_config.contains("outlet")) {
      for (const auto& block_name : junction_config["outlet"]) {
        connections.push_back({junction_name, block_name});
      }
    }
    DEBUG_MSG("Created junction " << junction_name);
  }
}

void create_closed_loop(
    Model& model,
    std::vector<std::tuple<std::string, std::string>>& connections,
    const svzero_json& config, const std::string& component,
    std::vector<std::string>& closed_loop_bcs) {
  ///< Flag to check if heart block is present (requires different handling)
  bool heartpulmonary_block_present = false;

  // Loop all closed loop blocks (dictionary format: name -> {type, values,
  // ...})
  for (auto& [block_name, closed_loop_config] : config[component].items()) {
    std::string closed_loop_type = closed_loop_config["type"];
    if (closed_loop_type == "ClosedLoopHeartAndPulmonary") {
      if (heartpulmonary_block_present == false) {
        heartpulmonary_block_present = true;
        std::string heartpulmonary_name = block_name;
        double cycle_period = closed_loop_config["cardiac_cycle_period"];
        if ((model.cardiac_cycle_period > 0.0) &&
            (cycle_period != model.cardiac_cycle_period)) {
          throw std::runtime_error(
              "Inconsistent cardiac cycle period defined in "
              "ClosedLoopHeartAndPulmonary.");
        } else {
          model.cardiac_cycle_period = cycle_period;
        }
        const auto& heart_params = closed_loop_config["values"];

        generate_block(model, heart_params, closed_loop_type,
                       heartpulmonary_name);

        // Junction at inlet to heart
        std::string heart_inlet_junction_name = "J_heart_inlet";
        connections.push_back({heart_inlet_junction_name, heartpulmonary_name});
        generate_block(model, {}, "NORMAL_JUNCTION", heart_inlet_junction_name);

        for (auto heart_inlet_elem : closed_loop_bcs) {
          connections.push_back({heart_inlet_elem, heart_inlet_junction_name});
        }

        // Junction at outlet from heart
        std::string heart_outlet_junction_name = "J_heart_outlet";
        connections.push_back(
            {heartpulmonary_name, heart_outlet_junction_name});
        generate_block(model, {}, "NORMAL_JUNCTION",
                       heart_outlet_junction_name);
        for (auto& outlet_block : closed_loop_config["outlet"]) {
          connections.push_back({heart_outlet_junction_name, outlet_block});
        }
      } else {
        throw std::runtime_error(
            "Error. Only one ClosedLoopHeartAndPulmonary can be included.");
      }
    }
    DEBUG_MSG("Created closed loop block " << block_name);
  }
}

void create_valves(
    Model& model,
    std::vector<std::tuple<std::string, std::string>>& connections,
    const svzero_json& config, const std::string& component) {
  // Loop all valves (dictionary format: name -> {type, values})
  for (auto& [valve_name, valve_config] : config[component].items()) {
    std::string valve_type = valve_config["type"];
    const auto& valve_values = valve_config["values"];
    generate_block(model, valve_values, valve_type, valve_name);
    connections.push_back({valve_values["upstream_block"], valve_name});
    connections.push_back({valve_name, valve_values["downstream_block"]});
    DEBUG_MSG("Created valve " << valve_name);
  }
}

void create_chambers(Model& model, const svzero_json& config,
                     const std::string& component) {
  // Loop all chambers (dictionary format: name -> {type, values})
  for (auto& [chamber_name, chamber_config] : config[component].items()) {
    std::string chamber_type = chamber_config["type"];
    generate_block(model, chamber_config["values"], chamber_type, chamber_name);
    DEBUG_MSG("Created chamber " << chamber_name);
  }
}

State load_initial_condition(const svzero_json& config, Model& model) {
  // Read initial condition
  auto initial_state = State::Zero(model.dofhandler.size());

  if (config.contains("initial_condition")) {
    const auto& initial_condition = config["initial_condition"];
    // Check for pressure_all or flow_all condition.
    // This will initialize all pressure:* and flow:* variables.
    double init_p, init_q;
    bool init_p_flag = initial_condition.contains("pressure_all");
    bool init_q_flag = initial_condition.contains("flow_all");
    if (init_p_flag) {
      init_p = initial_condition["pressure_all"];
    }
    if (init_q_flag) {
      init_q = initial_condition["flow_all"];
    }

    // Loop through variables and check for initial conditions.
    for (size_t i = 0; i < model.dofhandler.size(); i++) {
      std::string var_name = model.dofhandler.variables[i];
      double default_val = 0.0;
      // If initial condition is not specified for this variable,
      // check if pressure_all/flow_all are applicable
      if (!initial_condition.contains(var_name)) {
        if ((init_p_flag == true) && ((var_name.substr(0, 9) == "pressure:") ||
                                      (var_name.substr(0, 4) == "P_c:"))) {
          default_val = init_p;
          DEBUG_MSG("pressure_all initial condition for " << var_name);
        } else if ((init_q_flag == true) &&
                   (var_name.substr(0, 5) == "flow:")) {
          default_val = init_q;
          DEBUG_MSG("flow_all initial condition for " << var_name);
        } else {
          DEBUG_MSG("No initial condition found for "
                    << var_name << ". Using default value = 0.");
        }
      }
      initial_state.y[i] = initial_condition.value(var_name, default_val);
    }
  }
  if (config.contains("initial_condition_d")) {
    DEBUG_MSG("Reading initial condition derivative");
    const auto& initial_condition_d = config["initial_condition_d"];
    // Loop through variables and check for initial conditions.
    for (size_t i = 0; i < model.dofhandler.size(); i++) {
      std::string var_name = model.dofhandler.variables[i];
      if (!initial_condition_d.contains(var_name)) {
        DEBUG_MSG("No initial condition derivative found for " << var_name);
      }
      initial_state.ydot[i] = initial_condition_d.value(var_name, 0.0);
    }
  }
  return initial_state;
}
