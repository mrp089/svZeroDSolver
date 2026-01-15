// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "calibrate.h"

#include "LevenbergMarquardtOptimizer.h"
#include "SimulationParameters.h"

svzero_json calibrate(const svzero_json& config) {
  auto output_config = svzero_json(config);

  // Read calibration parameters
  DEBUG_MSG("Parse calibration parameters");
  auto const& calibration_parameters = config["calibration_parameters"];
  double gradient_tol =
      calibration_parameters.value("tolerance_gradient", 1e-5);
  double increment_tol =
      calibration_parameters.value("tolerance_increment", 1e-10);
  int max_iter = calibration_parameters.value("maximum_iterations", 100);
  bool calibrate_stenosis =
      calibration_parameters.value("calibrate_stenosis_coefficient", true);
  bool zero_capacitance =
      calibration_parameters.value("set_capacitance_to_zero", false);
  double lambda0 = calibration_parameters.value("initial_damping_factor", 1.0);

  int num_params = 3;
  if (calibrate_stenosis) {
    num_params = 4;
  }

  // Setup model
  auto model = Model();
  std::vector<std::tuple<std::string, std::string>> connections;
  std::vector<std::tuple<std::string, std::string>> inlet_connections;
  std::vector<std::tuple<std::string, std::string>> outlet_connections;

  // Create vessels (dictionary format: name -> {type, values})
  DEBUG_MSG("Load vessels");
  int param_counter = 0;
  for (auto const& [vessel_name, vessel_config] : config["vessels"].items()) {
    // Create parameter IDs
    std::vector<int> param_ids;
    for (size_t k = 0; k < num_params; k++)
      param_ids.push_back(param_counter++);
    std::string block_type = vessel_config["type"].get<std::string>();
    model.add_block(block_type, param_ids, vessel_name);
    DEBUG_MSG("Created vessel " << vessel_name);
  }

  // Create junctions (dictionary format: name -> {type, values, inlet, outlet})
  if (config.contains("junctions")) {
    for (auto const& [junction_name, junction_config] :
         config["junctions"].items()) {
      auto const& outlets = junction_config["outlet"];
      int num_outlets = outlets.size();

      if (num_outlets == 1) {
        model.add_block("NORMAL_JUNCTION", {}, junction_name);
      } else {
        std::vector<int> param_ids;
        for (size_t i = 0; i < (num_outlets * (num_params - 1)); i++)
          param_ids.push_back(param_counter++);
        model.add_block("BloodVesselJunction", param_ids, junction_name);
      }

      // Check for connections to inlet and outlet blocks
      for (auto const& block_name : junction_config["inlet"]) {
        connections.push_back({block_name, junction_name});
      }

      for (auto const& block_name : outlets) {
        connections.push_back({junction_name, block_name});
      }
      DEBUG_MSG("Created junction " << junction_name);
    }
  }

  // Process top-level connections array
  // Connections can be: [A, B] - simple connection between named blocks
  // For BC-to-vessel connections, the BC name is used directly
  if (config.contains("connections")) {
    for (const auto& conn : config["connections"]) {
      bool first_is_array = conn[0].is_array();
      bool second_is_array = conn[1].is_array();

      if (!first_is_array && !second_is_array) {
        // Simple connection: [A, B]
        std::string upstream = conn[0];
        std::string downstream = conn[1];
        // Check if upstream/downstream are vessels (have blocks) or BCs (no
        // blocks)
        bool upstream_is_vessel = model.has_block(upstream);
        bool downstream_is_vessel = model.has_block(downstream);

        if (upstream_is_vessel && downstream_is_vessel) {
          connections.push_back({upstream, downstream});
        } else if (!upstream_is_vessel && downstream_is_vessel) {
          // BC -> vessel (inlet)
          inlet_connections.push_back({upstream, downstream});
        } else if (upstream_is_vessel && !downstream_is_vessel) {
          // vessel -> BC (outlet)
          outlet_connections.push_back({upstream, downstream});
        }
      } else if (!first_is_array && second_is_array) {
        // Bifurcation: [A, [B, C, ...]]
        std::string upstream = conn[0];
        for (const auto& downstream : conn[1]) {
          connections.push_back({upstream, downstream});
        }
      } else if (first_is_array && !second_is_array) {
        // Confluence: [[A, B, ...], C]
        std::string downstream = conn[1];
        for (const auto& upstream : conn[0]) {
          connections.push_back({upstream, downstream});
        }
      }
    }
  }

  // Create Connections
  DEBUG_MSG("Create connections");
  for (auto& connection : connections) {
    auto ele1 = model.get_block(std::get<0>(connection));
    auto ele2 = model.get_block(std::get<1>(connection));
    model.add_node({ele1}, {ele2}, ele1->get_name() + ":" + ele2->get_name());
  }
  for (auto& connection : inlet_connections) {
    auto ele = model.get_block(std::get<1>(connection));
    model.add_node({}, {ele}, std::get<0>(connection) + ":" + ele->get_name());
  }
  for (auto& connection : outlet_connections) {
    auto ele = model.get_block(std::get<0>(connection));
    model.add_node({ele}, {}, ele->get_name() + ":" + std::get<1>(connection));
  }

  // Finalize model
  model.finalize();

  DEBUG_MSG("Number of parameters " << param_counter);

  // Read observations
  DEBUG_MSG("Reading observations");
  int num_obs = 0;
  std::vector<std::vector<double>> y_all;
  std::vector<std::vector<double>> dy_all;
  auto y_values = config["y"];
  auto dy_values = config["dy"];
  for (size_t i = 0; i < model.dofhandler.get_num_variables(); i++) {
    std::string var_name = model.dofhandler.variables[i];
    DEBUG_MSG("Reading observations for variable " << var_name);
    if (!y_values.contains(var_name)) {
      std::cout << "ERROR: Missing y observation for '" << var_name << "'"
                << std::endl;
      exit(1);
    }
    if (!dy_values.contains(var_name)) {
      std::cout << "ERROR: Missing dy observation for '" << var_name << "'"
                << std::endl;
      exit(1);
    }
    auto y_array = y_values[var_name].get<std::vector<double>>();
    auto dy_array = dy_values[var_name].get<std::vector<double>>();
    num_obs = y_array.size();
    if (i == 0) {
      y_all.resize(num_obs);
      dy_all.resize(num_obs);
    }
    for (size_t j = 0; j < num_obs; j++) {
      y_all[j].push_back(y_array[j]);
      dy_all[j].push_back(dy_array[j]);
    }
  }
  DEBUG_MSG("Number of observations: " << num_obs);

  // Setup start parameter vector
  Eigen::Matrix<double, Eigen::Dynamic, 1> alpha =
      Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(param_counter);
  DEBUG_MSG("Reading initial alpha");
  for (auto& [vessel_name, vessel_config] : output_config["vessels"].items()) {
    DEBUG_MSG("Reading initial alpha for " << vessel_name);
    auto block = model.get_block(vessel_name);
    const auto& values = vessel_config["values"];
    alpha[block->global_param_ids[0]] = values.value("resistance", 0.0);
    alpha[block->global_param_ids[1]] = values.value("capacitance", 0.0);
    alpha[block->global_param_ids[2]] = values.value("inductance", 0.0);
    if (num_params > 3) {
      alpha[block->global_param_ids[3]] = values.value("stenosis", 0.0);
    }
  }
  if (output_config.contains("junctions")) {
    for (auto& [junction_name, junction_config] :
         output_config["junctions"].items()) {
      DEBUG_MSG("Reading initial alpha for " << junction_name);
      auto block = model.get_block(junction_name);
      int num_outlets = block->outlet_nodes.size();

      if (num_outlets < 2) {
        continue;
      }

      for (size_t i = 0; i < num_outlets; i++) {
        alpha[block->global_param_ids[i]] = 0.0;
        alpha[block->global_param_ids[i + num_outlets]] = 0.0;
        if (num_params > 3) {
          alpha[block->global_param_ids[i + 2 * num_outlets]] = 0.0;
        }
      }
      if (junction_config["type"] == "BloodVesselJunction") {
        const auto& values = junction_config["values"];
        auto resistance = values["resistance"].get<std::vector<double>>();
        auto inductance = values["inductance"].get<std::vector<double>>();
        auto stenosis_coeff = values["stenosis"].get<std::vector<double>>();
        for (size_t i = 0; i < num_outlets; i++) {
          alpha[block->global_param_ids[i]] = resistance[i];
          alpha[block->global_param_ids[i + num_outlets]] = inductance[i];
          if (num_params > 3) {
            alpha[block->global_param_ids[i + 2 * num_outlets]] =
                stenosis_coeff[i];
          }
        }
      }
    }
  }

  // Run optimization
  DEBUG_MSG("Start optimization");
  auto lm_alg =
      LevenbergMarquardtOptimizer(&model, num_obs, param_counter, lambda0,
                                  gradient_tol, increment_tol, max_iter);

  alpha = lm_alg.run(alpha, y_all, dy_all);

  // Write optimized simulation config file
  for (auto& [vessel_name, vessel_config] : output_config["vessels"].items()) {
    auto block = model.get_block(vessel_name);
    double stenosis_coeff = 0.0;
    if (num_params > 3) {
      stenosis_coeff = alpha[block->global_param_ids[3]];
    }
    double c_value = 0.0;
    if (!zero_capacitance) {
      c_value = alpha[block->global_param_ids[1]];
    }
    vessel_config["values"] = {
        {"resistance", alpha[block->global_param_ids[0]]},
        {"capacitance", std::max(c_value, 0.0)},
        {"inductance", std::max(alpha[block->global_param_ids[2]], 0.0)},
        {"stenosis", stenosis_coeff}};
  }
  if (output_config.contains("junctions")) {
    for (auto& [junction_name, junction_config] :
         output_config["junctions"].items()) {
      auto block = model.get_block(junction_name);
      int num_outlets = block->outlet_nodes.size();

      if (num_outlets < 2) {
        continue;
      }

      std::vector<double> r_values;
      for (size_t i = 0; i < num_outlets; i++) {
        r_values.push_back(alpha[block->global_param_ids[i]]);
      }
      std::vector<double> l_values;
      for (size_t i = 0; i < num_outlets; i++) {
        l_values.push_back(
            std::max(alpha[block->global_param_ids[i + num_outlets]], 0.0));
      }

      std::vector<double> ste_values;

      if (num_params > 3) {
        for (size_t i = 0; i < num_outlets; i++) {
          ste_values.push_back(
              alpha[block->global_param_ids[i + 2 * num_outlets]]);
        }
      } else {
        for (size_t i = 0; i < num_outlets; i++) {
          ste_values.push_back(0.0);
        }
      }

      junction_config["type"] = "BloodVesselJunction";
      junction_config["values"] = {{"resistance", r_values},
                                   {"inductance", l_values},
                                   {"stenosis", ste_values}};
    }
  }

  output_config.erase("y");
  output_config.erase("dy");
  output_config.erase("calibration_parameters");

  return output_config;
}
