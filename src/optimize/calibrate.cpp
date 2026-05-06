// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "calibrate.h"

#include <cmath>

#include "Integrator.h"
#include "LevenbergMarquardtOptimizer.h"
#include "SimulationParameters.h"

nlohmann::json calibrate(const nlohmann::json& config) {
  auto output_config = nlohmann::json(config);

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

  // Create vessels
  DEBUG_MSG("Load vessels");
  std::map<std::int64_t, std::string> vessel_id_map;
  int param_counter = 0;
  for (auto const& vessel_config : config["vessels"]) {
    std::string vessel_name = vessel_config["vessel_name"];

    // Create parameter IDs
    std::vector<int> param_ids;
    for (size_t k = 0; k < num_params; k++)
      param_ids.push_back(param_counter++);
    std::string block_type =
        vessel_config["zero_d_element_type"].get<std::string>();
    model.add_block(block_type, param_ids, vessel_name);
    vessel_id_map.insert({vessel_config["vessel_id"], vessel_name});
    DEBUG_MSG("Created vessel " << vessel_name);

    // Read connected boundary conditions
    if (vessel_config.contains("boundary_conditions")) {
      auto const& vessel_bc_config = vessel_config["boundary_conditions"];
      if (vessel_bc_config.contains("inlet")) {
        inlet_connections.push_back({vessel_bc_config["inlet"], vessel_name});
      }
      if (vessel_bc_config.contains("outlet")) {
        outlet_connections.push_back({vessel_name, vessel_bc_config["outlet"]});
      }
    }
  }

  // Create junctions
  for (auto const& junction_config : config["junctions"]) {
    std::string junction_name = junction_config["junction_name"];
    auto const& outlet_vessels = junction_config["outlet_vessels"];
    int num_outlets = outlet_vessels.size();

    if (num_outlets == 1) {
      model.add_block("NORMAL_JUNCTION", {}, junction_name);

    } else {
      std::vector<int> param_ids;
      for (size_t i = 0; i < (num_outlets * (num_params - 1)); i++)
        param_ids.push_back(param_counter++);
      model.add_block("BloodVesselJunction", param_ids, junction_name);
    }

    // Check for connections to inlet and outlet vessels and append to
    // connections list
    for (auto vessel_id : junction_config["inlet_vessels"]) {
      connections.push_back({vessel_id_map[vessel_id], junction_name});
    }

    for (auto vessel_id : outlet_vessels) {
      connections.push_back({junction_name, vessel_id_map[vessel_id]});
    }
    DEBUG_MSG("Created junction " << junction_name);
  }

  // Create Connections
  DEBUG_MSG("Created connection");
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

  // Read observations of the state vector y from the forward solver
  DEBUG_MSG("Reading observations");
  int num_obs = 0;
  int num_vars = model.dofhandler.get_num_variables();
  std::vector<std::vector<double>> y_all;
  auto y_values = config["y"];
  for (size_t i = 0; i < num_vars; i++) {
    std::string var_name = model.dofhandler.variables[i];
    DEBUG_MSG("Reading observations for variable " << var_name);
    if (!y_values.contains(var_name)) {
      std::cout << "ERROR: Missing y observation for '" << var_name << "'"
                << std::endl;
      exit(1);
    }
    auto y_array = y_values[var_name].get<std::vector<double>>();
    num_obs = y_array.size();
    if (i == 0) {
      y_all.resize(num_obs);
    }
    for (size_t j = 0; j < num_obs; j++) {
      y_all[j].push_back(y_array[j]);
    }
  }
  DEBUG_MSG("Number of observations: " << num_obs);

  // Compute the time derivative dy from y consistent with the forward
  // solver's generalized-alpha integration. The forward solver produces
  // (y_n, ydot_n) pairs satisfying
  //     y_{n+1} = y_n + dt * ((1 - gamma) * ydot_n + gamma * ydot_{n+1}).
  // Solving for ydot_{n+1} gives the recurrence
  //     ydot_{n+1} = (y_{n+1} - y_n) / (dt*gamma) + (1 - 1/gamma) * ydot_n.
  // We close the recurrence with the periodic assumption ydot_0 = ydot_{N-1},
  // which holds for cycle-converged simulation output.
  DEBUG_MSG("Computing dy from y consistent with generalized-alpha");
  double rho_infty = 0.5;
  double cardiac_period = -1.0;
  if (config.contains("simulation_parameters")) {
    auto const& sp = config["simulation_parameters"];
    rho_infty = sp.value("rho_infty", 0.5);
    cardiac_period = sp.value("cardiac_period", -1.0);
  }
  if (cardiac_period <= 0.0 && config.contains("boundary_conditions")) {
    for (auto const& bc : config["boundary_conditions"]) {
      if (bc.contains("bc_values") && bc["bc_values"].contains("t")) {
        auto const& t_array = bc["bc_values"]["t"];
        if (t_array.is_array() && t_array.size() >= 2) {
          cardiac_period = t_array.back().get<double>() -
                           t_array.front().get<double>();
          break;
        }
      }
    }
  }
  if (cardiac_period <= 0.0) {
    cardiac_period = 1.0;
  }

  GenAlphaCoefficients ga(rho_infty);
  double dt = cardiac_period / static_cast<double>(num_obs - 1);
  double dt_gamma = dt * ga.gamma;

  std::vector<std::vector<double>> dy_all(num_obs,
                                          std::vector<double>(num_vars, 0.0));
  for (size_t v = 0; v < num_vars; v++) {
    // Forcing terms c_n = (y_n - y_{n-1}) / (dt*gamma) for n = 1, ..., N-1.
    std::vector<double> c(num_obs, 0.0);
    for (int n = 1; n < num_obs; n++) {
      c[n] = (y_all[n][v] - y_all[n - 1][v]) / dt_gamma;
    }
    // Solve ydot_0 = sum_{k=1}^{N-1} r^{N-1-k} c_k / (1 - r^{N-1})
    // by accumulating from k = N-1 down to k = 1.
    double sum = 0.0;
    double rk = 1.0;
    for (int k = num_obs - 1; k >= 1; k--) {
      sum += rk * c[k];
      rk *= ga.ydot_init_coeff;
    }
    double r_N1 = std::pow(ga.ydot_init_coeff, num_obs - 1);
    double ydot0 = sum / (1.0 - r_N1);

    // Forward propagate to get ydot at each observation.
    dy_all[0][v] = ydot0;
    for (int n = 1; n < num_obs; n++) {
      dy_all[n][v] = c[n] + ga.ydot_init_coeff * dy_all[n - 1][v];
    }
  }

  // Setup start parameter vector
  Eigen::Matrix<double, Eigen::Dynamic, 1> alpha =
      Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(param_counter);
  DEBUG_MSG("Reading initial alpha");
  for (auto& vessel_config : output_config["vessels"]) {
    std::string vessel_name = vessel_config["vessel_name"];
    DEBUG_MSG("Reading initial alpha for " << vessel_name);
    auto block = model.get_block(vessel_name);
    alpha[block->global_param_ids[0]] =
        vessel_config["zero_d_element_values"].value("R_poiseuille", 0.0);
    alpha[block->global_param_ids[1]] =
        vessel_config["zero_d_element_values"].value("C", 0.0);
    alpha[block->global_param_ids[2]] =
        vessel_config["zero_d_element_values"].value("L", 0.0);
    if (num_params > 3) {
      alpha[block->global_param_ids[3]] =
          vessel_config["zero_d_element_values"].value("stenosis_coefficient",
                                                       0.0);
    }
  }
  for (auto& junction_config : output_config["junctions"]) {
    std::string junction_name = junction_config["junction_name"];
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
    if (junction_config["junction_type"] == "BloodVesselJunction") {
      auto resistance = junction_config["junction_values"]["R_poiseuille"]
                            .get<std::vector<double>>();
      auto inductance =
          junction_config["junction_values"]["L"].get<std::vector<double>>();
      auto stenosis_coeff =
          junction_config["junction_values"]["stenosis_coefficient"]
              .get<std::vector<double>>();
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

  // Run optimization
  DEBUG_MSG("Start optimization");
  auto lm_alg =
      LevenbergMarquardtOptimizer(&model, num_obs, param_counter, lambda0,
                                  gradient_tol, increment_tol, max_iter);

  alpha = lm_alg.run(alpha, y_all, dy_all);

  // Write optimized simulation config file
  for (auto& vessel_config : output_config["vessels"]) {
    std::string vessel_name = vessel_config["vessel_name"];
    auto block = model.get_block(vessel_name);
    double stenosis_coeff = 0.0;
    if (num_params > 3) {
      stenosis_coeff = alpha[block->global_param_ids[3]];
    }
    double c_value = 0.0;
    if (!zero_capacitance) {
      c_value = alpha[block->global_param_ids[1]];
    }
    vessel_config["zero_d_element_values"] = {
        {"R_poiseuille", alpha[block->global_param_ids[0]]},
        {"C", std::max(c_value, 0.0)},
        {"L", std::max(alpha[block->global_param_ids[2]], 0.0)},
        {"stenosis_coefficient", stenosis_coeff}};
  }
  for (auto& junction_config : output_config["junctions"]) {
    std::string junction_name = junction_config["junction_name"];
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

    junction_config["junction_type"] = "BloodVesselJunction";
    junction_config["junction_values"] = {{"R_poiseuille", r_values},
                                          {"L", l_values},
                                          {"stenosis_coefficient", ste_values}};
  }

  output_config.erase("y");
  output_config.erase("dy");
  output_config.erase("calibration_parameters");

  return output_config;
}
