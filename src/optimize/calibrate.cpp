// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "calibrate.h"

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

  // Setup model using unified model loading
  auto model = Model();
  load_calibration_model_legacy(config, model);

  // Build mapping between optimization parameters and model parameters
  std::vector<int> opt_to_model_param_ids;

  // Map vessel parameters
  for (auto const& vessel_config : config["vessels"]) {
    std::string vessel_name = vessel_config["vessel_name"];
    auto block = model.get_block(vessel_name);

    // Add all vessel parameters to optimization vector
    for (size_t k = 0; k < num_params; k++) {
      if (k >= block->global_param_ids.size()) {
        throw std::runtime_error("Block " + vessel_name +
                                 " doesn't have enough parameters for calibration. Expected " +
                                 std::to_string(num_params) + ", got " +
                                 std::to_string(block->global_param_ids.size()));
      }
      opt_to_model_param_ids.push_back(block->global_param_ids[k]);
    }
  }

  // Map junction parameters
  for (auto const& junction_config : config["junctions"]) {
    std::string junction_name = junction_config["junction_name"];
    auto block = model.get_block(junction_name);
    int num_outlets = block->outlet_nodes.size();

    if (num_outlets > 1) {
      // BloodVesselJunction blocks have parameters for each outlet
      int expected_params = num_outlets * (num_params - 1);

      // Skip junctions that don't have parameters (e.g., NORMAL_JUNCTION without junction_values)
      // These will not be optimized and will keep their default behavior
      if (block->global_param_ids.size() == 0) {
        DEBUG_MSG("Skipping junction " << junction_name << " with no parameters");
        continue;
      }

      // For junctions with parameters, add them to optimization
      for (size_t i = 0; i < expected_params; i++) {
        if (i >= block->global_param_ids.size()) {
          throw std::runtime_error("Block " + junction_name +
                                   " doesn't have enough parameters for calibration. Expected " +
                                   std::to_string(expected_params) + ", got " +
                                   std::to_string(block->global_param_ids.size()));
        }
        opt_to_model_param_ids.push_back(block->global_param_ids[i]);
      }
    }
  }

  int param_counter = opt_to_model_param_ids.size();
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

  // Setup start parameter vector by reading from model parameters
  Eigen::Matrix<double, Eigen::Dynamic, 1> alpha =
      Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(param_counter);
  DEBUG_MSG("Reading initial alpha from model parameters");

  for (size_t opt_idx = 0; opt_idx < opt_to_model_param_ids.size(); opt_idx++) {
    int model_param_id = opt_to_model_param_ids[opt_idx];
    alpha[opt_idx] = model.get_parameter_value(model_param_id);
  }

  // Run optimization
  DEBUG_MSG("Start optimization");
  auto lm_alg =
      LevenbergMarquardtOptimizer(&model, num_obs, param_counter, lambda0,
                                  gradient_tol, increment_tol, max_iter);

  alpha = lm_alg.run(alpha, y_all, dy_all);

  // Build reverse mapping: model_param_id -> optimization index
  std::map<int, int> model_to_opt_param_ids;
  for (size_t opt_idx = 0; opt_idx < opt_to_model_param_ids.size(); opt_idx++) {
    model_to_opt_param_ids[opt_to_model_param_ids[opt_idx]] = opt_idx;
  }

  // Write optimized simulation config file
  for (auto& vessel_config : output_config["vessels"]) {
    std::string vessel_name = vessel_config["vessel_name"];
    auto block = model.get_block(vessel_name);

    double r_value = alpha[model_to_opt_param_ids[block->global_param_ids[0]]];
    double c_value = alpha[model_to_opt_param_ids[block->global_param_ids[1]]];
    double l_value = alpha[model_to_opt_param_ids[block->global_param_ids[2]]];
    double stenosis_coeff = 0.0;

    if (num_params > 3) {
      stenosis_coeff =
          alpha[model_to_opt_param_ids[block->global_param_ids[3]]];
    }
    if (zero_capacitance) {
      c_value = 0.0;
    }

    vessel_config["zero_d_element_values"] = {
        {"R_poiseuille", r_value},
        {"C", std::max(c_value, 0.0)},
        {"L", std::max(l_value, 0.0)},
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
      r_values.push_back(
          alpha[model_to_opt_param_ids[block->global_param_ids[i]]]);
    }
    std::vector<double> l_values;
    for (size_t i = 0; i < num_outlets; i++) {
      l_values.push_back(std::max(
          alpha[model_to_opt_param_ids[block->global_param_ids[i +
                                                               num_outlets]]],
          0.0));
    }

    std::vector<double> ste_values;

    if (num_params > 3) {
      for (size_t i = 0; i < num_outlets; i++) {
        ste_values.push_back(
            alpha[model_to_opt_param_ids
                      [block->global_param_ids[i + 2 * num_outlets]]]);
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
