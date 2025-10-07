// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "calibrate.h"

#include "LevenbergMarquardtOptimizer.h"
#include "Parameters.h"

nlohmann::json calibrate(const nlohmann::json& config) {
  auto output_config = nlohmann::json(config);

  // Read calibration parameters
  CalibrationParameters cali_params = load_calibration_params(config);

  // Setup model
  auto model = Model();
  load_calibration_model_legacy(config, model);

  // Build mapping between optimization parameters and model parameters
  // todo: have a list of optimization parameters for every block (default
  // empty) and add all calibrated parameters for each block
  std::vector<int> opt_to_model_param_ids;

  // Add all vessel parameters to optimization vector
  for (auto const& vessel_config : config["vessels"]) {
    std::string vessel_name = vessel_config["vessel_name"];
    auto block = model.get_block(vessel_name);
    for (size_t k = 0; k < 4; k++) {
      opt_to_model_param_ids.push_back(block->global_param_ids[k]);
    }
  }
  // Map junction parameters
  for (auto const& junction_config : config["junctions"]) {
    std::string junction_name = junction_config["junction_name"];
    auto block = model.get_block(junction_name);
    int num_outlets = block->outlet_nodes.size();
    if (num_outlets > 1) {
      for (size_t i = 0; i < num_outlets * 3; i++) {
        opt_to_model_param_ids.push_back(block->global_param_ids[i]);
      }
    }
  }

  const int param_counter = opt_to_model_param_ids.size();
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
    alpha[block->global_param_ids[3]] =
        vessel_config["zero_d_element_values"].value("stenosis_coefficient",
                                                     0.0);
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
      alpha[block->global_param_ids[i + 2 * num_outlets]] = 0.0;
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
        alpha[block->global_param_ids[i + 2 * num_outlets]] = stenosis_coeff[i];
      }
    }
  }

  // Run optimization
  DEBUG_MSG("Start optimization");
  auto lm_alg =
      LevenbergMarquardtOptimizer(&model, num_obs, param_counter, cali_params);

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
    double stenosis_coeff = alpha[block->global_param_ids[3]];
    double c_value = alpha[block->global_param_ids[1]];

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

    for (size_t i = 0; i < num_outlets; i++) {
      ste_values.push_back(alpha[block->global_param_ids[i + 2 * num_outlets]]);
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
