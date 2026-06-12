// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "ChamberSphere.h"

#include "Model.h"

void ChamberSphere::setup_dofs(DOFHandler& dofhandler) {
  Block::setup_dofs_(dofhandler, 7,
                     {"radius", "velo", "stress", "tau", "volume"});
}

void ChamberSphere::update_constant(SparseSystem& system,
                                    std::vector<double>& parameters) {
  const double rho = parameters[global_param_ids[ParamId::rho]];
  const double thick0 = parameters[global_param_ids[ParamId::thick0]];

  // balance of linear momentum
  system.E.coeffRef(global_eqn_ids[0], global_var_ids[5]) = rho * thick0;

  // spherical stress
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[6]) = -1;
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[7]) = 1;

  // volume change
  system.E.coeffRef(global_eqn_ids[2], global_var_ids[8]) = -1;

  // active stress
  system.E.coeffRef(global_eqn_ids[3], global_var_ids[7]) = 1;

  // acceleration
  system.E.coeffRef(global_eqn_ids[4], global_var_ids[4]) = 1;
  system.F.coeffRef(global_eqn_ids[4], global_var_ids[5]) = -1;

  // conservation of mass
  system.F.coeffRef(global_eqn_ids[5], global_var_ids[1]) = 1;
  system.F.coeffRef(global_eqn_ids[5], global_var_ids[3]) = -1;
  system.E.coeffRef(global_eqn_ids[5], global_var_ids[8]) = -1;

  // pressure equality
  system.F.coeffRef(global_eqn_ids[6], global_var_ids[0]) = 1;
  system.F.coeffRef(global_eqn_ids[6], global_var_ids[2]) = -1;
}

void ChamberSphere::update_time(SparseSystem& system,
                                std::vector<double>& parameters) {
  // active stress
  get_elastance_values(parameters);
  system.F.coeffRef(global_eqn_ids[3], global_var_ids[7]) = act;
}

void ChamberSphere::update_solution(
    SparseSystem& system, std::vector<double>& parameters,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& dy) {
  const double W1 = parameters[global_param_ids[ParamId::W1]];
  const double W2 = parameters[global_param_ids[ParamId::W2]];
  const double eta = parameters[global_param_ids[ParamId::eta]];
  const double thick0 = parameters[global_param_ids[ParamId::thick0]];
  const double sigma_max = parameters[global_param_ids[ParamId::sigma_max]];

  const double radius0 = parameters[global_param_ids[ParamId::radius0]];
  const double velo = y[global_var_ids[5]];
  const double Pout = y[global_var_ids[2]];
  const double stress = y[global_var_ids[6]];
  const double dradius_dt = dy[global_var_ids[4]];
  const double radius = y[global_var_ids[4]];

  // balance of momentum
  system.C.coeffRef(global_eqn_ids[0]) =
      (radius + radius0) * (-Pout * (radius + radius0) + stress * thick0) /
      pow(radius0, 2);
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[2]) =
      -pow(radius + radius0, 2) / pow(radius0, 2);
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[4]) =
      (-2 * Pout * (radius + radius0) + stress * thick0) / pow(radius0, 2);
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[6]) =
      thick0 * (radius + radius0) / pow(radius0, 2);

  // spherical stress
  system.C.coeffRef(global_eqn_ids[1]) =
      2 *
      (dradius_dt * eta * (2 * pow(radius0, 12) + pow(radius + radius0, 12)) +
       2 * pow(radius + radius0, 5) *
           (-pow(radius0, 6) + pow(radius + radius0, 6)) *
           (W1 * pow(radius0, 2) + W2 * pow(radius + radius0, 2))) /
      (pow(radius0, 2) * pow(radius + radius0, 11));
  system.dC_dy.coeffRef(global_eqn_ids[1], global_var_ids[4]) =
      24 * W1 * pow(radius0, 6) / pow(radius + radius0, 7) +
      8 * W2 * radius / pow(radius0, 2) +
      16 * W2 * pow(radius0, 4) / pow(radius + radius0, 5) + 8 * W2 / radius0 -
      44 * dradius_dt * eta * pow(radius0, 10) / pow(radius + radius0, 12) +
      2 * dradius_dt * eta / pow(radius0, 2);
  system.dC_dydot.coeffRef(global_eqn_ids[1], global_var_ids[4]) =
      2 * eta * (2 * pow(radius0, 12) + pow(radius + radius0, 12)) /
      (pow(radius0, 2) * pow(radius + radius0, 11));

  // volume change
  system.C.coeffRef(global_eqn_ids[2]) =
      4 * M_PI * velo * pow(radius + radius0, 2);
  system.dC_dy.coeffRef(global_eqn_ids[2], global_var_ids[4]) =
      8 * M_PI * velo * (radius + radius0);
  system.dC_dy.coeffRef(global_eqn_ids[2], global_var_ids[5]) =
      4 * M_PI * pow(radius + radius0, 2);

  // active stress
  system.C.coeffRef(global_eqn_ids[3]) = -act_plus * sigma_max;
}

void ChamberSphere::get_elastance_values(std::vector<double>& parameters) {
  const double alpha_max = parameters[global_param_ids[ParamId::alpha_max]];
  const double alpha_min = parameters[global_param_ids[ParamId::alpha_min]];
  const double tsys = parameters[global_param_ids[ParamId::tsys]];
  const double tdias = parameters[global_param_ids[ParamId::tdias]];
  const double steepness = parameters[global_param_ids[ParamId::steepness]];

  const double t = model->time;

  const auto T_cardiac = model->cardiac_cycle_period;
  const auto t_in_cycle = fmod(model->time, T_cardiac);

  const double S_plus = 0.5 * (1.0 + tanh((t_in_cycle - tsys) / steepness));
  const double S_minus = 0.5 * (1.0 - tanh((t_in_cycle - tdias) / steepness));

  // indicator function
  const double f = S_plus * S_minus;

  // activation rates
  const double act_t = alpha_max * f + alpha_min * (1 - f);

  act = std::abs(act_t);
  act_plus = std::max(act_t, 0.0);
}

void ChamberSphere::update_gradient(
    Eigen::SparseMatrix<double>& jacobian,
    Eigen::Matrix<double, Eigen::Dynamic, 1>& residual,
    Eigen::Matrix<double, Eigen::Dynamic, 1>& alpha, std::vector<double>& y,
    std::vector<double>& dy) {
  // The six time-independent parameters (rho, thick0, radius0, W1, W2, eta)
  // appear in the momentum (eq 0), spherical-stress (eq 1) and volume (eq 2)
  // equations, which are pure functions of the state. The active-stress
  // equation (eq 3) is time-dependent: its parameters (sigma_max and the
  // activation timing) are only calibrated when the optimizer provides the
  // observation time via ``model->time`` (i.e. when a ``t`` vector is supplied).
  const double rho = alpha[global_param_ids[ParamId::rho]];
  const double thick0 = alpha[global_param_ids[ParamId::thick0]];
  const double radius0 = alpha[global_param_ids[ParamId::radius0]];
  const double W1 = alpha[global_param_ids[ParamId::W1]];
  const double W2 = alpha[global_param_ids[ParamId::W2]];
  const double eta = alpha[global_param_ids[ParamId::eta]];
  const double sigma_max = alpha[global_param_ids[ParamId::sigma_max]];
  const double alpha_max = alpha[global_param_ids[ParamId::alpha_max]];
  const double alpha_min = alpha[global_param_ids[ParamId::alpha_min]];
  const double tsys = alpha[global_param_ids[ParamId::tsys]];
  const double tdias = alpha[global_param_ids[ParamId::tdias]];
  const double steepness = alpha[global_param_ids[ParamId::steepness]];

  const double Pin = y[global_var_ids[0]];
  const double Qin = y[global_var_ids[1]];
  const double Pout = y[global_var_ids[2]];
  const double Qout = y[global_var_ids[3]];
  const double radius = y[global_var_ids[4]];
  const double velo = y[global_var_ids[5]];
  const double stress = y[global_var_ids[6]];
  const double tau = y[global_var_ids[7]];
  const double dradius_dt = dy[global_var_ids[4]];
  const double dvelo_dt = dy[global_var_ids[5]];
  const double dtau_dt = dy[global_var_ids[7]];
  const double dvolume_dt = dy[global_var_ids[8]];
  const double t = model->cardiac_cycle_period > 0.0
                       ? fmod(model->time, model->cardiac_cycle_period)
                       : model->time;

  // Residuals (active-stress equation 3 is time-dependent and omitted).
  residual(global_eqn_ids[0]) =
      -Pout * pow(radius / radius0 + 1, 2) + dvelo_dt * rho * thick0 +
      stress * thick0 * (radius / radius0 + 1) / radius0;
  residual(global_eqn_ids[1]) =
      dradius_dt * eta * (1 + 2 / pow(radius / radius0 + 1, 12)) *
          (2 * radius / radius0 + 2) / radius0 -
      stress + tau +
      (4 - 4 / pow(radius / radius0 + 1, 6)) *
          (W1 + W2 * pow(radius / radius0 + 1, 2));
  residual(global_eqn_ids[2]) =
      -dvolume_dt + 4 * M_PI * pow(radius0, 2) * velo * pow(radius / radius0 + 1, 2);
  residual(global_eqn_ids[3]) = dtau_dt - sigma_max*fmax(0, alpha_max*(0.5 - 0.5*tanh((t - tdias)/steepness))*(0.5*tanh((t - tsys)/steepness) + 0.5) + alpha_min*(-(0.5 - 0.5*tanh((t - tdias)/steepness))*(0.5*tanh((t - tsys)/steepness) + 0.5) + 1)) + tau*fabs(alpha_max*(0.5*tanh((t - tdias)/steepness) - 0.5)*(0.5*tanh((t - tsys)/steepness) + 0.5) - alpha_min*((0.5*tanh((t - tdias)/steepness) - 0.5)*(0.5*tanh((t - tsys)/steepness) + 0.5) + 1));
  residual(global_eqn_ids[4]) = dradius_dt - velo;
  residual(global_eqn_ids[5]) = Qin - Qout - dvolume_dt;
  residual(global_eqn_ids[6]) = Pin - Pout;

  // Parameter Jacobian (only the six time-independent parameters).
  jacobian.coeffRef(global_eqn_ids[0], global_param_ids[ParamId::rho]) =
      dvelo_dt * thick0;
  jacobian.coeffRef(global_eqn_ids[0], global_param_ids[ParamId::thick0]) =
      dvelo_dt * rho + radius * stress / pow(radius0, 2) + stress / radius0;
  jacobian.coeffRef(global_eqn_ids[0], global_param_ids[ParamId::radius0]) =
      (2 * Pout * radius * (radius + radius0) - radius * stress * thick0 -
       stress * thick0 * (radius + radius0)) /
      pow(radius0, 3);
  jacobian.coeffRef(global_eqn_ids[1], global_param_ids[ParamId::radius0]) =
      -24 * W1 * radius * pow(radius0, 5) / pow(radius + radius0, 7) -
      8 * W2 * pow(radius, 2) / pow(radius0, 3) -
      16 * W2 * radius * pow(radius0, 3) / pow(radius + radius0, 5) -
      8 * W2 * radius / pow(radius0, 2) +
      40 * dradius_dt * eta * radius * pow(radius0, 9) / pow(radius + radius0, 12) -
      4 * dradius_dt * eta * radius / pow(radius0, 3) -
      4 * dradius_dt * eta * pow(radius0, 10) / pow(radius + radius0, 12) -
      2 * dradius_dt * eta / pow(radius0, 2);
  jacobian.coeffRef(global_eqn_ids[2], global_param_ids[ParamId::radius0]) =
      8 * M_PI * velo * (radius + radius0);
  jacobian.coeffRef(global_eqn_ids[1], global_param_ids[ParamId::W1]) =
      -4 * pow(radius0, 6) / pow(radius + radius0, 6) + 4;
  jacobian.coeffRef(global_eqn_ids[1], global_param_ids[ParamId::W2]) =
      4 * (-pow(radius0, 6) + pow(radius + radius0, 6)) /
      (pow(radius0, 2) * pow(radius + radius0, 4));
  jacobian.coeffRef(global_eqn_ids[1], global_param_ids[ParamId::eta]) =
      2 * dradius_dt * (2 * pow(radius0, 12) + pow(radius + radius0, 12)) /
      (pow(radius0, 2) * pow(radius + radius0, 11));

  // Active-stress (eq 3) parameter columns. These depend on the observation
  // time t and are only meaningful when a time vector is supplied; if the
  // activation parameters are not selected for calibration the columns are
  // unused.
  jacobian.coeffRef(global_eqn_ids[3], global_param_ids[ParamId::sigma_max]) = -fmax(0, alpha_max*(0.5 - 0.5*tanh((t - tdias)/steepness))*(0.5*tanh((t - tsys)/steepness) + 0.5) + alpha_min*(-(0.5 - 0.5*tanh((t - tdias)/steepness))*(0.5*tanh((t - tsys)/steepness) + 0.5) + 1));
  jacobian.coeffRef(global_eqn_ids[3], global_param_ids[ParamId::alpha_max]) = -sigma_max*(0.5 - 0.5*tanh((t - tdias)/steepness))*(0.5*tanh((t - tsys)/steepness) + 0.5)*(alpha_max*(0.5 - 0.5*tanh((t - tdias)/steepness))*(0.5*tanh((t - tsys)/steepness) + 0.5) + alpha_min*(-(0.5 - 0.5*tanh((t - tdias)/steepness))*(0.5*tanh((t - tsys)/steepness) + 0.5) + 1) > 0 ? 1.0 : 0.0) + tau*(0.5*tanh((t - tdias)/steepness) - 0.5)*(0.5*tanh((t - tsys)/steepness) + 0.5)*(((alpha_max*(0.5*tanh((t - tdias)/steepness) - 0.5)*(0.5*tanh((t - tsys)/steepness) + 0.5) - alpha_min*((0.5*tanh((t - tdias)/steepness) - 0.5)*(0.5*tanh((t - tsys)/steepness) + 0.5) + 1)) > 0) - ((alpha_max*(0.5*tanh((t - tdias)/steepness) - 0.5)*(0.5*tanh((t - tsys)/steepness) + 0.5) - alpha_min*((0.5*tanh((t - tdias)/steepness) - 0.5)*(0.5*tanh((t - tsys)/steepness) + 0.5) + 1)) < 0));
  jacobian.coeffRef(global_eqn_ids[3], global_param_ids[ParamId::alpha_min]) = -sigma_max*(-(0.5 - 0.5*tanh((t - tdias)/steepness))*(0.5*tanh((t - tsys)/steepness) + 0.5) + 1)*(alpha_max*(0.5 - 0.5*tanh((t - tdias)/steepness))*(0.5*tanh((t - tsys)/steepness) + 0.5) + alpha_min*(-(0.5 - 0.5*tanh((t - tdias)/steepness))*(0.5*tanh((t - tsys)/steepness) + 0.5) + 1) > 0 ? 1.0 : 0.0) + tau*(-(0.5*tanh((t - tdias)/steepness) - 0.5)*(0.5*tanh((t - tsys)/steepness) + 0.5) - 1)*(((alpha_max*(0.5*tanh((t - tdias)/steepness) - 0.5)*(0.5*tanh((t - tsys)/steepness) + 0.5) - alpha_min*((0.5*tanh((t - tdias)/steepness) - 0.5)*(0.5*tanh((t - tsys)/steepness) + 0.5) + 1)) > 0) - ((alpha_max*(0.5*tanh((t - tdias)/steepness) - 0.5)*(0.5*tanh((t - tsys)/steepness) + 0.5) - alpha_min*((0.5*tanh((t - tdias)/steepness) - 0.5)*(0.5*tanh((t - tsys)/steepness) + 0.5) + 1)) < 0));
  jacobian.coeffRef(global_eqn_ids[3], global_param_ids[ParamId::tsys]) = -sigma_max*(-0.5*alpha_max*(0.5 - 0.5*tanh((t - tdias)/steepness))*(1 - pow(tanh((t - tsys)/steepness), 2))/steepness + 0.5*alpha_min*(0.5 - 0.5*tanh((t - tdias)/steepness))*(1 - pow(tanh((t - tsys)/steepness), 2))/steepness)*(alpha_max*(0.5 - 0.5*tanh((t - tdias)/steepness))*(0.5*tanh((t - tsys)/steepness) + 0.5) + alpha_min*(-(0.5 - 0.5*tanh((t - tdias)/steepness))*(0.5*tanh((t - tsys)/steepness) + 0.5) + 1) > 0 ? 1.0 : 0.0) + tau*(-0.5*alpha_max*(1 - pow(tanh((t - tsys)/steepness), 2))*(0.5*tanh((t - tdias)/steepness) - 0.5)/steepness + 0.5*alpha_min*(1 - pow(tanh((t - tsys)/steepness), 2))*(0.5*tanh((t - tdias)/steepness) - 0.5)/steepness)*(((alpha_max*(0.5*tanh((t - tdias)/steepness) - 0.5)*(0.5*tanh((t - tsys)/steepness) + 0.5) - alpha_min*((0.5*tanh((t - tdias)/steepness) - 0.5)*(0.5*tanh((t - tsys)/steepness) + 0.5) + 1)) > 0) - ((alpha_max*(0.5*tanh((t - tdias)/steepness) - 0.5)*(0.5*tanh((t - tsys)/steepness) + 0.5) - alpha_min*((0.5*tanh((t - tdias)/steepness) - 0.5)*(0.5*tanh((t - tsys)/steepness) + 0.5) + 1)) < 0));
  jacobian.coeffRef(global_eqn_ids[3], global_param_ids[ParamId::tdias]) = -sigma_max*(0.5*alpha_max*(1 - pow(tanh((t - tdias)/steepness), 2))*(0.5*tanh((t - tsys)/steepness) + 0.5)/steepness - 0.5*alpha_min*(1 - pow(tanh((t - tdias)/steepness), 2))*(0.5*tanh((t - tsys)/steepness) + 0.5)/steepness)*(alpha_max*(0.5 - 0.5*tanh((t - tdias)/steepness))*(0.5*tanh((t - tsys)/steepness) + 0.5) + alpha_min*(-(0.5 - 0.5*tanh((t - tdias)/steepness))*(0.5*tanh((t - tsys)/steepness) + 0.5) + 1) > 0 ? 1.0 : 0.0) + tau*(-0.5*alpha_max*(1 - pow(tanh((t - tdias)/steepness), 2))*(0.5*tanh((t - tsys)/steepness) + 0.5)/steepness + 0.5*alpha_min*(1 - pow(tanh((t - tdias)/steepness), 2))*(0.5*tanh((t - tsys)/steepness) + 0.5)/steepness)*(((alpha_max*(0.5*tanh((t - tdias)/steepness) - 0.5)*(0.5*tanh((t - tsys)/steepness) + 0.5) - alpha_min*((0.5*tanh((t - tdias)/steepness) - 0.5)*(0.5*tanh((t - tsys)/steepness) + 0.5) + 1)) > 0) - ((alpha_max*(0.5*tanh((t - tdias)/steepness) - 0.5)*(0.5*tanh((t - tsys)/steepness) + 0.5) - alpha_min*((0.5*tanh((t - tdias)/steepness) - 0.5)*(0.5*tanh((t - tsys)/steepness) + 0.5) + 1)) < 0));
  jacobian.coeffRef(global_eqn_ids[3], global_param_ids[ParamId::steepness]) = -sigma_max*(-0.5*alpha_max*(0.5 - 0.5*tanh((t - tdias)/steepness))*(1 - pow(tanh((t - tsys)/steepness), 2))*(t - tsys)/pow(steepness, 2) + 0.5*alpha_max*(1 - pow(tanh((t - tdias)/steepness), 2))*(t - tdias)*(0.5*tanh((t - tsys)/steepness) + 0.5)/pow(steepness, 2) + alpha_min*(0.5*(0.5 - 0.5*tanh((t - tdias)/steepness))*(1 - pow(tanh((t - tsys)/steepness), 2))*(t - tsys)/pow(steepness, 2) - 0.5*(1 - pow(tanh((t - tdias)/steepness), 2))*(t - tdias)*(0.5*tanh((t - tsys)/steepness) + 0.5)/pow(steepness, 2)))*(alpha_max*(0.5 - 0.5*tanh((t - tdias)/steepness))*(0.5*tanh((t - tsys)/steepness) + 0.5) + alpha_min*(-(0.5 - 0.5*tanh((t - tdias)/steepness))*(0.5*tanh((t - tsys)/steepness) + 0.5) + 1) > 0 ? 1.0 : 0.0) + tau*(-0.5*alpha_max*(1 - pow(tanh((t - tdias)/steepness), 2))*(t - tdias)*(0.5*tanh((t - tsys)/steepness) + 0.5)/pow(steepness, 2) - 0.5*alpha_max*(1 - pow(tanh((t - tsys)/steepness), 2))*(t - tsys)*(0.5*tanh((t - tdias)/steepness) - 0.5)/pow(steepness, 2) - alpha_min*(-0.5*(1 - pow(tanh((t - tdias)/steepness), 2))*(t - tdias)*(0.5*tanh((t - tsys)/steepness) + 0.5)/pow(steepness, 2) - 0.5*(1 - pow(tanh((t - tsys)/steepness), 2))*(t - tsys)*(0.5*tanh((t - tdias)/steepness) - 0.5)/pow(steepness, 2)))*(((alpha_max*(0.5*tanh((t - tdias)/steepness) - 0.5)*(0.5*tanh((t - tsys)/steepness) + 0.5) - alpha_min*((0.5*tanh((t - tdias)/steepness) - 0.5)*(0.5*tanh((t - tsys)/steepness) + 0.5) + 1)) > 0) - ((alpha_max*(0.5*tanh((t - tdias)/steepness) - 0.5)*(0.5*tanh((t - tsys)/steepness) + 0.5) - alpha_min*((0.5*tanh((t - tdias)/steepness) - 0.5)*(0.5*tanh((t - tsys)/steepness) + 0.5) + 1)) < 0));
}