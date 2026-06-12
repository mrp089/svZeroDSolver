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

  auto warp_signed = [T_cardiac](double dt) {
    return fmod(dt + 1.5 * T_cardiac, T_cardiac) - 0.5 * T_cardiac;
  };

  const double phase_tsys = warp_signed(t_in_cycle - tsys);
  const double phase_tdias = warp_signed(t_in_cycle - tdias);

  const double S_plus = 0.5 * (1.0 + tanh(phase_tsys / steepness));
  const double S_minus = 0.5 * (1.0 - tanh(phase_tdias / steepness));

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
  // Calibrates the chamber parameters from a full-state observation set. The
  // six time-independent parameters appear in the momentum, spherical-stress
  // and volume equations; the active-stress equation (eq 3) and its parameters
  // (sigma_max and the activation timing) depend on the observation time, which
  // the optimizer supplies via model->time when a time vector is given. Shared
  // subexpressions are hoisted into temporaries (xN). Generated by
  // scripts/jacobian.py from scripts/ChamberSphere.yaml.
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
  const double t = model->cardiac_cycle_period > 0.0 ? fmod(model->time, model->cardiac_cycle_period) : model->time;
  const double T_cardiac = model->cardiac_cycle_period;
  auto warp_signed = [T_cardiac](double dt) {
    return T_cardiac > 0.0 ? fmod(dt + 1.5 * T_cardiac, T_cardiac) - 0.5 * T_cardiac : dt;
  };

  const double x0 = dvelo_dt*thick0;
  const double x1 = 1.0/radius0;
  const double x2 = radius*x1;
  const double x3 = x2 + 1;
  const double x4 = pow(x3, 2);
  const double x5 = stress*x1;
  const double x6 = dradius_dt*eta;
  const double x7 = pow(radius0, 2);
  const double x8 = 1.0/steepness;
  const double x9 = warp_signed(t - tdias);
  const double x10 = tanh(x8*x9);
  const double x11 = 0.5*x10 - 0.5;
  const double x12 = warp_signed(t - tsys);
  const double x13 = tanh(x12*x8);
  const double x14 = 0.5*x13 + 0.5;
  const double x15 = x11*x14;
  const double x16 = x15 + 1;
  const double x17 = alpha_max*x15 - alpha_min*x16;
  const double x18 = -x11;
  const double x19 = x14*x18;
  const double x20 = 1 - x19;
  const double x21 = alpha_max*x19 + alpha_min*x20;
  const double x22 = fmax(0, x21);
  const double x23 = 1.0/x7;
  const double x24 = radius*stress;
  const double x25 = pow(radius0, 3);
  const double x26 = 1.0/x25;
  const double x27 = radius + radius0;
  const double x28 = 8*W2;
  const double x29 = 2*x23;
  const double x30 = 4*x6;
  const double x31 = pow(x27, 12);
  const double x32 = 1.0/x31;
  const double x33 = pow(radius0, 6);
  const double x34 = pow(x27, 6);
  const double x35 = sigma_max*(x21 > 0 ? 1.0 : 0.0);
  const double x36 = (((x17) > 0) - ((x17) < 0));
  const double x37 = 1 - pow(x13, 2);
  const double x38 = 0.5*x8;
  const double x39 = x37*x38;
  const double x40 = 1 - pow(x10, 2);
  const double x41 = x14*x40;
  const double x42 = alpha_max*x38*x41 - 0.5*alpha_min*x14*x40*x8;
  const double x43 = 0.5/pow(steepness, 2);
  const double x44 = x12*x37*x43;
  const double x45 = x18*x44;
  const double x46 = x41*x43*x9;
  const double x47 = alpha_max*x46;
  const double x48 = x11*x44;

  residual(global_eqn_ids[0]) = -Pout*x4 + rho*x0 + thick0*x3*x5;
  residual(global_eqn_ids[1]) = -stress + tau + x1*x6*(1 + 2/pow(x3, 12))*(2*x2 + 2) + (4 - 4/pow(x3, 6))*(W1 + W2*x4);
  residual(global_eqn_ids[2]) = -dvolume_dt + 4*M_PI*velo*x4*x7;
  residual(global_eqn_ids[3]) = dtau_dt - sigma_max*x22 + tau*fabs(x17);
  residual(global_eqn_ids[4]) = dradius_dt - velo;
  residual(global_eqn_ids[5]) = Qin - Qout - dvolume_dt;
  residual(global_eqn_ids[6]) = Pin - Pout;
  jacobian.coeffRef(global_eqn_ids[0], global_param_ids[ParamId::rho]) = x0;
  jacobian.coeffRef(global_eqn_ids[0], global_param_ids[ParamId::thick0]) = dvelo_dt*rho + x23*x24 + x5;
  jacobian.coeffRef(global_eqn_ids[0], global_param_ids[ParamId::radius0]) = x26*(2*Pout*radius*x27 - stress*thick0*x27 - thick0*x24);
  jacobian.coeffRef(global_eqn_ids[1], global_param_ids[ParamId::radius0]) = -24*W1*radius*pow(radius0, 5)/pow(x27, 7) - 16*W2*radius*x25/pow(x27, 5) + 40*dradius_dt*eta*radius*pow(radius0, 9)*x32 - pow(radius, 2)*x26*x28 - radius*x23*x28 - radius*x26*x30 - pow(radius0, 10)*x30*x32 - x29*x6;
  jacobian.coeffRef(global_eqn_ids[2], global_param_ids[ParamId::radius0]) = 8*M_PI*velo*x27;
  jacobian.coeffRef(global_eqn_ids[1], global_param_ids[ParamId::W1]) = -4*x33/x34 + 4;
  jacobian.coeffRef(global_eqn_ids[1], global_param_ids[ParamId::W2]) = 4*x23*(-x33 + x34)/pow(x27, 4);
  jacobian.coeffRef(global_eqn_ids[1], global_param_ids[ParamId::eta]) = dradius_dt*x29*(2*pow(radius0, 12) + x31)/pow(x27, 11);
  jacobian.coeffRef(global_eqn_ids[3], global_param_ids[ParamId::sigma_max]) = -x22;
  jacobian.coeffRef(global_eqn_ids[3], global_param_ids[ParamId::alpha_max]) = tau*x11*x14*x36 - x19*x35;
  jacobian.coeffRef(global_eqn_ids[3], global_param_ids[ParamId::alpha_min]) = -tau*x16*x36 - x20*x35;
  jacobian.coeffRef(global_eqn_ids[3], global_param_ids[ParamId::tsys]) = tau*x36*(-alpha_max*x11*x39 + 0.5*alpha_min*x11*x37*x8) - x35*(-alpha_max*x18*x39 + 0.5*alpha_min*x18*x37*x8);
  jacobian.coeffRef(global_eqn_ids[3], global_param_ids[ParamId::tdias]) = -tau*x36*x42 - x35*x42;
  jacobian.coeffRef(global_eqn_ids[3], global_param_ids[ParamId::steepness]) = tau*x36*(-alpha_max*x48 - alpha_min*(-x46 - x48) - x47) - x35*(-alpha_max*x45 + alpha_min*(x45 - x46) + x47);
}
