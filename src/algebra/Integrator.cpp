// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "Integrator.h"

#include <cmath>
#include <stdexcept>

Integrator::Integrator(Model* model, double time_step_size, double rho,
                       double atol, int max_iter) {
  this->model = model;
  alpha_m = 0.5 * (3.0 - rho) / (1.0 + rho);
  alpha_f = 1.0 / (1.0 + rho);
  gamma = 0.5 + alpha_m - alpha_f;
  ydot_init_coeff = 1.0 - 1.0 / gamma;

  y_coeff = gamma * time_step_size;
  y_coeff_jacobian = alpha_f * y_coeff;

  size = model->dofhandler.size();
  system = SparseSystem(size);
  this->time_step_size = time_step_size;
  this->atol = atol;
  this->max_iter = max_iter;

  y_af = Eigen::Matrix<double, Eigen::Dynamic, 1>(size);
  ydot_am = Eigen::Matrix<double, Eigen::Dynamic, 1>(size);

  // Make some memory reservations
  system.reserve(model);
}

// Must declare default constructord and dedtructor
// because of Eigen.
Integrator::Integrator() {}
Integrator::~Integrator() {}

void Integrator::clean() {
  // Cannot be in destructor because dynamically allocated pointers will be lost
  // when objects are assigned from temporary objects.
  system.clean();
}

void Integrator::update_params(double time_step_size) {
  this->time_step_size = time_step_size;
  y_coeff = gamma * time_step_size;
  y_coeff_jacobian = alpha_f * y_coeff;
  model->update_constant(system);
  model->update_time(system, 0.0);
}

State Integrator::step(const State& old_state, double time) {
  // Predictor: Constant y, consistent ydot
  State new_state = State::Zero(size);
  new_state.ydot += old_state.ydot * ydot_init_coeff;
  new_state.y += old_state.y;

  // Determine new time (evaluate terms at generalized mid-point)
  double new_time = time + alpha_f * time_step_size;

  // Evaluate time-dependent element contributions in system
  model->update_time(system, new_time);

  // Count total number of step calls
  n_iter++;

  // Non-linear Newton-Raphson iterations
  for (size_t i = 0; i < max_iter; i++) {
    // Initiator: Evaluate the iterates at the intermediate time levels
    ydot_am.setZero();
    y_af.setZero();
    ydot_am += old_state.ydot + (new_state.ydot - old_state.ydot) * alpha_m;
    y_af += old_state.y + (new_state.y - old_state.y) * alpha_f;

    // Update solution-dependent element contribitions
    model->update_solution(system, y_af, ydot_am);

    // Evaluate residual
    system.update_residual(y_af, ydot_am);

    // Check termination criterium
    if (system.residual.cwiseAbs().maxCoeff() < atol) {
      break;
    }

    // Abort if maximum number of non-linear iterations is reached
    else if (i == max_iter - 1) {
      throw std::runtime_error(
          "Maximum number of non-linear iterations reached at time " +
          std::to_string(time));
    }

    // Evaluate Jacobian
    system.update_jacobian(alpha_m, y_coeff_jacobian);

    // Solve system for increment in ydot
    system.solve();

    // Perform post-solve actions on blocks
    model->post_solve(new_state.y);

    // Update the solution
    new_state.ydot += system.dydot;
    new_state.y += system.dydot * y_coeff;

    // Count total number of nonlinear iterations
    n_nonlin_iter++;
  }

  return new_state;
}

double Integrator::avg_nonlin_iter() {
  return (double)n_nonlin_iter / (double)n_iter;
}

std::vector<std::vector<double>> Integrator::compute_ydot(
    const std::vector<std::vector<double>>& y_series,
    const std::vector<double>& times, double rho) {
  const size_t num_steps = y_series.size();
  if (num_steps != times.size()) {
    throw std::runtime_error(
        "[Integrator::compute_ydot] Number of time points (" +
        std::to_string(times.size()) +
        ") does not match number of state observations (" +
        std::to_string(num_steps) + ").");
  }

  std::vector<std::vector<double>> ydot;
  if (num_steps == 0) return ydot;
  const size_t num_vars = y_series[0].size();
  ydot.assign(num_steps, std::vector<double>(num_vars, 0.0));
  if (num_steps == 1) return ydot;

  // Generalized-alpha coefficients (identical to the constructor).
  const double alpha_m = 0.5 * (3.0 - rho) / (1.0 + rho);
  const double alpha_f = 1.0 / (1.0 + rho);
  const double gamma = 0.5 + alpha_m - alpha_f;
  const double ydot_init_coeff = 1.0 - 1.0 / gamma;

  // Particular solution: forward recursion seeded with ydot[0] = 0. This
  // leaves the (unknown) initial derivative as a free homogeneous mode that
  // scales as ydot_init_coeff^n.
  for (size_t n = 1; n < num_steps; n++) {
    const double dt = times[n] - times[n - 1];
    if (dt <= 0.0) {
      throw std::runtime_error(
          "[Integrator::compute_ydot] Time vector must be strictly "
          "increasing.");
    }
    const double inv = 1.0 / (gamma * dt);
    for (size_t v = 0; v < num_vars; v++) {
      ydot[n][v] = (y_series[n][v] - y_series[n - 1][v]) * inv +
                   ydot_init_coeff * ydot[n - 1][v];
    }
  }

  // Close the periodic loop. With ydot[0] free, the full solution is
  // ydot[n] = ydot_init_coeff^n * ydot[0] + p[n], where p[n] is the particular
  // solution above. Enforcing periodicity ydot[N-1] == ydot[0] fixes
  //   ydot[0] = p[N-1] / (1 - ydot_init_coeff^(N-1)).
  const double homog = std::pow(ydot_init_coeff,
                                static_cast<double>(num_steps - 1));
  const double denom = 1.0 - homog;
  if (std::abs(denom) > 1e-12) {
    for (size_t v = 0; v < num_vars; v++) {
      const double ydot0 = ydot[num_steps - 1][v] / denom;
      double coeff = 1.0;
      for (size_t n = 0; n < num_steps; n++) {
        ydot[n][v] += coeff * ydot0;
        coeff *= ydot_init_coeff;
      }
    }
  }

  return ydot;
}

void Integrator::compute_collocation(
    const std::vector<std::vector<double>>& y_series,
    const std::vector<double>& times, double rho,
    std::vector<std::vector<double>>& y_collocation,
    std::vector<std::vector<double>>& ydot_collocation) {
  const double alpha_m = 0.5 * (3.0 - rho) / (1.0 + rho);
  const double alpha_f = 1.0 / (1.0 + rho);

  // Reconstruct the nodal time derivatives with the same generalized-alpha
  // relation, then interpolate to the alpha-collocation points.
  const std::vector<std::vector<double>> ydot =
      compute_ydot(y_series, times, rho);

  y_collocation.clear();
  ydot_collocation.clear();
  const size_t num_steps = y_series.size();
  if (num_steps < 2) return;
  const size_t num_vars = y_series[0].size();
  y_collocation.assign(num_steps - 1, std::vector<double>(num_vars));
  ydot_collocation.assign(num_steps - 1, std::vector<double>(num_vars));
  for (size_t n = 0; n + 1 < num_steps; n++) {
    for (size_t v = 0; v < num_vars; v++) {
      y_collocation[n][v] =
          y_series[n][v] + alpha_f * (y_series[n + 1][v] - y_series[n][v]);
      ydot_collocation[n][v] =
          ydot[n][v] + alpha_m * (ydot[n + 1][v] - ydot[n][v]);
    }
  }
}
