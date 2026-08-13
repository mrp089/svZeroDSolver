// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file Integrator.h
 * @brief Integrator source file
 */
#ifndef SVZERODSOLVER_ALGEBRA_INTEGRATOR_HPP_
#define SVZERODSOLVER_ALGEBRA_INTEGRATOR_HPP_

#include <Eigen/Dense>
#include <memory>
#include <string>

#include "Model.h"
#include "State.h"

/**
 * @brief Available time-integration schemes.
 */
enum class IntegratorType {
  generalized_alpha,  ///< Generalized-alpha (default)
  consistent_stiff    ///< Stiffly-stable damped-Newton with adaptive substeps
};

/**
 * @brief Abstract time integrator (abstraction layer over integration schemes).
 *
 * Concrete schemes (GeneralizedAlpha, ConsistentStiffIntegrator) share the
 * generalized-\f$\alpha\f$ machinery held here (system matrices, intermediate
 * states, coefficients) and differ only in their step / nonlinear-solve
 * strategy. See \cite JANSEN2000305 and the <a
 * href="https://simvascular.github.io/documentation/rom_simulation.html#0d-solver-theory">SimVascular
 * documentation</a>.
 */
class TimeIntegrator {
 protected:
  double alpha_m{0.0};
  double alpha_f{0.0};
  double gamma{0.0};
  double time_step_size{0.0};
  double ydot_init_coeff{0.0};
  double y_coeff{0.0};
  double y_coeff_jacobian{0.0};
  double atol{0.0};
  int max_iter{0};
  int size{0};
  int n_iter{0};
  int n_nonlin_iter{0};
  Eigen::Matrix<double, Eigen::Dynamic, 1> y_af;
  Eigen::Matrix<double, Eigen::Dynamic, 1> ydot_am;
  SparseSystem system;
  Model* model{nullptr};

  /**
   * @brief Set up coefficients and reserve system memory (shared by schemes).
   */
  void init(Model* model, double time_step_size, double rho, double atol,
            int max_iter);

  /**
   * @brief Evaluate the residual at the generalized mid-point for a candidate
   * new state and return its infinity norm. Also refreshes the solution-
   * dependent system contributions (needed before update_jacobian).
   */
  double assemble_residual(const State& old_state, const State& new_state);

 public:
  TimeIntegrator() {}
  virtual ~TimeIntegrator() {}

  /**
   * @brief Perform a time step from `state` at `time` and return the new state.
   */
  virtual State step(const State& state, double time) = 0;

  /**
   * @brief Update the time-step size and refresh constant/time contributions.
   */
  virtual void update_params(double time_step_size);

  /**
   * @brief Delete dynamically allocated memory in the sparse system.
   */
  void clean();

  /**
   * @brief Average number of nonlinear iterations over all steps.
   */
  double avg_nonlin_iter();
};

/**
 * @brief Generalized-alpha integrator (default scheme).
 */
class GeneralizedAlpha : public TimeIntegrator {
 public:
  GeneralizedAlpha() {}
  GeneralizedAlpha(Model* model, double time_step_size, double rho, double atol,
                   int max_iter);
  State step(const State& state, double time) override;
};

/**
 * @brief Stiffly-stable integrator for stiff active mechanics.
 *
 * Uses the generalized-alpha update with a maximally dissipative spectral
 * radius, augmented with (i) a backtracking line search on the residual norm
 * (damped Newton) and (ii) adaptive time sub-stepping (a step that fails to
 * converge is retried with a halved time step, recursively). This provides the
 * robustness needed to integrate the Bestel-Clement-Sorine active law with the
 * velocity-dependent (force-velocity) term, which plain generalized-alpha
 * cannot carry (the reference \cite genet23 relies on a bespoke energy-
 * consistent scheme for the same reason). It reduces to generalized-alpha when
 * a single sub-step with a full Newton step succeeds.
 */
class ConsistentStiffIntegrator : public TimeIntegrator {
  int max_subdiv{9};      ///< Max number of times the step size is halved
  int max_backtrack{12};  ///< Max backtracking line-search iterations

  /**
   * @brief One damped-Newton step of size `dt` starting at `time`.
   * @return true on convergence (new_state filled), false otherwise.
   */
  bool damped_step(const State& old_state, State& new_state, double dt,
                   double time);

 public:
  ConsistentStiffIntegrator() {}
  ConsistentStiffIntegrator(Model* model, double time_step_size, double rho,
                            double atol, int max_iter);
  State step(const State& state, double time) override;
};

/**
 * @brief Create a time integrator of the requested type.
 */
std::unique_ptr<TimeIntegrator> make_integrator(IntegratorType type,
                                                Model* model,
                                                double time_step_size,
                                                double rho, double atol,
                                                int max_iter);

#endif  // SVZERODSOLVER_ALGEBRA_INTEGRATOR_HPP_
