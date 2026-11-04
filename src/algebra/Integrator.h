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
  generalized_alpha,  ///< Generalized-alpha (default, \cite JANSEN2000305)
  consistent_stiff    ///< Stiffly-stable damped-Newton with adaptive substeps,
                      ///< for the stiff active law of \cite genet23
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
 *
 * The coefficients stored here depend only on the spectral radius. The
 * step-size-dependent coefficients (gamma*dt and its Jacobian factor) are
 * computed locally in each step, so that a sub-stepping scheme can vary dt
 * without saving and restoring shared state.
 */
class TimeIntegrator {
 protected:
  double alpha_m{0.0};          ///< Generalized-alpha coefficient alpha_m
  double alpha_f{0.0};          ///< Generalized-alpha coefficient alpha_f
  double gamma{0.0};            ///< Generalized-alpha coefficient gamma
  double time_step_size{0.0};   ///< Nominal time step size
  double ydot_init_coeff{0.0};  ///< Predictor coefficient for ydot
  double atol{0.0};             ///< Absolute nonlinear residual tolerance
  int max_iter{0};              ///< Maximum number of nonlinear iterations
  bool max_iter_error_to_warning{false};  ///< Warn instead of throwing when
                                          ///< max_iter is reached
  int size{0};                            ///< Size of the system
  int n_iter{0};                          ///< Number of step calls
  int n_nonlin_iter{0};  ///< Number of nonlinear iterations over all steps
  Eigen::Matrix<double, Eigen::Dynamic, 1> y_af;     ///< y at the mid-point
  Eigen::Matrix<double, Eigen::Dynamic, 1> ydot_am;  ///< ydot at the mid-point
  SparseSystem system;    ///< System matrices and residual
  Model* model{nullptr};  ///< Model to integrate

  /**
   * @brief Set up coefficients and reserve system memory (shared by schemes).
   *
   * @param model The model to simulate
   * @param time_step_size Time step size
   * @param rho Spectral radius of the generalized-alpha method
   * @param atol Absolute tolerance for the nonlinear iterations
   * @param max_iter Maximum number of nonlinear iterations
   * @param max_iter_error_to_warning If true, print a warning instead of
   * throwing an error when the maximum nonlinear iterations is reached.
   */
  void init(Model* model, double time_step_size, double rho, double atol,
            int max_iter, bool max_iter_error_to_warning = false);

  /**
   * @brief Evaluate the residual at the generalized mid-point for a candidate
   * new state. Also refreshes the solution-dependent system contributions
   * (needed before update_jacobian).
   *
   * @param old_state State at the beginning of the step
   * @param new_state Candidate state at the end of the step
   * @return Infinity norm of the residual
   */
  double assemble_residual(const State& old_state, const State& new_state);

 public:
  /**
   * @brief Construct an empty time integrator
   */
  TimeIntegrator() {}

  /**
   * @brief Destroy the time integrator
   */
  virtual ~TimeIntegrator() {}

  /**
   * @brief Perform a time step
   *
   * @param state Current state
   * @param time Current time
   * @return New state
   */
  virtual State step(const State& state, double time) = 0;

  /**
   * @brief Update the time-step size and refresh constant/time contributions.
   *
   * @param time_step_size Time step size for 0D model
   */
  virtual void update_params(double time_step_size);

  /**
   * @brief Delete dynamically allocated memory in the sparse system.
   */
  void clean();

  /**
   * @brief Get average number of nonlinear iterations in all step calls
   *
   * @return Average number of nonlinear iterations in all step calls
   */
  double avg_nonlin_iter();
};

/**
 * @brief Generalized-alpha integrator (default scheme).
 */
class GeneralizedAlpha : public TimeIntegrator {
 public:
  /**
   * @brief Construct an empty generalized-alpha integrator
   */
  GeneralizedAlpha() {}

  /**
   * @brief Construct a new generalized-alpha integrator
   *
   * @param model The model to simulate
   * @param time_step_size Time step size for generalized-alpha step
   * @param rho Spectral radius for generalized-alpha step
   * @param atol Absolute tolerance for non-linear iteration termination
   * @param max_iter Maximum number of non-linear iterations
   * @param max_iter_error_to_warning If true, print a warning instead of
   * throwing error when maximum iterations is reached
   */
  GeneralizedAlpha(Model* model, double time_step_size, double rho, double atol,
                   int max_iter, bool max_iter_error_to_warning = false);

  /**
   * @brief Perform a generalized-alpha time step
   *
   * @param state Current state
   * @param time Current time
   * @return New state
   */
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
 * a single sub-step with a full Newton step succeeds. Blocks with persistent
 * memory (e.g. ImpedanceBC) see one accepted step of the nominal size.
 */
class ConsistentStiffIntegrator : public TimeIntegrator {
  int max_subdiv{9};      ///< Max number of times the step size is halved
  int max_backtrack{12};  ///< Max backtracking line-search iterations

  /**
   * @brief One damped-Newton step of size `dt` starting at `time`.
   *
   * @param old_state State at the beginning of the sub-step
   * @param new_state State at the end of the sub-step (output)
   * @param dt Sub-step size
   * @param time Time at the beginning of the sub-step
   * @return true on convergence (new_state filled), false otherwise.
   */
  bool damped_step(const State& old_state, State& new_state, double dt,
                   double time);

 public:
  /**
   * @brief Construct an empty stiff integrator
   */
  ConsistentStiffIntegrator() {}

  /**
   * @brief Construct a new stiff integrator
   *
   * @param model The model to simulate
   * @param time_step_size Nominal time step size
   * @param rho Ignored; the scheme always uses spectral radius 0
   * @param atol Absolute tolerance for non-linear iteration termination
   * @param max_iter Maximum number of non-linear iterations per sub-step
   * @param max_iter_error_to_warning Kept for interface compatibility; a
   * sub-step that does not converge is retried with a smaller step instead
   */
  ConsistentStiffIntegrator(Model* model, double time_step_size, double rho,
                            double atol, int max_iter,
                            bool max_iter_error_to_warning = false);

  /**
   * @brief Perform a time step, sub-stepping if needed
   *
   * @param state Current state
   * @param time Current time
   * @return New state
   */
  State step(const State& state, double time) override;
};

/**
 * @brief Create a time integrator of the requested type.
 *
 * @param type Integration scheme
 * @param model The model to simulate
 * @param time_step_size Time step size
 * @param rho Spectral radius of the generalized-alpha method
 * @param atol Absolute tolerance for non-linear iteration termination
 * @param max_iter Maximum number of non-linear iterations
 * @param max_iter_error_to_warning If true, print a warning instead of
 * throwing error when maximum iterations is reached
 * @return The time integrator
 */
std::unique_ptr<TimeIntegrator> make_integrator(
    IntegratorType type, Model* model, double time_step_size, double rho,
    double atol, int max_iter, bool max_iter_error_to_warning = false);

#endif  // SVZERODSOLVER_ALGEBRA_INTEGRATOR_HPP_
