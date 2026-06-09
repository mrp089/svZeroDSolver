// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file Integrator.h
 * @brief Integrator source file
 */
#ifndef SVZERODSOLVER_ALGEBRA_INTEGRATOR_HPP_
#define SVZERODSOLVER_ALGEBRA_INTEGRATOR_HPP_

#include <Eigen/Dense>
#include <vector>

#include "Model.h"
#include "State.h"

/**
 * @brief Generalized-alpha integrator
 *
 * This class handles the time integration scheme for solving 0D blood
 * flow system using the generalized-\f$\alpha\f$ method \cite JANSEN2000305.
 *
 * Mathematical details are available on the <a
 * href="https://simvascular.github.io/documentation/rom_simulation.html#0d-solver-theory">SimVascular
 * documentation</a>.
 */

class Integrator {
 private:
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

 public:
  /**
   * @brief Construct a new Integrator object
   *
   * @param model The model to simulate
   * @param time_step_size Time step size for generalized-alpha step
   * @param rho Spectral radius for generalized-alpha step
   * @param atol Absolut tolerance for non-linear iteration termination
   * @param max_iter Maximum number of non-linear iterations
   */
  Integrator(Model* model, double time_step_size, double rho, double atol,
             int max_iter);

  /**
   * @brief Construct a new Integrator object
   *
   */
  Integrator();

  /**
   * @brief Destroy the Integrator object
   *
   */
  ~Integrator();

  /**
   * @brief Delete dynamically allocated memory (in class member
   * SparseSystem<double> system).
   */
  void clean();

  /**
   * @brief Update integrator parameter and system matrices with model parameter
   * updates.
   *
   * @param time_step_size Time step size for 0D model
   */
  void update_params(double time_step_size);

  /**
   * @brief Perform a time step
   *
   * @param state Current state
   * @param time Current time
   * @return New state
   */
  State step(const State& state, double time);

  /**
   * @brief Get average number of nonlinear iterations in all step calls
   *
   * @return Average number of nonlinear iterations in all step calls
   *
   */
  double avg_nonlin_iter();

  /**
   * @brief Reconstruct consistent time derivatives from a state time series
   *
   * Computes the time derivatives \f$\dot{\boldsymbol{y}}\f$ that are consistent
   * with a given time series of states \f$\boldsymbol{y}\f$ under the same
   * generalized-\f$\alpha\f$ relation used by step(). Over a single step the
   * scheme enforces
   * \f[
   *   \dot{\boldsymbol{y}}_{n+1} =
   *     \frac{\boldsymbol{y}_{n+1} - \boldsymbol{y}_n}{\gamma \, \Delta t_n}
   *     + \left(1 - \frac{1}{\gamma}\right) \dot{\boldsymbol{y}}_n,
   * \f]
   * with \f$\gamma\f$ derived from the spectral radius \f$\rho\f$. The free
   * initial derivative is fixed by assuming the series is periodic, i.e.
   * \f$\dot{\boldsymbol{y}}_0 = \dot{\boldsymbol{y}}_{N-1}\f$, which is the
   * natural choice for calibration data spanning full cardiac cycles. This
   * lets the calibrator take only the states and a time vector as input and
   * derive consistent derivatives instead of requiring them as a separate
   * input.
   *
   * @param y_series State time series indexed as [time point][variable]
   * @param times Time vector, one entry per time point
   * @param rho Spectral radius of generalized-alpha
   * @return Time derivatives with the same shape as y_series
   */
  static std::vector<std::vector<double>> compute_ydot(
      const std::vector<std::vector<double>>& y_series,
      const std::vector<double>& times, double rho);

  /**
   * @brief Build generalized-alpha collocation states from a state time series
   *
   * The generalized-alpha integrator does not enforce the 0D residual at the
   * time nodes, but at the intermediate collocation points
   * \f$\boldsymbol{y}_{n+\alpha_f}\f$ and \f$\dot{\boldsymbol{y}}_{n+\alpha_m}\f$.
   * Evaluating the calibration residual at the nodes therefore leaves an
   * \f$O(\Delta t)\f$ inconsistency that limits parameter recovery at coarse
   * sampling. This helper reconstructs the nodal derivatives with
   * compute_ydot() and then forms, for each of the \f$N-1\f$ intervals, the
   * collocation states
   * \f[
   *   \boldsymbol{y}_{n+\alpha_f} = \boldsymbol{y}_n
   *     + \alpha_f (\boldsymbol{y}_{n+1} - \boldsymbol{y}_n), \qquad
   *   \dot{\boldsymbol{y}}_{n+\alpha_m} = \dot{\boldsymbol{y}}_n
   *     + \alpha_m (\dot{\boldsymbol{y}}_{n+1} - \dot{\boldsymbol{y}}_n).
   * \f]
   * Evaluating the residual at these points makes it exactly consistent with
   * the integrator, so the calibrator recovers the generating parameters at any
   * sampling resolution.
   *
   * @param y_series State time series indexed as [time point][variable]
   * @param times Time vector, one entry per time point
   * @param rho Spectral radius of generalized-alpha
   * @param[out] y_collocation Collocation states, [interval][variable], size N-1
   * @param[out] ydot_collocation Collocation derivatives, same shape
   */
  static void compute_collocation(
      const std::vector<std::vector<double>>& y_series,
      const std::vector<double>& times, double rho,
      std::vector<std::vector<double>>& y_collocation,
      std::vector<std::vector<double>>& ydot_collocation);
};

#endif  // SVZERODSOLVER_ALGEBRA_INTEGRATOR_HPP_
