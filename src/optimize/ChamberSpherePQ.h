// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file ChamberSpherePQ.h
 * @brief Pressure/flow-only calibration of a ChamberSphere block
 *
 * Unlike BloodVessel, a ChamberSphere has a hidden storage state (its
 * volume/radius), equal to the time-integral of net flow. Its first-order DAE
 * residuals depend on the internal variables (radius, velo, stress, tau) and
 * their derivatives, which instantaneous port (P, Q) data cannot provide. The
 * point-wise calibration path (Block::update_gradient) therefore cannot be used
 * with P/Q-only observations.
 *
 * This module implements the alternative: reconstruct the full internal
 * trajectory globally from the port pressure and net flow (integrate flow to get
 * volume/radius, differentiate to get velocities, then evaluate the algebraic
 * stress/active-stress relations), and fit the one remaining dynamic equation
 * (active stress, equation 3) over a user-selected subset of parameters. The
 * absolute cavity volume is fixed to the reference state (radius r = 0 at the
 * first time point, i.e. V0 = 4/3 * pi * radius0^3).
 */

#ifndef SVZERODSOLVER_OPTIMIZE_CHAMBERSPHEREPQ_HPP_
#define SVZERODSOLVER_OPTIMIZE_CHAMBERSPHEREPQ_HPP_

#include <nlohmann/json.hpp>

/**
 * @brief Detect whether the configuration calibrates a ChamberSphere.
 *
 * A ChamberSphere is always calibrated from port pressure/flow only (the data
 * are assumed to be P/Q), so this returns true whenever the config contains a
 * ChamberSphere vessel.
 *
 * @param config JSON configuration for the calibration
 * @return true if the P/Q ChamberSphere path should be used
 */
bool is_chamber_pq_mode(const nlohmann::json& config);

/**
 * @brief Calibrate a ChamberSphere from port pressure and flow only.
 *
 * @param config JSON configuration (a single ChamberSphere vessel with a
 * ``calibrate`` list, ``y``/``dy``/``time`` holding only the four port
 * variables, and ``calibration_parameters``)
 * @return Calibrated JSON configuration
 */
nlohmann::json calibrate_chamber_pq(const nlohmann::json& config);

#endif  // SVZERODSOLVER_OPTIMIZE_CHAMBERSPHEREPQ_HPP_
