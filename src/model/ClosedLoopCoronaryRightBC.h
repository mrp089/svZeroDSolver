// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file ClosedLoopCoronaryRightBC.h
 * @brief Right side of ClosedLoopCoronaryBC
 */
#ifndef SVZERODSOLVER_MODEL_CLOSEDLOOPCORONARYRIGHTBC_HPP_
#define SVZERODSOLVER_MODEL_CLOSEDLOOPCORONARYRIGHTBC_HPP_

#include "ClosedLoopCoronaryBC.h"

/**
 * @brief Right side of closed loop coronary boundary condition
 * ClosedLoopCoronaryBC.
 *
 * ### Usage in json configuration file
 *
 *     "boundary_conditions": {
 *         "RCA": {
 *             "type": "ClosedLoopCoronaryRight",
 *             "values": {
 *                 "resistance_artery": 5.757416349,
 *                 "resistance_micro": 9.355801566,
 *                 "resistance_vein": 20.580381989,
 *                 "capacitance_im": 0.003463134,
 *                 "capacitance_artery": 0.001055957
 *             }
 *         }
 *     }
 */
class ClosedLoopCoronaryRightBC : public ClosedLoopCoronaryBC {
 public:
  /**
   * @brief Construct a new ClosedLoopCoronaryRightBC object
   *
   * @param id Global ID of the block
   * @param model The model to which the block belongs
   */
  ClosedLoopCoronaryRightBC(int id, Model* model)
      : ClosedLoopCoronaryBC(id, model,
                             BlockType::closed_loop_coronary_right_bc) {}

  /**
   * @brief Setup parameters that depend on the model
   *
   */
  void setup_model_dependent_params();
};

#endif  // SVZERODSOLVER_MODEL_CLOSEDLOOPCORONARYRIGHTBC_HPP_
