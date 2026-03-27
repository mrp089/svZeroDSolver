// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "ClosedLoopCoronaryRightBC.h"

#include "ClosedLoopHeartPulmonary.h"
#include "Model.h"

void ClosedLoopCoronaryRightBC::setup_model_dependent_params() {
  if (!ventricle_block_name.empty()) {
    // Decomposed mode: read im from own parameters and look up ventricle by
    // name
    im_param_id = global_param_ids[ParamId::IM];
    auto ventricle_block = model->get_block(ventricle_block_name);
    ventricle_var_id = ventricle_block->global_var_ids[0];  // P_in of chamber
  } else {
    // Monolithic mode: reach into CLH block
    auto heart_block = model->get_block("CLH");
    im_param_id =
        heart_block->global_param_ids[ClosedLoopHeartPulmonary::ParamId::IMR];
    ventricle_var_id =
        heart_block->global_var_ids[6];  // Solution ID for RV pressure
  }
}
