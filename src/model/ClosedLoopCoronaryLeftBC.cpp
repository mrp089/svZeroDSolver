// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "ClosedLoopCoronaryLeftBC.h"

#include "Model.h"

void ClosedLoopCoronaryLeftBC::setup_model_dependent_params() {
  if (!ventricle_block_name.empty()) {
    // Decomposed mode: im from own parameters, ventricle P from named block
    im_param_id = global_param_ids[ParamId::IM];
    auto ventricle_block = model->get_block(ventricle_block_name);
    ventricle_var_id = ventricle_block->global_var_ids[0];  // P_in of chamber
  }
}
