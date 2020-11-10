#!/usr/bin/env python3
#
# Copyright 2019 PSB
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""ec_model_data_set_up_model.py"""

import cobra
from ec_model_2019_06_25_data_standard_exchange_scenario import ec_model_shut_down_reactions


def set_up_ec_model(model: cobra.Model, prot_pool: float):
    try:
        model.reactions.get_by_id("ER_pool_TG_").upper_bound = prot_pool
    except Exception:
        print("INFO: Model has no protein pool reaction")
    model.reactions.EX_glc__D_e.lower_bound = .0

    for shut_down_reaction in ec_model_shut_down_reactions:
        model.reactions.get_by_id(shut_down_reaction).knock_out()

    return model


def set_up_ec_model_with_sbml(sbml_path: cobra.Model, prot_pool: float):
    model = cobra.io.read_sbml_model(sbml_path)
    model = set_up_ec_model(model, prot_pool)
    return model
