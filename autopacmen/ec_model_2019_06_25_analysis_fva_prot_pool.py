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
"""ec_model_analysis_fva_prot_pool.py


"""

from submodules.fva_prot_pool import fva_prot_pool
from ec_model_data_set_up_model import set_up_ec_model_with_sbml

model = set_up_ec_model_with_sbml("ec_model_output/iJO1366_GECKO.xml", .225)
model.reactions.EX_glc__D_e.lower_bound = -1000
prot_pools = [.3]
prot_pool_metabolite = model.metabolites.prot_pool

fva_prot_pool(model, prot_pools, "")
