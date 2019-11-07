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
"""ec_model_original_generate_comparable_model.py"""

import cobra
from ec_model_data_set_up_model import set_up_ec_model_with_sbml
from submodules.helper_create_model import get_irreversible_model

model = set_up_ec_model_with_sbml("ec_model_2019_06_25_input/iJO1366.xml", .225)
model = get_irreversible_model(model)
cobra.io.write_sbml_model(model, "ec_model_2019_06_25_input/iJO1366_saved_by_cobrapy_and_separated_reversible_reactions.xml")
model = set_up_ec_model_with_sbml("ec_model_2019_06_25_input/iJO1366_saved_by_cobrapy_and_separated_reversible_reactions.xml", .225)
cobra.io.write_sbml_model(model, "ec_model_2019_06_25_input/iJO1366_saved_by_cobrapy_and_separated_reversible_reactions_SHUT_DOWN_SCENARIO.xml")
