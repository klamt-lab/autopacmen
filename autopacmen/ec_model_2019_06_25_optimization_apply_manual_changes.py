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
"""ec_model_2019_06_25_optimization_apply_manual_changes.py"""

import cobra
from submodules.apply_manual_changes import apply_manual_changes


KCAT_CHANGE_FACTORS = {
    "PFL": ("", 10),
    "ACALD": ("reverse", 5),
    "ALCD2x": ("reverse", 1/10),
    "LDH_D": ("reverse", 1/10),
}

model: cobra.Model = cobra.io.read_sbml_model("ec_model_2019_06_25_output/iJO1366_sMOMENT_2019_06_25_STANDARD_EXCHANGE_SCENARIO.xml")
model = apply_manual_changes(model, KCAT_CHANGE_FACTORS)
cobra.io.write_sbml_model(model, "ec_model_2019_06_25_output_optimization/iJO1366_sMOMENT_2019_06_25_STANDARD_EXCHANGE_SCENARIO_MANUAL_CHANGES.xml")
