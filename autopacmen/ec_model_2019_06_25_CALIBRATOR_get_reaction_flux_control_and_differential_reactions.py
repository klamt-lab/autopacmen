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
"""ec_model_2019_06_25_CALIBRATOR_get_reaction_flux_control_and_differential_reactions.py

This script constains an exemplary use of the AutoPACMEN Model Calibrator
Python functions as described in AutoPACMEN's manual.
"""

# External modules
import cobra
# Internal modules
from ec_model_data_scenarios_for_optimization import ec_model_scenarios_for_optimization
from ec_model_data_set_up_model import set_up_ec_model_with_sbml
from submodules.reaction_flux_control_by_scenario import reaction_flux_control_by_scenario
from submodules.get_differential_reactions import get_differential_reactions
from submodules.helper_general import json_write


# Set-up of project
flux_control_folder = "ec_model_2019_06_25_output_optimization/flux_control_data_2019_06_25_manual_changes/"
project_name = "psb_orth"

# Read SBML model
print("Reading SBML model...")
original_thermogecko_sbml_path: str = "ec_model_2019_06_25_output_optimization/iJO1366_sMOMENT_2019_06_25_STANDARD_EXCHANGE_SCENARIO_MANUAL_CHANGES.xml"
model: cobra.Model = set_up_ec_model_with_sbml(original_thermogecko_sbml_path, .095)

# Set protein bound
model.reactions.get_by_id("ER_pool_TG_").upper_bound = .095

# Get flux controlling proteins
print("Getting flux control files...")
reaction_flux_control_by_scenario(model, flux_control_folder, project_name, ec_model_scenarios_for_optimization)

# Get differential proteins
print("Getting differential reactions (Growth)...")
unique_differential_reactions_of_scenarios, _ = \
    get_differential_reactions(list(ec_model_scenarios_for_optimization.keys()), flux_control_folder, project_name,
                               ec_model_scenarios_for_optimization,
                               threshold=(.1) / 1000, print_result=True)

# Get unique reactions in MATLAB style
for scenario_key in unique_differential_reactions_of_scenarios.keys():
    print(f"% {scenario_key}")
    unique_reactions = unique_differential_reactions_of_scenarios[scenario_key]
    for unique_reaction in unique_reactions:
        print(f'"R_{unique_reaction}",')
json_write("ec_model_2019_06_25_output_optimization/iJO1366_sMOMENT_2019_06_25_STANDARD_EXCHANGE_SCENARIO_MANUAL_CHANGES_unique_differential_reactions_of_scenarios.json", unique_differential_reactions_of_scenarios)
