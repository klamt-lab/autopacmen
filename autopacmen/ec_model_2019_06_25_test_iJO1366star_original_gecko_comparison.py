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
"""ec_model_2019_06_25_test_iJO*_gecko_analogon_original_gecko_comparison.py

Compares iJO1366* with a GECKOed version (Sánchez et al., 2017) of iJO1366.
"""

from ec_model_2019_06_25_data_set_up_model import set_up_ec_model_with_sbml


def fba_with_glucose_levels(model_gecko, model_smoment, glucose_levels, name):
    for glucose_level in glucose_levels:
        print("***")
        print(name+f", {glucose_level}:")
        print("~1. GECKO :D...")
        with model_gecko:
            model_gecko.reactions.EX_glc__D_e.lower_bound = -glucose_level
            fba_solution = model_gecko.optimize()
            print("Max prot pool?:", fba_solution.fluxes.ER_pool_TG_ == model_gecko.reactions.ER_pool_TG_.upper_bound)
            print("Max glucose?:", fba_solution.fluxes.EX_glc__D_e == model_gecko.reactions.EX_glc__D_e.lower_bound)
            model_gecko.summary(fva=1.0)
            print("Used protein pool:", fba_solution.fluxes.ER_pool_TG_)
            print("Used glucose uptake:", fba_solution.fluxes.EX_glc__D_e)
        print("\n~2. sMOMENT :D...")
        with model_smoment:
            model_smoment.reactions.EX_glc__D_e.lower_bound = -glucose_level
            fba_solution = model_smoment.optimize()
            print("Max prot pool?:", fba_solution.fluxes.ER_pool_TG_ == model_smoment.reactions.ER_pool_TG_.upper_bound)
            print("Max glucose?:", fba_solution.fluxes.EX_glc__D_e == model_smoment.reactions.EX_glc__D_e.lower_bound)
            model_smoment.summary(fva=1.0)
            print("Used protein pool:", fba_solution.fluxes.ER_pool_TG_)
            print("Used glucose uptake:", fba_solution.fluxes.EX_glc__D_e)
        print("=======")
        print("")
        print("")
        print("")


model_gecko = set_up_ec_model_with_sbml("ec_model_2019_06_25_output/iJO1366_2019_06_25_GECKO.xml", 0.095)
model_smoment = set_up_ec_model_with_sbml("ec_model_2019_06_25_output/iJO1366_sMOMENT_2019_06_25_GECKO_ANALOGON.xml", 0.095)

print(len(model_gecko.reactions))
print(len(model_smoment.reactions))
print("***")
print(len(model_gecko.metabolites))
print(len(model_smoment.metabolites))

fba_with_glucose_levels(model_gecko, model_smoment, [1000, 50, 13.9, 9.53, 8.5, 5, 2.5][::-1], "Aerobe")
print("===")
model_gecko.reactions.EX_o2_e.lower_bound = 0
model_smoment.reactions.EX_o2_e.lower_bound = 0
fba_with_glucose_levels(model_gecko, model_smoment, [1000, 50, 20, 16.69, 10, 5, 2.5][::-1], "Anaerobe")
