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
"""ec_model_2019_06_25_test_comparison_with_other_models_adadi_et_al_values.py


"""

import cobra
import copy
import matplotlib.pyplot as plt
from scipy.stats import linregress, spearmanr
from ec_model_2019_06_25_data_scenarios_for_moment_comparison import exchange_reactions_by_c_source
from ec_model_2019_06_25_data_standard_exchange_scenario import ec_model_shut_down_reactions
from ec_model_2019_06_25_data_set_up_model import set_up_ec_model_with_sbml

with open("ec_model_2019_06_25_input/c_sources_S3_Adadi_2012.txt", "r") as f:
    lines = f.readlines()
lines = [x.replace("\n", "") for x in lines][1:]
c_sources = []
measured = []
moment = []
fbawmc = []
for line in lines:
    line = line.replace("Â", "±")
    split_lines = line.split("\t")
    c_sources += [split_lines[0]]
    measured += [float(split_lines[1].split("±")[0])]
    moment += [float(split_lines[2])]
    fbawmc += [float(split_lines[3])]

model = set_up_ec_model_with_sbml("ec_model_2019_06_25_output_optimization/iJO1366_sMOMENT_2019_06_25_STANDARD_EXCHANGE_SCENARIO_MANUAL_CHANGES_FMINCON_CHANGE_FACTOR_50.xml", 0.095)

prot_bounds = [0.095]
with model:
    model.reactions.get_by_id("ER_pool_TG_").upper_bound = prot_bounds[0]
    model.reactions.EX_glc__D_e.lower_bound = -1000
    model.optimize()
    model.summary()

thermogecko_prot_pool = []
for prot_bound in prot_bounds:
    model.reactions.get_by_id("ER_pool_TG_").upper_bound = prot_bound
    results = []
    for c_source in c_sources:
        with model:
            for shut_down_reaction in ec_model_shut_down_reactions:
                if "EX_ac_" not in shut_down_reaction:
                    model.reactions.get_by_id(shut_down_reaction).knock_out()
            exchange_id = exchange_reactions_by_c_source[c_source]
            if exchange_id == "":
                results += [0.0]
                continue
            model.reactions.get_by_id(exchange_id).lower_bound = -1000
            solution = model.optimize()
            print("***sMOMENT-enhanced***")
            model.summary()
            results += [solution.fluxes.BIOMASS_Ec_iJO1366_core_53p95M]
    thermogecko_prot_pool.append(copy.deepcopy(results))

original_model = cobra.io.read_sbml_model("ec_model_2019_06_25_input/iJO1366.xml")
original_model.reactions.EX_glc__D_e.lower_bound = .0
normal_fba = []
for c_source in c_sources:
    with original_model:
        for shut_down_reaction in ec_model_shut_down_reactions:
            if "EX_ac_" not in shut_down_reaction:
                original_model.reactions.get_by_id(shut_down_reaction).knock_out()
        exchange_id = exchange_reactions_by_c_source[c_source]
        if exchange_id == "":
            results += [0.0]
            continue
        original_model.reactions.get_by_id(exchange_id).lower_bound = -10
        solution = original_model.optimize()
        # print("~~~Original~~~")
        # original_model.summary()
        normal_fba += [solution.fluxes.BIOMASS_Ec_iJO1366_core_53p95M]

print("")
print(" ================== ")

print("Normal FBA R²:")
slope, intercept, r_value, p_value, std_err = linregress(normal_fba, measured)
print(r_value**2)
print("_Spearman (correlation coefficient, p-value):")
print(spearmanr(normal_fba, measured))
print("_Mean difference:")
differences = [abs(i - j) for i, j in zip(normal_fba, measured)]
print(sum(differences) / len(differences))
plt.scatter(measured, normal_fba)
plt.ylabel('Normal FBA')
plt.xlabel('in vivo')
plt.xlim(left=0)
plt.xlim(right=max(max(normal_fba), max(measured)))
plt.ylim(bottom=0)
plt.ylim(top=max(max(normal_fba), max(measured)))
plt.show()
print("")

print("FBAwMC R²:")
slope, intercept, r_value, p_value, std_err = linregress(fbawmc, measured)
print(r_value**2)
print("_Spearman (correlation coefficient, p-value):")
print(spearmanr(fbawmc, measured))
print("_Mean difference:")
differences = [abs(i - j) for i, j in zip(fbawmc, measured)]
print(sum(differences) / len(differences))
plt.scatter(measured, fbawmc)
plt.ylabel('FBAwMC')
plt.xlabel('in vivo')
plt.xlim(left=0)
plt.xlim(right=max(max(fbawmc), max(measured)))
plt.ylim(bottom=0)
plt.ylim(top=max(max(fbawmc), max(measured)))
plt.show()
print("")

print("MOMENT R²:")
slope, intercept, r_value, p_value, std_err = linregress(moment, measured)
print(r_value**2)
print("_Spearman (correlation coefficient, p-value):")
print(spearmanr(moment, measured))
print("_Mean difference:")
differences = [abs(i - j) for i, j in zip(moment, measured)]
print(sum(differences) / len(differences))
plt.scatter(measured, moment)
plt.ylabel('MOMENT')
plt.xlabel('in vivo')
plt.xlim(left=0)
plt.xlim(right=max(max(moment), max(measured)))
plt.ylim(bottom=0)
plt.ylim(top=max(max(moment), max(measured)))
plt.show()
print("")

i = 0
for result in thermogecko_prot_pool:
    print(f"Thermogecko R² with prot_pool {prot_bounds[i]}:")
    slope, intercept, r_value, p_value, std_err = linregress(result, measured)
    print(r_value**2)
    print("_Spearman (correlation coefficient, p-value):")
    print(spearmanr(result, measured))
    print("_Mean difference:")
    differences = [abs(i - j) for i, j in zip(result, measured)]
    print(sum(differences) / len(differences))

    plt.scatter(measured, result)
    plt.ylabel(f'Thermogecko, {prot_bounds[i]}')
    plt.xlabel('in vivo')
    plt.xlim(left=0)
    plt.xlim(right=max(max(result), max(measured)))
    plt.ylim(bottom=0)
    plt.ylim(top=max(max(result), max(measured)))
    plt.show()
    print("")
    i += 1
