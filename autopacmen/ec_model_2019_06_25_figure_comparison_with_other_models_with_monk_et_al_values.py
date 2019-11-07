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
"""ec_model_2019_06_25_figure_comparison_with_other_models_with_monk_et_al_values.py

This is the script for the generation of the figure and the calculation of the corresponding
statistical data for the comparison of measured growth rates and predicted growth rates of iJO1366*.
"""

import cobra
import copy
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress, pearsonr, spearmanr
from ec_model_2019_06_25_data_scenarios_for_moment_comparison import exchange_reactions_by_c_source
from ec_model_2019_06_25_data_standard_exchange_scenario import ec_model_shut_down_reactions
from ec_model_2019_06_25_data_set_up_model import set_up_ec_model_with_sbml

with open("ec_model_2019_06_25_input/c_sources_S3_Adadi_2012_with_monk_values.txt", "r") as f:
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
thermogecko_prot_pool = []
for prot_bound in prot_bounds:
    model.reactions.get_by_id("ER_pool_TG_").upper_bound = prot_bound
    results = []
    i = 0
    for c_source in c_sources:
        with model:
            if "anaerobe" in c_source:
                model.reactions.get_by_id("EX_o2_e").lower_bound = 0
                print("O2 shutdown")

            for shut_down_reaction in ec_model_shut_down_reactions:
                if "EX_ac_" not in shut_down_reaction:
                    model.reactions.get_by_id(shut_down_reaction).knock_out()
            exchange_id = exchange_reactions_by_c_source[c_source]
            model.reactions.get_by_id(exchange_id).lower_bound = -1000
            solution = model.optimize()
            print(c_source, ";", solution.fluxes.BIOMASS_Ec_iJO1366_core_53p95M)
            # print(measured[i])
            # print()
            # print("***")
            # model.summary()
            results += [solution.fluxes.BIOMASS_Ec_iJO1366_core_53p95M]

            i += 1
    thermogecko_prot_pool.append(copy.deepcopy(results))

font = {"family": "normal",
        "size": 14}
plt.rc("font", **font)

print("")
print(" ================== ")
i = 0
for result in thermogecko_prot_pool:
    print(f"Thermogecko R² with prot_pool {prot_bounds[i]}:")
    slope, intercept, r_value, p_value, std_err = linregress(result, measured)
    print(r_value**2)
    print("_Spearman (correlation coefficient, p-value):")
    print(spearmanr(result, measured))
    print("_Pearson (correlation coefficient, p-value):")
    print(pearsonr(result, measured))
    print("_Mean difference:")
    differences = [abs(i - j) for i, j in zip(result, measured)]
    print(sum(differences) / len(differences))

    identity = np.linspace(0, max(max(result), max(measured)), 100)

    plt.plot(identity, identity, "k-", linewidth=.5)
    plt.grid()
    plt.rc("axes", axisbelow=True)
    plt.scatter(measured, result, color="k")
    plt.ylabel(f"Predicted growth rate of iJO1366*")
    plt.xlabel("Growth rates measured in vivo")
    plt.xlim(left=0)
    plt.xlim(right=max(max(result), max(measured)))
    plt.ylim(bottom=0)
    plt.ylim(top=max(max(result), max(measured)))
    plt.show()
    print("")
    i += 1
