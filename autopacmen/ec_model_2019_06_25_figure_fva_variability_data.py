#!/usr/bin/env python3
#
# Copyright 2019-2020 PSB
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
"""ec_model_2019_06_25_test_fva_variability_data.py

This is the script for the generation of the comparative cumulative distributions
of flux variabilities of iJO1366 vs. iJO1366* at a given maximal glucose uptake
rate. The box-form in the legends was edited later to lines for aesthetic reasons (no
actual data was changed :-).
"""

import cobra
import copy
import statistics
import xlsxwriter

import matplotlib.pyplot as plt


def _get_fva_statistics(model):
    fva_dict = {}
    minimums = []
    maximums = []
    variabilities = []
    for reaction in model.reactions:
        if reaction.id == "ER_pool_TG_":
            continue

        # Minimum
        with model:
            if reaction.id.endswith("reverse"):
                forward_reaction_name = reaction.id.replace("reverse", "forward")
                model.reactions.get_by_id(forward_reaction_name).knock_out()
            if reaction.id.endswith("forward"):
                reverse_reaction_name = reaction.id.replace("forward", "reverse")
                model.reactions.get_by_id(reverse_reaction_name).knock_out()

            model.objective = model_original.problem.Objective(-1 * model.reactions.get_by_id(reaction.id).flux_expression)
            minimization_solution = model.optimize()
            minimum = minimization_solution.fluxes[reaction.id]

        # Maximum
        with model:
            if reaction.id.endswith("reverse"):
                forward_reaction_name = reaction.id.replace("reverse", "forward")
                model.reactions.get_by_id(forward_reaction_name).knock_out()
            if reaction.id.endswith("forward"):
                reverse_reaction_name = reaction.id.replace("forward", "reverse")
                model.reactions.get_by_id(reverse_reaction_name).knock_out()

            model.objective = model_original.problem.Objective(1 * model.reactions.get_by_id(reaction.id).flux_expression)

            maximization_solution = model.optimize()
            maximum = maximization_solution.fluxes[reaction.id]

        variability = abs(maximum - minimum)
        variabilities.append(variability)
        minimums.append(minimum)
        maximums.append(maximum)

        fva_dict[reaction.id] = {}
        fva_dict[reaction.id]["minimum"] = minimum
        fva_dict[reaction.id]["maximum"] = maximum
        fva_dict[reaction.id]["variability"] = variability

    mean_minimum = statistics.mean(minimums)
    mean_maximum = statistics.mean(maximums)
    num_zeroes = len([x for x in variabilities if x < 10e-10])
    mean_variability = statistics.mean([x for x in variabilities if x > 10e-10])
    median_variability = statistics.median([x for x in variabilities if x > 10e-10])

    fva_dict["MEAN"] = {}
    fva_dict["MEAN"]["minimum"] = mean_minimum
    fva_dict["MEAN"]["maximum"] = mean_maximum
    fva_dict["MEAN"]["variability"] = mean_variability
    fva_dict["MEAN"]["variability_median"] = median_variability
    fva_dict["MEAN"]["num_zeroes"] = num_zeroes

    return fva_dict, variabilities


def split_all_reversibles(model: cobra.Model) -> cobra.Model:
    model_reaction_ids = [x.id for x in model.reactions]
    for reaction_id in model_reaction_ids:
        reaction = model.reactions.get_by_id(reaction_id)

        if reaction.lower_bound >= 0:
            continue

        forward_reaction = copy.deepcopy(reaction)
        forward_reaction.upper_bound = reaction.upper_bound
        forward_reaction.lower_bound = 0
        forward_reaction.id += "_FWD"
        model.add_reactions([forward_reaction])

        reverse_reaction = copy.deepcopy(reaction)
        reverse_reaction.id += "_REV"
        reverse_reaction.upper_bound = -reaction.lower_bound
        reverse_reaction.lower_bound = 0
        reverse_reaction_metabolites_copy = copy.deepcopy(reverse_reaction.metabolites)
        for key in list(reverse_reaction_metabolites_copy.keys()):
            reverse_reaction_metabolites_copy[key] *= -2
        reverse_reaction.add_metabolites(reverse_reaction_metabolites_copy)
        model.add_reactions([reverse_reaction])

        model.remove_reactions([reaction])
    return model


def _add_worksheet(workbook, model_name, fva_dict):
    worksheet = workbook.add_worksheet(model_name)

    worksheet.write(0, 0, "Reaction")
    worksheet.write(0, 1, "Minimal flux")
    worksheet.write(0, 2, "Maximal flux")
    worksheet.write(0, 3, "Variability")

    line = 1
    for reaction_id in fva_dict.keys():
        if reaction_id == "MEAN":
            continue
        fva_entry = fva_dict[reaction_id]
        worksheet.write(line, 0, reaction_id)
        worksheet.write(line, 1, fva_entry["minimum"])
        worksheet.write(line, 2, fva_entry["maximum"])
        worksheet.write(line, 3, fva_entry["variability"])
        line += 1

    return workbook


name_original = "iJO1366_shut_down_scenario"
name_smoment = "sMOMENT_2019_06_25_sds_mc_fm_50"
model_smoment = cobra.io.read_sbml_model("ec_model_2019_06_25_output_optimization/iJO1366_sMOMENT_2019_06_25_STANDARD_EXCHANGE_SCENARIO_MANUAL_CHANGES_FMINCON_CHANGE_FACTOR_50.xml")
model_original = cobra.io.read_sbml_model("ec_model_2019_06_25_input/iJO1366_saved_by_cobrapy_and_separated_reversible_reactions_STANDARD_EXCHANGE_SCENARIO.xml")

model_smoment.reactions.get_by_id("EX_glc__D_e").lower_bound = -9.53
model_smoment.reactions.get_by_id("ER_pool_TG_").upper_bound = .095
with model_smoment:
    max_growth_smoment = model_smoment.slim_optimize()

# model_original.reactions.get_by_id("BIOMASS_Ec_iJO1366_core_53p95M").upper_bound = max_growth_smoment
model_original.reactions.get_by_id("EX_glc__D_e").lower_bound = -9.53

model_original = split_all_reversibles(model_original)
model_smoment = split_all_reversibles(model_smoment)

print("First FVA...")
fva_dict_original, variabilites_original = _get_fva_statistics(model_original)
print("Second FVA...")
fva_dict_smoment, variabilities_smoment = _get_fva_statistics(model_smoment)

fig, ax = plt.subplots(figsize=(8, 4))

# From: https://matplotlib.org/3.1.0/gallery/statistics/histogram_cumulative.html
# plot the cumulative histogram
variabilites_original = [x for x in variabilites_original if x > 10e-10]
variabilities_smoment = [x for x in variabilities_smoment if x > 10e-10]
n_bins = 250
n, bins, patches = ax.hist(variabilites_original, n_bins, density=True, histtype='step',
                           cumulative=True, label=f'iJO1366 (n={len(variabilites_original)})')

# Overlay a reversed cumulative histogram.
ax.hist(variabilities_smoment, bins=bins, density=True, histtype='step', cumulative=True,
        label=f'iJO1366* (n={len(variabilities_smoment)})')

# Add titles
ax.grid(True)
ax.legend(loc='lower center')
ax.set_xlabel('Variability (mmol/(gDW*h))')
ax.set_ylabel('Cumulative distribution')
plt.show()
fig.savefig('ec_model_2019_06_25_test_fva/fva_comparison_2019_06_25_fitted_mc_fm_50_vs_original_iJO1366_both_with_standard_exchange_scenario.png')


# Create XLSX
print("Write XLSX \\o/")
workbook = xlsxwriter.Workbook(f"ec_model_2019_06_25_test_fva/fva_comparison_2019_06_25_fitted_mc_fm_50_vs_original_iJO1366_both_with_standard_exchange_scenario.xlsx")
workbook = _add_worksheet(workbook, name_original, fva_dict_original)
workbook = _add_worksheet(workbook, name_smoment, fva_dict_smoment)

worksheet = workbook.add_worksheet("Statistics")

worksheet.write(0, 0, "Model name")
worksheet.write(0, 1, "Mean variability (variabilities > 10⁻10)")
worksheet.write(0, 2, "Median variability (variabilities > 10⁻10)")
worksheet.write(0, 3, "Number reactions with 0 variability")

line = 1
names = [name_original, name_smoment]
for fva_dict in (fva_dict_original, fva_dict_smoment):
    worksheet.write(line, 0, names[line-1])
    worksheet.write(line, 1, fva_dict["MEAN"]["variability"])
    worksheet.write(line, 2, fva_dict["MEAN"]["variability_median"])
    worksheet.write(line, 3, fva_dict["MEAN"]["num_zeroes"])

    line += 1

workbook.close()
