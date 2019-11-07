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
"""ec_model_2019_06_25_test_secretion_behavior_on_different_maximal_glucose_uptake_rates.py

Shows the secretion and general FVA results of a protein-constrained model under different maximal
glucose uptake rates, and returns them as comma-separated values which were later imported into a
spreadsheet program in order to generate the associated figures.
"""

import cobra
import numpy as np
from ec_model_2019_06_25_data_set_up_model import set_up_ec_model_with_sbml


def fba_with_glucose_levels(model, glucose_levels, name):
    values = {}
    values["growth"] = []
    values["growth_var"] = []
    values["prot_pool"] = []
    values["prot_pool_var"] = []
    values["glucose"] = []
    values["glucose_var"] = []
    values["ethanol"] = []
    values["ethanol_var"] = []
    values["formate"] = []
    values["formate_var"] = []
    values["acetate"] = []
    values["acetate_var"] = []
    values["lactate"] = []
    values["lactate_var"] = []
    values["succinate"] = []
    values["succinate_var"] = []
    values["co2"] = []
    values["co2_var"] = []
    values["o2"] = []
    values["o2_var"] = []
    values["o2_min"] = []
    values["o2_max"] = []
    for glucose_level in glucose_levels:
        print("***")
        print(name+f", {glucose_level}:")
        with model:
            model.reactions.EX_glc__D_e.lower_bound = -glucose_level
            fba_solution = model.optimize()
            print("Max prot pool?:", fba_solution.fluxes.ER_pool_TG_ == model.reactions.ER_pool_TG_.upper_bound)
            print("Max glucose?:", fba_solution.fluxes.EX_glc__D_e == model.reactions.EX_glc__D_e.lower_bound)
            # model.summary(fva=1.0)
            print("Used protein pool:", fba_solution.fluxes.ER_pool_TG_)
            print("Used glucose uptake:", fba_solution.fluxes.EX_glc__D_e)
            # model.metabolites.prot_pool.summary()

            fva_reactions = [
                             model.reactions.BIOMASS_Ec_iJO1366_core_53p95M,
                             model.reactions.ER_pool_TG_,
                             model.reactions.EX_glc__D_e,
                             model.reactions.EX_etoh_e,
                             model.reactions.EX_for_e,
                             model.reactions.EX_ac_e,
                             model.reactions.EX_lac__D_e,
                             model.reactions.EX_succ_e,
                             model.reactions.EX_co2_e,
                             model.reactions.EX_o2_e,
                            ]
            fva_solution = cobra.flux_analysis.flux_variability_analysis(model, reaction_list=fva_reactions, loopless=False, fraction_of_optimum=1.0)
            values["growth"].append((fva_solution.loc["BIOMASS_Ec_iJO1366_core_53p95M", "minimum"] + fva_solution.loc["BIOMASS_Ec_iJO1366_core_53p95M", "maximum"]) / 2)
            values["growth_var"].append(abs(fva_solution.loc["BIOMASS_Ec_iJO1366_core_53p95M", "maximum"] - fva_solution.loc["BIOMASS_Ec_iJO1366_core_53p95M", "minimum"]))

            values["prot_pool"].append((fva_solution.loc["ER_pool_TG_", "minimum"] + fva_solution.loc["ER_pool_TG_", "maximum"]) / 2)
            values["prot_pool_var"].append(abs(fva_solution.loc["ER_pool_TG_", "maximum"] - fva_solution.loc["ER_pool_TG_", "minimum"]))

            values["glucose"].append((fva_solution.loc["EX_glc__D_e", "minimum"] + fva_solution.loc["EX_glc__D_e", "maximum"]) / 2)
            values["glucose_var"].append(abs(fva_solution.loc["EX_glc__D_e", "maximum"] - fva_solution.loc["EX_glc__D_e", "minimum"]))

            values["ethanol"].append((fva_solution.loc["EX_etoh_e", "minimum"] + fva_solution.loc["EX_etoh_e", "maximum"]) / 2)
            values["ethanol_var"].append(abs(fva_solution.loc["EX_etoh_e", "maximum"] - fva_solution.loc["EX_etoh_e", "minimum"]))

            values["formate"].append((fva_solution.loc["EX_for_e", "minimum"] + fva_solution.loc["EX_for_e", "maximum"]) / 2)
            values["formate_var"].append(abs(fva_solution.loc["EX_for_e", "maximum"] - fva_solution.loc["EX_for_e", "minimum"]))

            values["acetate"].append((fva_solution.loc["EX_ac_e", "minimum"] + fva_solution.loc["EX_ac_e", "maximum"]) / 2)
            values["acetate_var"].append(abs(fva_solution.loc["EX_ac_e", "maximum"] - fva_solution.loc["EX_ac_e", "minimum"]))

            values["lactate"].append((fva_solution.loc["EX_lac__D_e", "minimum"] + fva_solution.loc["EX_lac__D_e", "maximum"]) / 2)
            values["lactate_var"].append(abs(fva_solution.loc["EX_lac__D_e", "maximum"] - fva_solution.loc["EX_lac__D_e", "minimum"]))

            values["succinate"].append((fva_solution.loc["EX_succ_e", "minimum"] + fva_solution.loc["EX_succ_e", "maximum"]) / 2)
            values["succinate_var"].append(abs(fva_solution.loc["EX_succ_e", "maximum"] - fva_solution.loc["EX_succ_e", "minimum"]))

            values["co2"].append((fva_solution.loc["EX_co2_e", "minimum"] + fva_solution.loc["EX_co2_e", "maximum"]) / 2)
            values["co2_var"].append(abs(fva_solution.loc["EX_co2_e", "maximum"] - fva_solution.loc["EX_co2_e", "minimum"]))

            # Only maximum as the O2 variability at very low C source uptakes can be insanely high D:
            values["o2"].append(-fva_solution.loc["EX_o2_e", "maximum"])
            values["o2_var"].append(abs(fva_solution.loc["EX_o2_e", "maximum"] - fva_solution.loc["EX_o2_e", "minimum"]))
            values["o2_min"].append(abs(fva_solution.loc["EX_o2_e", "minimum"]))
            values["o2_max"].append(abs(fva_solution.loc["EX_o2_e", "maximum"]))

    glucose_level_strs = [str(x) for x in glucose_levels]

    print("Max glucose uptake" + ";" + ";".join(glucose_level_strs))
    print("Actual glucose uptake" + ";" + ";".join([str(-x) for x in values["glucose"]]))
    print("Acetate" + ";" + ";".join([str(x) for x in values["acetate"]]))
    print("Ethanol" + ";" + ";".join([str(x) for x in values["ethanol"]]))
    print("Formate" + ";" + ";".join([str(x) for x in values["formate"]]))
    print("Lactate" + ";" + ";".join([str(x) for x in values["lactate"]]))
    print("Succinate" + ";" + ";".join([str(x) for x in values["succinate"]]))
    print("CO2" + ";" + ";".join([str(x) for x in values["co2"]]))
    print("O2" + ";" + ";".join([str(x) for x in values["o2"]]))
    print("Used protein pool" + ";" + ";".join([str(x) for x in values["prot_pool"]]))
    print("Growth rate" + ";" + ";".join([str(x) for x in values["growth"]]))

    print("Variability of glucose" + ";" + ";".join([str(x) for x in values["glucose_var"]]))
    print("Variability of acetate" + ";" + ";".join([str(x) for x in values["acetate_var"]]))
    print("Variability of ethanol" + ";" + ";".join([str(x) for x in values["ethanol_var"]]))
    print("Variability of formate" + ";" + ";".join([str(x) for x in values["formate_var"]]))
    print("Variability of lactate" + ";" + ";".join([str(x) for x in values["lactate_var"]]))
    print("Variability of succinate" + ";" + ";".join([str(x) for x in values["succinate_var"]]))
    print("Variability of CO2" + ";" + ";".join([str(x) for x in values["co2_var"]]))
    print("Variability of O2" + ";" + ";".join([str(x) for x in values["o2_var"]]))
    print("Variability of protein pool" + ";" + ";".join([str(x) for x in values["prot_pool_var"]]))

    print("Minimal O2" + ";" + ";".join([str(x) for x in values["o2_min"]]))
    print("Maximal O2" + ";" + ";".join([str(x) for x in values["o2_max"]]))


model = set_up_ec_model_with_sbml("ec_model_2019_06_25_output_optimization/iJO1366_sMOMENT_2019_06_25_STANDARD_EXCHANGE_SCENARIO_MANUAL_CHANGES_FMINCON_CHANGE_FACTOR_50.xml", 0.095)

aerobic_space = list(np.linspace(.14, 9.53, 25)) + list(np.linspace(9.54, 13.829, 25))
fba_with_glucose_levels(model, aerobic_space, "Aerobe")
"""
print("===")
model.reactions.EX_o2_e.lower_bound = 0
anaerobic_space = list(np.linspace(1.26, 16.69, 25)) + list(np.linspace(16.7, 24.987, 25))
fba_with_glucose_levels(model, anaerobic_space, "Anaerobe")
"""
