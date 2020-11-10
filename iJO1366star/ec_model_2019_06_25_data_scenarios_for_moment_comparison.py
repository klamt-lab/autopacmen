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
"""ec_model_data_scenarios_for_moment_comparison.py

Source:
*Adadi et al., 2012: Prediction of Microbial Growth Rate versus Biomass Yield
 by a Metabolic Network with Kinetic Parameters (PLoS Computational Biology
 e1002575), Supplementary Table 3
"""

exchange_reactions_by_c_source = {
    "Acetate": "EX_ac_e",
    "N-acetylglucosamine": "EX_acgam_e",
    "Glycerol": "EX_glyc_e",
    "Oxoglutarate": "EX_akg_e",
    "L-Alanine": "EX_ala__L_e",
    "Pyruvate": "EX_pyr_e",
    "Fructose": "EX_fru_e",
    "Guanosine": "EX_gsn_e",
    "Fumarate": "EX_fum_e",
    "Ribose": "EX_rib__D_e",
    "Galactose": "EX_gal_e",
    "L-Lactate": "EX_lac__L_e",
    "Gluconate": "EX_glcn_e",
    "Sorbitol": "EX_sbt__D_e",
    "Glucosamine": "EX_gam_e",
    "L-Malate": "EX_mal__L_e",
    "Succinate": "EX_succ_e",
    "Glucose": "EX_glc__D_e",
    "Glucose_anaerobe": "EX_glc__D_e",
    "Glucose_aerobe": "EX_glc__D_e",
    "Maltose": "EX_malt_e",
    "Glucose 6-Phosphate": "EX_g6p_e",
    "Mannitol": "EX_mnl_e",
    "Trehalose": "EX_tre_e",
    "Xylose": "EX_xyl__D_e",
    "Mannose": "EX_man_e",
}

in_vivo_results_by_c_source = {
    "Acetate": 0.29,
    "Glycerol": 0.47,
    "Oxoglutarate": 0.24,
    "L-Alanine": 0.24,
    "Pyruvate": 0.41,
    "Fructose": 0.54,
    "Guanosine": 0.37,
    "N-acetylglucosamine": 0.61,
    "Fumarate": 0.47,
    "Ribose": 0.41,
    "L-Lactate": 0.41,
    "Gluconate": 0.68,
    "Sorbitol": 0.48,
    "Glucosamine": 0.4,
    "L-Malate": 0.55,
    "Succinate": 0.5,
    "Glucose": 0.66,
    "Maltose": 0.52,
    "Glucose 6-Phosphate": 0.78,
    "Mannitol": 0.61,
    "Trehalose": 0.48,
    "Xylose": 0.51,
    "Mannose": 0.35,
    "Galactose": 0.24,
}


def get_ec_model_optimization_scenarios():
    scenarios = {}
    c_sources = in_vivo_results_by_c_source.keys()
    for c_source in c_sources:
        scenario = {}
        exchange_reaction_id = exchange_reactions_by_c_source[c_source]
        scenario["setup"] = {
            exchange_reaction_id: {
                "lower_bound": -1000,
            }
        }
        in_vivo_result = in_vivo_results_by_c_source[c_source]
        scenario["target"] = {
            "reaction": "BIOMASS_Ec_iJO1366_core_53p95M",
            "value": in_vivo_result,
        }
        scenarios[c_source] = scenario

    return scenarios


ec_model_optimization_scenarios = get_ec_model_optimization_scenarios()
