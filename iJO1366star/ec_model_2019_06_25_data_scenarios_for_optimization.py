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
"""ec_model_data_scenarios_for_optimization.py

This file contains the scenarios for the selection of kcat values for the
automatic optimization of iJO1366*.

Sources:
*Monk et al., 2016: Multi-omics Quantification of Species Variation of
 Escherichia coli Links Molecular Features with Strain Phenotypes
 (Cell Systems 3, 238–251), Table 1
*Adadi et al., 2012: Prediction of Microbial Growth Rate versus Biomass Yield
 by a Metabolic Network with Kinetic Parameters (PLoS Computational Biology
 e1002575), Supplementary Table 3
"""

ec_model_scenarios_for_optimization = {
    # Glucose as C source (aerobic) from Monk 2016 paper, target is growth
    "Glucose_Aerobic_Growth": {
        "setup": {
            "EX_glc__D_e": {
                "lower_bound": -1000,
            },
        },
        "target": {
            "reaction": "BIOMASS_Ec_iJO1366_core_53p95M",
            "value": 0.73,
        },
        "weight": 3,
        "substitution_name": "Glucose (aerobic and anaerobic)"
    },
    # Glucose as C source (aerobic) from Monk 2016 paper, target is growth
    "Glucose_Anaerobic_Growth": {
        "setup": {
            "EX_glc__D_e": {
                "lower_bound": -1000,
            },
            "EX_o2_e": {
                "lower_bound": 0,
                "upper_bound": 0,
            },
            "EX_succ_e": {
                "upper_bound": 1000
            },
            "EX_for_e": {
                "upper_bound": 1000
            },
            "EX_glyclt_e": {
                "upper_bound": 1000
            },
            "EX_lac__D_e": {
                "upper_bound": 1000
            },
        },
        "target": {
            "reaction": "BIOMASS_Ec_iJO1366_core_53p95M",
            "value": 0.46,
        },
        "substitution_name": "Glucose (aerobic and anaerobic)"
    },
    # Acetate as C source (aerobic) from Adadi 2012 paper
    "Acetate": {
        "setup": {
            "EX_ac_e": {
                "lower_bound": -1000,
            },
        },
        "target": {
            "reaction": "BIOMASS_Ec_iJO1366_core_53p95M",
            "value": 0.29,
        },
    },
    # Glycerol as C source (aerobic) from Adadi 2012 paper
    "Glycerol": {
        "setup": {
            "EX_glyc_e": {
                "lower_bound": -1000,
            },
        },
        "target": {
            "reaction": "BIOMASS_Ec_iJO1366_core_53p95M",
            "value": 0.47,
        },
    },
    # Oxoglutarate as C source (aerobic)  from Adadi 2012 paper
    "Oxoglutarate": {
        "setup": {
            "EX_akg_e": {
                "lower_bound": -1000,
            },
        },
        "target": {
            "reaction": "BIOMASS_Ec_iJO1366_core_53p95M",
            "value": 0.24,
        },
    },
    # L-Alanine as C source (aerobic) from Adadi 2012 paper
    "L-Alanine": {
        "setup": {
            "EX_ala__L_e": {
                "lower_bound": -1000,
            },
        },
        "target": {
            "reaction": "BIOMASS_Ec_iJO1366_core_53p95M",
            "value": 0.24,
        },
    },
    # Pyruvate as C source (aerobic) from Adadi 2012 paper
    "Pyruvate": {
        "setup": {
            "EX_pyr_e": {
                "lower_bound": -1000,
            },
        },
        "target": {
            "reaction": "BIOMASS_Ec_iJO1366_core_53p95M",
            "value": 0.41,
        },
    },
    # Fructose as C source (aerobic) from Adadi 2012 paper
    "Fructose": {
        "setup": {
            "EX_fru_e": {
                "lower_bound": -1000,
            },
        },
        "target": {
            "reaction": "BIOMASS_Ec_iJO1366_core_53p95M",
            "value": 0.54,
        },
    },
    # Guanosine as C source (aerobic) from Adadi 2012 paper
    "Guanosine": {
        "setup": {
            "EX_gsn_e": {
                "lower_bound": -1000,
            },
        },
        "target": {
            "reaction": "BIOMASS_Ec_iJO1366_core_53p95M",
            "value": 0.37,
        },
    },
    # N-acetylglucosamine as C source (aerobic) from Adadi 2012 paper
    "N-acetylglucosamine": {
        "setup": {
            "EX_acgam_e": {
                "lower_bound": -1000,
            },
        },
        "target": {
            "reaction": "BIOMASS_Ec_iJO1366_core_53p95M",
            "value": 0.61,
        },
    },
    # Fumarate as C source (aerobic) from Adadi 2012 paper
    "Fumarate": {
        "setup": {
            "EX_fum_e": {
                "lower_bound": -1000,
            },
        },
        "target": {
            "reaction": "BIOMASS_Ec_iJO1366_core_53p95M",
            "value": 0.47,
        },
    },
    # Ribose as C source (aerobic) from Adadi 2012 paper
    "Ribose": {
        "setup": {
            "EX_rib__D_e": {
                "lower_bound": -1000,
            },
        },
        "target": {
            "reaction": "BIOMASS_Ec_iJO1366_core_53p95M",
            "value": 0.41,
        },
    },
    # L-Lactate as C source (aerobic) from Adadi 2012 paper
    "L-Lactate": {
        "setup": {
            "EX_lac__L_e": {
                "lower_bound": -1000,
            },
        },
        "target": {
            "reaction": "BIOMASS_Ec_iJO1366_core_53p95M",
            "value": 0.41,
        },
    },
    # Gluconate as C source (aerobic) from Adadi 2012 paper
    "Gluconate": {
        "setup": {
            "EX_glcn_e": {
                "lower_bound": -1000,
            },
        },
        "target": {
            "reaction": "BIOMASS_Ec_iJO1366_core_53p95M",
            "value": 0.68,
        },
    },
    # Sorbitol as C source (aerobic) from Adadi 2012 paper
    "Sorbitol": {
        "setup": {
            "EX_sbt__D_e": {
                "lower_bound": -1000,
            },
        },
        "target": {
            "reaction": "BIOMASS_Ec_iJO1366_core_53p95M",
            "value": 0.48,
        },
    },
    # Glucosamine as C source (aerobic) from Adadi 2012 paper
    "Glucosamine": {
        "setup": {
            "EX_gam_e": {
                "lower_bound": -1000,
            },
        },
        "target": {
            "reaction": "BIOMASS_Ec_iJO1366_core_53p95M",
            "value": 0.4,
        },
    },
    # L-Malate as C source (aerobic) from Adadi 2012 paper
    "L-Malate": {
        "setup": {
            "EX_mal__L_e": {
                "lower_bound": -1000,
            },
        },
        "target": {
            "reaction": "BIOMASS_Ec_iJO1366_core_53p95M",
            "value": 0.55,
        },
    },
    # Succinate as C source (aerobic) from Adadi 2012 paper
    "Succinate": {
        "setup": {
            "EX_succ_e": {
                "lower_bound": -1000,
            },
        },
        "target": {
            "reaction": "BIOMASS_Ec_iJO1366_core_53p95M",
            "value": 0.5,
        },
    },
    # Maltose as C source (aerobic) from Adadi 2012 paper
    "Maltose": {
        "setup": {
            "EX_malt_e": {
                "lower_bound": -1000,
            },
        },
        "target": {
            "reaction": "BIOMASS_Ec_iJO1366_core_53p95M",
            "value": 0.52,
        },
    },
    # Glucose 6-Phosphate as C source (aerobic) from Adadi 2012 paper
    "Glucose 6-Phosphate": {
        "setup": {
            "EX_g6p_e": {
                "lower_bound": -1000,
            },
        },
        "target": {
            "reaction": "BIOMASS_Ec_iJO1366_core_53p95M",
            "value": 0.78,
        },
    },
    # Mannitol as C source (aerobic) from Adadi 2012 paper
    "Mannitol": {
        "setup": {
            "EX_mnl_e": {
                "lower_bound": -1000,
            },
        },
        "target": {
            "reaction": "BIOMASS_Ec_iJO1366_core_53p95M",
            "value": 0.61,
        },
    },
    # Trehalose as C source (aerobic) from Adadi 2012 paper
    "Trehalose": {
        "setup": {
            "EX_tre_e": {
                "lower_bound": -1000,
            },
        },
        "target": {
            "reaction": "BIOMASS_Ec_iJO1366_core_53p95M",
            "value": 0.48,
        },
    },
    # Xylose as C source (aerobic) from Adadi 2012 paper
    "Xylose": {
        "setup": {
            "EX_xyl__D_e": {
                "lower_bound": -1000,
            },
        },
        "target": {
            "reaction": "BIOMASS_Ec_iJO1366_core_53p95M",
            "value": 0.51,
        },
    },
    # Mannose as C source (aerobic) from Adadi 2012 paper
    "Mannose": {
        "setup": {
            "EX_man_e": {
                "lower_bound": -1000,
            },
        },
        "target": {
            "reaction": "BIOMASS_Ec_iJO1366_core_53p95M",
            "value": 0.35,
        },
    },
    # Galactose as C source (aerobic) from Adadi 2012 paper
    "Galactose": {
        "setup": {
            "EX_gal_e": {
                "lower_bound": -1000,
            },
        },
        "target": {
            "reaction": "BIOMASS_Ec_iJO1366_core_53p95M",
            "value": 0.24,
        },
    },
}

ec_model_scenarios_for_optimization_uptake = {
    # Glucose as C source (aerobic) from Monk 2016 paper, target is growth
    "Glucose_Anaerobic_Uptake": {
        "setup": {
            "EX_glc__D_e": {
                "lower_bound": -1000,
            },
            "EX_o2_e": {
                "lower_bound": 0,
                "upper_bound": 0,
            },
            "EX_succ_e": {
                "upper_bound": 1000
            },
            "EX_for_e": {
                "upper_bound": 1000
            },
            "EX_glyclt_e": {
                "upper_bound": 1000
            },
            "EX_lac__D_e": {
                "upper_bound": 1000
            },
        },
        "target": {
            "reaction": "EX_glc__D_e",
            "value": -16.69,
        },
    },
    # Glucose as C source (aerobic) from Monk 2016 paper, target is growth
    "Glucose_Aerobic_Uptake": {
        "setup": {
            "EX_glc__D_e": {
                "lower_bound": -1000,
            },
        },
        "target": {
            "reaction": "EX_glc__D_e",
            "value": -9.53,
        },
    },
    # Glucose as C source (aerobic) from Monk 2016 paper, target is growth
    "Glucose_Aerobic_Acetate_Secretion": {
        "setup": {
            "EX_glc__D_e": {
                "lower_bound": -1000,
            },
        },
        "target": {
            "reaction": "EX_ac_e",
            "value": 3.49,
        },
    },

    "Glucose_Anaerobic_Acetate_Secretion": {
        "setup": {
            "EX_glc__D_e": {
                "lower_bound": -1000,
            },
            "EX_o2_e": {
                "lower_bound": 0,
                "upper_bound": 0,
            },
            "EX_succ_e": {
                "upper_bound": 1000
            },
            "EX_for_e": {
                "upper_bound": 1000
            },
            "EX_glyclt_e": {
                "upper_bound": 1000
            },
            "EX_lac__D_e": {
                "upper_bound": 1000
            },
        },
        "target": {
            "reaction": "EX_ac_e",
            "value": 11.71,
        },
    },
}

ec_model_scenario_glucose_uptake = {
    "Glucose_Aerobic_Uptake": {
        "setup": {
            "EX_glc__D_e": {
                "lower_bound": -1000,
            },
        },
        "target": {
            "reaction": "EX_glc__D_e",
            "value": -9.53,
        },
    },
}
