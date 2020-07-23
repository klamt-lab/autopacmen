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
"""This script consists of all steps needed for the genration of iJO1366-sMOMENT.

The single steps are explained in the supplementary file 1 of sMOMENT's publication (Bekiaris & Klamt, 2020).
"""

# IMPORTS
# External modules
import cobra
# Internal modules
from ec_model_2019_06_25_data_set_up_model import set_up_ec_model_with_sbml
from submodules.create_combined_kcat_database import create_combined_kcat_database
from submodules.create_smoment_model_reaction_wise import create_smoment_model_reaction_wise_with_sbml
from submodules.get_initial_spreadsheets import get_initial_spreadsheets_with_sbml
from submodules.get_protein_mass_mapping import get_protein_mass_mapping_with_sbml
from submodules.get_reactions_kcat_mapping import get_reactions_kcat_mapping
from submodules.parse_bigg_metabolites_file import parse_bigg_metabolites_file
from submodules.parse_brenda_textfile import parse_brenda_textfile
from submodules.parse_brenda_json_for_model import parse_brenda_json_for_model
from submodules.parse_sabio_rk_for_model import parse_sabio_rk_for_model_with_sbml


# Step 1
bigg_metabolites_file_path = "ec_model_2019_06_25_input/bigg_models_metabolites.txt"
json_output_folder = "./ec_model_2019_06_25_output/"
parse_bigg_metabolites_file(bigg_metabolites_file_path, json_output_folder)

# Step 2
brenda_textfile_path = "./ec_model_2019_06_25_input/brenda_download.txt"
bigg_metabolites_json_folder = "./ec_model_2019_06_25_output/"
json_output_path = "./ec_model_2019_06_25_output/kcat_database_brenda.json"
parse_brenda_textfile(brenda_textfile_path, bigg_metabolites_json_folder, json_output_path)

# Step 3
sbml_path = "./ec_model_2019_06_25_input/iJO1366.xml"
brenda_json_path = "./ec_model_2019_06_25_output/kcat_database_brenda.json"
output_json_path = "./ec_model_2019_06_25_output/kcat_database_brenda_for_model.json"
parse_brenda_json_for_model(sbml_path, brenda_json_path, output_json_path)

# Step 4
sbml_path = "./ec_model_2019_06_25_input/iJO1366.xml"
json_output_path = "./ec_model_2019_06_25_output/kcat_database_sabio_rk.json"
bigg_id_name_mapping_path: str = "./ec_model_2019_06_25_output/bigg_id_name_mapping.json"
parse_sabio_rk_for_model_with_sbml(sbml_path, json_output_path, bigg_id_name_mapping_path)

# Step 5
sabio_rk_kcat_database_path = "./ec_model_2019_06_25_output/kcat_database_sabio_rk.json"
brenda_kcat_database_path = "./ec_model_2019_06_25_output/kcat_database_brenda_for_model.json"
output_path = "./ec_model_2019_06_25_output/kcat_database_combined.json"
create_combined_kcat_database(sabio_rk_kcat_database_path, brenda_kcat_database_path, output_path)

# Step 6
sbml_path = "ec_model_2019_06_25_input/iJO1366.xml"
project_folder = "./ec_model_2019_06_25_output/"
project_name = "psb_orth"
organism = "Escherichia coli"
kcat_database_path = "./ec_model_2019_06_25_output/kcat_database_combined.json"
protein_kcat_database_path = "ec_model_2019_06_25_input_keff_paper/gene_id_data_mapping.json"
get_reactions_kcat_mapping(sbml_path, project_folder, project_name, organism, kcat_database_path, protein_kcat_database_path)

# Step 7
input_sbml = "./ec_model_2019_06_25_input/iJO1366.xml"
project_folder = "./ec_model_2019_06_25_output/"
project_name = "psb_orth"
get_initial_spreadsheets_with_sbml(input_sbml, project_folder, project_name)

# Step 8
input_sbml = "./ec_model_2019_06_25_input/iJO1366.xml"
project_folder = "./ec_model_2019_06_25_output/"
project_name = "psb_orth"
get_protein_mass_mapping_with_sbml(input_sbml, project_folder, project_name)

# Step 9
input_sbml = "./ec_model_2019_06_25_input/iJO1366.xml"
output_sbml = "iJO1366_sMOMENT_2019_06_25_GECKO_ANALOGON.xml"
project_folder = "./ec_model_2019_06_25_output/"
project_name = "psb_orth"
excluded_reactions = ["CO2tex", "O2tex", "H2tex"]
create_smoment_model_reaction_wise_with_sbml(input_sbml, output_sbml, project_folder, project_name, excluded_reactions)

# Step 10 - Generate model with standard exchanges scenario
model = set_up_ec_model_with_sbml("ec_model_2019_06_25_output/iJO1366_sMOMENT_2019_06_25.xml", .25)
cobra.io.write_sbml_model(model, "ec_model_2019_06_25_output/iJO1366_sMOMENT_2019_06_25_STANDARD_EXCHANGE_SCENARIO.xml")
