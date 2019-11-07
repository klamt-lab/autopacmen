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
"""ec_model_2019_06_25_data_read_supplementary_table_of_keff_paper.py

Converts the data from the keff paper into a more machine-readable JSON.
"""

import copy
import cobra
import openpyxl
from typing import Any, Dict, List

from submodules.helper_general import json_write

# Get kcats
workbook = openpyxl.load_workbook("ec_model_2019_06_25_input_keff_paper/c3mb70119k-1.xlsx", read_only=True)
worksheet = workbook["Table S2"]

gene_id_data_mapping: Dict[str, Dict[str, Any]] = {}
current_line = 1
for line in worksheet.rows:
    current_cell = 1

    if current_line < 13:
        current_line += 1
        continue

    current_gene_kcats: List[float] = []
    for cell in line:
        if current_cell == 1:
            current_gene_id = cell.value
        elif current_cell == 2:
            current_gene_names = cell.value
        elif current_cell == 3:
            current_gene_description = cell.value
        elif current_cell in [x for x in range(49, 54)]:
            current_gene_kcats.append(cell.value)
        current_cell += 1

    if None in current_gene_kcats:
        continue

    gene_id_data: Dict[str, Any] = {}
    gene_id_data["names"] = current_gene_names
    gene_id_data["description"] = current_gene_description
    gene_id_data["kcats"] = current_gene_kcats

    gene_id_data_mapping[current_gene_id] = copy.deepcopy(gene_id_data)


# Get directionality data
model = cobra.io.read_sbml_model("ec_model_2019_06_25_input/iJO1366.xml")
pfba_solution = cobra.flux_analysis.pfba(model)
for reaction in model.reactions:
    gene_reaction_rule = reaction.gene_reaction_rule
    if gene_reaction_rule == "":
        continue

    gene_reaction_rule = gene_reaction_rule.replace(" or ", "\t")
    gene_reaction_rule = gene_reaction_rule.replace(" and ", "\t")
    gene_names = gene_reaction_rule.split("\t")

    has_negative_flux = pfba_solution.fluxes[reaction.id] < 0
    if has_negative_flux:
        direction = "reverse"
    else:
        direction = "forward"
    for gene_name in gene_names:
        if gene_name not in gene_id_data_mapping.keys():
            continue

        if "direction" not in gene_id_data_mapping[gene_name].keys():
            gene_id_data_mapping[gene_name]["direction"] = {}

        gene_id_data_mapping[gene_name]["direction"][reaction.id] = direction

# Write JSON :D
json_write("ec_model_2019_06_25_input_keff_paper/gene_id_data_mapping.json", gene_id_data_mapping)
