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
"""get_initial_spreadsheets.py

This module contains functions which allow to genreate 'initial spreadsheets', i.e.
spreadsheets in which the user may enter necessary information for the generation of
AutoPACMEN models.
This is typically the first script that is used in the course of a AutoPACMEN
workflow, after the wanted SBML model was selected and generated.
"""

# IMPORTS
# External modules
import cobra
import xlsxwriter
from typing import Any, Dict, List, Union
# Internal modules
from .kegg import kegg_rest_get_batch
from .helper_general import standardize_folder


# INTERNAL FUNCTIONS SECTION
def _gene_rule_as_list(gene_rule: str) -> List[Any]:
    """Returns a given string gene rule in list form.

    I.e. (b0001 or b0002) and b0003 is returned as
    [["b0001", "b0002"], "b0003"]

    Arguments:
    *gene_rule: str ~ The gene rule which shall be converted into the list form.
    """
    # Gene rules: Only ) or (, (in blocks only and); No ) and (
    gene_rule_blocks = gene_rule.split(" ) or ( ")
    gene_rule_blocks = [x.replace("(", "").replace(")", "") for x in gene_rule_blocks]
    gene_rules_array: List[Union[str, List[str]]] = []
    for block in gene_rule_blocks:
        if " or " in block:
            block_list = block.split(" or ")
            block_list = [x.lstrip().rstrip() for x in block_list]
            gene_rules_array += block_list
        elif " and " in block:
            block_list = block.split(" and ")
            block_list = [x.lstrip().rstrip() for x in block_list]
            gene_rules_array.append(block_list)
        else:  # single enzyme
            gene_rules_array.append(block)
    return gene_rules_array


# PUBLIC FUNCTIONS SECTION
def get_initial_spreadsheets(model: cobra.Model, project_folder: str, project_name: str) -> None:
    """Creates a number of initially needed XLSX spreadsheets in the given folder.

    Output
    ----------
    The following spreadsheets are going to be created (all file names start with project_name+"_"):
    * reactions.xlsx ~ A list of all KEGG IDs for each reaction of the model, the user can then select
      one of the given KEGG IDs. It is checked (using the KEGG REST API) if a reaction does not occur
      in any KEGG pathway, and such reactions are automatically marked as not included.
    * metabolites.xlsx ~ A list of all KEGG IDs for each metabolite of the model, the user can then select
      one of the given KEGG IDs. It is checked (using the KEGG REST API) if a metabolite does not occur
      in any KEGG reaction, and such metabolites are automatically marked as not included.
    * compartments.xlsx ~ A list of the model's compartments, which allows the user to set
      a pH and ionic strength value for each of them. THis data is used for the calculation of
      thermodaynamic data.
    * metabolite_concentration.xlsx ~ Allows to set default maximal and minimal metabolite
      concentrations (used in THermoFBA or OptMDFPathway in CellNetAnalyzer) for each of the
      model's metabolites.
    * protein_data.xlsx ~ Allows to set the total protein pool as well as single protein concentrations
      of the model's proteins.
    * protein_data.xlsx ~ Allows to set the total protein pool as well as single protein concentrations
      of the model's proteins.
    * enzyme_stoichiometries.xlsx ~ Allows to set the internal stoichiometries for each enzyme of a given reaction.

    Arguments
    ----------
    * model: cobra.Model ~ The cobrapy model for which the initial spreadsheets will be created.
    * project_folder: str ~ In this folder, all XLSX files will be stored.
    * project_name: str ~ This project name is added before each XLSX file name.
    """
    # Standardize project folder
    project_folder = standardize_folder(project_folder)

    # Get basepath for all spreadsheets
    basepath = project_folder + project_name

    # GET REACTION KEGG IDS AND AMBIGUITIES
    reaction_id_kegg_id_mapping: Dict[str, str] = {}
    reaction_id_reaction_name_mapping: Dict[str, str] = {}
    reaction_id_eligible_ids_mapping: Dict[str, str] = {}
    for reaction in model.reactions:
        if "kegg.reaction" not in reaction.annotation.keys():
            print(f"INFO: Reaction {reaction.id} does not have a KEGG ID annotation")
            continue
        kegg_ids = reaction.annotation["kegg.reaction"]
        if type(kegg_ids) is str:  # Single ID :O
            kegg_ids = [kegg_ids]
        else:  # Multiple IDs :O
            entries = kegg_rest_get_batch(kegg_ids, batch_size=len(kegg_ids))
            i = 0
            eligible_ids = []
            for entry in entries:
                entry_as_str = "".join(entry)
                if "PATHWAY" in entry_as_str:
                    eligible_ids.append(kegg_ids[i])
                i += 1
            if len(eligible_ids) == 1:
                reaction_id_eligible_ids_mapping[reaction.id] = eligible_ids[0]

        reaction_id_kegg_id_mapping[reaction.id] = kegg_ids
        reaction_id_reaction_name_mapping[reaction.id] = reaction.name

    # GET METABOLITE KEGG IDS AND AMBIGUITIES
    metabolite_id_kegg_id_mapping: Dict[str, str] = {}
    metabolite_id_metabolite_name_mapping: Dict[str, str] = {}
    metabolite_name_eligible_ids_mapping: Dict[str, str] = {}
    searched_metabolites: List[str] = []
    for metabolite in model.metabolites:
        if "kegg.compound" not in metabolite.annotation.keys():
            print(f"INFO: Metabolite {metabolite.id} does not have a KEGG ID annotation")
            continue
        kegg_ids = metabolite.annotation["kegg.compound"]
        if type(kegg_ids) is str:  # Single ID :O
            kegg_ids = [kegg_ids]
        else:  # Multiple IDs :O
            if metabolite.name not in searched_metabolites:
                entries = kegg_rest_get_batch(kegg_ids, batch_size=len(kegg_ids))
                i = 0
                eligible_ids = []
                for entry in entries:
                    entry_as_str = "".join(entry)
                    if "REACTION" in entry_as_str:
                        eligible_ids.append(kegg_ids[i])
                    i += 1
                if len(eligible_ids) == 1:
                    metabolite_name_eligible_ids_mapping[metabolite.name] = eligible_ids[0]
                searched_metabolites.append(metabolite.name)
        metabolite_id_kegg_id_mapping[metabolite.id] = kegg_ids
        metabolite_id_metabolite_name_mapping[metabolite.id] = metabolite.name

    # Reactions <-> KEGG ID mapping XLSX :D
    workbook = xlsxwriter.Workbook(basepath+"_reactions.xlsx")
    worksheet = workbook.add_worksheet("Reaction IDs")

    yellow = workbook.add_format()
    yellow.set_bg_color("#FFFF00")
    blue = workbook.add_format()
    blue.set_bg_color("#FF00FF")

    row = 0
    for key in reaction_id_reaction_name_mapping.keys():
        worksheet.write(row, 0, key)
        worksheet.write(row, 1, reaction_id_reaction_name_mapping[key])
        column = 2
        kegg_ids = reaction_id_kegg_id_mapping[key]
        for kegg_id in kegg_ids:
            if len(kegg_ids) == 1:
                worksheet.write_url(row, column, "https://www.kegg.jp/dbget-bin/www_bget?rn:"+kegg_id)
                worksheet.write(row, column+1, kegg_id)
                worksheet.write(row, column+2, "Yes")
            else:
                if key in reaction_id_eligible_ids_mapping.keys():
                    if kegg_id == reaction_id_eligible_ids_mapping[key]:
                        worksheet.write_url(row, column, "https://www.kegg.jp/dbget-bin/www_bget?rn:"+kegg_id, blue)
                        worksheet.write(row, column+1, kegg_id, blue)
                        worksheet.write(row, column+2, "Yes")
                    else:
                        worksheet.write_url(row, column, "https://www.kegg.jp/dbget-bin/www_bget?rn:"+kegg_id, blue)
                        worksheet.write(row, column+1, kegg_id, blue)
                else:
                    worksheet.write_url(row, column, "https://www.kegg.jp/dbget-bin/www_bget?rn:"+kegg_id, yellow)
                    worksheet.write(row, column+1, kegg_id, yellow)
            column += 3
        row += 1
    workbook.close()

    # Metabolites<->KEGG ID mapping XLSX :D
    workbook = xlsxwriter.Workbook(basepath+"_metabolites.xlsx")
    worksheet = workbook.add_worksheet("Metabolite IDs")

    yellow = workbook.add_format()
    yellow.set_bg_color("#FFFF00")
    blue = workbook.add_format()
    blue.set_bg_color("#FF00FF")

    row = 0
    for key in metabolite_id_metabolite_name_mapping.keys():
        metabolite_name = metabolite_id_metabolite_name_mapping[key]
        worksheet.write(row, 0, key)
        worksheet.write(row, 1, metabolite_name)
        column = 2
        kegg_ids = metabolite_id_kegg_id_mapping[key]
        for kegg_id in kegg_ids:
            if len(kegg_ids) == 1:
                worksheet.write_url(row, column, "https://www.kegg.jp/dbget-bin/www_bget?cpd:"+kegg_id)
                worksheet.write(row, column+1, kegg_id)
                worksheet.write(row, column+2, "Yes")
            else:
                if metabolite_name in metabolite_name_eligible_ids_mapping.keys():
                    if kegg_id == metabolite_name_eligible_ids_mapping[metabolite_name]:
                        worksheet.write_url(row, column, "https://www.kegg.jp/dbget-bin/www_bget?cpd:"+kegg_id, blue)
                        worksheet.write(row, column+1, kegg_id, blue)
                        worksheet.write(row, column+2, "Yes")
                    else:
                        worksheet.write_url(row, column, "https://www.kegg.jp/dbget-bin/www_bget?cpd:"+kegg_id, blue)
                        worksheet.write(row, column+1, kegg_id, blue)
                else:
                    worksheet.write_url(row, column, "https://www.kegg.jp/dbget-bin/www_bget?cpd:"+kegg_id, yellow)
                    worksheet.write(row, column+1, kegg_id, yellow)
            column += 3
        row += 1
    workbook.close()

    # Compartment data XLSX :D
    workbook = xlsxwriter.Workbook(basepath+"_compartments.xlsx")
    worksheet = workbook.add_worksheet("Compartment pH")

    row = 0
    for compartment_key in model.compartments.keys():
        compartment_name = model.compartments[compartment_key]
        worksheet.write(row, 0, compartment_key)
        worksheet.write(row, 1, compartment_name)
        worksheet.write(row, 2, "NA")
        worksheet.write(row, 3, "NA")
        row += 1
    workbook.close()

    # Protein data XLSX :D
    print("NOTE: "+project_name+"_protein_data.xlsx has as default value for the enzyme pool P 0.095 mmol/gDW.")
    print("Please adjust the value accordingly for your model!")
    workbook = xlsxwriter.Workbook(basepath+"_protein_data.xlsx")
    worksheet = workbook.add_worksheet("Total protein data")
    worksheet.write(0, 0, "Total protein content [g/gDW]:")
    worksheet.write(0, 1, .095)
    worksheet.write(1, 0, "Fraction of masses of model-included enzymes in comparison to all enzymes (0.0 to 1.0):")
    worksheet.write(1, 1, 1.0)
    worksheet.write(2, 0, "Average saturation level (0.0 to 1.0):")
    worksheet.write(2, 1, 1.0)
    worksheet2 = workbook.add_worksheet("Single protein data")
    worksheet2.write(0, 0, "Protein ID (as in SBML model)")
    worksheet2.write(0, 1, "Protein concentration [mmol/gDW]")
    workbook.close()

    # Metabolite concentrations XLSX :D
    workbook = xlsxwriter.Workbook(basepath+"_metabolite_concentrations.xlsx")
    worksheet = workbook.add_worksheet("Default data")
    worksheet.write(0, 0, "Default minimal metabolite concentration [M]:")
    worksheet.write(1, 0, "Default maximal metabolite concentration [M]:")
    worksheet2 = workbook.add_worksheet("Single metabolite data")
    worksheet2.write(0, 0, "Metabolite ID (as in SBML model)")
    worksheet2.write(0, 1, "Minimal metabolite concentration [M]")
    worksheet2.write(0, 2, "Maximal metabolite concentration [M]")
    workbook.close()

    # Enzyme stoichiometry XLSX :D
    # Get gene rule <-> Reaction ID mapping
    reaction_id_gene_rules_mapping = {}
    for reaction in model.reactions:
        listed_gene_rules = _gene_rule_as_list(reaction.gene_reaction_rule)
        if listed_gene_rules != ['']:
            reaction_id_gene_rules_mapping[reaction.id] = listed_gene_rules

    # Write XLSX
    workbook = xlsxwriter.Workbook(basepath+"_enzyme_stoichiometries.xlsx")
    # Gene stoichiometry worksheets
    worksheet_stoichiometry = workbook.add_worksheet("Stoichiometries of complexes")
    line = 0
    for reaction_id in reaction_id_gene_rules_mapping.keys():
        gene_rule = reaction_id_gene_rules_mapping[reaction_id]
        if gene_rule == [""]:
            continue
        worksheet_stoichiometry.write(line, 0, reaction_id)
        row = 1
        for or_part in gene_rule:
            worksheet_stoichiometry.write(line, row, str(or_part))
            if type(or_part) is str:
                default_stoichiometry = "1"
            else:
                default_stoichiometry = ";".join(["1" for _ in range(len(or_part))])
            worksheet_stoichiometry.write(line, row+1, default_stoichiometry)
            row += 2
        line += 1

    workbook.close()


def get_initial_spreadsheets_with_sbml(sbml_path: str, project_folder: str, project_name: str) -> None:
    """Allows to call get_initial_spreadsheets with an SBML.

    Arguments
    ----------
    * sbml_path: str ~ The SBML of the analyzed model.
    * project_folder: str ~ In this folder, all XLSX files will be stored.
    * project_name: str ~ This project name is added before each XLSX file name.
    """
    # LOAD SBML MODEL
    model: cobra.Model = cobra.io.read_sbml_model(sbml_path)
    get_initial_spreadsheets(model, project_folder, project_name)
