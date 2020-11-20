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
"""parse_brenda_for_model.py

Contains functions which allow to vreate a model-specific and BRENDA-depending
kcat database.
"""

# IMPORTS
# External modules
import cobra
import copy
# Internal modules
from .helper_general import is_fitting_ec_numbers, json_load, json_write
from typing import Any, Dict, List


# PRIVATE FUNCTIONS
def _get_transfer_ec_number_entry(ec_number_entry_key: str, brenda_kcat_database_original: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    """Returns the new EC number to which the given EC number was transferred.

    This is indicated in the given dictionary by the 'TRANSFER' key.
    Especially since the EC class 7's (translocases) introduction in 2018, many EC numbers are being transferred to new ones.

    Arguments
    ----------
    *ec_number_entry_key: str ~ The EC number for which its newly assigned one shall be searched.
    *brenda_kcat_database_original: Dict[str, Dict[str, Any]] ~ The BRENDA database dictionary with TRANSFER data.

    Output
    ----------
    Either...
    * {'ERROR': None} if the transferred new EC number is transferred to the old one (a futile cycle :O), or if the new EC number
      is not in BRENDA's database :O, or if the new EC number does not have any kcat entries :O
    * the dictionary containing the substrates, organisms and associated kcats of the new EC number
    """
    ec_number_entry = brenda_kcat_database_original[ec_number_entry_key]
    while "TRANSFER" in ec_number_entry.keys():
        new_ec_number = ec_number_entry["TRANSFER"]
        if new_ec_number == ec_number_entry_key:
            return {"ERROR": None}
        if new_ec_number not in brenda_kcat_database_original.keys():
            return {"ERROR": None}
        ec_number_entry = brenda_kcat_database_original[new_ec_number]

    if ec_number_entry == {}:
        return {"ERROR": None}

    return copy.deepcopy(ec_number_entry)


# PUBLIC FUNCTIONS
def parse_brenda_json_for_model(sbml_path: str, brenda_json_path: str, output_json_path: str) -> None:
    """Reads out a BRENDA JSON file created with parse_brenda_textfile and creates a model-specific JSON.

    Arguments
    ----------
    * sbml_path: str ~ The path of the SBML model of which a specific BRENDA JSON kcat database
      shall be created
    * brenda_json_path: str ~ The full path to the BRENDA JSON created with parse_brenda_textfile.
    * output_json_path: str ~ The full path to the newly created JSON.

    Output
    ----------
    A JSON in the given folder and the name 'kcat_database_brenda.json', and with the following structure:
    <pre>
    {
        '$EC_NUMBER': {
            '$BIGG_ID_METABOLITE': {
                '$ORGANISM': [
                    kcat_list: float
                ],
                (...)
            },
            (...)
        },
        (...)
    }
    </pre>
    """
    model: cobra.Model = cobra.io.read_sbml_model(sbml_path)

    # Get EC numbers of the model's reactions
    ec_numbers_of_model: List[str] = []
    for reaction in model.reactions:
        if "ec-code" not in reaction.annotation.keys():
            continue

        ec_numbers_of_reaction = reaction.annotation["ec-code"]
        if type(ec_numbers_of_reaction) is str:
            ec_numbers_of_reaction = [ec_numbers_of_reaction]
        ec_numbers_of_model += ec_numbers_of_reaction
    ec_numbers_of_model = list(set(ec_numbers_of_model))

    # Get EC number entries for each EC number of the model
    brenda_kcat_database_original = json_load(brenda_json_path)
    brenda_kcat_database_for_model = {}
    for ec_number in ec_numbers_of_model:
        entry_error = False
        if ec_number in brenda_kcat_database_original.keys():
            ec_number_entry = _get_transfer_ec_number_entry(
                ec_number, brenda_kcat_database_original)
            if "ERROR" in ec_number_entry.keys():
                entry_error = True
            else:
                ec_number_entry["WILDCARD"] = False
                brenda_kcat_database_for_model[ec_number] = ec_number_entry

        if (ec_number not in brenda_kcat_database_original.keys()) or entry_error:
            eligible_ec_number_entries: List[Dict[str, Any]] = []
            for wildcard_level in range(1, 5):
                for database_ec_number in list(brenda_kcat_database_original.keys()):
                    if is_fitting_ec_numbers(ec_number, database_ec_number, wildcard_level):
                        database_ec_number_entry = _get_transfer_ec_number_entry(
                            database_ec_number, brenda_kcat_database_original)
                        if "ERROR" not in database_ec_number_entry.keys():
                            eligible_ec_number_entries.append(
                                database_ec_number_entry)
                if len(eligible_ec_number_entries) > 0:
                    break
            ec_number_entry = {}
            for eligible_ec_number_entry in eligible_ec_number_entries:
                for metabolite_key in eligible_ec_number_entry.keys():
                    metabolite_entry = eligible_ec_number_entry[metabolite_key]
                    if metabolite_key not in ec_number_entry.keys():
                        ec_number_entry[metabolite_key] = metabolite_entry
                    else:
                        ec_number_entry[metabolite_key] = {
                            **ec_number_entry[metabolite_key], **metabolite_entry}
            ec_number_entry["WILDCARD"] = True
            brenda_kcat_database_for_model[ec_number] = ec_number_entry

    json_write(output_json_path, brenda_kcat_database_for_model)
