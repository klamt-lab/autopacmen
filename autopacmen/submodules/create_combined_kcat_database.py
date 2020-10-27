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
"""create_combined_kcat_database.py

This module contains functions which return a combined kcat database from
SABIO-RK and BRENDA data.
"""

# IMPORTS
# Internal modules
from .helper_general import json_load, json_write
# External modules
from typing import Any, Dict, List


# PUBLIC FUNCTIONS
def create_combined_kcat_database(sabio_rk_kcat_database_path: str, brenda_kcat_database_path: str, output_path: str) -> None:
    """Creates a combined JSON of the given SABIO-K and BRENDA kcat databases with non-wildcard entries only.

    Arguments
    ----------
    * sabio_rk_kcat_database_path: str ~ The path to the SABIO-RK kcat database JSON
    * brenda_kcat_database_path: str ~ The path to the BRENDA kcat database JSON
    * output_path: str ~ The outputh path (with filename) of the genreated combined kcat database JSON

    Output:
    A JSON with the following format:
    <pre>
    {
        '$EC_NUMBER': {
            '$BIGG_IDS_OF_SUBSTRATES': {
                '$ORGANISM': {
                    kcat: float
                },
                (...)
            },
            (...),
            'SOURCE': 'SABIO_RK' or 'BRENDA' or 'BRENDA and SABIO-RK',
            'WILDCARD': false
        },
        (...)
    }
    </pre>
    """
    # Load the two given databases as JSONs
    sabio_rk_database = json_load(sabio_rk_kcat_database_path)
    brenda_database = json_load(brenda_kcat_database_path)

    # Get all EC number keys (BRENDA contains all relevant EC numbers)
    ec_number_keys: List[str] = list(brenda_database.keys())
    # Set-up combined kcat database dictionary
    combined_database: Dict[str, Dict[str, Any]] = {}
    # Go through each EC number :D...
    for ec_number_key in ec_number_keys:
        # Get the wildcard status (i.e., found with a * wildcard?) and check if the EC number occurs anywhere...
        # ...for SABIO-RK
        if ec_number_key not in sabio_rk_database.keys():
            print(f"WARNING: EC number {ec_number_key} could not be found in SABIO-RK, even with wildcards")
            print("Possible reasons: The EC number format is invalid or there was an SABIO-RK API error")
            is_sabio_rk_from_wildcard = True  # Let the combined database ignore this entry
        else:
            is_sabio_rk_from_wildcard = sabio_rk_database[ec_number_key]["WILDCARD"]
        # ...and for BRENDA
        if ec_number_key not in brenda_database.keys():
            print(f"WARNING: EC number {ec_number_key} could not be found in SABIO-RK, even with wildcards")
            print("Possible reason: The EC number format is invalid")
            is_brenda_from_wildcard = True  # Let the combined database ignore this entry
        else:
            is_brenda_from_wildcard = brenda_database[ec_number_key]["WILDCARD"]

        # If both are from wildcards, ignore them :3
        if (is_sabio_rk_from_wildcard) and (is_brenda_from_wildcard):
            continue

        # Set-up dictionary for the EC number since at least one of the two databases
        # is not from a wildcarded search :D
        combined_database[ec_number_key] = {}
        # If both are not from wildcards, combine them :D...
        if (not is_sabio_rk_from_wildcard) and (not is_brenda_from_wildcard):
            # ...by reading their metabolites...
            sabio_rk_metabolite_keys = list(sabio_rk_database[ec_number_key].keys())
            brenda_metabolite_keys = list(brenda_database[ec_number_key].keys())
            metabolite_keys = list(set(sabio_rk_metabolite_keys + brenda_metabolite_keys))
            # ...going through them...
            for metabolite_key in metabolite_keys:
                # ...excluding the WILDCARD key...
                if metabolite_key == "WILDCARD":
                    continue
                # ...and adding the metabolites according to their presence in the databases :D
                is_metabolite_in_brenda: bool = metabolite_key in brenda_metabolite_keys
                is_metabolite_in_sabio_rk: bool = metabolite_key in sabio_rk_metabolite_keys
                if is_metabolite_in_brenda and is_metabolite_in_sabio_rk:
                    sabio_rk_entry = sabio_rk_database[ec_number_key][metabolite_key]
                    brenda_entry = brenda_database[ec_number_key][metabolite_key]
                    combined_database[ec_number_key][metabolite_key] = {**sabio_rk_entry, **brenda_entry}
                elif is_metabolite_in_brenda:
                    brenda_entry = brenda_database[ec_number_key][metabolite_key]
                    combined_database[ec_number_key][metabolite_key] = brenda_entry
                else:
                    sabio_rk_entry = sabio_rk_database[ec_number_key][metabolite_key]
                    combined_database[ec_number_key][metabolite_key] = sabio_rk_entry
            combined_database[ec_number_key]["WILDCARD"] = is_sabio_rk_from_wildcard
            combined_database[ec_number_key]["SOURCE"] = "BRENDA and SABIO-RK"
        # If only the SABIO-RK entry does not come from a wildcard, use it :D
        elif not is_sabio_rk_from_wildcard:
            combined_database[ec_number_key] = sabio_rk_database[ec_number_key]
            combined_database[ec_number_key]["WILDCARD"] = False
            combined_database[ec_number_key]["SOURCE"] = "SABIO-RK"
        # If only the BRENDA entry does not come from a wildcard, use it :-)
        elif not is_brenda_from_wildcard:
            combined_database[ec_number_key] = brenda_database[ec_number_key]
            combined_database[ec_number_key]["WILDCARD"] = False
            combined_database[ec_number_key]["SOURCE"] = "BRENDA"
    json_write(output_path, combined_database)
