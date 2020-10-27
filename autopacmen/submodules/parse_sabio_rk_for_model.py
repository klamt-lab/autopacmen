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
"""parse_sabio_rk_for_model.py

The function of this module returns SABIO-RK kcat entries for each reaction in the
given metabolic model.
"""

# IMPORTS
# External modules
import cobra
from typing import List
# Internal modules
from .parse_sabio_rk import get_ec_number_kcats_wildcard_search
from .helper_general import json_write


# PUBLIC FUNCTIONS SECTION
def parse_sabio_rk_for_model(model: cobra.Model, json_output_path: str, bigg_id_name_mapping_path: str) -> None:
    """Retrieves kcats from SABIO-RK for the given model and stores it in a JSON for the given model in the given path.

    Algorithm
    ----------
    Using the SABIO-RK REST API (as of 2019/30/04, it is explained under
    http://sabiork.h-its.org/layouts/content/docuRESTfulWeb/RESTWebserviceIntro.gsp),


    Arguments
    ----------
    * model: cobra.Model ~ The model for which kcats shall be retrieved from SABIO-RK.
    * json_output_path: str ~ The path of the JSON that shall be created

    Output
    ----------
    * A JSON in the given project folder with the following structure:
    <pre>
        {
            "$EC_NUMBER_OR_KEGG_REACTION_ID": {
                "$SUBSTRATE_WITH_BIGG_ID_1": {
                    "$ORGANISM_1": [
                        $kcat_1,
                        (...)
                        $kcat_n,
                    ]
                },
                (...),
                "REST": {
                    "$ORGANISM_1": [
                        $kcat_1,
                        (...)
                        $kcat_n,
                    ]
                }
            }
            (...),
        }
    </pre>
    'REST' stands for a substrate without found BIGG ID.
    """
    # GET LIST OF EC NUMBERS
    ec_numbers_list: List[str] = []
    for reaction in model.reactions:
        if "ec-code" not in reaction.annotation.keys():
            continue
        ec_codes = reaction.annotation["ec-code"]
        if type(ec_codes) is str:
            ec_codes = [ec_codes]
        ec_numbers_list += ec_codes
    ec_numbers_list = list(set(ec_numbers_list))

    # GET KCATS FOR EC NUMBERS
    ec_number_kcat_mapping = get_ec_number_kcats_wildcard_search(
        ec_numbers_list, bigg_id_name_mapping_path)

    json_write(json_output_path, ec_number_kcat_mapping)


def parse_sabio_rk_for_model_with_sbml(sbml_path: str, json_output_path: str, bigg_id_name_mapping_path: str) -> None:
    """See this module's parse_sabio_rk_for_model() documentation. This function uses an SBML path.

    Arguments
    ----------
    * sbml_path: str ~ The model's SBML path.
    * json_output_path: str ~ The path of the JSON that shall be created
    """
    # LOAD SBML MODEL
    model: cobra.Model = cobra.io.read_sbml_model(sbml_path)
    parse_sabio_rk_for_model(model, json_output_path,
                             bigg_id_name_mapping_path)
