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
"""parse_bigg_metabolites_file.py

This module contains a function which transforms a BIGG metabolites .txt list
into an AutoPACMEN-readable JSON.
"""

# IMPORTS
# Internal modules
from .helper_general import json_write, standardize_folder


# PUBLIC FUNCTIONS SECTIONS
def parse_bigg_metabolites_file(bigg_metabolites_file_path: str, json_output_folder: str) -> None:
    """Parses a BIGG metabolites text file and returns a dictionary for this file.

    As of 29/04/2019, a BIGG metabolites list of all BIGG-included metabolites
    is retrievable under http://bigg.ucsd.edu/data_access

    Arguments
    ----------
    * bigg_metabolites_file_path: str ~ The file path to the BIGG metabolites file.
      The usual file name (which has to be included too in this argument) is
      bigg_models_metabolites.txt
    * output_folder: str ~ The folder in which the JSON including the parsed BIGG
      metabolites file data is stored with the name 'bigg_id_name_mapping.json'

    Output
    ----------
    * A JSON file with the name 'bigg_id_name_mapping.json' in the given output folder,
      with the following structure:
    <pre>
     {
         "$BIGG_ID": "$CHEMICAL_OR_USUAL_NAME",
         (...),
         "$BIGG_ID": "$BIGG_ID",
         (...),
     }
    </pre>
    The BIGG ID <-> BIGG ID mapping is done for models which already use the BIGG IDs.
    """
    # Standardize output folder
    json_output_folder = standardize_folder(json_output_folder)

    # Open the BIGG metabolites file as string list, and remove all newlines
    with open(bigg_metabolites_file_path, "r") as f:
        lines = f.readlines()
    lines = [x.replace("\n", "") for x in lines if len(x) > 0]

    # Mapping variable which will store the BIGG ID<->
    bigg_id_name_mapping = {}
    # Go through each BIGG metabolites file line (which is a tab-separated file)
    # and retrieve the BIGG ID and the name (if there is a name for the given BIGG
    # ID)
    for line in lines:
        bigg_id = line.split("\t")[1]
        # Exception to check if there is no name :O
        try:
            name = line.split("\t")[2].lower()
        except Exception:
            continue

        bigg_id_name_mapping[name] = bigg_id
        bigg_id_name_mapping[bigg_id] = bigg_id

    # Write the JSON in the given folder :D
    json_write(json_output_folder+"bigg_id_name_mapping.json",
               bigg_id_name_mapping)
