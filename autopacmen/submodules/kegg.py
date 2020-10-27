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
"""kegg.py

This module contains functions which acan access the multi-omics
KEGG database.
"""

# IMPORTS
# External modules
import copy
import time
from Bio.KEGG import REST
from typing import List


def get_full_organism_name_from_kegg_id(organism_kegg_id):
    """Get the given organism's full name using KEGG info. Returned as str.

    Arguments
    ----------
    * organism_kegg_id: str ~ The organism's full name, e.g. "Xanthomonas
      campestris pv. campesris B100"
    """
    # Use the KEGG API's "info" function
    kegg_organism_info = kegg_rest("info", organism_kegg_id)
    # The name is after many spaces and without unneccesary spaces at the beginning and the end
    organism_full_name = kegg_organism_info[0].split(
        "    ")[-1].replace("KEGG Genes Database", "").lstrip().rstrip()
    # Return the retrieved name :D
    return organism_full_name


def kegg_rest(type: str, argument: str, optional_argument: str = "", sleep_time: float = .5) -> List[str]:
    """This function calls Biopython's KEGG REST function and returns the lines as a string list.

    All empty lines are deleted from the list as they do not contain any information.

    Arguments
    ----------
    * type: str ~ The KEGG REST action. Can be either 'info', 'get', 'link' or 'list.
    * argument: str ~ The argument for the KEGG order.
    * optional_argument: str="" ~ The second argument which is necessary for 'link' and 'list' actions
      to work correctly.
    * sleep_time: float=10.0 ~ The time that shall be waited after a REST action is performed.
      Its default value of 10.0 seconds is in accordance with the NCBI
      rule that its servers shall not be contacted more often than every
      10 seconds. KEGG might have lower or higher required sleep times,
      but I did not find any specified time step.
    """
    # Execute correct Biotpython KEGG REST function.
    if type == "info":
        kegg_data = REST.kegg_info(argument)
    elif type == "get":
        kegg_data = REST.kegg_get(argument)
    elif type == "link":
        kegg_data = REST.kegg_link(argument, optional_argument)
    elif type == "list":
        kegg_data = REST.kegg_list(argument, optional_argument)
    elif type == "find":
        kegg_data = REST.kegg_find(argument, optional_argument)

    # Wait the sleep time doing nothing.
    time.sleep(sleep_time)

    # Get one string per line of the KEGG REST result.
    lines: List[str] = kegg_data.read().split("\n")

    # Delete empty lines.
    not_empty_lines: List[str] = [i for i in lines if len(i) > 0]

    return not_empty_lines


def kegg_rest_get_batch(input_ids: List[str], batch_size: int = 4) -> List[List[str]]:
    """Using Biopython's KEGG REST function, this function returns a list of lists of KEGG REST get results.

    See this module's kegg_rest() comments for more.

    Arguments
    ----------
    * input_ids: List[str] ~ The list of searched KEGG GET IDs.
    * batch_size: int = 4 ~ The size of batches, i.e. the maximal number of looked up IDs (batching
      is done in order to lower the number of necessary KEGG API calls).
    """

    batched_input_ids: List[str] = []
    current_batch_start = 0
    while current_batch_start < len(input_ids):
        combined_input_ids = "+".join(
            input_ids[current_batch_start:current_batch_start + batch_size])
        batched_input_ids.append(combined_input_ids)
        current_batch_start += batch_size

    output: List[List[str]] = []
    for batch in batched_input_ids:
        kegg_output = kegg_rest("get", batch)

        single_id_data: List[str] = []
        for line in kegg_output:
            if "///" in line:
                output.append(copy.deepcopy(single_id_data))
                single_id_data = []
            else:
                single_id_data.append(line)

    return output
