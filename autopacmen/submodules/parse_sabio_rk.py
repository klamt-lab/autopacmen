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
"""sabio_rk.py

This module contains functions which acan access the SABIO-RK kinetic data
platform.
"""

# IMPORTS
# External libraries
import copy
import csv
import io
import requests
import time
from typing import Any, Dict, List
# Internal libraries
from .helper_general import ensure_folder_existence, get_files, is_fitting_ec_numbers, json_load, json_write


# SCRIPT-WIDE CONSTANTS
# URL for SABIO-RK's kcat REST API
QUERY_URL = "http://sabiork.h-its.org/sabioRestWebServices/kineticlawsExportTsv"
# Time in seconds to wait between SABIO-RK API calls
WAIT_TIME = 1.5
# Constant unit multipliers since SABIO-RK contains kcats in different units :O
# The multiplication numbers show with which number a kcat of the corresponding unit
# has to be multiplied in order to get the kcat in s^(-1), the kcat unit that is
# used all over AutoPACMEN.
UNIT_MULTIPLIER: Dict[str, float] = {
    "s^(-1)": 1.0,
    "min^(-1)": 1/60,
    "h^(-1)": 1/(60*60),
}


# PRIVATE FUNCTIONS
def _add_wildcard_to_ec_number(ec_number: str, level: int) -> str:
    """Adds asterisk wildcards to the given EC number.

    Input
    ----------
    * ec_number: str ~ The EC number which shall be wildcarded
    * level: int ~ The wildcard level. It has to be at least 0 and
      maximally 5


    Example
    ----------
    <pre>
    >>> _add_wildcard_to_ec_number("1.1.1.1", 1)
    1.1.1.*
    </pre>

    Output
    ----------
    * The wildcarded EC number as string
    """
    ec_number_list = ec_number.split(".")[::-1]
    i = 0
    while i <= level:
        if i < level:
            ec_number_list[i] = "*"
        i += 1
    wildcarded_ec_number = ".".join(ec_number_list[::-1])
    return wildcarded_ec_number


def _extract_kcat_lines(result: csv.DictReader) -> List[str]:
    """Converts a SABIO-RK csv.DictReader line into a list of strings with the content of these lines.

    Input
    ----------
    * result: csv.DictReader ~ The csv.DictReader instance of a SABIO-RK kcat saerch output line

    Output
    ----------
    A list of all kcat lines as strings
    """
    kcat_lines: List[Any] = []
    for row in result:
        if (row["parameter.type"] == "kcat") and (row["parameter.startValue"] != ""):
            kcat_lines.append(row)
    return kcat_lines


def _get_species_results(result: csv.DictReader) -> List[str]:
    """Returns the organism data from the given SABIO-RK API csv.DictReader

    Arguments
    ----------
    * result: csv.DictReader ~ The csv.DictReader instance of a SABIO-RK kcat saerch output line

    Output
    ----------
    The organism data from the SABIO-RK API result
    """
    species_results: List[str] = []
    for row in result:
        species = row["Organism"]
        if species not in species_results:
            species_results.append(species)
    return species_results


# PUBLIC FUNCTIONS
def sabio_rk_query_with_string(query_string: str) -> str:
    """Call SABIO-RK API with the given query string.

    The query string is the actual term that is searched

    Arguments
    ----------
    * query_string: str ~ The query string

    Output
    ----------
    SABIO-RK's request result as str:

    If results could be found, it is a csv-like structure that can
    be futher processed e.g. with the Python standard library function csv.DictReader.
    The returned result fields are the EC number, the KEGG reaction ID, the organism
    and the 'parameter', i.e. the category of the given result. A possible category
    is 'kcat'.
    If no results could be found, it is the str "NO_RESULT".
    """
    # Build-up SABIO-RK query with a dictionary format that can be understood by the
    # requests library (this is also the way to retrieve potential huge amounts of
    # information as shown in SABIO-RK's API documentation).
    query = {"fields[]": ["ECNumber", "KeggReactionID", "Organism", "Parameter", "Substrate"], "q": query_string}

    # Send the request to SABIO-RK :D
    request = requests.post(QUERY_URL, params=query)

    # Error check whether the API call was successful or
    # not. 'Not successful' means that no search result with the given query
    # could be found. In this case, "NO_RESULT" is returned.
    try:
        request.raise_for_status()
    except Exception:
        print("SABIO-RK API error with query:")
        print(query_string)
        time.sleep(WAIT_TIME)
        return "NO_RESULT"

    # Wait time in order to not overload SABIO_RK's server
    time.sleep(WAIT_TIME)

    # Return successful search result in the given csv-like structure format.
    return request.text


def sabio_rk_query_with_query_dicts(query_dicts: List[Any]) -> str:
    """Performs a SABIO-RK query with given query dicts and converts them into a string

    Input
    ----------
    * query_dicts: List[Any] ~ A list of query dicts in the form of {id_name: ec_number, "Parametertype": "kcat", "EnzymeType": "wildtype"}

    Output
    ----------
    * As in sabio_rk_query_with_string()
    """
    i = 0
    # Add 'AND' for all variables of a query dict
    for i in range(len(query_dicts)):
        query_dicts[i] = " AND ".join([f"{k}:{v}" for k, v in query_dicts[i].items()])
    # Add 'OR' between each individual query dict
    query_string = "(" + ") OR (".join(query_dicts) + ")"
    # Add parentheses for URL conformity
    query_string = "(" + query_string + ")"

    # Perform the API call :D
    return sabio_rk_query_with_string(query_string)


def sabio_rk_query_get_csv_lines(query_dicts: List[Any]) -> Any:
    """Returns the result of a SABIo-RK API call with the given query dicts as csvReader lines

    Input
    ----------
    * query_dicts: List[Any] ~ A list of query dicts in the form of {id_name: ec_number, "Parametertype": "kcat", "EnzymeType": "wildtype"}

    Output
    ----------
    * The sabio_rk_query_with_string() result in a different form:
      Either "NO_RESULT" (a str) if no search esult was found,
      or a list of strings which includes the lines of the CSV that
      is returned by SABIO-RK
    """
    result = sabio_rk_query_with_query_dicts(query_dicts)
    if result == "NO_RESULT":
        return "NO_RESULT"
    else:
        return list(csv.DictReader(io.StringIO(result), delimiter="\t"))


def get_id_associated_kcats(searched_ids: List[str], id_type: str,
                            bigg_id_name_mapping_path: str, batch_size: int = 5) -> Dict[str, Any]:
    """Returns a dictionary with SABIO-RK kcat data for the given EC numbers or KEGG IDs.

    This function calls the SABIO-RK API.

    Input
    ----------
    * searched_ids: List[str] ~ The list of searched IDs
    * id_type: str ~ Must be either 'EC' or 'KEGG', depending on whether you are looking for kcats for EC numbers
      or KEGG IDs.
    * batch_size: int = 5 ~ The SABIO-RK API search batching number (i.e., with satch_size=5 five IDs are searched at once)

    Output
    ----------
    A dictionary with the following content:
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
    # Set-up the cache if it does not exist yet \o/
    cache_basepath = "./_cache/sabio_rk_total/"
    ensure_folder_existence("./_cache/")
    ensure_folder_existence(cache_basepath)
    cache_files = get_files(cache_basepath)
    # Load the given BIGG ID<->metabolite common name mapping
    bigg_id_name_mapping = json_load(bigg_id_name_mapping_path)
    # In order to save search time, use the seat (i.e., a list where
    # every member occurs only once) of the given searched IDs
    searched_ids = list(set(searched_ids))

    # Set the given ID name to the name which SABIO-RK uses for them
    if id_type == "EC":
        id_name = "ECNumber"
    elif id_type == "KEGG":
        id_name = "KeggReactionID"

    # Depending on the wildcard level which is serched, either
    # the output or the wildcard output will be used as output
    # These central dictionaries will contain the ID<->kcat mapping
    output = {}
    wildcard_output = {}
    # We use batched searched in order to save search time :D
    batch_start = 0
    # Loop while not all IDs were searched \o/
    while batch_start < len(searched_ids):
        # Get the batch for the search :-)
        batch = searched_ids[batch_start: batch_start + batch_size]
        # The query dicts contain a list of dictionaries which contain
        # the data for a SABIO-RK search entry
        query_dicts: List[Dict[str, str]] = []
        # Go through each single EC number in the search bath
        for ec_number in batch:
            # Create the cache filename
            cache_filename = ec_number.replace(".", "_").replace("*", "W") + ".json"
            # If the EC number is already searched, i.e. it can be found in the cache,
            # take the results from there in order to save much search time :D
            if cache_filename in cache_files:
                cache_filepath = cache_basepath + cache_filename
                output[ec_number] = json_load(cache_filepath)
                print(f"Loading {cache_filename}...")
            # Otherwise, create an actual SABIO-RK API search query
            else:
                query_dicts.append({id_name: ec_number, "Parametertype": "kcat", "EnzymeType": "wildtype"})
        # If not all of the searched IDs are present in the cache...
        if len(query_dicts) > 0:
            # ...use SABIO-RK's API :D
            print(f"Performing query {query_dicts}...")
            result = sabio_rk_query_get_csv_lines(query_dicts)

            # If there was an error with the SABIO-RK result (i.e., no result found or an invalid given ID),
            # continue with the next batch
            if result == "NO_RESULT":
                batch_start += batch_size
                continue
        # ...otherwise set the query result to nothing
        else:
            result = []

        # Loop through every SABIO-RK API query call result :D
        temp_ec_numbers_found_in_search = []
        result = _extract_kcat_lines(result)
        for row in result:
            # Get the unit of the parameter
            unit = row["parameter.unit"]
            # If it is a weird unusable unit, do not use this result and continue with the next result \o/
            if unit not in list(UNIT_MULTIPLIER.keys()):  # e.g. (s^-1)*(mg^-1)
                continue

            # Get the serached ID
            ec_number = row[id_name]
            # Generate a lowercarse and semicolon seperated list of substrates
            substrates_names = row["Substrate"]
            substrates_list = [x.lower() for x in substrates_names.replace("+", "").split(";")]
            substrates_list = sorted(substrates_list)
            # Convert the substrates name list into a BIGG ID list (only works
            # if there is a name<->BIGG ID mapping present for each substrate)
            bigg_ig_substrates_list = []
            for substrate in substrates_list:
                if substrate in bigg_id_name_mapping.keys():
                    bigg_id = bigg_id_name_mapping[substrate]
                    bigg_ig_substrates_list.append(bigg_id)
                # If one of the substrates cannot be found, use the pseudometabolite "REST"
                # and break :O
                else:
                    bigg_ig_substrates_list = ["REST"]
                    break
            # Set the substrate list to a semicolon-connected string
            substrate = ";".join(bigg_ig_substrates_list)
            # Get the result's organism :D
            species = row["Organism"]
            # Get the kcat and set
            # it to 1/s for consistent behaviour :D
            raw_kcat = float(row["parameter.startValue"])  # Without unit correction
            kcat = raw_kcat * UNIT_MULTIPLIER[unit]  # With unit correction 🎉

            # Add the result to the output for the given EC number, sustrate and species
            if ec_number not in output.keys():
                output[ec_number] = {}
            if substrate not in output[ec_number].keys():
                output[ec_number][substrate] = {}
            if species not in output[ec_number][substrate].keys():
                output[ec_number][substrate][species] = []
            output[ec_number][substrate][species].append(kcat)

            # Since we found a result, add the EC number :D
            temp_ec_numbers_found_in_search.append(ec_number)

        # Create cache files for all newly found EC numbers which were not present
        # in the cache
        temp_ec_numbers_found_in_search = list(set(temp_ec_numbers_found_in_search))
        for ec_number in temp_ec_numbers_found_in_search:
            cache_filename = ec_number.replace(".", "_") + ".json"
            if cache_filename not in cache_files:
                json_write(cache_basepath + cache_filename, output[ec_number])

        # Get all wildcarded searched EC numbers...
        wildcarded_searched_ec_numbers = [x for x in batch if "*" in x]
        # ...and loop through them in order to create a result for the EC numbers
        # which fit into the wildcard (i.e 1.1.1.123 in 1.1.1.*) :D
        for wildcarded_ec_number in wildcarded_searched_ec_numbers:
            # Ste the cache name for the wildcarded EC number
            cache_filename = wildcarded_ec_number.replace(".", "_").replace("*", "W") + ".json"
            # If the wildcarded EC number cannot be found in the cache, search for
            # fitting EC numbers, and combine their entries into a huge entry for the
            # wildcarded EC number
            if cache_filename not in cache_files:
                fitting_ec_numbers = []
                for found_ec_number in temp_ec_numbers_found_in_search:
                    if is_fitting_ec_numbers(wildcarded_ec_number, found_ec_number, wildcarded_ec_number.count("*")):
                        fitting_ec_numbers.append(found_ec_number)

                # Combine the EC number entries of fitting EC numbers :D
                wildcarded_ec_number_dict: Dict[str, Any] = {}
                for fitting_ec_number in fitting_ec_numbers:
                    fitting_ec_number_result = output[fitting_ec_number]
                    for metabolite_key in fitting_ec_number_result.keys():
                        if metabolite_key not in wildcarded_ec_number_dict.keys():
                            wildcarded_ec_number_dict[metabolite_key] = fitting_ec_number_result[metabolite_key]
                        else:
                            for organism_key in fitting_ec_number_result[metabolite_key].keys():
                                if organism_key not in wildcarded_ec_number_dict[metabolite_key].keys():
                                    wildcarded_ec_number_dict[metabolite_key][organism_key] =\
                                        copy.deepcopy(fitting_ec_number_result[metabolite_key][organism_key])
                                else:
                                    wildcarded_ec_number_dict[metabolite_key][organism_key] +=\
                                        copy.deepcopy(fitting_ec_number_result[metabolite_key][organism_key])
                                wildcarded_ec_number_dict[metabolite_key][organism_key] =\
                                    list(set(wildcarded_ec_number_dict[metabolite_key][organism_key]))
                # Create cache files for the searched wildcarded EC numbers \o/
                if wildcarded_ec_number_dict != {}:
                    json_write(cache_basepath + cache_filename, wildcarded_ec_number_dict)
                    wildcard_output[wildcarded_ec_number] = wildcarded_ec_number_dict
            # If the wildcarded EC number is in the cache, load the cache file :D
            else:
                wildcard_output[wildcarded_ec_number] = json_load(cache_basepath + cache_filename)
                print(f"Loading {cache_filename}...")

        # Continue with the next searched ID batch :D
        batch_start += batch_size

    # If the wildcard level is greater than 0, set the wildcard output as output
    if len(wildcard_output.keys()) > 0:
        output = wildcard_output

    return output


def get_ec_number_kcats_wildcard_search(ec_numbers: List[str],
                                        bigg_id_name_mapping_path: str,
                                        batch_size: int = 5) -> Dict[str, Any]:
    """Returns EC number-dependent kcats using an incremental wildcard level until kcarts were found for all ECs.

    Arguments
    ----------
    *ec_numbers: List[str] ~ The list of checked EC numbers.
    *bigg_id_name_mapping_path: str ~ The full path to a BIGG ID<->Name mapping
    *batch_size: int = 5 ~ How many EC numbers shall be looked up in parallel in SABIO-RK. Do not set this too high!

    Output
    ----------
    A dictionary with the following content:
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
    print("Starting EC numbers kcat search in SABIO-RK...")

    # We will look-up all EC numbers in SABIO-RK. If an EC number does not have an entry in
    # SABIO-RK, the wildcard level will become higher, e.g. 1.1.1.1 would become 1.1.1.*, and
    # if no entry is found in SABIO-EK, it would become 1.1.*.*, and so on
    # We start with no wildcards :D
    wildcard_level = 0
    # A list containing all found EC numbers, which won't be looked up using a higher wildcard level
    all_found_ec_numbers: List[str] = []
    # The central returned dictionary which will contain all combined kcat entries, divided
    # for organisms and metabolites
    ec_number_kcat_mapping: Dict[str, Any] = {}
    # Since an EC number has a maximum of 5 numbers and there is at least one EC number entry for
    # all EC number major category (e.g. 1.*.*.*), wildcards levels from 0 to 4 are reasonable,
    # and we loop throigh them :D
    for wildcard_level in range(5):
        # Get the list of all EC numbers which we want to search, i.e. all EC numbers which
        # are not already found.
        ec_numbers_to_analyze = list(set(ec_numbers) ^ set(all_found_ec_numbers))  # Difference
        # Add the current wildcard level to the searched EC numbers
        searched_ec_numbers = [_add_wildcard_to_ec_number(x, wildcard_level) for x in ec_numbers_to_analyze]
        # If no searched EC numbers are left with the current wildcard level, quit the for loop since we
        # are done :D
        if searched_ec_numbers == []:
            break
        # If there are searched EC numbers left, get the EC-number associated kcat entries
        print(f"Wildcard level {wildcard_level}...")
        print(searched_ec_numbers)
        if wildcard_level < 3:
            # With a low wildcard level, the default batch size of 5 is acceptable and
            # helps to search quicker
            resulting_ec_number_kcat_mapping = get_id_associated_kcats(searched_ec_numbers, "EC",
                                                                       bigg_id_name_mapping_path)
        else:
            # With a high wildcard level, no batching of searched EC numbers should be done since
            # the high number of results would make the search much slower D:
            resulting_ec_number_kcat_mapping = get_id_associated_kcats(searched_ec_numbers, "EC",
                                                                       bigg_id_name_mapping_path,
                                                                       batch_size=1)

        # In the following loops, the resulting dictionaries for wildcarded EC numbers are mapped to
        # the searched EC numbers, e.g. the if 1.1.1.1233 and 1.1.1.2143 are searched and no results
        # were found with wildcard level 0, the results from wildcard level 1.1.1.* are given
        # to them with the entry that it comes from a wildcard
        resulting_found_ec_numbers = resulting_ec_number_kcat_mapping.keys()
        temp_all_found_ec_numbers: List[str] = []
        for ec_number in ec_numbers_to_analyze:
            if ec_number in all_found_ec_numbers:
                continue
            for found_ec_number in resulting_found_ec_numbers:
                if not is_fitting_ec_numbers(ec_number, found_ec_number, wildcard_level):
                    continue
                kcat = resulting_ec_number_kcat_mapping[found_ec_number]
                if ec_number not in list(ec_number_kcat_mapping.keys()):
                    ec_number_kcat_mapping[ec_number] = kcat
                    temp_all_found_ec_numbers += [ec_number]
                else:
                    ec_number_kcat_mapping[ec_number] = {**ec_number_kcat_mapping[ec_number], **kcat}
                if wildcard_level == 0:
                    ec_number_kcat_mapping[ec_number]["WILDCARD"] = False
                else:
                    ec_number_kcat_mapping[ec_number]["WILDCARD"] = True

        # Continue with the new list of found EC numbers and the next wildcard level :D
        all_found_ec_numbers += temp_all_found_ec_numbers
        wildcard_level += 1

    # Return the found major EC-number<->kcat mapping \o/
    return ec_number_kcat_mapping


"""
Exemplary usage:
kegg_ids = ["R00006", "R00286", "R00086"]
result = get_id_associated_kcats(kegg_ids, "KEGG")
print(result)

ec_numbers = ["1.1.1.213441", "2.7.1.124243", "1.234243.123213.123213", "2.12323.123213.12213"]
result = get_ec_number_kcats_wildcard_search(ec_numbers)
print(result)
"""
