#!/usr/bin/env python3
#
# Copyright 2018-2019 PSB
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
"""helper_general.py

This module contains functions which are useful for a multitude of programs,
and which are hard to be categorized.
"""

# IMPORTS
# External modules
import json
import openpyxl
import os
import pickle
import requests
import sys
import time
from pebble import concurrent
from typing import Any, Dict, List


# CONSTANT SECTION
# As no useful UniProt API access libraries were found, the programs use the UniProt URL API directly.
# UNIPROT_UPLOADLISTS_URL represents the UniProt ID conversion search URL. The paremeters for the
# conversion are appendet after this string.
UNIPROT_UPLOADLISTS_URL = "https://www.uniprot.org/uploadlists"
# UNIPROT_URL represents the UniProtKB database search URL. After this string, the parameters are appended
# (e.g. a UniProt ID to refer to the protein's entry).
UNIPROT_URL = "https://www.uniprot.org/uniprot"


# INTERNAL FUNCTIONS SECTION
@concurrent.process(timeout=300)
def _get(url: str) -> str:
    """Internal function which tries to get the text of the given internet site.

    Using pebble as concurrency library (see the decorator @concurrent.process),
    this function stops with a timeout error after 300 seconds = 5 minutes.
    This value was chosen as in the programs which use this internal function,
    no request to biological databases should take longer than this time.

    Argument
    ----------
    * url: str ~ The URL of the internet site.
    """
    # Retrieve the text from the GET request representation directly.
    text = requests.get(url).text

    return text


# PUBLIC FUNCTIONS SECTION
def check_argument(argument: str, command: str) -> str:
    """Checks wheter the given mandatory argument variable is None or not. If none, the program stops.

    It is checked wheter the argument is None as argparser sets an argument variable to None
    if the argument is not set.

    Arguments
    ----------
    * argument: str ~ The argument that shall be checked and which shall not be None.
    * command: str ~ The argument long name (i.e. the name after --) which is displayed
      in an error message if the argument is None.
    """
    if argument is None:
        print(f"Mandatory argument --{command} not given. Stopping program.")

        # Exit program with error code (not 0).
        sys.exit(-1)

    return argument


def ensure_folder_existence(folder: str) -> None:
    """Checks if the given folder exists. If not, the folder is created.

    Argument
    ----------
    * folder: str ~ The folder whose existence shall be enforced.
    """
    if os.path.isdir(folder):
        return
    os.makedirs(folder)


def get(url: str, num_tries: int = 3, sleep_time: float = 10.0) -> List[str]:
    """Returns an internet site's content in the form of a list with one string per line.

    Arguments
    ----------
    * url:str ~ The URL of the internet site.
    * num_tries: int=3 ~ The number of maximal tries to retrieve the site.
      Its default value of 3 resulted of empirical tests.
    * sleep_time: float=10.0 ~ The number of seconds that shall be waited
      after an error or after the site's text was retrieved.
      Its default value of 10.0 seconds is in accordance with the NCBI
      rule that its servers shall not be contacted more often than every
      10 seconds.
    """
    # Set text variable which will contain the site's text if no error happens.
    text = ""

    # Try it until num_tries is reached or until the site's content was
    # successfully retrieved.
    for i in range(num_tries):
        # Wait the given time in doing nothing.
        time.sleep(sleep_time)

        # Try connection. May fail due to timeout (see _get()'s comments)
        # or due to a client- or server-side connection problem.
        try:
            # _get() is the internal function which actually sends the
            # request to the site. It is a separate function in order to
            # be able to start it in concurrent mode so that it is possible
            # to await a timeout.
            # Otherwise, if _get's content would be in
            # this function, this function and the whole program would be
            # blocked without the possibilty to stop after a set timeout.
            # A 'future' is a variable controlling a concurrently running
            # process.
            future = _get(url)
            text = future.result()

            # If no error/timeout occured while calling _get, we reach the
            # break statement which causes that no other try will be performed
            # as we got the site's data successfully.
            break
        except Exception as error:
            # Print error. The error may be a timeout or a connection problem.
            print(f"Error with {url}")
            print(f"Error: {error}")
            print(f"In try {i+1}")

            # Wait the time until the next try if the maximal number of tries
            # is not reached yet.
            time.sleep(sleep_time)

    # Split the internet site's text into a string list containing each line.
    text_lines: List[str] = text.split("\n")

    # Return the line list.
    return text_lines


def get_and_save(url: str, filepath: str) -> List[str]:
    """Retrieves the text content from the given URL and saves it in the given text file path.

    It also returns the content of the URL as a list of strings (one list element per newline-separated line).

    Arguments
    ----------
    * url: str ~ The URL of which the text content is read out
    * filepath: str ~ The path pof the text file which will include the URL text content
    """
    text = get(url)

    with open(filepath, "w") as f:
        f.write("".join([i + "\n" for i in text]))

    return text


def get_entry_id_kegg_id_mapping(xlsx_path: str, worksheet_name: str) -> Dict[str, str]:
    """Reads a worksheet in which the user can select the KEGG IDs of e.g. metabolites or reactions.

    This is used for the KEGG ID XLSX for reactions and metabolites in AutoPACMEN.
    THe functions returns a dictionary of IDs as keys and selected KEGG IDs as values.

    Arguments
    ----------
    * xlsx_path: str ~ The path to the analyzed XLSX
    * worksheet_name: str ~ THe name of the analyzed worksheet of the analyzed XLSX
    """
    workbook = openpyxl.load_workbook(filename=xlsx_path, read_only=True)
    worksheet = workbook[worksheet_name]

    entry_id_kegg_id_mapping: Dict[str, str] = {}
    for row in worksheet.rows:
        current_cell = 1
        assessment_found = False
        for cell in row:
            if current_cell == 1:
                reaction_id = cell.value
            if current_cell > 3:
                if ((current_cell - 1) % 3) == 0:
                    current_kegg_id = cell.value

                if ((current_cell - 2) % 3) == 0:
                    assessment = cell.value
                    if assessment is not None:
                        if assessment == "Yes":
                            assessment_found = True
                            break
            current_cell += 1
        if assessment_found:
            entry_id_kegg_id_mapping[reaction_id] = current_kegg_id
        else:
            print(f"ERROR: No assessment for reaction {reaction_id}!")
            print(f"Spreadsheet: {xlsx_path}")
            print(f"Terminating program...")
            sys.exit(-1)
    return entry_id_kegg_id_mapping


def get_files(path: str) -> List[str]:
    """Returns the names of the files in the given folder as a list of strings.

    Arguments
    ----------
    * path: str ~ The path to the folder of which the file names shall be returned
    """
    files: List[str] = []
    for (_, _, filenames) in os.walk(path):
        files.extend(filenames)
    return files


def get_float_cell_value(cell_value) -> float:
    """Returns the value of an openpyxl cell value.

    This function is used in oder to operate with spreadsheets
    which use the comma as the decimal mark instead of a period.

    Arguments
    ----------
    * cell_value ~ The openpyxl cell value
    """
    if type(cell_value) not in (int, float):
        cell_value = cell_value.replace(",", ".")
        cell_value = float(cell_value)
    return cell_value


def is_fitting_ec_numbers(ec_number_one: str, ec_number_two: str, wildcard_level: int) -> bool:
    """Check whether the EC numbers are the same under the used wildcard level.

    Arguments
    ----------
    * ec_number_one: str ~ The first given EC number.
    * ec_number_two: str ~ The second given EC number.
    * wildcard_level: int ~ The wildcard level.
    """
    if wildcard_level == 0:
        ec_number_one_full_numbers = ec_number_one.split(".")
        ec_number_two_full_numbers = ec_number_two.split(".")
    else:
        ec_number_one_full_numbers = ec_number_one.split(".")[:-wildcard_level]
        ec_number_two_full_numbers = ec_number_two.split(".")[:-wildcard_level]

    if ec_number_one_full_numbers == ec_number_two_full_numbers:
        return True
    else:
        return False


def json_load(path: str) -> Dict[Any, Any]:
    """Loads the given JSON file and returns it as dictionary.

    Arguments
    ----------
    * path: str ~ The path of the JSON file
    """
    with open(path) as f:
        dictionary = json.load(f)
    return dictionary


def json_write(path: str, dictionary: Dict[Any, Any]) -> None:
    """Writes a JSON file at the given path with the given dictionary as content.

    Arguments
    ----------
    * path: str ~  The path of the JSON file that shall be written
    * dictionary: Dict[Any, Any] ~ The dictionary which shalll be the content of
      the created JSON file
    """
    json_output = json.dumps(dictionary, indent=4)
    with open(path, "w", encoding="utf-8") as f:
        f.write(json_output)


def pickle_load(path: str) -> Any:
    """Returns the value of the given pickle file.

    Arguments
    ----------
    * path: str ~ The path to the pickle file.
    """
    pickle_file = open(path, 'rb')
    pickled_object = pickle.load(pickle_file)
    pickle_file.close()
    return pickled_object


def pickle_write(path: str, pickled_object: Any) -> None:
    """Writes the given object as pickled file with the given path

    Arguments
    ----------
    * path: str ~ The path of the pickled file that shall be created
    * pickled_object: Any ~ The object which shall be saved in the pickle file
    """
    pickle_file = open(path, 'wb')
    pickle.dump(pickled_object, pickle_file)
    pickle_file.close()


def mkdir(path: str) -> None:
    """Creates a directory if the given directory does not exist yet.

    Uses os.mkdir and os.path.exists internally.

    Argument
    ----------
    * path: str ~ The directory's path.
    """
    if not os.path.exists(path):
        os.mkdir(path)


def resolve_pathway_ids(pathway_ids: str, pathways: List[Any]) -> List[int]:
    """Returns a list index integer list of the given pathway_ids for the given pathway list.

    It gets 'all' and returns all indices of pathways or a comma-separated list of pathway numbers.

    Example
    ----------
    If an organism has 5 pathways in total, this function would return for 'all' as pathway_ids
    [0, 1, 2, 3, 4], and for '1,3,4' it would return [0, 2, 3]

    Arguments
    ----------
    * pathway_ids ~ The IDs of the pathways in a string. Can be 'all' to get all indices, or
      a comma-separated list of numbers.
    * pathways ~ The list of pathways of which an index integer list shall be created.
    """
    # Check if 'all' is given...
    if pathway_ids == "all":
        return [i for i in range(len(pathways))]
    # ...if not, it should be a comma-separated list of numbers.
    else:
        # Return a list with each number as integer and substracted by 1.
        # The substraction is done as e.g. the 1st pathway of an organism
        # has the index 0 in a Python pathway list.
        pathway_ids_list = pathway_ids.split(",")
        return [int(i) - 1 for i in pathway_ids_list]


def sanitize_path(text: str) -> str:
    """Replaces all invalid characters for a path with valid ones.

    E.g. useful for creating files with the name of KEGG pathways,
    as these names may contain invalid characters.

    Argument
    ----------
    * text: str ~ The string that may contain invalid characters.
    """
    return text.replace("\\", "_").replace("/", "_").\
        replace(":", "_").replace("*", "_").\
        replace("<", "_").replace(">", "_").\
        replace("|", "_")


def standardize_folder(folder: str) -> str:
    """Returns for the given folder path is returned in a more standardized way.

    I.e., folder paths with potential \\ are replaced with /. In addition, if
    a path does not end with / will get an added /.

    Argument
    ----------
    * folder: str ~ The folder path that shall be standardized.
    """
    # Standardize for \ or / as path separator character.
    folder = folder.replace("\\", "/")

    # If the last character is not a path separator, it is
    # added so that all standardized folder path strings
    # contain it.
    if folder[-1] != "/":
        folder += "/"

    return folder
