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
"""parse_brenda_textfile.py

Includes a function which converts a BRENDA textfile into a machine-readable JSON :D
"""

# IMPORTS
# External modules
from typing import Any, Dict, List
# Internal modules
from .helper_general import json_write, json_load, standardize_folder


# PUBLIC FUNCTIONS SECTION
def parse_brenda_textfile(brenda_textfile_path: str, bigg_metabolites_json_folder: str,
                          json_output_path: str) -> None:
    """Goes through a BRENDA database textfile and converts it into a machine-readable JSON.

    The JSON includes kcats for found organisms and substrates.
    As of 29/04/2019, the BRENDA database can be downloaded as textfile under
    https://www.brenda-enzymes.org/download_brenda_without_registration.php

    The BRENDA database is not in a completely standardized format, so that this functions
    contains many convoluted checks and circumventions of non-standardized data.

    kcats from mutated enzymes are excluded.

    Arguments
    ----------
    * brenda_textfile_path: str ~ The BRENDA database text file path
    * bigg_metabolites_json_folder: str ~ The folder in which the BIGG metabolites
      database is stored (it has to have the name 'bigg_id_name_mapping.json').
    * json_output_path: str ~ The path of the JSON that shall be created

    Output
    ----------
    * A JSON containing the BRENDA textfile kcat data in a machine-readable format:
    <pre>
        {
            "$EC_NUMBER": {
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
    # Standardize output folder
    bigg_metabolites_json_folder = standardize_folder(
        bigg_metabolites_json_folder)

    # Load BIGG ID <-> metabolite name mapping :D
    bigg_id_name_mapping: Dict[str, str] = json_load(
        bigg_metabolites_json_folder+"bigg_id_name_mapping.json")

    # Load BRENDA textfile as list of strings without newlines :D
    with open(brenda_textfile_path, "r", encoding="utf-8") as f:
        lines = f.readlines()
    lines = [x.replace("\n", "") for x in lines]

    # Go through each line and collect the organism lines and kcat lines for each EC number
    in_turnover_numbers = False
    in_organism_reference = False
    ec_number_kcat_lines_mapping: Dict[str, List[str]] = {}
    ec_number_organsism_lines_mapping: Dict[str, List[str]] = {}
    current_ec_number = ""
    organism_lines: List[str] = []
    kcat_lines: List[str] = []
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("ID\t"):
            if current_ec_number != "":
                ec_number_organsism_lines_mapping[current_ec_number] = organism_lines
                ec_number_kcat_lines_mapping[current_ec_number] = kcat_lines
            current_ec_number = line.replace("ID\t", "").replace(" ()", "")
            organism_lines = []
            kcat_lines = []

        if len(line) == 0:
            in_turnover_numbers = False
            in_organism_reference = False
        elif line.startswith("PROTEIN"):
            in_organism_reference = True
            i += 1
            line = lines[i]
        elif line.startswith("TURNOVER_NUMBER"):
            in_turnover_numbers = True
            i += 1
            line = lines[i]

        if in_organism_reference:
            if line.startswith("PR"):
                organism_lines.append("")
            if len(organism_lines[-1]) > 0:
                organism_lines[-1] += " "
            organism_lines[-1] += " " + line

        elif in_turnover_numbers:
            if line.startswith("TN"):
                kcat_lines.append("")
            if len(kcat_lines[-1]) > 0:
                kcat_lines[-1] += " "
            kcat_lines[-1] += line

        if len(line) == 0:
            in_turnover_numbers = False
            in_organism_reference = False

        i += 1

    # Create the BRENDA database dictionary using the collected kcat and organism lines
    # of each EC number :D
    ec_numbers = list(ec_number_kcat_lines_mapping.keys())
    brenda_kcat_database: Dict[str, Any] = {}
    for ec_number in ec_numbers:
        if "(transferred to " in ec_number:
            actual_ec_number = ec_number.split(" (transferred")[0]
            try:
                brenda_kcat_database[actual_ec_number] = {}
                brenda_kcat_database[actual_ec_number]["TRANSFER"] = \
                    ec_number.lower().replace("  ", " ").split(
                        "(transferred to ec")[1].replace(")", "").lstrip()
            except Exception:
                # Some transfers go to general subgroups instead of single EC numbers so that
                # no kcat database can be built from it D:
                print("WARNING: BRENDA text file line " +
                      ec_number + " is not interpretable!")
            continue

        brenda_kcat_database[ec_number] = {}

        reference_number_organism_mapping = {}
        organism_lines = ec_number_organsism_lines_mapping[ec_number]
        for organism_line in organism_lines:
            reference_number = organism_line.split("#")[1]
            organism_line_split_first_part = organism_line.split("# ")[1]
            organism_line_split = organism_line_split_first_part.split(" ")
            organism_line_split = [
                x for x in organism_line_split if len(x) > 0]

            end = 1
            for part in organism_line_split:
                # Some organism names contain their SwissProt or UniProt ID,
                # since we don't nned them they are excluded
                if ("swissprot" in part.lower()) or \
                    (part.lower() == "and") or \
                    ("uniprot" in part.lower()) or \
                    ("genbank" in part.lower()) or \
                        ("trembl" in part.lower()):
                    end -= 2
                    break

                if ("<" in part) or ("(" in part):
                    end -= 1
                    break

                end += 1
            organism_name = " ".join(organism_line_split[:end])
            reference_number_organism_mapping[reference_number] = organism_name

        kcat_lines = ec_number_kcat_lines_mapping[ec_number]
        for kcat_line in kcat_lines:
            kcat_line = kcat_line
            # Exclude kcats of mutated/changed proteins since
            # they may not have a biological relevance
            if ("mutant" in kcat_line.lower()) or ("mutated" in kcat_line.lower()):
                continue
            reference_number = kcat_line.split("#")[1].split(",")[0]
            organism = reference_number_organism_mapping[reference_number]
            kcat_str = "".join(kcat_line.split("#")[2]).split("{")[
                0].lstrip().rstrip()
            kcat = max([float(x) for x in kcat_str.split("-") if len(x) > 0])
            substrate = "".join(kcat_line.split("{")[1]).split("}")[0]

            substrate = substrate.lower()
            if substrate in bigg_id_name_mapping.keys():
                substrate = bigg_id_name_mapping[substrate]
            else:
                substrate = "REST"

            if substrate not in brenda_kcat_database[ec_number].keys():
                brenda_kcat_database[ec_number][substrate] = {}
            if organism not in brenda_kcat_database[ec_number][substrate].keys():
                brenda_kcat_database[ec_number][substrate][organism] = []
            brenda_kcat_database[ec_number][substrate][organism].append(kcat)

    # Write final BRENDA kcat database :D
    json_write(json_output_path, brenda_kcat_database)
