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
"""get_protein_mass_mapping.py

Functions for the generation of a model's mapping of its proteins and their masses.
"""

# IMPORTS
# External modules
import cobra
import requests
import time
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from typing import Dict, List
# Internal modules
from .helper_general import ensure_folder_existence, get_files, json_write, pickle_write, pickle_load, standardize_folder


# PUBLIC FUNCTIONS SECTION
def get_protein_mass_mapping(model: cobra.Model, project_folder: str, project_name: str) -> None:
    """Returns a JSON with a mapping of protein IDs as keys, and as values the protein mass in kDa.

    The protein masses are calculated using the amino acid sequence from UniProt (retrieved using
    UniProt's REST API).

    Arguments
    ----------
    * model: cobra.Model ~ The model in the cobrapy format
    * project_folder: str ~ The folder in which the JSON shall be created
    * project_name: str ~ The beginning of the JSON's file name

    Output
    ----------
    A JSON file with the path project_folder+project_name+'_protein_id_mass_mapping.json'
    and the following structure:
    <pre>
    {
        "$PROTEIN_ID": $PROTEIN_MASS_IN_KDA,
        (...),
    }
    </pre>
    """
    # Standardize project folder
    project_folder = standardize_folder(project_folder)

    # The beginning of the created JSON's path :D
    basepath: str = project_folder + project_name

    # GET UNIPROT ID - PROTEIN MAPPING
    uniprot_id_protein_id_mapping: Dict[str, List[str]] = {}
    for gene in model.genes:
        # Without a UniProt ID, no mass mapping can be found
        if "uniprot" not in gene.annotation:
            continue
        uniprot_id = gene.annotation["uniprot"]
        if uniprot_id in uniprot_id_protein_id_mapping.keys():
            uniprot_id_protein_id_mapping[uniprot_id].append(gene.id)
        else:
            uniprot_id_protein_id_mapping[uniprot_id] = [gene.id]

    # GET UNIPROT ID<->PROTEIN MASS MAPPING
    uniprot_id_protein_mass_mapping: Dict[str, float] = {}
    # The cache stored UniProt masses for already searched
    # UniProt IDs (each file in the cache folder has the name
    # of the corresponding UniProt ID). This prevents searching
    # UniProt for already found protein masses. :-)
    cache_basepath = "./_cache/uniprot/"
    ensure_folder_existence("./_cache/")
    ensure_folder_existence(cache_basepath)
    cache_files = get_files(cache_basepath)
    # Go through each batch of UniProt IDs (multiple UniProt IDs
    # are searched at once in order to save an amount of UniProt API calls)
    # and retrieve the amino acid sequences and using these sequences, their
    # masses.
    print("Starting UniProt ID<->Protein mass search using UniProt API...")
    uniprot_ids = list(uniprot_id_protein_id_mapping.keys())
    batch_size = 5
    batch_start = 0
    while batch_start < len(uniprot_ids):
        # Create the batch with all UniProt IDs
        prebatch = uniprot_ids[batch_start:batch_start+batch_size]
        batch = []
        # Remove all IDs which are present in the cache (i.e.,
        # which were searched for already).
        # The cache consists of pickled protein mass floats, each
        # onein a file with the name of the associated protein.
        for uniprot_id in prebatch:
            if uniprot_id not in cache_files:
                batch.append(uniprot_id)
            else:
                cache_filepath = cache_basepath + uniprot_id
                uniprot_id_protein_mass_mapping[uniprot_id] = pickle_load(cache_filepath)
                print(uniprot_id+":", uniprot_id_protein_mass_mapping[uniprot_id])

        # If all IDs could be found in the cache, continue with the next batch.
        if len(batch) == 0:
            batch_start += batch_size
            continue

        # Create the UniProt query for the batch
        # With 'OR', all given IDs are searched, and subsequently in this script,
        # the right associated masses are being picked.
        query = " OR ".join(batch)
        uniprot_query_url = f"https://www.uniprot.org/uniprot/?query={query}&format=tab&columns=id,sequence"
        print(f"UniProt batch search for: {query}")

        # Call UniProt's API :-)
        uniprot_data = requests.get(uniprot_query_url).text.split("\n")
        # Wait in order to cool down their server :-)
        time.sleep(2.0)

        # Read out the API-returned lines
        for line in uniprot_data[1:]:
            if line == "":
                continue
            uniprot_id = line.split("\t")[0]
            sequence = line.split("\t")[1]
            # Get the protein mass using biopython's associated function for amino acid sequences
            try:
                mass = ProteinAnalysis(sequence, monoisotopic=False).molecular_weight()
            except ValueError:  # e.g. if an "X" is in a sequence
                continue
            uniprot_id_protein_mass_mapping[uniprot_id] = float(mass)

        # Create the pickled cache files for the searched protein masses
        for uniprot_id in batch:
            cache_filepath = cache_basepath + uniprot_id
            pickle_write(cache_filepath, uniprot_id_protein_mass_mapping[uniprot_id])

        # Continue with the next batch :D
        batch_start += batch_size

    # Create the final protein ID <-> mass mapping
    protein_id_mass_mapping: Dict[str, float] = {}
    for uniprot_id in list(uniprot_id_protein_mass_mapping.keys()):
        try:
            protein_ids = uniprot_id_protein_id_mapping[uniprot_id]
        except Exception:
            print(f"No mass found for {uniprot_id}!")
            continue
        for protein_id in protein_ids:
            protein_id_mass_mapping[protein_id] = uniprot_id_protein_mass_mapping[uniprot_id]

    # Write protein mass list JSON :D
    print("Protein ID<->Mass mapping done!")
    json_write(basepath+"_protein_id_mass_mapping.json", protein_id_mass_mapping)


def get_protein_mass_mapping_with_sbml(sbml_path: str, project_folder: str, project_name: str) -> None:
    """This module's get_protein_mass_mapping() with SBML instead of a cobrapy module as argument.

    Arguments
    ----------
    * sbml_path: str ~ The path to the model's SBML
    * project_folder: str ~ The folder in which the JSON shall be created
    * project_name: str ~ The beginning of the JSON's file name
    """
    model: cobra.Model = cobra.io.read_sbml_model(sbml_path)
    get_protein_mass_mapping(model, project_folder, project_name)
