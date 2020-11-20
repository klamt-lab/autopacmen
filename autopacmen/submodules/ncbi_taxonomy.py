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
"""ncbi_taxonomy.py

This module contains functions which can access NCBI TAXONOMY.
"""

# IMPORTS
# External modules
import time
from Bio import Entrez
from typing import Dict, List


# SCRIPT-WIDE CONSTANTS
WAIT_TIME = .5  # Time to wait for each API call
NCBI_BATCH_SIZE = 20


# PUBLIC FUNCTIONS
def get_entrez_id_from_organism_full_name(organism_full_name):
    """Get organism's Entrez numeric identifier.

    This numeric identifier is neccessary for BLAST and NCBI TAXONOMY
    searches.
    This function uses Biopython functions. Returns BLAST-compatible ID as
    txid + NCBI ID + [ORGN].

    Arguments:
    >organism_kegg_id: str ~ The organism's full name, e.g. "Xanthomonas
     campestris pv. campesris B100"
    """
    # An e-mail has to be set, you may change it to yours if you want to
    # be notified if any problems occur.
    Entrez.email = "x@x.x"
    # Set the Entrez search to the NCBI TAXONOMY database.
    handle = Entrez.esearch(db="Taxonomy", term=organism_full_name)
    # Wait in order to not overload the NCBI's server
    time.sleep(WAIT_TIME)
    # Reformat the Entrez search result in order to extract the Entrez ID
    record = Entrez.read(handle)
    organism_ncbi_id = record["IdList"][0]
    # txid+NUMBER+[ORGN] is the form that is used for NCBI BLASTP searches to restrict a search
    # to an organism using the Entrez query constraint input.
    organism_ncbi_id = "txid"+organism_ncbi_id+"[ORGN]"
    # Return the retrieved ID :D
    return organism_ncbi_id


def get_taxonomy_from_organism_ncbi_id(organism_ncbi_id):
    """Get organism's taxonomy from NCBI Taxonomy using Biopython functions.

    The taxonomy is returned as list, starting with the nearest and
    ending with the highest taxonomic level above the organism.

    Arguments:
    >ncbi_organism_id: str ~ The organism's NCBI ID, e.g. retrieved by
     this module's "get_entrez_id_from_organism_full_name" function, in
     the format txid + NCBI ID + [ORGN]
    """
    Entrez.email = "x@x.x"
    handle = Entrez.efetch(db="Taxonomy", id=organism_ncbi_id, retmode="xml")
    records = Entrez.read(handle)
    taxonomy = records[0]["Lineage"].split(";")[::-1]
    taxonomy = [i.lstrip() for i in taxonomy]
    return taxonomy


def get_entrez_id_from_organism_full_name_batch(organism_full_names: List[str]) -> List[str]:
    """Retrieves the Entrez numeric ID of the given organisms.

    This numeric identifier is neccessary for BLAST and NCBI TAXONOMY
    searches.
    This function uses Biopython functions. Returns BLAST-compatible ID as
    txid + NCBI ID + [ORGN].

    Arguments:
    >organism_full_names: List[str] ~ A list of full names of organisms, e.g. "Xanthomonas
     campestris pv. campesris B100"
    """
    batch_start = 0
    organism_ncbi_ids_result: List[str] = []
    # Go through each organism :D
    while batch_start < len(organism_full_names):
        organism_full_names_slice = organism_full_names[batch_start:batch_start+NCBI_BATCH_SIZE]
        query_names = " OR ".join(organism_full_names_slice)
        # An e-mail has to be set, you may change it to yours if you want to
        # be notified if any problems occur.
        Entrez.email = "x@x.x"
        # Set the Entrez search to the NCBI TAXONOMY database.
        handle = Entrez.esearch(db="Taxonomy", term=query_names)
        # Wait in order to not overload the NCBI's server
        time.sleep(WAIT_TIME)
        # Reformat the Entrez search result in order to extract the Entrez ID
        record = Entrez.read(handle)
        organism_ncbi_ids = record["IdList"][::-1]
        # txid+NUMBER+[ORGN] is the form that is used for NCBI BLASTP searches to restrict a search
        # to an organism using the Entrez query constraint input.
        organism_ncbi_ids_result += ["txid"+x +
                                     "[ORGN]" for x in organism_ncbi_ids]

        batch_start += NCBI_BATCH_SIZE
        time.sleep(WAIT_TIME)
    # Return the retrieved IDs :D
    return organism_ncbi_ids_result


def get_taxonomy_from_organism_ncbi_id_batch(organism_ncbi_ids: List[str]) -> Dict[str, List[str]]:
    """Get the taxonomy from NCBI Taxonomy of the given organisms using Biopython functions.

    The taxonomy is returned as Dictionary (Dict[str, List[str]) for each organism,
    where each value is a string list starting with the nearest and
    ending with the highest taxonomic level above the organism.

    Arguments:
    >organism_ncbi_ids: List[str] ~ The list of the NCBI IDs of the organisms,
     e.g. retrieved by this module's "get_entrez_id_from_organism_full_name"
     function, in the format txid + NCBI ID + [ORGN]
    """
    taxonomies: Dict[str, List[str]] = {}
    batch_start = 0
    while batch_start < len(organism_ncbi_ids):
        organism_ncbi_ids_slice = organism_ncbi_ids[batch_start:batch_start+NCBI_BATCH_SIZE]
        query_ids = " OR ".join(organism_ncbi_ids_slice)
        Entrez.email = "x@x.x"
        handle = Entrez.efetch(db="Taxonomy", id=query_ids, retmode="xml")
        records = Entrez.read(handle)
        for record in records:
            taxonomy = record["Lineage"].split(";")[::-1]
            taxonomy = [i.lstrip() for i in taxonomy]
            taxonomies[record["ScientificName"]] = taxonomy
        batch_start += NCBI_BATCH_SIZE
    return taxonomies


def most_taxonomic_similar(base_species: str, taxonomy_dict: Dict[str, List[str]]) -> Dict[str, int]:
    """Returns a dictionary with a score of taxonomic distance from the given organism.

    e.g. if base_species is "Escherichia coli" and taxonomy_dict is
    <pre>
    {
        "Escherichia coli": ["Escherichia", "Bacteria", "Organism"],
        "Pseudomonas": ["Pseudomonas", "Bacteria", "Organism"],
        "Homo sapiens": ["Homo", "Mammalia", "Animalia", "Organism"],
    }
    </pre>
    this function would return
    <pre>
    {
        "Escherichia coli": 0,
        "Pseudomonas": 1,
        "Homo sapiens": 2,
    }
    </pre>

    Arguments
    ----------
    * base_species: str ~ The species to which a relation is made.
    * taxonomy_dict: Dict[str, List[str]] ~ A dictionary with organism names as keys and
      their taxonomic levels (sorted from nearest to farthest) as string list.
    """
    base_taxonomy = taxonomy_dict[base_species]
    level: int = 0
    level_dict: Dict[str, int] = {}
    for taxonomic_level in base_taxonomy:
        level_dict[taxonomic_level] = level
        level += 1

    score_dict: Dict[str, int] = {}
    for species in taxonomy_dict.keys():
        for taxonomic_level in taxonomy_dict[species]:
            if taxonomic_level in list(level_dict.keys()):
                score_dict[species] = level_dict[taxonomic_level]
                break

    return score_dict


"""
# Example:
organism_ncbi_ids = get_entrez_id_from_organism_full_name_batch(["Escherichia coli", "Escherichia fergusonii", "Vibrio natriegens", "Mus musculus"])
taxonomies = get_taxonomy_from_organism_ncbi_id_batch(organism_ncbi_ids)
print(taxonomies)
print(most_taxonomic_similar("Escherichia coli", taxonomies))
"""
