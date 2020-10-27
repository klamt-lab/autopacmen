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
"""get_reactions_kcat_mapping.py

This module contains functions which, in total, return a mapping of
reactions to forward and reverse reaction direction kcats.
"""

# IMPORTS
# External modules
import cobra
import math
import random
import sys
import statistics
from typing import Any, Dict, List
# Internal modules
from .ncbi_taxonomy import get_entrez_id_from_organism_full_name_batch, get_taxonomy_from_organism_ncbi_id_batch, most_taxonomic_similar
from .helper_general import ensure_folder_existence, get_files, json_load, json_write, pickle_load, pickle_write, standardize_folder


# PRIVATE FUNCTIONS
def _get_kcat_from_protein_kcat_database(searched_direction: str, reaction: cobra.Reaction, protein_kcat_database):
    """Returns the kcat from the given protein<->kcat database for the given reaction, if there is one.

    The kcat is minimum of the maximal kcats for each protein of the reaction which can be found in the database.

    Arguments:
    *searched_direction: str ~ The direction in which a kcat is searched
    *reaction: cobra.Reaction ~ The reaction for which the kcat is searched
    *protein_kcat_database ~ The protein<->kcat database

    Output:
    Either the protein database kcat of the reaction in the given direction, or math.nan if there is no
    kcat for this reaction.
    """
    # Get the reaction's gene names
    gene_reaction_rule = reaction.gene_reaction_rule
    gene_reaction_rule = gene_reaction_rule.replace(" or ", "\t")
    gene_reaction_rule = gene_reaction_rule.replace(" and ", "\t")
    gene_names = gene_reaction_rule.split("\t")

    # Get the maximal kcats for each gene name in the given reaction direction
    max_kcats = []
    for gene_name in gene_names:
        if gene_name not in protein_kcat_database.keys():
            continue
        kcat_direction = protein_kcat_database[gene_name]["direction"][reaction.id]
        max_kcat = max(protein_kcat_database[gene_name]["kcats"])

        if kcat_direction == searched_direction == "forward":
            max_kcats.append(max_kcat)
        else:
            max_kcats.append(max_kcat)

    # Get the minimal maximal kcat and return it :D
    if len(max_kcats) > 0:
        min_max_kcat = min(max_kcats)
    else:
        min_max_kcat = math.nan

    return min_max_kcat


def _get_kcat_list(searched_metabolites: List[str], complete_entry: Dict[str, Any], organism: str,
                   searched_direction: str, reaction: cobra.Reaction,
                   protein_kcat_database) -> List[float]:
    """Returns a list of kcats for the given reaction, created with a taxonomic search and usingg the protein kcat database.

    Algorithm
    ----------
    The given kcat entries are changed so that the kcats are now ordered and associated for each organism and
    the searched metabolies. Using NCBI Taxonomy, the taxonomic distance of the kcat tnry organisms to the given
    organisms is calculated, and - until the highest taxonomic distance or the desired minimal kcat entry length is
    reached - the kcats from the taxonomically nearest organisms are added. If there is a protein database kcat,
    it is added to the list too.

    Arguments
    ----------
    * searched_metabolites: List[str] ~ The list of fitting metabolites for the reaction
    * complete_entry: Dict[str, Any] ~ The complete entry of kcats for this entry
    * organism: str ~ The organism of the model which shall become protein-constraint-enhanced
    * searched_direction: str ~ The affected reaction's direction, e.g. 'forward'
    * reaction: cobra.Reaction ~ The affected reaction for which a kcat shall be calculated
    * protein_kcat_database ~ The protein<->kcat database

    Output
    ----------
    The list of kcats which results from the taxonomic search
    """
    # Create a dictionary with kcat entries for each organism, determined from the given original kcat entry
    species_kcat_mapping: Dict[str, List[float]] = {}
    for searched_metabolite in searched_metabolites:
        species_entries = complete_entry[searched_metabolite]
        for species in list(species_entries.keys()):
            if species not in species_kcat_mapping.keys():
                species_kcat_mapping[species] = []
            species_kcat_mapping[species] += species_entries[species]

    # Create the list of all species
    all_species = list(species_kcat_mapping.keys())
    if organism not in all_species:
        all_species.append(organism)
        organism_added = True
    else:
        organism_added = False

    # Get the taxonomy of the kcat entry's organisms, either by calling
    # NCBI Taxonomy or by reading the cache of previously searched organisms
    cache_basepath = "./_cache/ncbi_taxonomy/"
    ensure_folder_existence("./_cache/")
    ensure_folder_existence(cache_basepath)
    cache_files = get_files(cache_basepath)
    species_to_search: List[str] = []
    taxonomy_dict_cache: Dict[str, List[str]] = {}
    for species in all_species:
        cache_filename = species.replace("/", "") + "_taxonomy"
        if cache_filename in cache_files:
            cache_filepath = cache_basepath + cache_filename
            taxonomy_dict_cache[species] = pickle_load(cache_filepath)
        elif cache_filename+"_NA" in cache_files:
            cache_filepath = cache_basepath + cache_filename+"_NA"
            taxonomy_dict_cache[species] = pickle_load(cache_filepath)
        else:
            species_to_search.append(species)

    # If there are species which could be searched in NCBI Taxonomy, create a full
    # taxonomy dict which includes 'NOT FOUND' for all non-found organisms and the
    # taxonomies for each found organism.
    if len(species_to_search) > 0:
        ncbi_ids = get_entrez_id_from_organism_full_name_batch(species_to_search)
        taxonomy_dict_search = get_taxonomy_from_organism_ncbi_id_batch(ncbi_ids)

        for searched_species in species_to_search:
            if searched_species not in taxonomy_dict_search.keys():
                taxonomy_dict_search[searched_species] = ["NOT FOUND"]

        for species in list(taxonomy_dict_search.keys()):
            cache_filename = species.replace("/", "") + "_taxonomy"
            if taxonomy_dict_search[species] == ["NOT FOUND"]:
                cache_filename += "_NA"
            cache_filepath = cache_basepath + cache_filename
            pickle_write(cache_filepath, taxonomy_dict_search[species])
        full_taxonomy_dict = {**taxonomy_dict_search, **taxonomy_dict_cache}
    else:
        full_taxonomy_dict = taxonomy_dict_cache

    # Process the taxonomies in order to find the taxonomic distances of the given organism to the organisms
    # which are present in the kcat entries.
    score_dict = most_taxonomic_similar(organism, full_taxonomy_dict)
    for species in full_taxonomy_dict.keys():
        if species not in score_dict:
            score_dict[species] = max(score_dict.values())+1
    # If we added the organism without kcat entries for it, we delete its distance
    # since there is no kcat which can be retrieved from it
    if organism_added:
        del(score_dict[organism])

    # Loop through the given organisms taxonomically and start with the lowest distance
    # Keep looping to higher taxonomic distances until the highest distance is reached
    # or the desired minimal number of kcat entries is reached
    minimal_distance = min(score_dict.values())
    maximal_distance = max(score_dict.values())
    current_distance = minimal_distance
    num_min_kcat_entries = 10
    kcat_list: List[float] = []
    while (len(kcat_list) < num_min_kcat_entries) and (current_distance <= maximal_distance):
        for species in score_dict.keys():
            if species not in species_kcat_mapping.keys():  # e.g. soil bacterium -> bacterium
                continue
            if score_dict[species] == current_distance:
                kcat_list += species_kcat_mapping[species]
        current_distance += 1

    # Get the protein database kcat for this reaction (if there is one, otherwise it returns math.nan)
    if protein_kcat_database != {}:
        protein_database_kcat = _get_kcat_from_protein_kcat_database(searched_direction, reaction, protein_kcat_database)
    else:
        protein_database_kcat = math.nan

    # Add the protein database kcat if there is one, it will influence the resulting kcat since a mean is used
    if protein_database_kcat is not math.nan:
        kcat_list.append(protein_database_kcat)

    return kcat_list


def _get_kcat(searched_metabolites, complete_entry, organism: str, searched_direction: str, reaction: cobra.Reaction,
              protein_kcat_database, type_of_kcat_selection: str = "mean") -> float:
    """Returns the kcat for the reaction with the searched metabolites, the given organism and the given complete entry.

    Arguments
    ----------
    * searched_metabolites ~ The metabolites which were selected as substrates.
    * complete_entry ~ The complete entry of the eligible EC numbers for this reaction.
    * organism: str ~ The analyzed organism's name
    * searched_direction: str ~ The direction in which the
    * reaction: cobra.Reaction ~ The viewed reaction.
    * protein_kcat_database ~ The protein-dependent kcat database that may be added if one exists.
    * type_of_kcat_selection: str ~ The type of kcat selection. Can be "mean", "median" or "random". Is "mean" by default.

    Algorithm
    ----------
    _get_kcat_list is called with the given arguments, which returns a list of eligible kcats. If the length of the eligible
    kcat's list is shorter than 10, a new kcat search without metabolite constriants is performed in order to potentially get
    more kcats.
    With the final kcat list, the mean of the kcats is taken and returned as output.
    """
    # Get the list of all eligible kcats
    kcat_list = _get_kcat_list(searched_metabolites, complete_entry, organism, searched_direction, reaction, protein_kcat_database)

    # If the list is shorter than 10, search without any metabolite constraint in order to potentially get more kcats
    if len(kcat_list) < 10:
        searched_metabolites = ["ALL"]
        kcat_list = _get_kcat_list(searched_metabolites, complete_entry, organism, searched_direction, reaction, protein_kcat_database)

    # Take the eligible kcats using the given selection method
    if type_of_kcat_selection == "mean":
        kcat = statistics.mean(kcat_list)
    elif type_of_kcat_selection == "median":
        kcat = statistics.median(kcat_list)
    elif type_of_kcat_selection == "random":
        kcat = random.choice(kcat_list)
    else:
        print("Wrong type_of_kcat_selection set! Must be 'mean', 'median' or 'random'.")
        sys.exit(-1)

    # Return the determined kcat :D
    return kcat


def _get_searched_metabolites(complete_entry, reaction_part_bigg_ids: List[str]) -> List[str]:
    """Returns which metabolites have a valid BIGG ID that ca be found in the complete entry.

    Arguments
    ----------
    * complete_entry ~ The complete kcat entry of all eligible EC numbers of the analyzed reaction
    * reaction_part_bigg_ids: List[str] ~ A list of BIGG IDs of either the products or educts of a reaction

    Output
    ----------
    A list of BIGG IDs which can be used as identifiers, or "ALL" if no fitting metabolite BIGG ID could be found
    """
    # Go through every metabolite in the complete entry
    eligible_metabolites = []
    for complete_entry_key in complete_entry.keys():
        # Check for identical names
        for reaction_part_bigg_id in reaction_part_bigg_ids:
            if reaction_part_bigg_id == complete_entry_key:
                eligible_metabolites.append(complete_entry_key)

        # Check for names in SABIO-RK style, e.g. "etoh_c;atp;h2o"
        num_found_metabolites = 0
        for reaction_part_bigg_id in reaction_part_bigg_ids:
            if (f";{reaction_part_bigg_id};") in (f";{complete_entry_key};"):
                num_found_metabolites += 1
        if num_found_metabolites == len(complete_entry_key.split(";")):
            eligible_metabolites.append(complete_entry_key)

    # Get the unique set of eligible metabolite IDs
    eligible_metabolites = list(set(eligible_metabolites))

    # If no eligible metabolite ID could be found, set "ALL" as pseudo-metabolite
    # indicating that all kcats are usable
    if len(eligible_metabolites) == 0:
        eligible_metabolites = ["ALL"]

    return eligible_metabolites


def _print_assigned_kcats(reaction_id: str, forward_kcat: float, reverse_kcat: float) -> None:
    """Prints the assigned kcats in the terminal.

    Example
    ----------
    <pre>
    >>>_print_assigned_kcats('Test', 1, 2)
    ***
    Reaction: Test
    Forward kcat: 1
    Reverse kcat: 2
    </pre>

    Arguments
    ----------
    * reaction_id: str ~ The printed reaction ID
    * forward_kcat: float ~ The assigned forward kcat
    * reverse_kcat: float ~ The assigned reverse kcat
    """
    print("***")
    print("Reaction:", reaction_id)
    print("Forward kcat:", forward_kcat)
    print("Reverse kcat:", reverse_kcat)
    print("")


# PUBLIC FUNCTIONS
def get_reactions_kcat_mapping(sbml_path: str, project_folder: str, project_name: str,
                               organism: str, kcat_database_path: str, protein_kcat_database_path: str,
                               type_of_kcat_selection: str = "mean") -> None:
    """Returns a reaction<->kcat mapping for the given model :D

    The selection of kcats is depending on the affected metabolites of the reaction direction (one
    kcat is given for each the forward and reverse direction), and on the organism (the kcats
    from the taxonomically nearest organism is prefered).

    Arguments
    ----------
    *sbml_path: str ~ Te SBML path to the model
    *project_folder: str ~ The folder in which the model data files are sored
    *project_name: str ~ The name of the used project
    *organism: str ~ The organism's name
    *kcat_database_path: str ~ A path to an already created EC number<->kcats database
    *protein_kcat_database_path: str ~ A path to the custom protein<->kcat database
    *type_of_kcat_selection: str ~ Can be "mean", "median" or "random". Refers to the selection of found kcats of a reaction.
                                   Is "mean" by default.

    Output
    ----------
    A JSON in the given project folder with the name $project_name+'_reactions_kcat_mapping_combined.json' and
    the following structure:
    <pre>
    {
        "$REACTION_NAME": {
            "forward": $forward_kcat,
            "reverse": $reverse_kcat
        },
        (...)
    }
    </pre>
    """
    # Standardize project folder
    project_folder = standardize_folder(project_folder)
    # Set the path for the output JSON
    basepath = project_folder + project_name
    # Load the combined, EC-number-dependent kcat database :D
    kcat_database = json_load(kcat_database_path)
    # If given, load the protein-dependent kcat database :D
    if protein_kcat_database_path != "":
        protein_kcat_database = json_load(protein_kcat_database_path)
    else:
        protein_kcat_database = {}

    # Load the given stoichiometric model
    model = cobra.io.read_sbml_model(sbml_path)

    # Set-up dictionary which will be the content of the output JSON
    reactions_kcat_mapping: Dict[str, Dict[str, float]] = {}
    # Go through each reaction in order to assign kcats for it :D
    for reaction in model.reactions:
        # If no EC number is given in the reaction's annotations,
        # the protein-dependent database is read out in order to
        # find a kcat. This only works if at least one of the assigned
        # enzymes of the reaction's gene rule has a kcat in the
        # protein-dependent database.
        if "ec-code" not in reaction.annotation.keys():
            # 0 means that no kcat can be assigned
            forward_kcat: Any = 0
            reverse_kcat: Any = 0
            # Retrieve the kcats from the protein-dependent database :D
            forward_kcat = _get_kcat_from_protein_kcat_database("forward", reaction, protein_kcat_database)
            reverse_kcat = _get_kcat_from_protein_kcat_database("reverse", reaction, protein_kcat_database)

            # If no kcat could be assigned, set the kcat to math.nan
            # which indicates this case
            if forward_kcat == 0.0:
                forward_kcat = math.nan
            if reverse_kcat == 0.0:
                reverse_kcat = math.nan

            # Add the retrieved forward and reverse kcats to the reaction<->kcat mapping dictionary :D
            reactions_kcat_mapping[reaction.id] = {}
            reactions_kcat_mapping[reaction.id]["forward"] = forward_kcat
            reactions_kcat_mapping[reaction.id]["reverse"] = reverse_kcat

            # Print the assigned kcats
            _print_assigned_kcats(reaction.id, forward_kcat, reverse_kcat)
            continue

        # Retrieve the reaction's associated EC numbers
        reaction_ids = reaction.annotation["ec-code"]
        # If only one EC number is given, set the EC number string to
        # a list in order to make it work with the following code lines
        if type(reaction_ids) is str:
            reaction_ids = [reaction_ids]
        # Get all EC numbers which do not contain a - wildcard, such as
        # in 2.1.1.-
        # These wildcarded EC numbers are in general too permissive in order
        # to get useful kcats
        eligible_reaction_ids = [x for x in reaction_ids if "-" not in x]
        if len(eligible_reaction_ids) == 0:
            eligible_reaction_ids = [x for x in reaction_ids]

        # Create a 'complete entry' from all eligible (i.e., non-wildcarded)
        # EC numbers. This complete entry contains - for every organism
        # and substrate given in the EC number kcat entries - all kcats
        # of all eligible EC numbers. In addition, the pseudo-substrate
        # "ALL" is added which contains all organisms. "ALL" is used
        # later if no fitting substrate can be found.
        complete_entry: Dict[str, Any] = {}
        complete_entry["ALL"] = {}
        # Go through each reaction ID :D
        for reaction_id in eligible_reaction_ids:
            # If the EC number could not be found in the given EC number<->kcat
            # database, print it and proceed with the next eligible EC number
            if reaction_id not in kcat_database.keys():
                print(f"INFO: No entry for EC number {reaction_id}")
                print("")
                continue
            # Otherwise, get the reaction ID entry from the given database :D
            reaction_id_entry = kcat_database[reaction_id]
            # Exclude all kcat entries which come from a wildcard search
            # with *
            if reaction_id_entry["WILDCARD"]:
                continue
            # Go trough each metabolite in the EC number<->kcat database entries
            for metabolite_key in reaction_id_entry.keys():
                # Ignore the keys which show additional information
                # about the nature of the kcat data
                if metabolite_key in ("WILDCARD", "SOURCE", "TRANSFER"):
                    continue
                # Add the metabolite to the complete entry if it does not already occur
                if metabolite_key not in complete_entry:
                    complete_entry[metabolite_key] = {}
                # Go throudh each species in the currently analyzed EC number
                for species_key in reaction_id_entry[metabolite_key]:
                    # Add the species to the metabolite entry if it does not already occur
                    if species_key not in complete_entry[metabolite_key]:
                        complete_entry[metabolite_key][species_key] = []
                    # ...and do the same for the pseudo-metabolite "ALL"
                    if species_key not in complete_entry["ALL"].keys():
                        complete_entry["ALL"][species_key] = []
                    # Add the list of kcats of the currently analyzed EC number to the current species
                    # and the current metabolite, and for "ALL"
                    complete_entry[metabolite_key][species_key] += reaction_id_entry[metabolite_key][species_key]
                    complete_entry["ALL"][species_key] += reaction_id_entry[metabolite_key][species_key]

        # If no entries with kcats could be found for any of the eligible EC numbers, continue with the next reaction.
        if complete_entry["ALL"] == {}:
            continue

        # Get the BIGG IDs of the educts and products uusing the SBML's BIGG ID annotation
        educt_bigg_ids: List[str] = []
        for reactant in reaction.reactants:
            if "bigg.metabolite" in reactant.annotation.keys():
                educt_bigg_ids.append(reactant.annotation["bigg.metabolite"])
        product_bigg_ids: List[str] = []
        for product in reaction.products:
            if "bigg.metabolite" in product.annotation.keys():
                product_bigg_ids.append(product.annotation["bigg.metabolite"])
        # If no bigg IDs could be found in the SBML, add the pseudo-metabolite "X"
        # which indicated that "ALL" should be used later.
        if len(educt_bigg_ids) == 0:
            educt_bigg_ids = ["X"]
        if len(product_bigg_ids) == 0:
            product_bigg_ids = ["X"]

        # Get the metabolites which are used in the subsequent forward kcat search
        searched_educts = _get_searched_metabolites(complete_entry, educt_bigg_ids)
        # Get the forward kcat depending on the educts and the organism
        forward_kcat = _get_kcat(searched_educts, complete_entry, organism, "forward", reaction, protein_kcat_database, type_of_kcat_selection)

        # Get the metabolites which are used in the subsequent forward kcat search
        searched_products = _get_searched_metabolites(complete_entry, product_bigg_ids)
        # Get the reverse kcat depending on the products and the organism
        reverse_kcat = _get_kcat(searched_products, complete_entry, organism, "reverse", reaction, protein_kcat_database, type_of_kcat_selection)

        # Set the found out kcats in the reactions<->kcat mapping :D
        reactions_kcat_mapping[reaction.id] = {}
        reactions_kcat_mapping[reaction.id]["forward"] = forward_kcat
        reactions_kcat_mapping[reaction.id]["reverse"] = reverse_kcat

        # display the found out kcats for this reaction \o/
        _print_assigned_kcats(reaction.id, forward_kcat, reverse_kcat)

    # Export the kcat mapping results as JSON :D
    json_write(basepath+"_reactions_kcat_mapping_combined.json", reactions_kcat_mapping)
