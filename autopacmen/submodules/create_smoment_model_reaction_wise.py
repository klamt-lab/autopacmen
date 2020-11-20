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
"""create_smoment_model_reaction_wise.py

Includes the central function which allows to create a proteome-constraint-enhanced
stoichiometric model :D
"""

# IMPORTS
# External modules
import cobra
import math
import random
import statistics
import sys
from typing import Dict, List
# Internal modules
from .helper_general import json_load, standardize_folder
from .helper_create_model import add_prot_pool_reaction, get_irreversible_model, get_p_measured, \
    read_enzyme_stoichiometries_xlsx, read_protein_data_xlsx, \
    get_model_with_separated_measured_enzyme_reactions


# PUBLIC FUNCTIONS
def create_smoment_model_reaction_wise(model: cobra.Model, output_sbml_name: str,
                                       project_folder: str, project_name: str,
                                       excluded_reactions: List[str],
                                       type_of_default_kcat_selection: str = "median") -> None:
    """Adds proteomic constraints according to sMOMENT to the given stoichiometric model and stores it as SBML.

    Arguments
    ----------

    * model: cobra.Model ~ A cobra Model representation of the metabolic network. This model will
      be changed using cobrapy functions in order to add the proteomic constraints.
    * output_sbml_name: str ~ The base name of the created SBML.
    * project_folder: str ~ The folder in which the spreadsheets and JSONs with the model's supplemental
      data can be found.
    * project_name: str ~ The sMOMENTed model creation's name, which will be added at the beginning
      of the created SBML's name.
    * excluded_reactions: List[str] ~ A string list of reaction IDs (the 'reverse' and 'forward'
      name additions must not be added, i.e. for 'ACALD_forward' just 'ACALD' has to be given) to
      which no kcat shall be added. Typically used for gas exchange reactions such as 'CO2tex'.
    * type_of_default_kcat_selection: str ~ The type of selection of default kcat values. Can be "mean",
      "median" or "random". Is "median" by default.

    Output
    ----------
    An SBML in the given folder with the given name, which describes the given stoichiometric model
    enhanced by the protein constraint introduction with this function.
    """
    # Standardize project folder
    project_folder = standardize_folder(project_folder)

    # Set folder path for newly created SBML and name for the reaction ID addition (added at the end,
    # and used in order to have a programatically convinient way to separate additions such as 'reverse'
    # from the 'actual' reaction ID).
    basepath: str = project_folder + project_name
    id_addition: str = "_TG_"

    # READ REACTIONS<->KEGG ID XLSX
    protein_id_mass_mapping: Dict[str, float] = json_load(
        basepath + "_protein_id_mass_mapping.json")

    # Load protein data XLSX
    protein_id_concentration_mapping, p_total, unmeasured_protein_fraction, mean_saturation = \
        read_protein_data_xlsx(basepath)

    # Read enzyme stoichiometries xlsx
    reaction_id_gene_rules_mapping, reaction_id_gene_rules_protein_stoichiometry_mapping = \
        read_enzyme_stoichiometries_xlsx(basepath)

    # Calculate p_measured
    p_measured = get_p_measured(
        protein_id_concentration_mapping, protein_id_mass_mapping)

    # Split reactions with measured enzymes
    model, reaction_id_gene_rules_mapping, reaction_id_gene_rules_protein_stoichiometry_mapping = \
        get_model_with_separated_measured_enzyme_reactions(model,
                                                           protein_id_concentration_mapping,
                                                           reaction_id_gene_rules_mapping,
                                                           reaction_id_gene_rules_protein_stoichiometry_mapping,
                                                           excluded_reactions,
                                                           protein_id_mass_mapping)

    # Make model irreversible, separating all reversible reactions to which a gene rule is given
    # in order to save some reactions.
    model = get_irreversible_model(model, id_addition)

    # Add prot_pool reaction according to the given protein pool values
    model, prot_pool_metabolite = add_prot_pool_reaction(model, id_addition, p_total, p_measured,
                                                         unmeasured_protein_fraction, mean_saturation)

    # Read reaction <-> kcat mapping :-)
    reactions_kcat_mapping_database = json_load(
        basepath + "_reactions_kcat_mapping_combined.json")

    # sMOMENT :D
    # Get all kcats which are not math.nan and calculate the median of them, which will be used as default kcat
    all_kcats = [x["forward"] for x in reactions_kcat_mapping_database.values()] + \
                [x["reverse"] for x in reactions_kcat_mapping_database.values()]
    all_kcats = [x for x in all_kcats if not math.isnan(x)]

    if type_of_default_kcat_selection == "median":
        default_kcat = statistics.median(all_kcats)
    elif type_of_default_kcat_selection == "mean":
        default_kcat = statistics.mean(all_kcats)
    elif type_of_default_kcat_selection == "random":
        default_kcat = random.choice(all_kcats)
    else:
        print('ERROR: Argument type_of_default_kcat_selection must be either "median", "mean" or "random".')
        sys.exit(-1)

    print(f"Default kcat is: {default_kcat}")

    # Get all reaction IDs of the given model
    model_reaction_ids = [x.id for x in model.reactions]

    # Add measured enzyme pseudo-metabolites and pseudo-reactions
    for protein_id in protein_id_concentration_mapping.keys():
        new_metabolite = cobra.Metabolite(id="ENZYME_"+protein_id,
                                          name="Pseudo-metabolite of protein "+protein_id,
                                          compartment="sMOMENT")
        max_protein_concentration = protein_id_concentration_mapping[protein_id]
        new_reaction = cobra.Reaction(id="ENZYME_DELIVERY_"+protein_id,
                                      name="Delivery reaction of pseudo-metabolite "+protein_id,
                                      lower_bound=0,
                                      upper_bound=max_protein_concentration)
        new_reaction.add_metabolites({new_metabolite: 1})
        model.add_reactions([new_reaction])

    # Main loop :D, add enzyme constraints to reactions \o/
    for model_reaction_id in model_reaction_ids:
        # Get the reaction and split the ID at the ID addition
        reaction = model.reactions.get_by_id(model_reaction_id)
        splitted_id = reaction.id.split(id_addition)

        # If the reaction has no name, ignore it
        if splitted_id[0] == "":
            continue
        # Take the reaction ID from the first part of the split
        reaction_id = splitted_id[0]
        # Remove GPRSPLIT name addition from reactions with measured protein concentrations
        if "_GPRSPLIT_" in reaction_id:
            reaction_id = reaction_id.split("_GPRSPLIT_")[0]

        # If the reaction has no associated enzyme stoichiometries, ignore it
        if reaction_id not in list(reaction_id_gene_rules_mapping.keys()):
            continue
        # If the reaction has no gene rule, ignore it
        gene_rule = reaction_id_gene_rules_mapping[reaction_id]
        if gene_rule == [""]:
            continue
        # If the reaction is manually excluded, ignore it
        if reaction_id in excluded_reactions:
            continue

        # Check if all proteins in the reaction's gene rule have a found mass
        # This is not the case for e.g. spontaneous reactions which often get the pseudo-enzyme 's0001'
        all_available = True
        for enzyme in gene_rule:
            if type(enzyme) == str:
                if enzyme not in list(protein_id_mass_mapping.keys()):
                    print(enzyme)
                    all_available = False
                    break
            else:
                for enzyme_id in enzyme:
                    if enzyme_id not in list(protein_id_mass_mapping.keys()):
                        all_available = False
                        break
        # If not all of the mass-checked enzymes have a found mass, ignore this reaction
        if not all_available:
            continue

        # Retrieve the reaction's forward and reverse kcats from the given reaction<->kcat database
        if reaction_id in reactions_kcat_mapping_database.keys():
            forward_kcat = reactions_kcat_mapping_database[reaction_id]["forward"]
            reverse_kcat = reactions_kcat_mapping_database[reaction_id]["reverse"]
        # If the reaction is not in the database, set the default kcat
        else:
            forward_kcat = default_kcat
            reverse_kcat = default_kcat

        # If the given reaction<->kcat database contains math.nan as the reaction's kcat,
        # set the default kcat as math.nan means that no kcat could be found.
        if math.isnan(forward_kcat):
            forward_kcat = default_kcat
        if math.isnan(reverse_kcat):
            reverse_kcat = default_kcat

        # Add the given forward or reverse kcat is the reaction was
        # splitted due to its reversibility.
        # If the reaction is not splitted, add the forward kcat (this
        # is the only possible direction for non-splitted=non-reversible
        # reactions)
        if model_reaction_id.endswith(id_addition + "forward"):
            reaction_kcat = forward_kcat
        elif model_reaction_id.endswith(id_addition + "reverse"):
            reaction_kcat = reverse_kcat
        else:
            reaction_kcat = forward_kcat

        # Add protein pool pseudo-metabolite depending on isozyme complex presence
        # List of selectable MW/kcat stoichiometries (the most conservative constraint will be chosen)
        stoichiometries: List[float] = []
        # List of enzyme names and stoichiometries (semicolon-separated) for a console report
        stoichiometry_enzyme_name_list: List[str] = []
        for isozyme_id in gene_rule:
            # If it's not a complex :O...
            if type(isozyme_id) is str:
                # ...get the reaction ID without the additions...
                reaction_id = reaction_id.split("_TG_")[0]

                # ...get the number of units for this protein...
                number_units = reaction_id_gene_rules_protein_stoichiometry_mapping[
                    reaction_id][isozyme_id][isozyme_id]
                stoichiometry = number_units
                # ...and determine the protein pool stoichiometry by
                # 1) Multiplying the number of units for this protein with its mass (converted from kDa to mDa, since the reaction
                #    flux is defined for mmol/(gDW*h) and not mol/(gDW*h))
                stoichiometry *= (protein_id_mass_mapping[isozyme_id] / 1000)
                # 2) Dividing it with the reaction's kcat (converted from 1/s to 1/h)
                stoichiometry /= (reaction_kcat * 3600)
                # 3) Setting the right direction (educt)
                stoichiometry *= -1
                stoichiometries.append(stoichiometry)
                stoichiometry_enzyme_name_list.append(
                    isozyme_id + ";" + str(number_units))

                # Add proteomics constraints
                if isozyme_id in protein_id_concentration_mapping.keys():
                    enzyme_pseudo_metabolite = model.metabolites.get_by_id(
                        "ENZYME_"+isozyme_id)
                    stoichiometry = reaction_id_gene_rules_protein_stoichiometry_mapping[
                        reaction_id][isozyme_id][isozyme_id]
                    stoichiometry *= 1 / (reaction_kcat * 3600)
                    stoichiometry *= -1
                    reaction.add_metabolites(
                        {enzyme_pseudo_metabolite: stoichiometry})
            # If it is a complex :O...
            else:
                # ...convert the complex IDs to a hashable tuple (used for the stoichiometry selection)...
                isozyme_id = tuple(isozyme_id)
                stoichiometry = 0

                # ...go through each single ID of the complex...
                stoichiometry_enzyme_name_list.append("")
                for single_id in isozyme_id:
                    # ...get the reaction ID without additions...
                    reaction_id = reaction_id.split("_TG_")[0]

                    # ...get the number of units for this protein...
                    number_units = reaction_id_gene_rules_protein_stoichiometry_mapping[
                        reaction_id][isozyme_id][single_id]
                    single_stoichiometry = number_units
                    # ...and determine the protein pool stoichiometry addition by
                    # 1) Multiplying the number of units for this protein with its mass (converted from kDa to Da)
                    single_stoichiometry *= (
                        protein_id_mass_mapping[single_id] / 1000)
                    # 2) Dividing it with the reaction's kcat (converted from 1/s to 1/h)
                    single_stoichiometry /= (reaction_kcat * 3600)
                    # 3) Setting the right direction (educt)
                    single_stoichiometry *= -1
                    # 4) and add it to the complex's stoichiometry
                    stoichiometry += single_stoichiometry
                    # Add name of current single ID
                    stoichiometry_enzyme_name_list[-1] += single_id + \
                        ";" + str(number_units) + " "
                stoichiometry_enzyme_name_list[-1] = stoichiometry_enzyme_name_list[-1].rstrip()
                # Add to list of stoichiometries
                stoichiometries.append(stoichiometry)

                # Add proteomics constraints
                for single_id in isozyme_id:
                    if single_id in protein_id_concentration_mapping.keys():
                        enzyme_pseudo_metabolite = model.metabolites.get_by_id(
                            "ENZYME_"+single_id)
                        stoichiometry = reaction_id_gene_rules_protein_stoichiometry_mapping[
                            reaction_id][isozyme_id][single_id]
                        stoichiometry *= 1 / (reaction_kcat * 3600)
                        stoichiometry *= -1
                        reaction.add_metabolites(
                            {enzyme_pseudo_metabolite: stoichiometry})

        # Take the maximal stoichiometry (i.e., the one with the least cost since this one will usually be prefered
        # anyway in an FBA).
        metabolites = {}
        max_stoichiometry = max(stoichiometries)
        metabolites[prot_pool_metabolite] = max_stoichiometry
        reaction.add_metabolites(metabolites)
        selected_enzyme = stoichiometry_enzyme_name_list[stoichiometries.index(
            max_stoichiometry)]

        # Print report of selected kcat and molecular weight for this reaction
        print("Reaction: ", model_reaction_id)
        print("Selected kcat: ", reaction_kcat)
        print("Selected molecular weight (kDa): ", end="")
        if " " in selected_enzyme:  # Multiple enzymes
            mass_sum = .0
            for single_enzyme in selected_enzyme.split(" "):
                enzyme_name = single_enzyme.split(";")[0]
                enzyme_unit_number = float(single_enzyme.split(";")[1])
                mass_sum += protein_id_mass_mapping[enzyme_name] * \
                    enzyme_unit_number
            print(mass_sum)
        else:  # Single enzyme
            enzyme_name = selected_enzyme.split(";")[0]
            enzyme_unit_number = float(selected_enzyme.split(";")[1])
            print(protein_id_mass_mapping[enzyme_name] * enzyme_unit_number)

    # Output as SBML (without constraints due to cobrapy limitations)
    cobra.io.write_sbml_model(model, project_folder + output_sbml_name)


def create_smoment_model_reaction_wise_with_sbml(input_sbml_path: str, output_sbml_name: str,
                                                 project_folder: str, project_name: str,
                                                 excluded_reactions: List[str],
                                                 type_of_default_kcat_selection: str = "median") -> None:
    """See this module's create_smoment_model_reaction_wise()"""
    # Load SBML model
    model: cobra.Model = cobra.io.read_sbml_model(input_sbml_path)
    # Call gecko model creation function :D
    create_smoment_model_reaction_wise(model, output_sbml_name,
                                       project_folder, project_name,
                                       excluded_reactions,
                                       type_of_default_kcat_selection)
