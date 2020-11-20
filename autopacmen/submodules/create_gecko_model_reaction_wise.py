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
"""modeling_create_gecko_model.py
"""

# IMPORTS
# External modules
import cobra
import copy
import math
import statistics
from typing import List, Dict
# Internal modules
from .helper_general import json_load, standardize_folder
from .helper_create_model import add_prot_pool_reaction, get_irreversible_model, get_p_measured, \
    read_enzyme_stoichiometries_xlsx, read_protein_data_xlsx


# PUBLIC FUNCTIONS
def create_gecko_model_reaction_wise(model: cobra.Model, output_sbml_name: str,
                                     project_folder: str, project_name: str, excluded_reactions: List[str]) -> cobra.Model:
    """Creates a GECKO model as described in

    <i>
    Sánchez, B. J., Zhang, C., Nilsson, A., Lahtvee, P. J., Kerkhoven, E. J., & Nielsen, J. (2017).
    Improving the phenotype predictions of a yeast genome‐scale metabolic model by incorporating enzymatic
    constraints. Molecular systems biology, 13(8).
    </i>

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
    """
    # Standardize project folder
    project_folder = standardize_folder(project_folder)

    # This base path is the location were the generated files wil be stored
    basepath: str = project_folder + project_name

    # READ REACTIONS<->KEGG ID XLSX
    protein_id_mass_mapping: Dict[str, float] = json_load(
        basepath+"_protein_id_mass_mapping.json")

    # LOAD XLSX WITH PROTEIN DATA
    # Load protein data XLSX
    protein_id_concentration_mapping, p_total, unmeasured_protein_fraction, mean_saturation = \
        read_protein_data_xlsx(basepath)
    # Read enzyme kinetics xlsx
    reaction_id_gene_rules_mapping, reaction_id_gene_rules_protein_stoichiometry_mapping = \
        read_enzyme_stoichiometries_xlsx(basepath)

    # Read reaction <-> kcat mapping :D
    reactions_kcat_mapping_database = json_load(
        basepath + "_reactions_kcat_mapping_combined.json")
    all_kcats = [x["forward"] for x in reactions_kcat_mapping_database.values()] + \
                [x["reverse"] for x in reactions_kcat_mapping_database.values()]
    all_kcats = [x for x in all_kcats if not math.isnan(x)]
    default_kcat = statistics.median(all_kcats)
    print(f"Default kcat is: {default_kcat}")

    # GECKO :D #
    # This ID addition will be added to all reactions which are modified by this method
    id_addition = "_TG_"
    # Calculate p_measured
    p_measured = get_p_measured(
        protein_id_concentration_mapping, protein_id_mass_mapping)
    # Make model irreversible
    model = get_irreversible_model(model, id_addition)
    # Add prot_pool reaction
    model, prot_pool_metabolite = add_prot_pool_reaction(model, id_addition, p_total, p_measured,
                                                         unmeasured_protein_fraction, mean_saturation)

    # Add enzyme source reaction for every unmeasured protein
    for protein_id in list(protein_id_mass_mapping.keys()):
        if protein_id in list(protein_id_concentration_mapping.keys()):  # Measured
            eu = cobra.Reaction(id=id_addition+"EU_"+protein_id,
                                name=f"Enzyme usage reaction of measured protein {protein_id}",
                                subsystem="AutoPACMEN")
            enzyme = cobra.Metabolite(id=protein_id+"_met",
                                      name=f"Protein {protein_id}",
                                      compartment="AutoPACMEN")

            eu.add_metabolites({enzyme: 1.0})
            eu.lower_bound = 0
            eu.upper_bound = protein_id_concentration_mapping[protein_id]
            model.add_reactions([eu])
        else:  # Unmeasured
            er = cobra.Reaction(id=id_addition+"ER_"+protein_id,
                                name=f"Enzyme usage reaction of unmeasured protein {protein_id}",
                                subsystem="AutoPACMEN")
            enzyme = cobra.Metabolite(id=protein_id+"_met",
                                      name=f"Protein {protein_id}",
                                      compartment="AutoPACMEN")
            # Mapping is in Da, GECKO uses kDa (g/mmol)
            molecular_weight = protein_id_mass_mapping[protein_id] / 1000
            er.add_metabolites({prot_pool_metabolite: -molecular_weight,
                                enzyme: 1})
            er.lower_bound = 0
            er.upper_bound = 1000.0
            model.add_reactions([er])

    # Add enzymes to reactions
    current_arm_reaction = 1
    model_reaction_ids = [x.id for x in model.reactions]
    for model_reaction_id in model_reaction_ids:
        reaction = model.reactions.get_by_id(model_reaction_id)
        splitted_id = reaction.id.split(id_addition)

        # If the reaction has no name, ignore it
        if splitted_id[0] == "":
            continue
        # Take the reaction ID from the first part of the split
        reaction_id = splitted_id[0]
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

        all_available = True
        for enzyme in gene_rule:
            if type(enzyme) == str:
                try:
                    model.metabolites.get_by_id(enzyme+"_met")
                except Exception:
                    all_available = False
                    break
            else:
                for enzyme_id in enzyme:
                    try:
                        model.metabolites.get_by_id(enzyme_id+"_met")
                    except Exception:
                        all_available = False
                        break
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

        # Add arm reaction if isozymes occur
        if len(gene_rule) > 1:  # Isozymes occur :O
            arm_reaction_id = id_addition + \
                f"arm_reaction_{current_arm_reaction}"
            arm_reaction_name = f"Arm reaction no. {current_arm_reaction} for gene rule {str(gene_rule)}"
            arm_reaction = cobra.Reaction(id=arm_reaction_id,
                                          name=arm_reaction_name,
                                          subsystem="AutoPACMEN")

            arm_reaction_metabolites = {}
            for metabolite in list(reaction.metabolites.keys()):
                stoichiometry = reaction.metabolites[metabolite]
                if stoichiometry < 0:  # Educt
                    arm_reaction_metabolites[metabolite] = stoichiometry
                    reaction.add_metabolites({metabolite: -stoichiometry})

            im_id = f"im_{current_arm_reaction}"
            im_name = f"Intermediate metabolite of arm reaction {current_arm_reaction}"
            intermediate_metabolite = cobra.Metabolite(id=im_id,
                                                       name=im_name,
                                                       compartment="AutoPACMEN")
            arm_reaction_metabolites[intermediate_metabolite] = 1
            arm_reaction.add_metabolites(arm_reaction_metabolites)
            reaction.add_metabolites({intermediate_metabolite: -1})

            arm_reaction.lower_bound = 0
            arm_reaction.upper_bound = reaction.upper_bound

            model.add_reactions([arm_reaction])

            current_arm_reaction += 1

        # Add reactions depending on isozyme complex presence
        new_reactions = []
        i = 1
        for isozyme_id in gene_rule:
            new_reaction = copy.deepcopy(reaction)
            new_reaction.id = new_reaction.id + id_addition + str(i)
            protein_ids = []
            if type(isozyme_id) is str:  # No complex :O
                protein = model.metabolites.get_by_id(isozyme_id+"_met")
                reaction_id = reaction_id.split("_TG_")[0]

                stoichiometry = reaction_id_gene_rules_protein_stoichiometry_mapping[
                    reaction_id][isozyme_id][isozyme_id]
                stoichiometry /= (reaction_kcat * 3600)
                stoichiometry *= -1

                metabolites = {}
                metabolites[protein] = stoichiometry
                protein_ids.append(isozyme_id)
                new_reaction.add_metabolites(metabolites)
            else:  # Complex :O
                metabolites = {}
                isozyme_id = tuple(isozyme_id)

                for single_id in isozyme_id:
                    protein = model.metabolites.get_by_id(single_id+"_met")
                    reaction_id = reaction_id.split("_TG_")[0]

                    stoichiometry = reaction_id_gene_rules_protein_stoichiometry_mapping[
                        reaction_id][isozyme_id][single_id]
                    stoichiometry /= (reaction_kcat * 3600)
                    stoichiometry *= -1

                    metabolites[protein] = stoichiometry
                    protein_ids.append(single_id)
                new_reaction.add_metabolites(metabolites)

            gene_reaction_rule = " and ".join(protein_ids)
            new_reaction.gene_reaction_rule = gene_reaction_rule
            new_reactions.append(new_reaction)
            i += 1
        model.add_reactions(new_reactions)
        model.remove_reactions([reaction])

    cobra.io.write_sbml_model(model, project_folder+output_sbml_name)

    return model


def create_gecko_model_reaction_wise_with_sbml(input_sbml_path: str, output_sbml_name: str,
                                               project_folder: str, project_name: str, excluded_reactions: List[str]) -> cobra.Model:
    """See create_gecko_model_reaction_wise()"""
    # LOAD SBML MODEL
    model: cobra.Model = cobra.io.read_sbml_model(input_sbml_path)
    # Call gecko model creation function :D
    model = create_gecko_model_reaction_wise(model, output_sbml_name,
                                             project_folder, project_name,
                                             excluded_reactions)
    return model
