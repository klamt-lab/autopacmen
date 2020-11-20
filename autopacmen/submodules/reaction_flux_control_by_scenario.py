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
"""reaction_flux_control_by_scenario.py

Functions for the generation of reaction flux controls for given scenarios.
"""

# IMPORTS
# External modules
import cobra
# Internal modules
from .helper_create_model import apply_scenario_on_model
from .helper_general import standardize_folder


# PRIVATE FUNCTIONS
def _reaction_flux_control(model: cobra.Model, output_folder: str, project_name: str,
                           scenario_key: str, objective: str):
    """Calculate the reaction flux control for the given model.

    Arguments
    ----------
    * model: cobra.Model ~ The model for which the flux control analysis shall be performed.
    * output_folder: str ~ The folder in which the reaction flux control files shall be stored
    * project_name: str ~ The name of the current project
    * scenario_key: str ~ The names of the analyzed scenario
    * objective: str ~ The model's objective

    Output
    ----------
    Tab-separated file with the following content andthe file path
    output_folder+project_name+'reaction_flux_control'+scenario_key+'.txt':
    <pre>
    Line 1  Reaction ID Reaction name   Relative objective flux change
    Line 2  $REACTION_ID    $REACTION_NAME  $RELATIVE_OBJECTIVE_FLUX_CHANGE
    (...)
    </pre>
    """
    base_solution = model.optimize()

    separator = "\t"
    output = f"Reaction ID{separator}"\
             f"Reaction name{separator}"\
             f"Relative objective flux change\n"

    # Go through each reaction and calculate the influence of the deletion of the enzyme constraint:D
    fractions = {}
    prot_pool_metabolite = model.metabolites.get_by_id("prot_pool")
    for reaction in model.reactions:
        # Do not proceed with the protein pool reaction itself (it would be infeasible anyway)
        if reaction.id.startswith("ER_pool_TG_"):
            continue

        # Proceed only if the raction has a protein constraint
        is_with_protein_pool = "prot_pool" in [
            x.id for x in reaction.metabolites]
        if not is_with_protein_pool:
            continue

        # Change the model so that the protein constraint is removed,
        # get the objective solution, and calculate the fraction of this solution
        # to the original solution
        with model:
            prot_pool_stoichiometry = reaction.metabolites[prot_pool_metabolite]
            reaction.subtract_metabolites(
                {prot_pool_metabolite: prot_pool_stoichiometry})
            solution = model.optimize()
            base_solution_flux = base_solution.fluxes[objective]
            changed_solution_flux = solution.fluxes[objective]

            target_result = changed_solution_flux / base_solution_flux
            fractions[reaction.id] = target_result

            output += reaction.id
            output += separator
            output += str(reaction.name)
            output += separator
            output += str(fractions[reaction.id])
            output += "\n"

    # Write the file using the project name and the scenario name
    with open(output_folder+f"{project_name}_reaction_flux_control_{scenario_key}.txt", "w") as f:
        f.write(output)


# PUBLIC FUNCTIONS
def reaction_flux_control_by_scenario(model: cobra.Model, output_folder: str, project_name: str, scenarios):
    """Creates reaction flux control files for the given model in the given folder, and depending on the given scenarios.

    Definition of 'reaction flux control file':
    A reaction flux control file contains a list of reactions and the realative change of the objective
    solution if the protein constraint of the reaction is removed.

    Arguments
    ----------
    * model: cobra.Model ~ The metabolic model for which the reaction flux control file shall be created.
    * output_folder: str ~ The folder in which the reaction flux control files shall be created
    * project_name: str ~ The name of the current project
    * scenarios ~ The scenarios for which the reaction control shall be calculated; these scenarios are
      applied on the model, one by one

    Output
    ----------
    Reaction flux control files in the given folder (see _reaction_flux_control()'s comment for more)
    """
    # Standardize output folder
    output_folder = standardize_folder(output_folder)

    # Go through each given scenario :D
    for scenario_key in scenarios.keys():
        scenario = scenarios[scenario_key]
        objective = scenario["target"]["reaction"]
        with model:
            model = apply_scenario_on_model(model, scenario)
            _reaction_flux_control(model, output_folder,
                                   project_name, scenario_key, objective)


"""
Example:
from ec_model_data_set_up_model import set_up_ec_model_with_sbml
model: cobra.Model = set_up_ec_model_with_sbml("./iJO1366star/ec_model_2019_06_25_output/iJO1366_sMOMENT.xml")
reaction_flux_control_by_c_source(model, OUTPUT_FOLDER, PROJECT_NAME, CHANGE, exchange_reactions, OBJECTIVE)
"""
