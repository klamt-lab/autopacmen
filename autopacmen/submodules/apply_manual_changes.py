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
"""appy_manual_changes.py

Includes the function for the application of a manual kcat change
of a sMOMENT-enhanced model.
"""
# IMPORTS
# External modules
import cobra
from typing import Dict, Tuple, Union


# PUBLIC FUNCTIONS
def apply_manual_changes(model: cobra.Model, kcat_change_factors: Dict[str, Tuple[str, Union[float, int]]]) -> cobra.Model:
    """Applies the given kcat value changes.

    Arguments
    ----------
    * model: cobra.Model ~ The model on which the manual changes shall be applied.
    * kcat_change_factors: Dict[str, Tuple[str, Union[float, int]]] ~ See the following example :-)

    Example of kcat_change_factors
    ----------
    If the kcat of REAC1_TG_forward shall be lowered by the factor 10, and REAC2's kcat
    shall be 10 times higher (REAC2 is not splitted into a forward or reverse reaction):
    <pre>
    kcat_change_factors = {
        "REAC1": ("forward", 1/10),
        "REAC2": ("", 10)
    }
    </pre>

    Output
    ----------
    The given cobra.Model with the given manual changes.
    """
    # Get the protein pool metabolite, whose stoichiometry represents the MW/kcat value
    prot_pool_metabolite = model.metabolites.get_by_id("prot_pool")
    # Go through each reaction ID of the manuel changes :D
    for base_reaction_id in kcat_change_factors.keys():
        change_tuple = kcat_change_factors[base_reaction_id]

        # Get the right reaction if a reaction direction is given
        direction = change_tuple[0]
        if direction != "":
            addition = "_TG_" + direction
        else:
            addition = ""
        searched_reaction_id = base_reaction_id + addition

        # Get the kcat change factor
        # Since the given manual changes are seen with respect to the kcat
        # and the protein pool emtabolite representes MW/kcat, the
        # change factor is 1/kcat_change_factor
        kcat_change_factor = change_tuple[1]
        stoichiometry_change_factor = 1 / kcat_change_factor

        # Calculate the new stoichiometries
        reaction = model.reactions.get_by_id(searched_reaction_id)
        current_stoichiometry = reaction.metabolites[prot_pool_metabolite]
        new_stoichiometry = current_stoichiometry * stoichiometry_change_factor

        # Apply the new stoichimetries :D
        change_dict = {prot_pool_metabolite: -current_stoichiometry}
        reaction.add_metabolites(change_dict)
        change_dict = {prot_pool_metabolite: new_stoichiometry}
        reaction.add_metabolites(change_dict)

        # Display the changes
        print(f"Manual change of {searched_reaction_id}'s kcat, change factor {kcat_change_factor}")

    # Return the changed model :-)
    return model
