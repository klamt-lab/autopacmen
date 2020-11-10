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
"""modeling_get_protein_mass_mapping.py

Command-line interface for get_protein_mass_mapping
"""

# IMPORTS
# External modules
import click
import cobra
from typing import Dict, Tuple, Union
# Internal modules
from .submodules.apply_manual_changes import apply_manual_changes


# Set-up console arguments using click decorators
@click.command()
@click.option("--sbml_path_input",
              required=True,
              type=click.Path(exists=True, file_okay=True, dir_okay=True),
              prompt="Input SBML",
              help="Full path to the SBML of which the kcats shall be changed manually")
@click.option("--sbml_path_output",
              required=True,
              type=click.Path(exists=True, dir_okay=True),
              prompt="Output SBML",
              help="Full path to the newly created SBML with the manual kcat changes")
@click.option("--kcat_change_factors",
              required=True,
              type=str,
              prompt="Description of the manual kcat changes:",
              help="Textual description of the manual kcat changes. The format of this description is: . Basic format: basereaction_id,reaction_direction,change_factor;(...) "
                   "Example: if one wants to lower the kcat of the reverse reaction of ACALD by the factor 10, and to highen the kcat of the irreversible reaction CBD by the factor 5, then the textual description of the changes would be: "
                   "ACALD,reverse,1/10;CBD,,5")
# Actual CLI function
def apply_manual_changes_cli(sbml_path_input: str, sbml_path_output: str, kcat_change_factors: str) -> None:
    """Applies manually given kcat changes to the given sMOMENT model.

    This application is given by multiplying the stoichiometry of the protein pool pseudo-metabolite in a reaction with the
    inverse of the given manual change factor. The *inverse* is used as it is intended to change the kcat and as the protein
    pool stoichiometry is MW/kcat (MW is the molecular weight).

    Example
    ----------
    Suppose we want to lower the kcat of the reverse reaction of ACALD by the factor 10, and to highen the kcat of the irreversible reaction CBD by
    the factor 5, and the SBML of the model which shall be changed has the path 'C:\\models\\model.xml' and the newly created SBML with the manual
    changes has the path 'C:\\models\\model_new.xml', then the command would be:
    <pre>
    python optimization_apply_manual_changes --sbml_path_input C:\\models\\model.xml  --sbml_path_output C:\\models\\model_new.xml
    --kcat_change_factors ACALD,reverse,1/10;CBD,,5
    </pre>
    """
    model = cobra.io.read_sbml_model(sbml_path_input)
    kcat_change_factors_dict: Dict[str, Tuple[str, Union[float, int]]] = {}
    single_changes = kcat_change_factors.split(";")
    for single_change in single_changes:
        reaction = single_change.split(",")[0]
        direction = single_change.split(",")[1]
        change_factor = float(single_change.split(",")[2])
        kcat_change_factors_dict[reaction] = (direction, change_factor)
    model = apply_manual_changes(model, kcat_change_factors_dict)
    cobra.io.write_sbml_model(sbml_path_output, sbml_path_output)


# Start-up routine if script is called
if __name__ == '__main__':
    # Thanks to the click decorators, the command-line interface
    # function does not need to be called directly. The given
    # console arguments are added automatically.
    apply_manual_changes_cli()
    pass
