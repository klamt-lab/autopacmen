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
"""modeling_create_gecko_model.py"""

# IMPORTS
# External modules
import click
from typing import List
# Internal modules
from .submodules.create_gecko_model_reaction_wise import create_gecko_model_reaction_wise_with_sbml


# Set-up console arguments using click decorators
@click.command()
@click.option("--input_sbml_path",
              required=True,
              type=click.Path(exists=True, file_okay=True, dir_okay=True, readable=True),
              prompt="Path to input SBML",
              help="Path to the SBML model which describes the stoichiometric metabolic model "
                   "that shall be converted into a protein-constrained-enhanced model")
@click.option("--project_folder",
              required=True,
              type=click.Path(exists=True, dir_okay=True),
              prompt="Project folder",
              help="Path to project folder. This folder needs to include the reaction<->kcat mapping,"
                   "the protein<->mass mapping and the enzyme stoichiometry spreadsheet (all file names "
                   "have to start with the project name).")
@click.option("--project_name",
              required=True,
              type=str,
              prompt="Project name",
              help="Project name, which represents the start of the project folder's model files")
@click.option("--output_sbml_name",
              required=True,
              type=click.Path(dir_okay=True),
              prompt="Name of output SBML",
              help="File name of the output SBML which describes the protein-constraint-enhanced model. "
                   "This output SBML will be stored in the project folder.")
@click.option("--excluded_reactions",
              required=False,
              type=str,
              prompt="Name of output SBML",
              help="Excluded reactions for which the pseudo-metabolite of the protein pool shall not be introduced. Must be semicolon-separated.")
def create_gecko_model_cli(input_sbml_path: str, output_sbml_name: str,
                           project_folder: str, project_name: str,
                           excluded_reactions: str) -> None:
    """Applies the original GECKO method on the given SBML. All AutoPACMEN scripts starting with "data_" have had to be run first in order to get all necessary data.

    The GECKO method itself is described in:

    Sánchez, B. J., Zhang, C., Nilsson, A., Lahtvee, P. J., Kerkhoven, E. J., & Nielsen, J. (2017).
    Improving the phenotype predictions of a yeast genome‐scale metabolic model by incorporating enzymatic
    constraints. Molecular systems biology, 13(8).

    Example:

    Suppose we want to apply GECKO on the metabolic model described by the SBML 'C:\\models\\model.xml', the GECKO-enhanced model SBML shall have
    the path 'model_new.xml' (it will be stored in the project folder), the project folder containing additional information is 'C:\\project\\',
    the project name is 'example', and we want to exclude the introduction of the protein pool pseudo-metabolite in the reactions 'CBD' and 'ACALD',
    then our command would be:

    python modeling_create_gecko_model.py --input_sbml_path C:\\models\\model.xml --output_sbml_name model_new.xml --project_folder C:\\project\\ --project_name example --excluded_reactions CBD;ACALD
    """
    excluded_reactions_list: List[str] = excluded_reactions.split(";")
    create_gecko_model_reaction_wise_with_sbml(input_sbml_path, output_sbml_name, project_folder, project_name, excluded_reactions_list)


# Start-up routine if script is called
if __name__ == '__main__':
    # Thanks to the click decorators, the command-line interface
    # function does not need to be called directly. The given
    # console arguments are added automatically.
    create_gecko_model_cli()
    pass
