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
# Internal modules
from .submodules.get_protein_mass_mapping import get_protein_mass_mapping_with_sbml


# Set-up console arguments using click decorators
@click.command()
@click.option("--sbml_path",
              required=True,
              type=click.Path(exists=True, file_okay=True, dir_okay=True),
              prompt="Analyzed SBML",
              help="Path to the SBML representation of the metabolic model of whose gene rules the protein masses shall be read out")
@click.option("--project_folder",
              required=True,
              type=click.Path(exists=True, dir_okay=True),
              prompt="Project folder",
              help="Path to the folder in which the JSON with the protein mass data will be created")
@click.option("--project_name",
              required=True,
              type=str,
              prompt="Project name",
              help="Name of the project. This will be the prefix of the created JSON.")
# Actual CLI function
def get_protein_mass_list_cli(sbml_path: str, project_folder: str, project_name: str) -> None:
    """Creates a JSON with the protein masses for all proteins given in the gene rules of the given metabolic model.

    The protein masses are read out using the UniProt API. Therefore, the gene rules need protein names which can
    be found in UniProt. The newly created JSON will be stored in the given project folder and have the file name
    of the project plus '_protein_id_mass_mapping.json'

    Example
    ----------
    Suppose we want to get the protein masses of the SBML 'C:\\model.xml' and the project folder is 'C:\\folder\\' and the
    project name is 'exemplary', the command would be:
    <pre>
    python modeling_get_protein_mass_mapping.py --sbml_path 'C:\\model.xml' --project_folder 'C:\\folder\\' --project_name 'exemplary'
    </pre>
    """
    get_protein_mass_mapping_with_sbml(sbml_path, project_folder, project_name)


# Start-up routine if script is called
if __name__ == '__main__':
    # Thanks to the click decorators, the command-line interface
    # function does not need to be called directly. The given
    # console arguments are added automatically.
    get_protein_mass_list_cli()
    pass
