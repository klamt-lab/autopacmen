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
"""modeling_create_initial_spreadsheets.py

Command-line interface for AutoPACMEN's initial spreadsheets function
"""

import click
from .submodules.get_initial_spreadsheets import get_initial_spreadsheets_with_sbml


@click.command()
@click.option("--sbml_path",
              required=True,
              type=click.Path(exists=True, file_okay=True, dir_okay=True),
              prompt="Analyzed SBML",
              help="Path to SBML.")
@click.option("--project_folder",
              required=True,
              type=click.Path(exists=True, dir_okay=True),
              prompt="Project folder",
              help="Path to project.")
@click.option("--project_name",
              required=True,
              type=str,
              prompt="Project name",
              help="Name of project.")
def get_initial_spreadsheets_cli(sbml_path: str, project_folder: str, project_name: str) -> None:
    """Creates initial AutoPACMEN XLSX spreadsheets in which additional information can be entered by the use.

    For the application of sMOMENT, the following spreadsheets are relevant (§ stands for the given project name):
    * §_enzyme_stoichiometries.xlsx: here, the user can enter intra-complex stoichiometries of the proteins
      which constitute the complex
    * §_protein_data: here, the user has to enter the data for the calculation of the protein pool. Optionally,
      the user can also enter proteomic concentration data for each protein.

    Example
    ----------
    Suppose we want to create the initial spreadsheets for the SBML 'C:\\Test.xml' and have the project folder
    'C:\\folder' and the project name 'project':
    <pre>
    python modeling_get_initial_spreadsheets.py --sbml_path C:\\Test.xml --project_folder C:\\folder --project_name project
    </pre>
    """
    get_initial_spreadsheets_with_sbml(sbml_path, project_folder, project_name)


# Start-up routine if script is called
if __name__ == '__main__':
    # Thanks to the click decorators, the command-line interface
    # function does not need to be called directly. The given
    # console arguments are added automatically.
    get_initial_spreadsheets_cli()
    pass
