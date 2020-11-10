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
"""modeling_get_reactions_kcat_mapping.py

Command-line interface for the retrieval of a protein<->kcat mapping of a model.
"""

# IMPORTS
# External modules
import click
# Internal modules
from .submodules.get_reactions_kcat_mapping import get_reactions_kcat_mapping


# Set-up console arguments using click decorators
@click.command()
@click.option("--sbml_path",
              required=True,
              type=click.Path(exists=True, file_okay=True, dir_okay=True),
              prompt="Analyzed SBML",
              help="Full path to the SBML of the metabolic model that shall be sMOMENT-enhanced")
@click.option("--project_folder",
              required=True,
              type=click.Path(exists=True, dir_okay=True),
              prompt="Project folder",
              help="Path to the project folder in which the reactions<->kcat mapping JSON will be created.")
@click.option("--project_name",
              required=True,
              type=str,
              prompt="Project name",
              help="Name of the current project. The generated reactions<->kcat mapping JSON will have this name as prefix.")
@click.option("--organism",
              required=True,
              type=str,
              prompt="Organism name",
              help="The scientific name of the organism that the matabolic model shall represent, e.g. 'Escherichia coli'."
                   "This name is usied, together with NCBI TAXONOMY, for the taxonomy-dependent search of kcat values.")
@click.option("--kcat_database_path",
              required=True,
              type=click.Path(exists=True, file_okay=True, dir_okay=True),
              prompt="kcat database path",
              help="Full path to the SABIO-RK&BRENDA kcat<->reaction JSON which was creates using data_create_combined_kcat_database.py")
@click.option("--protein_kcat_database_path",
              required=True,
              type=click.Path(exists=True, file_okay=True, dir_okay=True),
              default="",
              prompt="protein kcat database path",
              help="Full path to the custom user-defined kcat<->protein JSON. It must be a dictionary containing the protein names"
                   "(as given in the metabolic network's gene rules) as keys and associated kcats as children. See this script's description for more."
                   "Default is '', which means that no user-defined database is given.")
@click.option("--type_of_kcat_selection",
              required=True,
              type=str,
              default="mean",
              prompt='Type of kcat selection, either "mean", "median" or "random".',
              help='Can be "mean", "median" or "random". Refers to the selection of found kcats of a reaction. Is "mean" by default.')
# Actual CLI function
def get_reactions_kcat_mapping_cli(sbml_path: str, project_folder: str, project_name: str,
                                   organism: str, kcat_database_path: str, protein_kcat_database_path: str,
                                   type_of_kcat_selection: str) -> None:
    """Generates the reaction<->kcat mapping for the given metabolic model.

    The algorithm of the kcat selection process is explained in this program package's publication.

    Details on the optional user-defined kcat database:
    It must be a JSON of the following format:
    <pre>
    {
        "$PROTEIN_IDENTIFIER": {
            "kcats": [
                $FLOAT_LIST_OF_ASSOCIATED_KCAT_VALUES
            ],
            "direction": {
                "$ASSOCIATED_REACTION_OF_WHICH_THE_KCATS_WERE_TAKEN": "$DIRECTION_OF_REACTION (either forward or reverse)"
            }
        }
    }
    </pre>
    In other words, this JSON contains protein-dependent measured kcat values, and the reactions and their directions for which
    these kcat values were measured.

    Example
    ----------
    Suppose we want to get the reactions<->kcat mapping of the metabolic model described by the SBML 'C:\\model\\test.xml', the additional
    information from AutoPACMEN's data_*.py scripts is in the project folder 'C:\\project\\', the project name is 'example', the metabolic
    model represents Escherichia coli, the combined reaction<->kcat database is in 'C:\\database\\database.json' and an optional
    user-defined protein database is given and can be found under 'C:\\database\\user.json', and we want to select the kcats
    using the mean, then our command would be:
    <pre>
    python modeling_get_reactions_kcat_mapping.py --sbml_path C:\\model\\test.xml --project_folder C:\\project\\ --project_name example --organism 'Escherichia coli' --kcat_database_path C:\\database\\database.json --protein_kcat_database_path C:\\database\\user.json --type_of_kcat_selection 'mean'
    </pre>
    """
    get_reactions_kcat_mapping(sbml_path, project_folder, project_name, organism,
                               kcat_database_path, protein_kcat_database_path, type_of_kcat_selection)


# Start-up routine if script is called
if __name__ == '__main__':
    # Thanks to the click decorators, the command-line interface
    # function does not need to be called directly. The given
    # console arguments are added automatically.
    get_reactions_kcat_mapping_cli()
    pass
