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
"""data_parse_bigg_metabolites_file.py

CLI for the conversion of a BIGG metabolites textfile into a
machine-readable JSON.
"""

# IMPORTS
# External modules
import click
# Internal modules
from .submodules.parse_bigg_metabolites_file import parse_bigg_metabolites_file


# Set-up command-line parameters using click decorators
@click.command()
@click.option("--bigg_metabolites_file_path",
              required=True,
              type=click.Path(exists=True, file_okay=True, dir_okay=True, readable=True),
              prompt="Full path to the BIGG metabolites",
              help="BIGG metabolites file path")
@click.option("--json_output_folder",
              required=True,
              type=click.Path(exists=True, dir_okay=True),
              prompt="JSON output folder",
              help="Path to the folder in which the newly generated JSON will be created")
def parse_brenda_textfile_cli(bigg_metabolites_file_path: str, json_output_folder: str) -> None:
    """Converts the given BIGG metabolites text file into a machine-readable JSON file.

    The BIGG metabolites text file can be downloaded as text file from http://bigg.ucsd.edu/data_access
    (accessed on August 22, 2019). This text file contains all metabolite ID data of the BIGG project.
    The newly created JSON will have the name 'bigg_id_name_mapping.json'.

    Example
    ----------
    Load the BIGG metabolites file 'C:\\bigg_database.txt' and store the newly generated
    JSON in the folder 'C:\\database\\'
    <pre>
    python data_parse_bigg_metabolites_file.py --bigg_metabolites_file_path C:\\bigg_database.txt  --json_output_folder C:\\database
    </pre>
    """
    parse_bigg_metabolites_file(bigg_metabolites_file_path, json_output_folder)


# Start-up routine if script is called
if __name__ == '__main__':
    # Thanks to the click decorators, the command-line interface
    # function does not need to be called directly. The given
    # console arguments are added automatically.
    parse_brenda_textfile_cli()
    pass
