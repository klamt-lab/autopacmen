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
"""data_parse_brenda_textfile_for_model.py

CLI for the conversion of a BRENDA database JSON into a model-specific
BRENDA database JSON.
"""

# IMPORTS
# External modules
import click
# Internal modules
from .submodules.parse_brenda_json_for_model import parse_brenda_json_for_model


# Set-up command-line parameters using click decorators
@click.command()
@click.option("--sbml_path",
              required=True,
              type=click.Path(exists=True, file_okay=True,
                              dir_okay=True, readable=True),
              prompt="Path to SBML model",
              help="Full path to the SBML with the model of which the BRENDA JSON will be derived.")
@click.option("--brenda_json_path",
              required=True,
              type=click.Path(exists=True, file_okay=True,
                              dir_okay=True, readable=True),
              prompt="BRENDA JSON path",
              help="Full path to the BRENDA JSON created with data_parse_brenda_textfile.py")
@click.option("--json_output_path",
              required=True,
              type=click.Path(file_okay=True, dir_okay=True, writable=True),
              prompt="JSON output path",
              help="The full path for the model-specific JSON file of the "
                   "BRENDA JSON created with data_parse_brenda_textfile.py")
# Command-line interface function
def parse_brenda_json_for_model_cli(sbml_path: str, brenda_json_path: str, json_output_path: str) -> None:
    """Converts the given BRENDA JSON created with data_parse_brenda_textfile.py into a even more easily readable model-specific JSON.

    This conversion is needed for all subsequent AutoPACMEN steps. The model-specific JSON contains all of the model's EC number entries
    contains a polished version of the BRENDA JSON without transferred EC numbers (i.e., a reassigned EC number) and other possible problems.

    Example
    ----------
    Create a model-specific JSON called 'C:\\folder\\brenda_model_specific.json' from the BRENDA JSON
    'C:\\folder\\brenda.json' using the model 'C:\\folder\\model.xml'
    <pre>
    python data_parse_bigg_metabolites_file.py --sbml_path C:\\folder\\model.xml --brenda_json_path C:\\folder\\brenda.json --json_output_path C:\\folder\\brenda_model_specific.json
    </pre>
    """
    parse_brenda_json_for_model(sbml_path, brenda_json_path, json_output_path)


# Start-up routine if script is called
if __name__ == '__main__':
    # Thanks to the click decorators, the command-line interface
    # function does not need to be called directly. The given
    # console arguments are added automatically.
    parse_brenda_json_for_model_cli()
    pass
