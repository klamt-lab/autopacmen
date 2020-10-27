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
"""data_parse_brenda_textfile.py

CLI for the conversion of a BRENDA database textfile into
a machine-readable JSON.
"""

# IMPORTS
# External modules
import click
# Internal modules
from .submodules.parse_brenda_textfile import parse_brenda_textfile


# Set-up command-line parameters using click decorators
@click.command()
@click.option("--brenda_textfile_path",
              required=True,
              type=click.Path(exists=True, file_okay=True, dir_okay=True, readable=True),
              prompt="Brenda textfile path",
              help="Full path to the BRENDA database text file.")
@click.option("--bigg_metabolites_json_folder",
              required=True,
              type=click.Path(exists=True, dir_okay=True),
              prompt="BIGG metabolites JSON folder",
              help="The folder of the BIGG metabolites JSON (created using "
                   "data_parse_bigg_metabolites_file.py). The name of the JSON "
                   "must not have been changed.")
@click.option("--json_output_path",
              required=True,
              type=click.Path(file_okay=True, dir_okay=True, writable=True),
              prompt="JSON output path",
              help="The full path for the more machine-readable JSON file of the"
                   " BRENDA database text file shall be placed. The resulting JSON "
                   "will have the name. This JSON will be created with this script.")
# Command-line interface function
def parse_brenda_textfile_cli(brenda_textfile_path: str, bigg_metabolites_json_folder: str,
                              json_output_path: str) -> None:
    """Converts the given BRENDA text file into a more machine-readable JSON file.

    The BRENDA database can be downloaded as text file from https://www.brenda-enzymes.org/download_brenda_without_registration.php
    (accessed on August 6, 2019).
    This text file contains all data of BRENDA in an EC-number-sorted manner.
    Since the text file does not have an easily readable standardized format, this
    script converts it into a JSON which contains only the sMOMENT-relevant data in an ordered
    fashion.
    Additionally, substrate names of BRENDA are converted with BIGG metabolite
    names if the name of the BRENDA substrate is equal to a BIGG metabolite name.

    Example
    ----------
    Load the BRENDA textfile 'C:\\folder\\download.txt' while the BIGG metabolites file
    is in 'C:\\bigg\\' and the newly generated JSON shall be named 'C:\\folder\\brenda.json'
    <pre>
    python data_parse_brenda_textfile.py --brenda_textfile_path C:\\folder\\download.txt --bigg_metabolites_json_folder C:\\bigg\\ --json_output_path C:\\folder\\brenda.json
    </pre>
    """
    parse_brenda_textfile(brenda_textfile_path, bigg_metabolites_json_folder, json_output_path)


# Start-up routine if script is called
if __name__ == '__main__':
    # Thanks to the click decorators, the command-line interface
    # function does not need to be called directly. The given
    # console arguments are added automatically.
    parse_brenda_textfile_cli()
    pass
