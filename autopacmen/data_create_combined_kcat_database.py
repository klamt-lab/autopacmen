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
"""data_create_combined_kcat_database.py

CLI for the construction of a combined (BRENDA+SABIO-RK) kcat database.
"""

# IMPORTS
# External modules
import click
# Internal modules
from .submodules.create_combined_kcat_database import create_combined_kcat_database


# Set-up command-line parameters using click decorators
@click.command()
@click.option("--sabio_rk_kcat_database_path",
              required=True,
              type=click.Path(exists=True, file_okay=True, dir_okay=True, readable=True),
              prompt="SABIO-RK JSON path",
              help="Full path SABIO-RK JSON created with data_parse_brenda_json_for_model.py")
@click.option("--brenda_kcat_database_path",
              required=True,
              type=click.Path(exists=True, file_okay=True, dir_okay=True, readable=True),
              prompt="BRENDA JSON path",
              help="")
@click.option("--output_path",
              required=True,
              type=click.Path(exists=True, file_okay=True, dir_okay=True),
              prompt="Output path",
              help="Full path to the newly created combined JSON")
def parse_create_combined_kcat_database(sabio_rk_kcat_database_path: str, brenda_kcat_database_path: str, output_path: str) -> None:
    """Combines the BRENDA and SABIO-RK JSONs into one big JSON which can be used by modeling_get_reactions_kcat_mapping.py

    The BRENDA JSON is to have been created with data_parse_brenda_json_for_model.py, the SABIO-RK JSON with
    data_parse_sabio_rk_for_model.py.
    Only entries which did not result from wildcard searches will be used for the combined database.

    Example
    ----------
    Combine the BRENDA JSON 'C:\\JSONS\\brenda.json' and the SABIO-RK JSON 'C:\\JSONS\\sabio.json' into
    the JSON 'C:\\JSONS\\combined.json:
    <pre>
    python data_create_combined_kcat_database.py --sabio_rk_kcat_database_path C:\\JSONS\\brenda.json --brenda_kcat_database_path C:\\JSONS\\sabio.json --output_path C:\\JSONS\\sabio.json
    </pre>
    """
    create_combined_kcat_database(sabio_rk_kcat_database_path, brenda_kcat_database_path, output_path)


# Start-up routine if script is called
if __name__ == '__main__':
    # Thanks to the click decorators, the command-line interface
    # function does not need to be called directly. The given
    # console arguments are added automatically.
    parse_create_combined_kcat_database()
    pass
