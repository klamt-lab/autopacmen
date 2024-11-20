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
"""analysis_fva_prot_pool.py

Command-line interface for the running of a Flux Variability Analysis
(FVA) with given protein bounds on a sMOMENT-enhanced model.
"""

from typing import List

# IMPORTS
# External modules
import click

# Internal modules
from .submodules.fva_prot_pool import fva_prot_pool_with_sbml


# Set-up command-line parameters using click decorators
@click.command()
@click.option("--sbml_path",
              required=True,
              type=click.Path(exists=True, file_okay=True, dir_okay=True, readable=True),
              prompt="SBML path",
              help="Full path to the sMOMENT model SBML")
@click.option("--protein_pool_bounds",
              required=True,
              type=str,
              prompt="Protein pool bounds",
              help="Protein pool bounds for which the FVAs are run, semicolon-separated.")
@click.option("--objective",
              required=False,
              default="",
              type=str,
              prompt="Objective",
              help="Objective of the FVA.")
# Command-line interface function
def fva_prot_pool_with_sbml_cli(sbml_path: str, objective: str, protein_pool_bounds: str) -> None:
    """Runs an FVA (Flux Variability Analysis) and prints the result with the sMOMENT-enhanced model and the given protein pool bounds.

    The FVA results are generated by cobrapy can give a hint how to fit the protein pool.

    Example
    ----------
    Run FVAs with the sMOMENT model 'C:\\model.xml' with the protein bounds 0.5, 0.1 and 0.05 and the objective reaction ACALD:
    <pre>
    python analysis_fva_prot_pool.py --sbml_path C:\\model.xml --protein_pool_bounds 0.5;0.1;0.05 --objective ACALD
    </pre>
    """
    # Parse protein pool bounds given by user as string with each protein pool separated with a semicolon
    protein_pool_bounds_list: List[str] = protein_pool_bounds.split(";")
    protein_pool_bounds_float_list: List[float] = [float(x) for x in protein_pool_bounds_list]
    # Run FVAs :D
    fva_prot_pool_with_sbml(sbml_path, protein_pool_bounds_float_list, objective)


# Start-up routine if script is called
if __name__ == '__main__':
    # Thanks to the click decorators, the command-line interface
    # function does not need to be called directly. The given
    # console arguments are added automatically.
    fva_prot_pool_with_sbml_cli()
    pass
