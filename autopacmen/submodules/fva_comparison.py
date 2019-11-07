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
"""fba_comparison.py

This module contains the functions which perform FVA for two cobrapy models.
"""

# IMPORTS
# External modules
import cobra


def fva_comparison(model_original: cobra.Model, model_smoment: cobra.Model, objective: str = "") -> None:
    """Compares the original model with the AutoPACMEN model using FVA for 100% of the objective value.

    The results are printed in the console, and are using cobrapy's FVA summary function.

    Arguments
    ----------
    * model_original: cobra.Model ~ The original SBML model for which the FVA shall be performed.
    * model_smoment: cobra.Model ~ The AutoPACMENed model for which an FVA shall be performed.
    * objective: str = "" ~ An alternative objective for both models. Keep in mind that due to the separation
      of reversible reactions in the sMOMENTed model, there may be differential reaction names for the raction
      that you want to be the objective.
    """
    # Set objective if given for the original model :D
    if objective != "":
        model_original.objective = objective
    # Perform and print FVA for the original model :D
    model_original.optimize()
    print("Original model, FVA solution (for 100% of the maximal objective value):")
    model_original.summary(fva=1.0)

    # Set objective if given for the sMOMENTed model :D
    if objective != "":
        model_smoment.objective = objective
    # Perform and print FVA for the sMOMENTed model :D
    model_smoment.optimize()
    print("\sMOMENT-enhanced model, FVA solution (for 100% of the maximal objective value):")
    model_smoment.summary()


def fva_comparison_with_sbml(sbml_original_path: str, sbml_smomented_path: str, objective: str) -> None:
    """Loads the two given SBMLs using cobrapy, and uses fva_comparison.

    Arguments
    ----------
    * sbml_original_path: str ~ The path to the non-sMOMENTed SBML model.
    * sbml_smomented_path: str ~ The path to the sMOMENT-enhanced SBML model.
    * objective: str = "" ~ An alternative objective for both models. Keep in mind that due to the separation
      of reversible reactions in the sMOMENTed model, there may be differential reaction names for the raction
      that you want to be the objective.
    """
    model_original: cobra.Model = cobra.io.read_sbml_model(sbml_original_path)
    model_smoment: cobra.Model = cobra.io.read_sbml_model(sbml_smomented_path)
    cobra.manipulation.delete_model_genes(model_smoment, ["b1091"])
    fva_comparison(model_original, model_smoment, objective)
