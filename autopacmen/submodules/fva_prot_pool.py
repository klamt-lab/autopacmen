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

This module contains functions which allow to perform an FVA of a given
protein-constraint-enhanced model with a given protein pool.
"""

# IMPORTS
# External modules
import cobra
from typing import List


# PUBLIC FUNCTIONS
def fva_prot_pool(model: cobra.Model, pp_upper_bounds: List[float], objective: str = "") -> None:
    """Performs an FVA for the given protein-constraint-enhanced model with the given protein pools.

    Output
    ----------
    Standard cobrapy FVA result summaries for each given protein pool.

    Arguments
    ----------
    * model: cobra.Model ~ The FVAed model as cobrapy Model instance
    * pp_upper_bounds: List[float] ~ The list of upper protein pool bounds.
    * objective: str = "" ~ If set, the objecive will be changed to the given term. Otherwise, the standard
      objective is used.
    """
    for pp_upper_bound in pp_upper_bounds:
        with model:
            if objective != "":
                model.objective = objective
            model.reactions.get_by_id(
                "ER_pool_TG_").upper_bound = pp_upper_bound
            solution = model.optimize()
            print(
                f"\sMOMENT-enhanced model, FVA solution for prot_pool upper bound of {pp_upper_bound}:")
            model.summary(fva=1.0)
            print(abs(solution.fluxes.EX_glc__D_e) /
                  abs(solution.fluxes.EX_ac_e))
            model.metabolites.prot_pool.summary()


def fva_prot_pool_with_sbml(sbml_path: str, pp_upper_bounds: List[float], objective: str) -> None:
    """SBML-loading version of fva_prot_pool(), see its comments for more."""
    model = cobra.io.read_sbml_model(sbml_path)
    fva_prot_pool(model, pp_upper_bounds, objective)
