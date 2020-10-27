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
"""get_differential_reactions.py

This module contains the functions for getting reactions whose protein-constraint
deletion has an influence on the result, depending on the given model scenarios.
"""

# IMPORTS
# External modules
import copy
from typing import Dict, List


#  PRIVATE FUNCTIONS
def _get_differential_reactions_from_flux_control_file(filepath: str, threshold: float) -> List[str]:
    """Returns a list of differential reactions of the flux control file in accordance with the given threshold.

    Arguments
    ----------
    * filepath: str ~ The path to the flux control file which shall be analyzed
    * threshold: float ~ The differential reaction threshold
    """
    # Get the file's lines as string list without newlines
    with open(filepath, "r") as f:
        lines = f.readlines()
    lines = [x.replace("\n", "") for x in lines[1:] if len(x) > 0]

    # Get the single differential reactions
    differential_reactions = []
    for line in lines:
        line_split = line.split("\t")
        change = float(line_split[2])

        if change >= (1.0 + threshold):
            protein_id = line_split[0]
            differential_reactions.append(protein_id)

    return differential_reactions


# PUBLIC FUNCTIONS
def get_all_differential_reactions(scenario_names: List[str], flux_control_files_path: str, project_name: str,
                                   threshold: float = (.1)/100):
    """Returns the set of all differential reactions (with the given threshold) of all given scenarios.

    Arguments
    ----------
    * scenario_names: str ~ The names of the scenarios for which the differential reactions shall be collected.
    * flux_control_files_path: str ~ The path to the flux control files.
    * project_name: str ~ The name of the project, which is at the beginning of the file name of the flux files.
    * threshold: float = (.1)/100 ~ The threshold for a change in the solution of a reaction.
    """
    for c_source in scenario_names:
        filepath = f"{flux_control_files_path}/{project_name}_reaction_flux_control_{c_source}.txt"
        differential_reactions = _get_differential_reactions_from_flux_control_file(
            filepath, threshold)
    differential_reactions = list(set(differential_reactions))
    return differential_reactions


def get_differential_reactions(scenario_names: List[str], flux_control_files_path: str, project_name: str,
                               scenarios, threshold: float = (.1)/100, print_result: bool = True):
    """Returns the differential reactions in the given flux files.

    Definition of 'differential reaction'
    ----------
    In a protein-constraint-enhanced metabolic network, the deletion of the protein constraint in a reaction
    can have an influence on the objective solution value, or not. A 'differential reaction' is a reaction in
    which the constraint's deletion has an influence over the given threshold. I.e., if the original objective solution
    is 1.0 and - after the deletion of the protein constraint for the reaction - again 1.0, the reaction is not
    differential. If the deletion of the protein constraint leads e.g. to the solution 1.001 or .999, and the threshold
    is smaller or equal to .001, the reaction is 'differential'.

    Arguments
    ----------
    * scenario_names: str ~ The names of the scenarios for which flux files were created earlier.
    * flux_control_files_path: str ~ The path in which the flux control files can be found.
    * project_name: str ~ The name of the current project (the flux control file names start with it).
    * threshold: float = (.1)/100 ~ The differential reaction threshold
    * print_result: bool = True ~ Whether or not a verbose text ouput of the found differential reactions
      shall be printed.

    Output
    ----------
    2 values:
    * unique_differential_proteins_of_scenario: Dict[str, List[str]] ~ A list of all differential
      reactions which occur in *only one* given scenario.
    * differential_reactions_of_all_scenarios: List[str] ~ A list of all differential reactions which
      occur in *all* given scenarios.
    """
    # Get the differential reactions of each single scenario
    differential_reactions_by_scenario: Dict[str, List[str]] = {}
    for scenario_name in scenario_names:
        filepath = f"{flux_control_files_path}/{project_name}_reaction_flux_control_{scenario_name}.txt"
        differential_reactions = _get_differential_reactions_from_flux_control_file(
            filepath, threshold)
        differential_reactions_by_scenario[scenario_name] = differential_reactions

    # Check for substitution name
    # A 'substitution name' is used in order to combine the differential reactions of multiple scenarios
    # If the substitution_name for two scenarios is the same, they are combined
    # A scenario does not need a substitution name
    substitutions = []
    # Get the substitution rules
    for scenario_name in scenario_names:
        if "substitution_name" not in scenarios[scenario_name].keys():
            continue
        substitution_name = scenarios[scenario_name]["substitution_name"]
        substitutions.append((scenario_name, substitution_name))
    # Combine the differential reactions of the scenarios with the same substitution names
    for substitution in substitutions:
        old_name = substitution[0]
        new_name = substitution[1]

        old_index = scenario_names.index(old_name)
        del(scenario_names[old_index])

        if new_name in differential_reactions_by_scenario.keys():
            combined_list = list(set(
                differential_reactions_by_scenario[new_name] + differential_reactions_by_scenario[old_name]))
            differential_reactions_by_scenario[new_name] = copy.deepcopy(
                combined_list)
        else:
            scenario_names.append(new_name)
            differential_reactions_by_scenario[new_name] = copy.deepcopy(
                differential_reactions_by_scenario[old_name])
        del(differential_reactions_by_scenario[old_name])

    # Get the unique differential reactions by checking each differential reaction of a scenario
    # with all other ones
    all_differential_proteins = []
    unique_differential_reactions_of_scenario: Dict[str, List[str]] = {}
    for scenario_name in scenario_names:
        differential_reactions_of_scenario = differential_reactions_by_scenario[scenario_name]
        all_differential_proteins.append(
            set(differential_reactions_of_scenario))

        other_differential_reactions: List[str] = []
        for other_c_source in scenario_names:
            if other_c_source != scenario_name:
                other_differential_reactions += differential_reactions_by_scenario[other_c_source]
        other_differential_reactions = list(set(other_differential_reactions))
        unique_differential_reactions_of_single_scenario = [x for x in differential_reactions_of_scenario
                                                            if x not in other_differential_reactions]
        unique_differential_reactions_of_scenario[scenario_name] = unique_differential_reactions_of_single_scenario

        if print_result:
            print("Unique differential reactions of "+scenario_name+":")
            print(unique_differential_reactions_of_single_scenario)

    differential_reactions_of_all_scenarios = set.intersection(
        *all_differential_proteins)
    if print_result:
        print("Differential reactions of all C sources:")
        print(differential_reactions_of_all_scenarios)

    return unique_differential_reactions_of_scenario, differential_reactions_of_all_scenarios
