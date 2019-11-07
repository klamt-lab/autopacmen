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
"""ec_model_2019_06_25_modeling_create_gecko_model.py

This script generates the GECKOed model version of iJO1366. The GECKO
method was followed as described in Sánchez et al., 2017.
"""

from submodules.create_gecko_model_reaction_wise import create_gecko_model_reaction_wise_with_sbml

INPUT_SBML: str = "./ec_model_2019_06_25_input/iJO1366.xml"
OUTPUT_SBML: str = "iJO1366_2019_06_25_GECKO.xml"
PROJECT_FOLDER: str = "./ec_model_2019_06_25_output/"
PROJECT_NAME: str = "psb_orth"
EXCLUDED_REACTIONS = ["CO2tex", "O2tex", "H2tex"]

create_gecko_model_reaction_wise_with_sbml(INPUT_SBML, OUTPUT_SBML, PROJECT_FOLDER, PROJECT_NAME, EXCLUDED_REACTIONS)
