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

import cobra

gecko = cobra.io.read_sbml_model("ec_model_2019_06_25_output/iJO1366_2019_06_25_GECKO.xml")
analogon = cobra.io.read_sbml_model("ec_model_2019_06_25_output/iJO1366_sMOMENT_2019_06_25_GECKO_ANALOGON.xml")

gecko_reactions_with_arm = []
for reaction in gecko.reactions:
    if "for gene rule" in reaction.name:
        continue
    metabolite_ids = [x.id for x in reaction.metabolites.keys() if x.id.startswith("im_")]
    if len(metabolite_ids) == 1:
        gecko_reactions_with_arm.append(reaction.id.split("_TG_")[0])

analogon_reactions_with_arm = []
for reaction in analogon.reactions:
    if "Arm reaction" in reaction.name:
        continue
    metabolite_ids = [x.id for x in reaction.metabolites.keys() if x.id.startswith("armm_")]
    if len(metabolite_ids) == 1:
        analogon_reactions_with_arm.append(reaction.id.split("_GPRSPLIT")[0])

gecko_reactions_with_arm = set(gecko_reactions_with_arm)
analogon_reactions_with_arm = set(analogon_reactions_with_arm)

print("===STRUCTURAL COMPARISON OF ORIGINAL GECKO AND SMOMENT-BASED GECKO-ANALOGOUS MODEL===")
print("Number of arm reactions - original GECKO:", len(gecko_reactions_with_arm))
print("Number of arm reactions - sMOMENT GECKO analogon:", len(analogon_reactions_with_arm))
difference = analogon_reactions_with_arm - gecko_reactions_with_arm
print("---")
print("Number of reactions - original GECKO: ", len(gecko.reactions))
print("Number of reactions - sMOMENT GECKO analogon: ", len(analogon.reactions))
print("---")
print("Number of metabolites - original GECKO: ", len(gecko.metabolites))
print("Number of metabolites - sMOMENT GECKO analogon: ", len(analogon.metabolites))
