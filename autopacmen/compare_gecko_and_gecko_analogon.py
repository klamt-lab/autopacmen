import cobra

gecko = cobra.io.read_sbml_model("ec_model_2019_06_25_output/iJO1366_sMOMENT_2019_06_25_GECKO.xml")
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

print(len(gecko_reactions_with_arm))
print(len(analogon_reactions_with_arm))
difference = analogon_reactions_with_arm - gecko_reactions_with_arm
print("A")
