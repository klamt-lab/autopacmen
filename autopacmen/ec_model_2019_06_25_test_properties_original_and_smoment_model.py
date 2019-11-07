import cobra

model_original = cobra.io.read_sbml_model("ec_model_2019_06_25_input/iJO1366.xml")
model_smoment = cobra.io.read_sbml_model("ec_model_2019_06_25_output/iJO1366_sMOMENT_2019_06_25.xml")

print("Original model:")
print(f"Number of reactions is {len(model_original.reactions)}")
print(f"Number of metabolites is {len(model_original.metabolites)}")

print("sMOMENT model:")
print(f"Number of reactions is {len(model_smoment.reactions)}")
print(f"Number of metabolites is {len(model_smoment.metabolites)}")
