import cobra

biomass_reaction_id = "BIOMASS_Ec_iJO1366_core_53p95M"
substrate_reaction_id = "EX_glc__D_e"
product_reaction_id = "EX_leu__L_e"
min_yield = .2
min_growth = .1

model = cobra.io.read_sbml_model("ec_model_2019_06_25_output_optimization/iJO1366star.xml")

biomass_reaction = model.reactions.get_by_id(biomass_reaction_id)
substrate_reaction = model.reactions.get_by_id(substrate_reaction_id)
product_reaction = model.reactions.get_by_id(product_reaction_id)

model.reactions.get_by_id("EX_o2_e").lower_bound = 0
model.reactions.get_by_id("EX_o2_e").lower_bound = 0
biomass_reaction.lower_bound = min_growth
substrate_reaction.lower_bound = -1000
product_reaction.upper_bound = 1000


# minYield*Substrate - Product <= 0
yield_constraint = model.problem.Constraint(
    -min_yield*substrate_reaction.flux_expression - product_reaction.flux_expression,
    ub=0,
    lb=-1000)
model.add_cons_vars(yield_constraint)

solution = model.optimize()
print(model.summary())
print(-solution.fluxes[substrate_reaction_id])
print(solution.fluxes[product_reaction_id])
print(solution.fluxes[product_reaction_id] / (-solution.fluxes[substrate_reaction_id]))
