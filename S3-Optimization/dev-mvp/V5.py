# %%
# !pip install pyomo
# !pip install gurobipy
# !pip install xlsxwriter

# %%
import pandas as pd
import numpy as np
import itertools
import pyomo.environ as pe
import pyomo.gdp as gdp
import gurobipy
import pyomo.opt as po

import modelfunctions as mfunc

# %%
sets = mfunc.create_sets()
sets.keys()

# %%
compat_dfs = mfunc.create_compatability_excel(sets, "..\Inputs\CompatabilityMatrices-V5.xlsx")
dfs = mfunc.create_parametric_excel(sets, "..\Inputs\ParametricTables-V5.xlsx", compat="..\Inputs\CompatabilityMatricesInputs-V5.xlsx")
dfs = mfunc.read_parametric_excel(sets, "..\Inputs\ParametricTablesInputs-V5.xlsx", compat="..\Inputs\CompatabilityMatricesInputs-V5.xlsx")
dfs.keys()

# %%
##################################################################################
# MODEL
##################################################################################


model = pe.ConcreteModel()
bigM = 10**5

####################
# SETS
####################

# Initialize the sets based on the 'sets' dictionary
model.T = pe.Set(initialize=sets['T'])  # Indices for periods
model.F = pe.Set(initialize=sets['F'])  # Indices for feedstocks

model.C_Solid = pe.Set(initialize=sets['C_Solid'])  # Indices for solid components
model.C_Liquid = pe.Set(initialize=sets['C_Liquid'])  # Indices for liquid components
model.C = model.C_Solid | model.C_Liquid # Indices for all components

model.BP = pe.Set(initialize=sets['BP'])  # Indices for by-products
model.U = pe.Set(initialize=sets['U'])  # Indices for utilities

model.K_prep = pe.Set(initialize=sets['K_prep'])  # Indices for preprocessing options
model.K_pret = pe.Set(initialize=sets['K_pret'])  # Indices for pretreatment options
model.K_conv = pe.Set(initialize=sets['K_conv'])  # Indices for fermentation options
model.K_pur = pe.Set(initialize=sets['K_pur'])  # Indices for purification options
model.K = model.K_prep | model.K_pret | model.K_conv | model.K_pur # Indices for all process options

model.J = pe.Set(initialize=sets['J'])  # Indices for fermentation strains
model.P = pe.Set(initialize=sets['P'])  # Indices for final products
model.Q = pe.Set(initialize=sets['Q'])  # Indices for technology capacity levels

# Auxiliary Sets
model.Flow = pe.Set(initialize=["IN", "OUT"])
model.MinMax = pe.Set(initialize=["Min", "Max"])
# Total flow of utility u consumed in process k at period t
model.F_util = pe.Var(model.U, model.K, model.T, domain=pe.NonNegativeReals)
# Total flow of byproduct bp generated in process k at period t
model.F_bypr = pe.Var(model.BP, model.K, model.T, domain=pe.NonNegativeReals)



# %%
####################
# Parameters
####################
# Compatibility_Purification
parameter_dict = dfs['Compatibility_Purification'].set_index(['P', 'K_pur']).Value.to_dict()
model.Compatibility_Purification = pe.Param(model.P, model.K_pur, initialize=parameter_dict, default=0)

# Compatibility_Conversion
parameter_dict = dfs['Compatibility_Conversion'].set_index(['C', 'K_pret', 'K_conv', 'J', 'P']).Value.to_dict()
model.Compatibility_Conversion = pe.Param(model.C, model.K_pret, model.K_conv, model.J, model.P, initialize=parameter_dict, default=0)

# Yield_Hydrolysis
parameter_dict = dfs['Yield_Hydrolysis'].set_index(['K_conv', 'K_pret', 'C', 'P']).Value.to_dict()
model.Yield_Hydrolysis = pe.Param(model.K_conv, model.K_pret, model.C, model.P, initialize=parameter_dict, default=0)

# Compatibility_Pretreatment
parameter_dict = dfs['Compatibility_Pretreatment'].set_index(['K_prep', 'K_pret']).Value.to_dict()
model.Compatibility_Pretreatment = pe.Param(model.K_prep, model.K_pret, initialize=parameter_dict, default=0)

# Feedstock_Composition
parameter_dict = dfs['Feedstock_Composition'].set_index(['C', 'F']).Value.to_dict()
model.Feedstock_Composition = pe.Param(model.C, model.F, initialize=parameter_dict, default=0)
# Yield_Pretreatment
parameter_dict = dfs['Yield_Pretreatment'].set_index(['C', 'K_pret']).Value.to_dict()
model.Yield_Pretreatment = pe.Param(model.C, model.K_pret, initialize=parameter_dict, default=0)
# Yield_Conversion
parameter_dict = dfs['Yield_Conversion'][['P', 'C', 'K_conv', 'J', "Value"]].set_index(['P', 'C', 'K_conv', 'J']).Value.to_dict()
model.Yield_Conversion = pe.Param(model.P, model.C, model.K_conv, model.J, initialize=parameter_dict, default=0)
# Yield_Purification 
parameter_dict = dfs['Yield_Purification'].set_index(['P', 'K_pur']).Value.to_dict()
model.Yield_Purification = pe.Param(model.P, model.K_pur, initialize=parameter_dict, default=0)
# Demand
parameter_dict = dfs['Demand'].set_index(['P', 'T']).Value.to_dict()
model.Demand = pe.Param(model.P, model.T, initialize=parameter_dict, default=0)
# Supply
parameter_dict = dfs['Supply'].set_index(['F', 'T']).Value.to_dict()
model.Supply = pe.Param(model.F, model.T, initialize=parameter_dict, default=0)
# Price
parameter_dict = dfs['Price_F'].set_index(['F', 'T']).Value.to_dict()
model.Price_F = pe.Param(model.F, model.T, initialize=parameter_dict, default=0)
parameter_dict = dfs['Price_U'].set_index(['U', 'T']).Value.to_dict()
model.Price_U = pe.Param(model.U, model.T, initialize=parameter_dict, default=0)
parameter_dict = dfs['Price_P'].set_index(['P', 'T']).Value.to_dict()
model.Price_P = pe.Param(model.P, model.T, initialize=parameter_dict, default=0)
# Technology_Cost
parameter_dict = dfs['Technology_Cost'].set_index(['K']).Value.to_dict()
model.Technology_Cost = pe.Param(model.K, initialize=parameter_dict, default=0)
parameter_dict = dfs['Technology_Cost'].set_index(['K']).ScaleFactor.to_dict()
model.Technology_ScaleFactor = pe.Param(model.K, initialize=parameter_dict, default=0)
parameter_dict = dfs['Technology_Cost'].set_index(['K']).InstalationFactor.to_dict()
model.Technology_InstalationFactor = pe.Param(model.K, initialize=parameter_dict, default=0)
parameter_dict = dfs['CEPCI'].set_index(['T']).Value.to_dict()
model.CEPCI = pe.Param(model.T, initialize=parameter_dict, default=0)
# Operating_Cost (same structure as Depreciation_Factor)
parameter_dict = dfs['Operating_Cost'].set_index(['T']).Value.to_dict()
model.Operating_Cost = pe.Param(model.T, initialize=parameter_dict, default=0)

# Maintenance_Cost (same structure as Depreciation_Factor)
parameter_dict = dfs['Maintenance_Cost'].set_index(['T']).Value.to_dict()
model.Maintenance_Cost = pe.Param(model.T, initialize=parameter_dict, default=0)

# Utility_Consumption
parameter_dict = dfs['Utility_Consumption'].set_index(['U', 'K']).Value.to_dict()
model.Utility_Consumption = pe.Param(model.U, model.K, initialize=parameter_dict, default=0)
# Byproduct_Production
parameter_dict = dfs['Byproduct_Production'].set_index(['BP', 'K']).Value.to_dict()
model.Byproduct_Production = pe.Param(model.BP, model.K, initialize=parameter_dict, default=0)



# %%
####################
# Variables
####################

# Indicates if technology k with capacity q is chosen to be installed at period t
model.Y_tech = pe.Var(model.K, model.Q, model.T, domain=pe.Binary)

# Indicates if preprocessing option k_prep is chosen for feedstock f at period t
model.Y_prep = pe.Var(model.F, model.K_prep, model.T, domain=pe.Binary)
# Indicates if pretreatment option k_pret is chosen for feedstock f at period t
model.Y_pret = pe.Var(model.F, model.K_pret, model.T, domain=pe.Binary)
# Indicates if fermentation option k_conv with strain j is chosen to produce final product p with component c at period t
model.Y_conv = pe.Var(model.C, model.P, model.K_conv, model.J, model.T, domain=pe.Binary)
# Indicates if purification option k_pur is chosen for product p at period t
model.Y_pur = pe.Var(model.P, model.K_pur, model.T, domain=pe.Binary)



# Quantity of feedstock f going into preprocessing k_prep at period t
model.F_in_prep = pe.Var(model.F, model.K_prep, model.T, domain=pe.NonNegativeReals)
# Quantity of feedstock f going out of preprocessing k_prep at period t
model.F_out_prep = pe.Var(model.F, model.K_prep, model.T, domain=pe.NonNegativeReals)

# Quantity of feedstock f going into pretreatment k_pret at period t
model.F_in_pret = pe.Var(model.F, model.K_pret, model.T, domain=pe.NonNegativeReals)
# Quantity of component c going out of pretreatment k_pret at period t
model.F_out_pret = pe.Var(model.C, model.K_pret, model.T, domain=pe.NonNegativeReals)

# Quantity of component c going into conversion k_conv with fermentation strain j at period t to produce product p
model.F_in_conv = pe.Var(model.C, model.P, model.K_conv, model.J, model.T, domain=pe.NonNegativeReals)
# Quantity of product p going out of conversion k_conv with fermentation strain j at period t from component c
model.F_out_conv = pe.Var(model.C, model.P, model.K_conv, model.J, model.T, domain=pe.NonNegativeReals)

# Quantity of product p going into purification k_conv at period t
model.F_in_pur = pe.Var(model.P, model.K_pur, model.T, domain=pe.NonNegativeReals)

# Quantity of product p going out of purification k_pur at period t
model.F_out_pur = pe.Var(model.P, model.K_pur, model.T, domain=pe.NonNegativeReals)
# Capital investment at period t
model.CAPEX = pe.Var(model.T, domain=pe.NonNegativeReals)



# %%
####################
# Constraints
####################

# If a technology is selected at time t, it remains until the end of T. Only one cpacity per technology
model.tech_continuity = pe.ConstraintList()
model.tech_capacity = pe.ConstraintList()
for t in model.T:
    for k in model.K:
        model.tech_capacity.add(sum(model.Y_tech[k, q, t] for q in model.Q) <= 1)
        if t < max(model.T):
            for q in model.Q:
                model.tech_continuity.add(model.Y_tech[k, q, t+1] >= model.Y_tech[k, q, t])

# At least one technology must exist at each time t per stage. Also, thy can only exist if a capacity is installed
model.tech_existance = pe.ConstraintList()
for t in model.T:
    model.tech_existance.add(sum(model.Y_prep[f, k, t] for f in model.F for k in model.K_prep) >= 1)
    model.tech_existance.add(sum(model.Y_pret[f, k, t] for f in model.F for k in model.K_pret) >= 1)
    model.tech_existance.add(sum(model.Y_conv[c, p, k, j, t] for c in model.C for p in model.P for k in model.K_conv for j in model.J) >= 1)
    model.tech_existance.add(sum(model.Y_pur[p, k, t] for p in model.P for k in model.K_pur) >= 1)

    for f in model.F:
        for k in model.K_prep:
            model.tech_existance.add(model.Y_prep[f, k, t] <= sum(model.Y_tech[k, q, t] for q in model.Q))
        for k in model.K_pret:
            model.tech_existance.add(model.Y_pret[f, k, t] <= sum(model.Y_tech[k, q, t] for q in model.Q))
    for p in model.P:
        for c in model.C:
            for j in model.J:
                for k in model.K_conv:
                    model.tech_existance.add(model.Y_conv[c, p, k, j, t] <= sum(model.Y_tech[k, q, t] for q in model.Q))
        for k in model.K_pur:
            model.tech_existance.add(model.Y_pur[p, k, t] <= sum(model.Y_tech[k, q, t] for q in model.Q))

# Calculate CAPEX
model.capex_calculation = pe.ConstraintList()
for t in model.T:
    if t > 1:
        model.capex_calculation.add(model.CAPEX[t] == sum((model.Technology_Cost[k] * q**model.Technology_ScaleFactor[k] * model.Technology_InstalationFactor[k] * model.CEPCI[t]/model.CEPCI[1] * (model.Y_tech[k, q, t] - model.Y_tech[k, q, t - 1])) for k in model.K for q in model.Q))
    else:
        model.capex_calculation.add(model.CAPEX[1] == sum((model.Technology_Cost[k] * q**model.Technology_ScaleFactor[k] * model.Technology_InstalationFactor[k] * model.Y_tech[k, q, 1]) for k in model.K for q in model.Q))

# Select only one process per transformation
model.single_technology = pe.ConstraintList()
model.single_prep_constraint = pe.Constraint(model.F, model.T, rule = lambda m, f, t: sum(m.Y_prep[f, k_prep, t] for k_prep in model.K_prep) <= 1)
model.single_pret_constraint = pe.Constraint(model.F, model.T,  rule = lambda m, f, t: sum(m.Y_pret[f, k_pret, t] for k_pret in model.K_pret) <= 1)
model.single_pur_constraint = pe.Constraint(model.P, model.T, rule = lambda m, p, t: sum(m.Y_pur[p, k_pur, t] for k_pur in model.K_pur) <= 1)
model.single_conv_constraint = pe.Constraint(model.C, model.P, model.K_conv, model.J, model.T, rule = lambda m, c, p, k_conv, j, t: sum(m.Y_conv[c, p, k_conv, j, t] for j in model.J) <= 1)


# Big-M: Choosing pathways
M = bigM 

model.bigM_pathways = pe.ConstraintList()
for t in model.T:
    for f in model.F:

        # Choosing preprocessing
        for k_prep in model.K_prep:
            # Corresponding to Y_prep[f, k_prep, t] = True in disjunctive formulation
            model.bigM_pathways.add(model.F_in_prep[f, k_prep, t]  <= model.Supply[f, t] + M * (1 - model.Y_prep[f, k_prep, t]))
            # Corresponding to Y_prep[f, k_prep, t] = False in  disjunctive formulation
            model.bigM_pathways.add(model.F_in_prep[f, k_prep, t] <= M * model.Y_prep[f, k_prep, t])
        
        # Choosing pretreatment
        for k_pret in model.K_pret:
                # Corresponding to Y_pret[f, k_pret, t] = True in disjunctive formulation
                model.bigM_pathways.add(
                    model.F_in_pret[f, k_pret, t] <= 
                    sum(model.Compatibility_Pretreatment[k_prep, k_pret] * model.F_out_prep[f, k_prep, t] for k_prep in model.K_prep) 
                    + M * (1 - model.Y_pret[f, k_pret, t]))
                # Corresponding to Y_pret[f, k_pret, t] = False in disjunctive formulation
                model.bigM_pathways.add(
                    model.F_in_pret[f, k_pret, t] <= M * model.Y_pret[f, k_pret, t])
    
    for p in model.P:
        # Choosing conversion
        for c in model.C:
            for k_conv in model.K_conv:
                for j in model.J:
                    # Corresponding to Y_conv[c, p, k_conv, j, t] = True in disjunctive formulation
                    model.bigM_pathways.add(
                        model.F_in_conv[c, p, k_conv, j, t] <= 
                        sum(model.Compatibility_Conversion[c, k_pret, k_conv, j, p] * model.F_out_pret[c, k_pret, t] for k_pret in model.K_pret) 
                        + M * (1 - model.Y_conv[c, p, k_conv, j, t]))
            
                    # Corresponding to Y_conv[c, p, k_conv, j, t] = False in disjunctive formulation
                    model.bigM_pathways.add(
                        model.F_in_conv[c, p, k_conv, j, t] <= M * model.Y_conv[c, p, k_conv, j, t])

        # Choosing purification
        for k_pur in model.K_pur:
                # Corresponding to Y_pur[p, k_pur, t] = True in disjunctive formulation
                model.bigM_pathways.add(
                    model.F_in_pur[p, k_pur, t] <= 
                    model.Compatibility_Purification[p, k_pur] * sum(model.F_out_conv[c, p, k, j, t] for c in model.C for k in model.K_conv for j in model.J)
                    + M * (1 - model.Y_pur[p, k_pur, t]))
                
                # Corresponding to Y_pur[p, k_pur, t] = False in disjunctive formulation
                model.bigM_pathways.add(
                    model.F_in_pur[p, k_pur, t] <= M * model.Y_pur[p, k_pur, t])


# Mass Balances
model.Mass_Balances = pe.ConstraintList()
for t in model.T:
    # Between Processes
    #for f in model.F:
        #model.Mass_Balances.add(sum(model.F_in_prep[f, k_prep, t] for k_prep in model.K_prep) <= model.Supply[f, t])
        #model.Mass_Balances.add(sum(model.F_in_pret[f, k, t] for k in model.K_pret) == sum(model.F_out_prep[f, k, t] for k in model.K_prep))
    for c in model.C:
        model.Mass_Balances.add(sum(model.F_in_conv[c, p, k, j, t] for p in model.P for k in model.K_conv for j in model.J) <= sum(model.F_out_pret[c, k, t] for k in model.K_pret))     
    
    for p in model.P:
        #model.Mass_Balances.add(sum(model.F_in_pur[p, k, t] for k in model.K_pur) == sum(model.F_out_conv[c, p, k, j, t] for c in model.C for k in model.K_conv for j in model.J))
        model.Mass_Balances.add(sum(model.F_out_pur[p, k, t] for k in model.K_pur) == model.Demand[p, t])

    # In Processes (Transformations)
    for f in model.F:
        for k_prep in model.K_prep:
            model.Mass_Balances.add(model.F_out_prep[f, k_prep, t]== model.F_in_prep[f, k_prep, t])
    for k_pret in model.K_pret:
        for c in model.C:
            model.Mass_Balances.add(model.F_out_pret[c, k_pret, t] == sum(model.F_in_pret[f, k_pret, t]  * model.Feedstock_Composition[c, f] for f in model.F) * model.Yield_Pretreatment[c, k_pret])
    for p in model.P:
        for c in model.C:
            for j in model.J:
                for k_conv in model.K_conv:
                    model.Mass_Balances.add(model.F_out_conv[c, p, k_conv, j, t] == model.F_in_conv[c, p, k_conv, j, t] * model.Yield_Conversion[p, c, k_conv, j])
        for k_pur in model.K_pur:
            model.Mass_Balances.add(model.F_out_pur[p, k_pur, t] == model.F_in_pur[p, k_pur, t] * model.Yield_Purification[p, k_pur])

# Utility Consumption

# Utility consumption per processed flow in technology k at t
model.utility_consumption = pe.ConstraintList()
for t in model.T:
    for u in model.U:
        for k in model.K_prep:
            model.utility_consumption.add(model.F_util[u, k, t] == sum(model.F_in_prep[f, k, t] for f in model.F) * model.Utility_Consumption[u, k])
        for k in model.K_pret:
            if model.Utility_Consumption[u, k] != "Enzymes":
                model.utility_consumption.add(model.F_util[u, k, t] == sum(model.F_in_pret[f, k, t] for f in model.F) * model.Utility_Consumption[u, k] )
            else:
                model.utility_consumption.add(model.F_util[u, k, t] == sum(model.F_out_pret[c, k, t] for c in model.C) * model.Utility_Consumption[u, k])
        for k in model.K_conv:
                model.utility_consumption.add(model.F_util[u, k, t] == sum(model.F_in_conv[c, p, k, j, t] for c in model.C for p in model.P for j in model.J)  * model.Utility_Consumption[u, k])
        for k in model.K_pur:
            model.utility_consumption.add(model.F_util[u, k, t] == sum(model.F_in_pur[p, k, t] for p in model.P) * model.Utility_Consumption[u, k])

# Byproduct Generation per processed flow in technology k at t
model.byproduct_generation = pe.ConstraintList()
for t in model.T:
    for bp in model.BP:
        for k in model.K_prep:
            model.byproduct_generation.add(model.F_bypr[bp, k, t] == sum(model.F_in_prep[f, k, t] for f in model.F) * model.Byproduct_Production[bp, k])
        for k in model.K_pret:
            model.byproduct_generation.add(model.F_bypr[bp, k, t] == sum(model.F_in_pret[f, k, t] for f in model.F) * model.Byproduct_Production[bp, k] )
        for k in model.K_conv:
                model.byproduct_generation.add(model.F_bypr[bp, k, t] == sum(model.F_in_conv[c, p, k, j, t] for c in model.C for p in model.P for j in model.J)  * model.Byproduct_Production[bp, k])
        for k in model.K_pur:
            model.byproduct_generation.add(model.F_bypr[bp, k, t] == sum(model.F_in_pur[p, k, t] for p in model.P) * model.Byproduct_Production[bp, k])



# %%

####################
# Objecitve
####################
model.EBIT = pe.Var(model.T, domain = pe.Reals)
model.EBIT_calc = pe.Constraint(model.T, 
                                rule=lambda m, t: 
                                m.EBIT[t] == sum(sum(m.F_out_pur[p, k_pur, t] for k_pur in m.K_pur) * m.Price_P[p, t] for p in model.P) \
                                - sum(sum(m.F_in_prep[f, k, t] for k in m.K_prep) * m.Price_F[f, t] for f in model.F) 
                                - sum(sum(m.F_util[u, k, t] for k in model.K) * m.Price_U[u, t] for u in model.U) 
                                #- (m.Operating_Cost[t] + m.Maintenance_Cost[t])* m.CAPEX[t]
                                )

def OP_rule(m):      
    return sum( m.EBIT[t] for t in m.T)
def ServiceLevel_rule(m):      
    return sum( (m.Demand[p, t] - sum(m.F_out_pur[p, k_pur, t] for k_pur in m.K_pur)) for p in m.P for t in m.T)
def NPV_rule(m):      
    return sum( (m.EBIT[t] - m.CAPEX[t]) for t in m.T)

model.Objective_function = pe.Objective(rule = OP_rule, sense = pe.maximize) 

####################
# Solve
####################

solver = po.SolverFactory("gurobi", solver_io="python")
solver_parameters = "ResultFile=iismodel.mps" # write an ILP file to print the IIS
results = solver.solve(model, 
                    options_string=solver_parameters
                    )
results.write()

# %%
mfunc.print_vars(model)


