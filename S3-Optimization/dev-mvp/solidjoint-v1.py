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
version = "7"
sets = mfunc.create_sets()
sets.keys()

# %%
compat_dfs = mfunc.create_compatability_excel(sets, f"Inputs\CompatabilityMatrices-V{version}.xlsx")
compat_dfs.keys()

# %%
dfs = mfunc.create_parametric_excel(sets, f"Inputs\ParametricTables-V{version}.xlsx", compat=f"Inputs\CompatabilityMatricesInputs-V{version}.xlsx")
dfs.keys()

# %%
dfs = mfunc.read_parametric_excel(sets, f"Inputs\ParametricTablesInputs-V{version}.xlsx", compat=f"Inputs\CompatabilityMatricesInputs-V{version}.xlsx")
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

model.S = pe.Set(initialize=sets['S'])  # Indices for suppliers
model.X = pe.Set(initialize=sets['X'])  # Indices for plant sites
model.M = pe.Set(initialize=sets['M'])  # Indices for markets

model.F = pe.Set(initialize=sets['F'])  # Indices for feedstocks

model.C_Solid = pe.Set(initialize=sets['C_Solid'])  # Indices for solid components
model.C_Liquid = pe.Set(initialize=sets['C_Liquid'])  # Indices for liquid components
model.C = model.C_Solid | model.C_Liquid # Indices for all components

model.P = pe.Set(initialize=sets['P'])  # Indices for final products

model.FP = model.F | model.P

model.BP = pe.Set(initialize=sets['BP'])  # Indices for by-products
model.U = pe.Set(initialize=sets['U'])  # Indices for utilities

model.K_prep = pe.Set(initialize=sets['K_prep'])  # Indices for preprocessing options
model.K_pret = pe.Set(initialize=sets['K_pret'])  # Indices for pretreatment options
model.K_conv = pe.Set(initialize=sets['K_conv'])  # Indices for fermentation options
model.K_pur = pe.Set(initialize=sets['K_pur'])  # Indices for purification options
model.K = model.K_prep | model.K_pret | model.K_conv | model.K_pur # Indices for all process options

model.J = pe.Set(initialize=sets['J'])  # Indices for fermentation strains
model.Q = pe.Set(initialize=sets['Q'])  # Indices for technology capacity levels

# Auxiliary Sets
model.Flow = pe.Set(initialize=["IN", "OUT"])
model.MinMax = pe.Set(initialize=["Min", "Max"])




# %%
####################
# Parameters
####################

## Compatibility_SupplierPlant
parameter_dict = dfs['Compatibility_SupplierPlant'].set_index(['S', 'X']).Value.to_dict()
model.Compatibility_SupplierPlant = pe.Param(model.S, model.X, initialize=parameter_dict, default=0)
# Compatibility_PlantMarket
parameter_dict = dfs['Compatibility_PlantMarket'].set_index(['X', 'M']).Value.to_dict()
model.Compatibility_PlantMarket = pe.Param(model.X, model.M, initialize=parameter_dict, default=0)
# Distances_SupplierPlant
parameter_dict = dfs['Distances_SupplierPlant'].set_index(['S', 'X']).Value.to_dict()
model.Distances_SupplierPlant = pe.Param(model.S, model.X, initialize=parameter_dict, default=0)
# Distances_PlantMarket
parameter_dict = dfs['Distances_PlantMarket'].set_index(['X', 'M']).Value.to_dict()
model.Distances_PlantMarket = pe.Param(model.X, model.M, initialize=parameter_dict, default=0)
# Transportation_Cost_Product
parameter_dict = dfs['Transportation_Cost_Product'].set_index(['P']).Value.to_dict()
model.Transportation_Cost_Product = pe.Param(model.P, initialize=parameter_dict, default=0)
# Transportation_Cost_Feedstock
parameter_dict = dfs['Transportation_Cost_Feedstock'].set_index(['F']).Value.to_dict()
model.Transportation_Cost_Feedstock = pe.Param(model.F, initialize=parameter_dict, default=0)

# Storage_Cost
parameter_dict = dfs['Storage_Cost'].set_index(['FP', 'X']).Value.to_dict()
model.Storage_Cost = pe.Param(model.FP, model.X, initialize=parameter_dict, default=0)

# Compatibility_Purification
parameter_dict = dfs['Compatibility_Purification'].set_index(['P', 'K_pur']).Value.to_dict()
model.Compatibility_Purification = pe.Param(model.P, model.K_pur, initialize=parameter_dict, default=0)

# Compatibility_Conversion
parameter_dict = dfs['Compatibility_Conversion'].set_index(['K_pret', 'K_conv', 'J']).Value.to_dict()
model.Compatibility_Conversion = pe.Param(model.K_pret, model.K_conv, model.J, initialize=parameter_dict, default=0)

# Yield_Hydrolysis
parameter_dict = dfs['Yield_Hydrolysis'].set_index(['K_conv', 'K_pret', 'C', 'P']).Value.to_dict()
model.Yield_Hydrolysis = pe.Param(model.K_conv, model.K_pret, model.C, model.P, initialize=parameter_dict, default=1)

# Solid_Loading
parameter_dict = dfs['Solid_Loading'].set_index(['K_conv']).Value.to_dict()
model.Solid_Loading = pe.Param(model.K_conv, initialize=parameter_dict, default=1)

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
parameter_dict = dfs['Demand'].set_index(['P', 'M', 'T']).Value.to_dict()
model.Demand = pe.Param(model.P, model.M, model.T, initialize=parameter_dict, default=0)
# Supply
parameter_dict = dfs['Supply'].set_index(['F', 'S', 'T']).Value.to_dict()
model.Supply = pe.Param(model.F, model.X, model.T, initialize=parameter_dict, default=0)
# Price
parameter_dict = dfs['Price_F'].set_index(['F', 'S', 'T']).Value.to_dict()
model.Price_F = pe.Param(model.F, model.S, model.T, initialize=parameter_dict, default=0)
parameter_dict = dfs['Price_U'].dropna().set_index(['U', 'X', 'T']).Value.to_dict()
model.Price_U = pe.Param(model.U, model.X, model.T, initialize=parameter_dict, default=0)
parameter_dict = dfs['Price_P'].set_index(['P', 'M', 'T']).Value.to_dict()
model.Price_P = pe.Param(model.P, model.M, model.T, initialize=parameter_dict, default=0)
# Technology_Cost
parameter_dict = dfs['Technology_Cost'].set_index(['K']).Value.to_dict()
model.Technology_Cost = pe.Param(model.K, initialize=parameter_dict, default=0)
parameter_dict = dfs['Technology_Cost'].set_index(['K']).BaseScale.to_dict()
model.Technology_BaseScale = pe.Param(model.K, initialize=parameter_dict, default=0)
parameter_dict = dfs['Technology_Cost'].set_index(['K']).Flow.to_dict()
model.Technology_Flow = pe.Param(model.K, initialize=parameter_dict, default="IN", within = model.Flow)
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

depreciationfactor = (1 - dfs['Economics'][dfs['Economics'].Indicator == 'SalvageValue'].reset_index().Value[0]) / dfs['Economics'][dfs['Economics'].Indicator == 'LifeCycle'].reset_index().Value[0]
model.DepreciatonFactor = pe.Param(initialize=depreciationfactor, default=0)

# Utility_Consumption
parameter_dict = dfs['Utility_Consumption'].set_index(['U', 'K']).Value.to_dict()
model.Utility_Consumption = pe.Param(model.U, model.K, initialize=parameter_dict, default=0)
# Byproduct_Production
parameter_dict = dfs['Byproduct_Production'].set_index(['BP', 'K']).Value.to_dict()
model.Byproduct_Production = pe.Param(model.BP, model.K, initialize=parameter_dict, default=0)



# %%
####################
# Binary Variables
####################

# Indicates if supplier S in chosen for feestock F in plant X i at period t
model.Y_supply = pe.Var(model.S, model.F, model.X,  model.T, initialize=0, domain=pe.Binary)
# Indicates if preprocessing option k_prep in plant x is chosen for feedstock f at period t
model.Y_prep = pe.Var(model.F, model.K_prep, model.X,  model.T, domain=pe.Binary)
# Indicates if pretreatment option k_pret is chosen for feedstock f at period t
model.Y_pret = pe.Var(model.F, model.K_pret, model.X, model.T, domain=pe.Binary)
# Indicates if fermentation option k_conv with strain j is chosen to produce final product p with component c at period t
model.Y_conv = pe.Var(model.F, model.P, model.K_conv, model.X, model.J, model.T, domain=pe.Binary)

# Indicates if purification option k_pur is chosen for product p at period t
model.Y_pur = pe.Var(model.P, model.K_pur, model.X, model.T, domain=pe.Binary)
# Indicates if market P in chosen for product P in plant X i at period t
model.Y_market = pe.Var(model.M, model.P, model.X,  model.T, domain=pe.Binary)

####################
# Continuous Variables
####################

# Inventory (Storage) fo fp at plant x at time t
model.ST = pe.Var(model.FP, model.X, model.T, domain = pe.NonNegativeReals)

# Supplied feedstock f from supplier s in plant x
model.F_in_feed = pe.Var(model.S, model.F, model.X, model.T, domain=pe.NonNegativeReals)

# Quantity of feedstock f going into preprocessing k_prep at period t
model.F_in_prep = pe.Var(model.F, model.K_prep, model.X, model.T, domain=pe.NonNegativeReals)

# Quantity of feedstock f going into pretreatment k_pret at period t
model.F_in_pret = pe.Var(model.F, model.K_pret, model.X, model.T, domain=pe.NonNegativeReals)
# Quantity of component c going out of pretreatment k_pret at period t
model.F_out_pret = pe.Var(model.C, model.F, model.K_pret, model.X, model.T, domain=pe.NonNegativeReals)


# Quantity of component c going into conversion k_conv with fermentation strain j at period t to produce product p
model.F_in_conv = pe.Var(model.C, model.F, model.P, model.K_conv, model.X, model.J, model.T, domain=pe.NonNegativeReals)
# Quantity of component c going into conversion k_conv with fermentation strain j at period t to produce product p
model.F_out_conv = pe.Var(model.P, model.F, model.K_conv, model.X, model.J, model.T, domain=pe.NonNegativeReals)

# Quantity of product p going into purification k_conv at period t
model.F_in_pur = pe.Var(model.P, model.K_pur, model.X, model.T, domain=pe.NonNegativeReals)
# Quantity of product p going out of purification k_pur at period t
model.F_out_pur = pe.Var(model.P, model.K_pur, model.X, model.T, domain=pe.NonNegativeReals)

# Quanitity of product leaving plant to market
model.F_out_prod = pe.Var(model.M, model.P, model.X,  model.T, domain=pe.NonNegativeReals)


# Total flow of utility u consumed in process k at period t
model.F_util = pe.Var(model.U, model.K, model.X, model.T, domain=pe.NonNegativeReals)
# Total flow of byproduct bp generated in process k at period t
model.F_bypr = pe.Var(model.BP, model.K, model.X, model.T, domain=pe.NonNegativeReals)

# %%
####################
# Constraints
####################


# Select only one process per transformation
model.single_technology = pe.ConstraintList()
model.single_prep_constraint = pe.Constraint(model.F, model.X, model.T, rule = lambda m, f, x, t: sum(m.Y_prep[f, k_prep, x, t] for k_prep in model.K_prep) <= 1)
model.single_pret_constraint = pe.Constraint(model.F, model.X, model.T,  rule = lambda m, f, x, t: sum(m.Y_pret[f, k_pret, x, t] for k_pret in model.K_pret) <= 1)
model.single_pur_constraint = pe.Constraint(model.P, model.X, model.T, rule = lambda m, p, x, t: sum(m.Y_pur[p, k_pur, x, t] for k_pur in model.K_pur) <= 1)
#model.single_conv_constraint = pe.Constraint(model.F, model.P, model.K_conv, model.X, model.T, rule = lambda m, f, p, k_conv, x, t: sum(m.Y_conv[f, p, k_conv, x, j, t] for j in model.J) <= 1)



# Choosing Processes

# Big-M: Choosing pathways
M = bigM 

model.supply_disjunction = pe.ConstraintList()
for t in model.T:
    for x in model.X:
        for s in model.S:
            for f in model.F:
                # Corresponding to Y_supply[s, f, x, t] = True in disjunctive formulation
                model.supply_disjunction.add(model.F_in_feed[s, f, x, t]  <= model.Supply[f, s, t] + M * (1 - model.Y_supply[s, f, x, t]))
                # Corresponding to Y_supply[s, f, x, t] = False in  disjunctive formulation
                model.supply_disjunction.add(model.F_in_feed[s, f, x, t] <= M * model.Y_supply[s, f, x, t])

model.bigM_pathways = pe.ConstraintList()
for t in model.T:
    for x in model.X:

        for f in model.F:

            # Choosing preprocessing
            for k_prep in model.K_prep:
                # Corresponding to Y_prep[f, k_prep, t] = True in disjunctive formulation
                if t > 1:
                    model.bigM_pathways.add(model.F_in_prep[f, k_prep, x, t]  <= 
                                            model.ST[f, x, t - 1] + sum(model.F_in_feed[s, f, x, t] for s in model.S)
                                            + M * (1 - model.Y_prep[f, k_prep, x, t]))
                else:
                    model.bigM_pathways.add(model.F_in_prep[f, k_prep, x, t]  <= sum(model.F_in_feed[s, f, x, t] for s in model.S) + M * (1 - model.Y_prep[f, k_prep, x, t]))
                # Corresponding to Y_prep[f, k_prep, t] = False in  disjunctive formulation
                model.bigM_pathways.add(model.F_in_prep[f, k_prep, x, t] <= M * model.Y_prep[f, k_prep, x, t])
            
            # Choosing pretreatment
            for k_pret in model.K_pret:
                    # Corresponding to Y_pret[f, k_pret, t] = True in disjunctive formulation
                    model.bigM_pathways.add(
                        model.F_in_pret[f, k_pret, x, t] <= 
                        sum(model.Compatibility_Pretreatment[k_prep, k_pret] * model.F_in_prep[f, k_prep, x, t] for k_prep in model.K_prep) 
                        + M * (1 - model.Y_pret[f, k_pret, x, t]))
                    # Corresponding to Y_pret[f, k_pret, t] = False in disjunctive formulation
                    model.bigM_pathways.add(
                        model.F_in_pret[f, k_pret, x, t] <= M * model.Y_pret[f, k_pret, x, t])
        

            # Choosing purification
            for k_pur in model.K_pur:
                for p in model.P:
                    # Corresponding to Y_pur[p, k_pur, t] = True in disjunctive formulation
                    model.bigM_pathways.add(
                        model.F_in_pur[p, k_pur, x, t] <= 
                        model.Compatibility_Purification[p, k_pur] * sum(model.F_out_conv[p, f, k, x, j, t] for k in model.K_conv for j in model.J for f in model.F)
                        + M * (1 - model.Y_pur[p, k_pur, x, t]))
                    
                    # Corresponding to Y_pur[p, k_pur, t] = False in disjunctive formulation
                    model.bigM_pathways.add(
                        model.F_in_pur[p, k_pur, x, t] <= M * model.Y_pur[p, k_pur, x, t])
                    
            for m in model.M:
                # Corresponding to Y_supply[s, f, x, t] = True in disjunctive formulation
                if t>1:
                    model.bigM_pathways.add(model.F_out_prod[m, p, x, t]  <= 
                                            sum(model.F_out_pur[p, k_pur, x, t] for k_pur in model.K_pur) + model.ST[p, x, t - 1]
                                            + M * (1 - model.Y_market[m, p, x, t]))
                else:
                    model.bigM_pathways.add(model.F_out_prod[m, p, x, t]  <= 
                                        sum(model.F_out_pur[p, k_pur, x, t] for k_pur in model.K_pur)
                                        + M * (1 - model.Y_market[m, p, x, t]))
                # Corresponding to Y_supply[s, f, x, t] = False in  disjunctive formulation
                model.bigM_pathways.add(model.F_out_prod[m, p, x, t] <= M * model.Y_market[m, p, x, t])


            # Choosing purification
            for k_conv in model.K_conv:
                for p in model.P:
                    for j in model.J:
                        for c in model.C:
                            for f in model.F:
                                # Corresponding to Y_pur[p, k_pur, t] = True in disjunctive formulation
                                model.bigM_pathways.add(
                                    model.F_in_conv[c, f, p, k_conv, x, j, t]  <= 
                                    sum(model.Compatibility_Conversion[k_pret, k_conv, j] * model.F_out_pret[c, f, k_pret, x, t] for k_pret in model.K_pret)
                                    + M * (1 - model.Y_conv[f, p, k_conv, x, j, t]))
                                
                                # Corresponding to Y_pur[p, k_pur, t] = False in disjunctive formulation
                                model.bigM_pathways.add(
                                    model.F_in_conv[c, f, p, k_conv, x, j, t]  <=  M * model.Y_conv[f, p, k_conv, x, j, t])

# Same flow disjunction
def rule_flow_disjunct(disjunct, f, p, k_pret, k_conv, x, j, t, flag):
    m = disjunct.model()
    
    if flag:
        disjunct.c1 = pe.ConstraintList() 
        disjunct.c1.add(m.Y_conv[f, p, k_conv, x, j, t] + m.Y_pret[f, k_pret, x, t] == 2)
        for c in model.C:
            mass_fraction = (m.Feedstock_Composition[c, f] * m.Yield_Pretreatment[c, k_pret]) / sum(m.Feedstock_Composition[ct, f] * m.Yield_Pretreatment[ct, k_pret] for ct in m.C)
            disjunct.c1.add(expr = 
                            m.F_in_conv[c, f, p, k_conv, x, j, t] 
                            == mass_fraction * sum(m.F_in_conv[ct, f, p, k_conv, x, j, t] for ct in model.C))
        
            
    else:
        disjunct.c0 = pe.ConstraintList() 
        disjunct.c0.add(m.Y_conv[f, p, k_conv, x, j, t] + m.Y_pret[f, k_pret, x, t] <= 1)
        #disjunct.c0.add(sum(m.F_in_conv[c, f, p, k_conv, x, j, t] for c in model.C) == 0)

model.flow_disjunct = gdp.Disjunct(model.F, model.P, model.K_pret, model.K_conv, model.X, model.J, model.T, [0, 1], rule = rule_flow_disjunct)

model.flow_disjunction = gdp.Disjunction(model.F, model.P, model.K_pret, model.K_conv, model.X, model.J, model.T, 
                                                rule = lambda m, f, p, k_pret, k_conv, x, j, t: [
                                                    m.flow_disjunct[f, p, k_pret, k_conv, x, j, t, i] for i in [0, 1]
                                                ])



# Mass Balances
model.Mass_Balances = pe.ConstraintList()
model.Mass_Balances_test = pe.ConstraintList()
for t in model.T:

    # Between Processes
    for s in model.S:
        for f in model.F:
            model.Mass_Balances.add(sum(model.F_in_feed[s, f, x, t] for x in model.X) <= M * model.Supply[f, s, t])
            if t> 1:
                model.Mass_Balances.add(sum(model.F_in_prep[f, k_prep, x, t] for k_prep in model.K_prep) <= 
                                        model.ST[f, x, t - 1] + sum(model.F_in_feed[s, f, x, t] for s in model.S))
            else:
                model.Mass_Balances.add(sum(model.F_in_prep[f, k_prep, x, t] for k_prep in model.K_prep) <= 
                                        sum(model.F_in_feed[s, f, x, t] for s in model.S))

    for m in model.M:
        for p in model.P:
            model.Mass_Balances.add(sum(model.F_out_prod[m, p, x, t] for x in model.X) == model.Demand[p, m, t])
            if t>1:
                model.Mass_Balances.add(sum(model.F_out_pur[p, k_pur, x, t] for k_pur in model.K_pur) >= 
                                        (model.ST[p, x, t - 1] + sum(model.F_out_prod[m, p, x, t] for m in model.M)))
            else:
                model.Mass_Balances.add(sum(model.F_out_pur[p, k_pur, x, t] for k_pur in model.K_pur) >= 
                                        sum(model.F_out_prod[m, p, x, t] for m in model.M))
    
    for x in model.X:

        # Storage
        if t > 1:
            for f in model.F:
                model.Mass_Balances.add(model.ST[f, x, t] == model.ST[f, x, t - 1] + sum(model.F_in_feed[s, f, x, t] for s in model.S) - sum(model.F_in_prep[f, k_prep, x, t] for k_prep in model.K_prep))
            for p in model.P:
                model.Mass_Balances.add(model.ST[p, x, t] == model.ST[p, x, t - 1] - sum(model.F_out_prod[m, p, x, t] for m in model.M) + sum(model.F_out_pur[p, k_pur, x, t] for k_pur in model.K_pur))
        else:
            for fp in model.FP:
                model.Mass_Balances.add(model.ST[fp, x, t] == 0)

        for c in model.C:
            for f in model.F:
                model.Mass_Balances.add(sum(model.F_in_conv[c, f, p, k, x, j, t] for p in model.P for k in model.K_conv for j in model.J) <= sum(model.F_out_pret[c, f, k, x, t] for k in model.K_pret))     
        

        # In Processes (Transformations)
        for f in model.F:
            """ for k_prep in model.K_prep:
                model.Mass_Balances.add(model.F_out_prep[f, k_prep, x, t] == model.F_in_prep[f, k_prep, x, t]) """
            for k_pret in model.K_pret:
                for c in model.C:
                    model.Mass_Balances.add(model.F_out_pret[c, f, k_pret, x, t] == model.F_in_pret[f, k_pret, x, t]  * model.Feedstock_Composition[c, f] * model.Yield_Pretreatment[c, k_pret])
                    #model.Mass_Balances.add(model.F_out_pret[c, k_pret, x, t] == sum(model.F_in_pret[f, k_pret, x, t]  * model.Feedstock_Composition[c, f] for f in model.F) * model.Yield_Pretreatment[c, k_pret])
            for p in model.P:
                for k_conv in model.K_conv:
                    for j in model.J:
                        for k_pret in model.K_pret:
                            model.Mass_Balances_test.add(model.F_out_conv[p, f, k_conv, x, j, t]  <= 
                                                         sum(model.F_in_conv[c, f, p, k_conv, x, j, t] *  model.Yield_Conversion[p, c, k_conv, j] * model.Yield_Hydrolysis[k_conv, k_pret, c, p] for c in model.C)
                                                         + M * (1 - model.Y_pret[f, k_pret, x, t])
                                                         )
                        model.Mass_Balances_test.add(model.F_out_conv[p, f, k_conv, x, j, t]  <= 
                                                        M * sum(model.Y_pret[f, k_pret, x, t] for k_pret in model.K_pret)
                                                        )
                for k_pur in model.K_pur:
                    model.Mass_Balances.add(model.F_out_pur[p, k_pur, x, t] == model.F_in_pur[p, k_pur, x, t] * model.Yield_Purification[p, k_pur])


# Utility consumption per processed flow in technology k at t
model.utility_consumption = pe.ConstraintList()
for t in model.T:
    for x in model.X:
        for u in model.U:
            for k in model.K_prep:
                model.utility_consumption.add(model.F_util[u, k, x, t] == sum(model.F_in_prep[f, k, x, t] for f in model.F) * model.Utility_Consumption[u, k])
            for k in model.K_pret:
                if u != "Enzymes":
                    model.utility_consumption.add(model.F_util[u, k, x, t] == sum(model.F_in_pret[f, k, x, t] for f in model.F) * model.Utility_Consumption[u, k] )
                else:
                    model.utility_consumption.add(model.F_util[u, k, x, t] == sum(model.F_out_pret[c, f, k, x, t] for c in ["Cellulose"] for f in model.F) * model.Utility_Consumption[u, k])
            for k in model.K_conv:
                    model.utility_consumption.add(model.F_util[u, k, x, t] == sum(model.F_in_conv[c, f, p, k, x, j, t] for f in model.F for c in model.C for p in model.P for j in model.J)  * model.Utility_Consumption[u, k])
            for k in model.K_pur:
                model.utility_consumption.add(model.F_util[u, k, x, t] == sum(model.F_out_pur[p, k, x, t] for p in model.P) * model.Utility_Consumption[u, k])

# Byproduct Generation per processed flow in technology k at t
model.byproduct_generation = pe.ConstraintList()
for t in model.T:
    for x in model.X:
        for bp in model.BP:
            for k in model.K_prep:
                model.byproduct_generation.add(model.F_bypr[bp, k, x, t] == sum(model.F_in_prep[f, k, x, t] for f in model.F) * model.Byproduct_Production[bp, k])
            for k in model.K_pret:
                model.byproduct_generation.add(model.F_bypr[bp, k, x, t] == sum(model.F_in_pret[f, k, x, t] for f in model.F) * model.Byproduct_Production[bp, k] )
            for k in model.K_conv:
                    model.byproduct_generation.add(model.F_bypr[bp, k, x, t] == sum(model.F_in_conv[c, f, p, k, x, j, t] for f in model.F for c in model.C for p in model.P for j in model.J)  * model.Byproduct_Production[bp, k])
            for k in model.K_pur:
                model.byproduct_generation.add(model.F_bypr[bp, k, x, t] == sum(model.F_in_pur[p, k, x, t] for p in model.P) * model.Byproduct_Production[bp, k])




# %%

####################
# Objecitve
####################

model.Util_Cost = pe.Var(domain = pe.Reals)
model.Util_Cost_calc = pe.Constraint(rule = lambda m: m.Util_Cost == sum(sum(m.F_util[u, k, x, t] for k in model.K) * m.Price_U[u, x, t] for u in model.U for x in m.X for t in m.T) )


model.EBIT_Plant = pe.Var(model.X, model.T, domain = pe.Reals)
model.EBIT_Plant_calc = pe.Constraint(model.X, model.T,
                                rule=lambda m, x, t: 
                                m.EBIT_Plant[x, t] == sum(sum(m.F_out_prod[mk, p, x, t] * m.Price_P[p, mk, t] for p in model.P) for mk in model.M) \
                                - sum(sum(m.F_in_feed[s, f, x, t] * m.Price_F[f, s, t] for f in model.F) for s in model.S) 
                                - sum(sum(m.F_util[u, k, x, t] for k in model.K) * m.Price_U[u, x, t] for u in model.U) 
                                # Transportation
                                #- sum(sum(m.F_in_feed[s, f, x, t] * m.Transportation_Cost_Feedstock[f] for f in model.F) * m.Distances_SupplierPlant[s, x] for s in model.S)
                                #- sum(sum(m.F_out_prod[mk, p, x, t] * m.Transportation_Cost_Product[p] for p in model.P) * m.Distances_PlantMarket[x, mk] for mk in model.M)
                                # Storage
                                #- sum(m.ST[fp, x, t] * m.Storage_Cost[fp, x] for fp in model.FP)
                                # Operating & Maintenance Cost
                                #- (m.Operating_Cost[t] + m.Maintenance_Cost[t])* m.CAPEX[x, t]
                                # Depreciaiton for Tax Shield at period T (all invested capex until t)
                                #- sum(m.CAPEX[x, t_] for t_ in range(1, t + 1)) * model.DepreciatonFactor 
                                )

model.EBIT = pe.Var(model.T, domain = pe.Reals)
model.EBIT_calc = pe.Constraint(model.T,
                                rule=lambda m, t: 
                                m.EBIT[t] == sum(m.EBIT_Plant[x, t] for x in m.X))


def OP_rule(m):      
    return sum( m.EBIT[t] for t in m.T)

model.Objective_function = pe.Objective(rule = OP_rule, sense = pe.maximize) 

####################
# Solve
####################
pe.TransformationFactory('gdp.bigm').apply_to(model, bigM=M)
solver = po.SolverFactory("gurobi", solver_io="python")
solver_parameters = "ResultFile=iismodel.mps" # write an ILP file to print the IIS
solver.options['DualReductions'] = 0
#solver.options['NonConvex'] = 2
results = solver.solve(model, 
                    options_string=solver_parameters
                    )
results.write()

# %%
mfunc.print_vars(model)

# %%
def results_to_df(model):    
    res = {}
    for v in model.component_objects(pe.Var, active=True):
        try:
            varobject = getattr(model, str(v))
            data = {}
            for index in varobject:
                if varobject[index].value is not None:
                    if varobject[index].value != 0:
                        data[index] = varobject[index].value
            res[v] = pd.DataFrame(data)
        except:
            pass
    return results
res = results_to_df(model)
res

# %%
def solve_model():

    return results.solver.wallclock_time

time_list = []
time_list.append(results.solver.wallclock_time)
for i in range(10):
    res = solve_model()
    time_list.append(res)

ds = pd.Series(time_list, name = "JulienRunTime")
ds.to_pickle("runtimes.pickle")

# %%
ds

# %%
ds = pd.read_pickle("runtimes.pickle")
ds

# %%



