# %%
import pandas as pd
import numpy as np
import itertools
import pyomo.environ as pe
import pyomo.gdp as gdp
import gurobipy
import pyomo.opt as po
import pickle
import cloudpickle

import modelfunctions as mfunc

# %%
sets = {
    "S"     : ["BE"],  # Indices for supply geogpraphies
    "M"     : ["DE"],  # Indices for markets
    "X"     : ["DE"],  # Indices for plant site locations
    "T"       : [1,3,5,7,9], #list(forecasts["Supply"]["T"].unique()),  # Indices for periods (10 years)
    "F"       : ["Wheat", "Maize", "Beet"],#list(forecasts["Supply"]["P"].unique()),  # Indices for feedstocks
    "C"       : [],  # Indices for components
    "C_Solid" : ["Cellulose", "Hemicellulose", "Lignin"],  # Indices for solid components
    "C_Liquid": ["Xylose", "Glucose", "Lignin_Liquid"],  # Indices for liquid components
    
    "BP"      : ['Succinic Acid','Glycolic Acid','Formic Acid','Acetic Acid','Phenolic','Furfural','HMF'
                ],  # Indices for by-products
    
    "RM"      : ['H2SO4','Water','Ammonia','Ionic Liquid mat.','Lime mat.'],
    "E"       : [
                    'Cooling water',
                    'Electricity',
                    'Heat',
                    'LP-Steam',
                    'Refrigeration -15,5°C',
                    'RNG',
                    'Cooling',
    ], 
    "UT"       : [
                    'Calcium hydroxide',
                    'Cyclohexane',
                    'Ethylene Glycol',
                    'Gypsum disposal',
                    'HCL',
                    'Iron (II) chloride',
                    'MeOH',
                    'Monoethanolamine',
                    'Steam 150 psig',
                    'Tri-n-octylamine (0,25 mol/kg) in 1-octanol',
                    'Waste to landfill',
                    'Olivine',
                    'Magnesium oxide',
                    'Enzymes',
                    'Transport'
                ],  # Indices for utilities
    
    "K"       : [],  # Indices for process options
    "K_prep"  : ['MILL3mm','MILL1.4mm','MILL10mm','MILL0.853mm','MILL0.15mm'],  # Indices for preprocessing options
    
    "K_pret"  : ['Dilute Acid','LHW','AFEX','Steam Explosion','Lime','Ionic Liquid'],  # Indices for pretreatment options
    
    "K_conv_F" : ['SHF','SSF'],
    "K_conv_CF": ['SHCF','SSCF'],
    "K_conv_O"  : [],  # Indices for fermentation options
    
    "K_pur"   : ["Reactive Extraction", 'Precipitation', 'Electrodialysis', 'Direct Crystalization', 'Azeotropic Dist', 'Extraction', 'Absorption', "Lignin Recovery"],  # Indices for purification options
    
    "J"       : ["Strain1","Strain2","Strain3","Strain4","Strain5","Strain6 ","Strain7","Strain8","Strain9" 
                #"8G-3", "SPT2OE", "Y2034", "Angle", "CRD51", "MD4", "TP50", "SR8/TFP3-XI76", "SR8/TFP3-XI76"
                ],  # Indices for fermentation strains
    
    "P"       : ['Ethanol', 'Lactic Acid', 'Succinic Acid', "Lignin"],  # Indices for final products
    
    "Q"       : [0.3, 0.5, 0.8, 1, 1.5, 2, 5, 10, 50, 100]  # Indices for technology capacity levels
}


sets["C"] = sets["C_Solid"]  + sets["C_Liquid"] 
sets["U"] = sets["RM"]  + sets["E"] + sets["UT"] 
sets["K_conv"] = sets["K_conv_F"]  + sets["K_conv_CF"] + sets["K_conv_O"]
sets["K"] = sets["K_prep"]  + sets["K_pret"] + sets["K_conv"] + sets["K_pur"]
sets["FP"] = sets["F"]  + sets["P"] 

# %%
version = "2027"
dfs = mfunc.import_paremeters(sets, f"Inputs/inputs-{version}.xlsx")
dfs.keys()

# %%
def create_n_solve(max_LCA):
    ##################################################################################
    # MODEL
    ##################################################################################


    model = pe.ConcreteModel()
    bigM = 10**6

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
    model.RM = pe.Set(initialize=sets['RM'])  # Indices for raw materials
    model.E = pe.Set(initialize=sets['E'])  # Indices for energy consupmtion
    model.UT = pe.Set(initialize=sets['UT'])  # Indices for other utilities
    model.U = model.RM | model.E | model.UT  # Indices for all utilities

    model.K_prep = pe.Set(initialize=sets['K_prep'])  # Indices for preprocessing options
    model.K_pret = pe.Set(initialize=sets['K_pret'])  # Indices for pretreatment options
    model.K_conv_F = pe.Set(initialize=sets['K_conv_F']) #Indices for fermentation options
    model.K_conv_CF = pe.Set(initialize=sets['K_conv_CF']) #Indices for co-fermentation options
    model.K_conv_O = pe.Set(initialize=sets['K_conv_O']) #Indices for other conversion options
    model.K_conv =  model.K_conv_F | model.K_conv_CF | model.K_conv_O  # Indices for all conversion options
    model.K_pur = pe.Set(initialize=sets['K_pur'])  # Indices for purification options
    model.K = model.K_prep | model.K_pret | model.K_conv | model.K_pur # Indices for all process options

    model.J = pe.Set(initialize=sets['J'])  # Indices for fermentation strains
    model.Q = pe.Set(initialize=sets['Q'])  # Indices for technology capacity levels

    # Auxiliary Sets
    model.Flow = pe.Set(initialize=["IN", "OUT"])
    model.MinMax = pe.Set(initialize=["Min", "Max"])

    ####################
    # Parameters
    ####################

    """
    ## Compatibility_SupplierPlant
    parameter_dict = dfs['Compatibility_SupplierPlant'].set_index(['S', 'X']).Value.dropna().to_dict()
    model.Compatibility_SupplierPlant = pe.Param(model.S, model.X, initialize=parameter_dict, default=0)
    # Compatibility_PlantMarket
    parameter_dict = dfs['Compatibility_PlantMarket'].set_index(['X', 'M']).Value.dropna().to_dict()
    model.Compatibility_PlantMarket = pe.Param(model.X, model.M, initialize=parameter_dict, default=0)
    """

    # Transportation_Cost
    parameter_dict = dfs['Economics'][dfs['Economics'].Indicator == 'TransportFactorKm'].reset_index().Value[0]
    model.TransportFactorKm = pe.Param(initialize=parameter_dict, default=0.4)
    parameter_dict = dfs['Economics'][dfs['Economics'].Indicator == 'TransportFactorQty'].reset_index().Value[0]
    model.TransportFactorQty = pe.Param(initialize=parameter_dict, default=0)

    # Transport Load
    parameter_dict = dfs['Economics'][dfs['Economics'].Indicator == 'TransportLoadMin'].reset_index().Value[0]
    model.Q_TR_min = pe.Param(initialize=parameter_dict, default=0.4)
    parameter_dict = dfs['Economics'][dfs['Economics'].Indicator == 'TransportLoadMax'].reset_index().Value[0]
    model.Q_TR_max = pe.Param(initialize=parameter_dict, default=0)

    # Distances
    parameter_dict = dfs['Distances_SupplierPlant'].set_index(['S', 'X']).Value.dropna().to_dict()
    model.Dist_S = pe.Param(model.S, model.X, initialize=parameter_dict, default=0)
    parameter_dict = dfs['Distances_PlantMarket'].set_index(['X', 'M']).Value.dropna().to_dict()
    model.Dist_M = pe.Param(model.X, model.M, initialize=parameter_dict, default=0)

    # Storage_Cost
    parameter_dict = dfs['Economics'][dfs['Economics'].Indicator == 'HoldingCostFactor'].reset_index().Value[0]
    model.HoldingCostFactor = pe.Param(initialize=parameter_dict, default=0.4)
    parameter_dict = dfs['Economics'][dfs['Economics'].Indicator == 'FeedstockStorageCost'].reset_index().Value[0]
    model.FeedstockStorageCost = pe.Param(initialize=parameter_dict, default=0)

    # Compatibility_Purification
    parameter_dict = dfs['Compatibility_Purification'].set_index(['P', 'K_pur']).Value.dropna().to_dict()
    model.Compatibility_Purification = pe.Param(model.P, model.K_pur, initialize=parameter_dict, default=0)

    # Compatibility_Conversion
    parameter_dict = dfs['Compatibility_Conversion'].set_index(['K_pret', 'K_conv', 'J']).Value.dropna().to_dict()
    model.Compatibility_Conversion = pe.Param(model.K_pret, model.K_conv, model.J, initialize=parameter_dict, default=0)

    # Yield_Hydrolysis
    parameter_dict = dfs['Yield_Hydrolysis'].set_index(['K_conv', 'K_pret', 'C', 'P']).Value.dropna().to_dict()
    model.Yield_Hydrolysis = pe.Param(model.K_conv, model.K_pret, model.C, model.P, initialize=parameter_dict, default=1)

    # Solid_Loading
    parameter_dict = dfs['Solid_Loading'].set_index(['K_conv']).Value.dropna().to_dict()
    model.Solid_Loading = pe.Param(model.K_conv, initialize=parameter_dict, default=1)

    # Compatibility_Pretreatment
    parameter_dict = dfs['Compatibility_Pretreatment'].set_index(['K_prep', 'K_pret']).Value.dropna().to_dict()
    model.Compatibility_Pretreatment = pe.Param(model.K_prep, model.K_pret, initialize=parameter_dict, default=0)

    # Feedstock_Composition
    parameter_dict = dfs['Feedstock_Composition'].set_index(['C', 'F']).Value.dropna().to_dict()
    model.Feedstock_Composition = pe.Param(model.C, model.F, initialize=parameter_dict, default=0)
    # Yield_Pretreatment
    parameter_dict = dfs['Yield_Pretreatment'].set_index(['C', 'K_pret']).Value.dropna().to_dict()
    model.Yield_Pretreatment = pe.Param(model.C, model.K_pret, initialize=parameter_dict, default=0)
    # Yield_Conversion
    parameter_dict = dfs['Yield_Conversion'][['P', 'C', 'K_conv', 'J', "Value"]].set_index(['P', 'C', 'K_conv', 'J']).Value.dropna().to_dict()
    model.Yield_Conversion = pe.Param(model.P, model.C, model.K_conv, model.J, initialize=parameter_dict, default=0)
    # Yield_Purification 
    parameter_dict = dfs['Yield_Purification'].set_index(['P', 'K_pur']).Value.dropna().to_dict()
    model.Yield_Purification = pe.Param(model.P, model.K_pur, initialize=parameter_dict, default=0)
    # Demand
    parameter_dict = dfs['Demand'].set_index(['P', 'M', 'T']).Value.dropna().to_dict()
    model.Demand = pe.Param(model.P, model.M, model.T, initialize=parameter_dict, default=0)
    # Supply
    parameter_dict = dfs['Supply'].set_index(['F', 'S', 'T']).Value.dropna().to_dict()
    model.Supply = pe.Param(model.F, model.S, model.T, initialize=parameter_dict, default=0)
    # Price
    parameter_dict = dfs['Price_F'].set_index(['F', 'S', 'T']).Value.dropna().to_dict()
    model.Price_F = pe.Param(model.F, model.S, model.T, initialize=parameter_dict, default=0)
    parameter_dict = dfs['Price_P'].set_index(['P', 'M', 'T']).Value.dropna().to_dict()
    model.Price_P = pe.Param(model.P, model.M, model.T, initialize=parameter_dict, default=0)

    parameter_dict = dfs['Price_U'].dropna().set_index(['U']).Value.dropna().to_dict()
    model.Price_U = pe.Param(model.U, initialize=parameter_dict, default=0)

    # Technology_Cost
    parameter_dict = dfs['Technology_Cost'].set_index(['K']).Value.dropna().to_dict()
    model.Technology_Cost = pe.Param(model.K, initialize=parameter_dict, default=0)
    parameter_dict = dfs['Technology_Cost'].set_index(['K']).BaseScale.dropna().to_dict()
    model.Technology_BaseScale = pe.Param(model.K, initialize=parameter_dict, default=0, within = pe.NonNegativeReals)
    parameter_dict = dfs['Technology_Cost'].set_index(['K']).Flow.dropna().to_dict()
    model.Technology_Flow = pe.Param(model.K, initialize=parameter_dict, default="IN", within = model.Flow)
    parameter_dict = dfs['Technology_Cost'].set_index(['K']).ScaleFactor.dropna().to_dict()
    model.Technology_ScaleFactor = pe.Param(model.K, initialize=parameter_dict, default=0, within = pe.NonNegativeReals)
    parameter_dict = dfs['Technology_Cost'].set_index(['K']).InstalationFactor.dropna().to_dict()
    model.Technology_InstalationFactor = pe.Param(model.K, initialize=parameter_dict, default=0, within = pe.NonNegativeReals)
    parameter_dict = dfs['CEPCI'].set_index(['T']).Value.dropna().to_dict()
    model.CEPCI = pe.Param(model.T, initialize=parameter_dict, default=0)

    # Operating_Cost (same structure as Depreciation_Factor)
    parameter_dict = dfs['Economics'][dfs['Economics'].Indicator == 'OperatingCostFactor'].reset_index().Value[0]
    model.Operating_Cost = pe.Param(initialize=parameter_dict, default=0)
    # Maintenance_Cost (same structure as Depreciation_Factor)
    parameter_dict = dfs['Economics'][dfs['Economics'].Indicator == 'MaintenanceCostFactor'].reset_index().Value[0]
    model.Maintenance_Cost = pe.Param(initialize=parameter_dict, default=0)
    # Depreciation
    depreciationfactor = (1 - dfs['Economics'][dfs['Economics'].Indicator == 'SalvageValue'].reset_index().Value[0]) / dfs['Economics'][dfs['Economics'].Indicator == 'LifeCycle'].reset_index().Value[0]
    model.DepreciatonFactor = pe.Param(initialize=depreciationfactor, default=0)

    # Number of oeprators
    N_OP = dfs['Economics'][dfs['Economics'].Indicator == 'Operators'].reset_index().Value[0]
    model.N_OP = pe.Param(initialize=N_OP, default=0)
    # Salaries
    Salary = dfs['Labor_Salary'].set_index(['X']).Salary.dropna().to_dict()
    model.Salary = pe.Param(model.X, initialize=Salary, default=0)

    # Utility_Consumption
    parameter_dict = dfs['Utility_Consumption'].set_index(['U', 'K']).Value.dropna().to_dict()
    model.Utility_Consumption = pe.Param(model.U, model.K, initialize=parameter_dict, default=0)
    parameter_dict = dfs['Utility_Consumption'].set_index(['U', 'K']).Flow.dropna().to_dict()
    model.Utility_Consumption_Flow = pe.Param(model.U, model.K, initialize=parameter_dict,  default="IN", within = model.Flow)
    # Byproduct_Production
    parameter_dict = dfs['Byproduct_Production'].set_index(['BP', 'K']).Value.dropna().to_dict()
    model.Byproduct_Production = pe.Param(model.BP, model.K, initialize=parameter_dict, default=0)

    #LCI
    parameter_dict = dfs['LCI_Factors'].set_index(['U'])['CO2'].dropna().to_dict()
    model.LCI_Factor = pe.Param(model.U, initialize=parameter_dict, default=0)

    ####################
    # Binary Variables
    ####################

    # Indicates if technology k with capacity q is chosen to be installed in location x, at period t
    model.Y_tech = pe.Var(model.K, model.Q, model.X, model.T, domain=pe.Binary)

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

    # Transport Trips
    model.NT_S = pe.Var(model.S, model.X, model.T, domain = pe.NonNegativeIntegers)
    model.NT_M = pe.Var(model.X, model.M, model.T, domain = pe.NonNegativeIntegers)

    # Supplied feedstock f from supplier s in plant x
    model.F_in_feed = pe.Var(model.S, model.F, model.X, model.T, domain=pe.NonNegativeReals)

    # Quantity of feedstock f going into preprocessing k_prep at period t
    model.F_in_prep = pe.Var(model.F, model.K_prep, model.X, model.T, domain=pe.NonNegativeReals)# Quantity of feedstock f going into preprocessing k_prep at period t
    model.F_out_prep = pe.Var(model.F, model.K_prep, model.X, model.T, domain=pe.NonNegativeReals)

    # Quantity of feedstock f going into pretreatment k_pret at period t
    model.F_in_pret = pe.Var(model.F, model.K_pret, model.X, model.T, domain=pe.NonNegativeReals)
    # Quantity of component c going out of pretreatment k_pret at period t
    model.F_out_pret = pe.Var(model.C, model.F, model.K_pret, model.X, model.T, domain=pe.NonNegativeReals)

    # Quantity of component c going into conversion k_conv with fermentation strain j at period t to produce product p
    model.F_in_conv = pe.Var(model.C, model.F, model.P, model.K_conv, model.X, model.J, model.T, domain=pe.NonNegativeReals)
    # Quantity of component c going into conversion k_conv with fermentation strain j at period t to produce product p
    model.F_out_conv = pe.Var(model.P, model.F, model.K_conv, model.X, model.J, model.T, domain=pe.NonNegativeReals)

    # Bioreactor Liquid
    model.F_liquid_conv = pe.Var(model.K_conv, model.X, model.T, domain = pe.Reals)
    model.F_water_conv = pe.Var(model.K_conv, model.X, model.T, domain = pe.NonNegativeReals)

    # Quantity of product p going into purification k_conv at period t
    model.F_in_pur = pe.Var(model.P, model.K_pur, model.X, model.T, domain=pe.NonNegativeReals)
    # Quantity of product p going out of purification k_pur at period t
    model.F_out_pur = pe.Var(model.P, model.K_pur, model.X, model.T, domain=pe.NonNegativeReals)

    # Quanitity of product leaving plant to market
    model.F_out_prod = pe.Var(model.M, model.P, model.X,  model.T, domain=pe.NonNegativeReals)


    # Total flow of utility u consumed in process k at period t
    model.F_util = pe.Var(model.U, model.K, model.X, model.T, domain=pe.Reals)
    # Total flow of byproduct bp generated in process k at period t
    model.F_bypr = pe.Var(model.BP, model.K, model.X, model.T, domain=pe.NonNegativeReals)


    # Capital investment
    model.PC = pe.Var(model.X, model.T, domain=pe.NonNegativeReals)
    model.CAPEX = pe.Var(model.X, model.T, domain=pe.NonNegativeReals)
    model.TCI = pe.Var(model.X, model.T, domain=pe.NonNegativeReals)
    model.COM = pe.Var(model.X, model.T, domain=pe.NonNegativeReals)

    # LCA
    model.LCA = pe.Var(domain=pe.NonNegativeReals)


    ####################
    # Constraints
    ####################

    # Technology and Capacity Allocation
    # If a technology is selected at time t, it remains until the end of T. Only one cpacity per technology
    model.tech_continuity = pe.ConstraintList()
    model.tech_capacity = pe.ConstraintList()
    for t in model.T:
        for x in model.X:
            for k in model.K:
                model.tech_capacity.add(sum(model.Y_tech[k, q, x, t] for q in model.Q) <= 1)
                if t < max(model.T):
                    for q in model.Q:
                        model.tech_continuity.add(model.Y_tech[k, q, x, t+2] >= model.Y_tech[k, q, x, t])

    # At least one technology must exist at each time t per stage. Also, it can only exist if a capacity is installed
    model.tech_existance = pe.ConstraintList()
    for t in model.T:
        model.tech_existance.add(sum(model.Y_prep[f, k, x, t] for f in model.F for k in model.K_prep for x in model.X) >= 1)
        model.tech_existance.add(sum(model.Y_pret[f, k, x, t] for f in model.F for k in model.K_pret for x in model.X) >= 1)
        model.tech_existance.add(sum(model.Y_conv[f, p, k, x, j, t] for f in model.F for p in model.P for k in model.K_conv for x in model.X for j in model.J) >= 1)
        model.tech_existance.add(sum(model.Y_pur[p, k, x, t] for p in model.P for k in model.K_pur for x in model.X) >= 1)

        for x in model.X:
            for f in model.F:
                for k in model.K_prep:
                    model.tech_existance.add(model.Y_prep[f, k, x, t] <= sum(model.Y_tech[k, q, x, t] for q in model.Q))
                for k in model.K_pret:
                    model.tech_existance.add(model.Y_pret[f, k, x, t] <= sum(model.Y_tech[k, q, x, t] for q in model.Q))
            for p in model.P:
                for f in model.F:
                    for j in model.J:
                        for k in model.K_conv:
                            model.tech_existance.add(model.Y_conv[f, p, k, x, j, t] <= sum(model.Y_tech[k, q, x, t] for q in model.Q))
                for k in model.K_pur:
                    model.tech_existance.add(model.Y_pur[p, k, x, t] <= sum(model.Y_tech[k, q, x, t] for q in model.Q))


    # Select only one process per transformation
    model.single_technology = pe.ConstraintList()
    model.single_prep_constraint = pe.Constraint(model.F, model.X, model.T, rule = lambda m, f, x, t: sum(m.Y_prep[f, k_prep, x, t] for k_prep in model.K_prep) <= 1)
    model.single_pret_constraint = pe.Constraint(model.F, model.X, model.T,  rule = lambda m, f, x, t: sum(m.Y_pret[f, k_pret, x, t] for k_pret in model.K_pret) <= 1)
    model.single_pur_constraint = pe.Constraint(model.P, model.X, model.T, rule = lambda m, p, x, t: sum(m.Y_pur[p, k_pur, x, t] for k_pur in model.K_pur) <= 1)



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
    
    
    model.Test = pe.ConstraintList()

    model.bigM_pathways = pe.ConstraintList()
    for t in model.T:
        for x in model.X:

            # Setting Capacity Constraints (choosing capacity)
            for q in model.Q:
                for k in model.K_prep:
                    if model.Technology_Flow[k] == "IN":
                        model.bigM_pathways.add(sum(model.F_in_prep[f, k, x, t] for f in model.F) <= q * model.Technology_BaseScale[k] + M * (1 - model.Y_tech[k, q, x, t]))
                    elif model.Technology_Flow[k] == "OUT":
                        model.bigM_pathways.add(sum(model.F_out_prep[f, k, x, t] for f in model.F) <= q * model.Technology_BaseScale[k] + M * (1 - model.Y_tech[k, q, x, t]))
                for k in model.K_pret:
                    if model.Technology_Flow[k] == "IN":
                        model.bigM_pathways.add(sum(model.F_in_pret[f, k, x, t] for f in model.F) <= q * model.Technology_BaseScale[k] + M * (1 - model.Y_tech[k, q, x, t]))
                    elif model.Technology_Flow[k] == "OUT":
                        model.bigM_pathways.add(sum(model.F_out_pret[c, f, k, x, t] for c in model.C for f in model.F) <= q * model.Technology_BaseScale[k] + M * (1 - model.Y_tech[k, q, x, t]))
                for k in model.K_conv:
                    if model.Technology_Flow[k] == "IN":
                        model.bigM_pathways.add(sum(model.F_in_conv[c, f, p, k, x, j, t] for c in model.C for f in model.F for p in model.P for j in model.J) <= q * model.Technology_BaseScale[k] + M * (1 - model.Y_tech[k, q, x, t]))
                    elif model.Technology_Flow[k] == "OUT":
                        model.bigM_pathways.add(sum(model.F_out_conv[p, f, k, x, j, t] for f in model.F for p in model.P for j in model.J) <= q * model.Technology_BaseScale[k] + M * (1 - model.Y_tech[k, q, x, t]))
                for k in model.K_pur:
                    if model.Technology_Flow[k] == "IN":
                        model.bigM_pathways.add(sum(model.F_in_pur[p, k, x, t] for p in model.P)  <= q * model.Technology_BaseScale[k] + M * (1 - model.Y_tech[k, q, x, t]))
                    elif model.Technology_Flow[k] == "OUT":
                        model.bigM_pathways.add(sum(model.F_out_pur[p, k, x, t] for p in model.P) <= q * model.Technology_BaseScale[k] + M * (1 - model.Y_tech[k, q, x, t]))

            for f in model.F:

                # Choosing preprocessing
                for k_prep in model.K_prep:
                    # Corresponding to Y_prep[f, k_prep, t] = True in disjunctive formulation
                    if t > 1:
                        model.bigM_pathways.add(model.F_in_prep[f, k_prep, x, t]  <= 
                                                model.ST[f, x, t - 2] + sum(model.F_in_feed[s, f, x, t] for s in model.S)
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
                            sum(model.Compatibility_Pretreatment[k_prep, k_pret] * model.F_out_prep[f, k_prep, x, t] for k_prep in model.K_prep) 
                            + M * (1 - model.Y_pret[f, k_pret, x, t]))
                        # Corresponding to Y_pret[f, k_pret, t] = False in disjunctive formulation
                        model.bigM_pathways.add(
                            model.F_in_pret[f, k_pret, x, t] <= M * model.Y_pret[f, k_pret, x, t])
            

                # Choosing purification
                for k_pur in model.K_pur:
                    for p in model.P:

                        # Corresponding to Y_pur[p, k_pur, t] = True in disjunctive formulation
                        model.Test.add(
                            model.F_in_pur[p, k_pur, x, t] <= 
                            model.Compatibility_Purification[p, k_pur] * sum(model.F_out_conv[p, f, k, x, j, t] for k in model.K_conv for j in model.J for f in model.F)
                            + M * (1 - model.Y_pur[p, k_pur, x, t]))
                        
                        # Corresponding to Y_pur[p, k_pur, t] = False in disjunctive formulation
                        model.Test.add(
                            model.F_in_pur[p, k_pur, x, t] <= M * model.Y_pur[p, k_pur, x, t])
                        
            for m in model.M:
                for p in model.P:
                    # Corresponding to Y_supply[s, f, x, t] = True in disjunctive formulation
                    if t>1:
                        model.bigM_pathways.add(model.F_out_prod[m, p, x, t]  <= 
                                                sum(model.F_out_pur[p, k_pur, x, t] for k_pur in model.K_pur) + model.ST[p, x, t - 2]
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
    def rule_flow_disjunct(disjunct, f, p, k_pret, k_conv, x, t, flag):
        m = disjunct.model()
        
        if flag:
            disjunct.c1 = pe.ConstraintList() 
            if k_conv in model.K_conv_F:
                for j in model.J:
                    disjunct.c1.add(m.Y_conv[f, p, k_conv, x, j, t] + m.Y_pret[f, k_pret, x, t] == 2)
                    for c in model.C_Solid:
                        total = sum(m.Feedstock_Composition[cs, f] * m.Yield_Pretreatment[cs, k_pret] for cs in m.C_Solid)
                        mass_fraction = (m.Feedstock_Composition[c, f] * m.Yield_Pretreatment[c, k_pret]) / total if total > 0 else 0
                        disjunct.c1.add(expr = 
                                        m.F_in_conv[c, f, p, k_conv, x, j, t] 
                                        == mass_fraction * sum(m.F_in_conv[cs, f, p, k_conv, x, j, t] for cs in model.C_Solid))
                    for c in model.C_Liquid:
                        total = sum(m.Feedstock_Composition[cs, f] * m.Yield_Pretreatment[cs, k_pret] for cs in m.C_Liquid)
                        mass_fraction = (m.Feedstock_Composition[c, f] * m.Yield_Pretreatment[c, k_pret]) / total if total > 0 else 0
                        disjunct.c1.add(expr = 
                                        m.F_in_conv[c, f, p, k_conv, x, j, t] 
                                        == mass_fraction * sum(m.F_in_conv[cs, f, p, k_conv, x, j, t] for cs in model.C_Liquid))
                for c in model.C:
                    total = sum(m.Feedstock_Composition[ct, f] * m.Yield_Pretreatment[ct, k_pret] for ct in m.C)
                    mass_fraction = (m.Feedstock_Composition[c, f] * m.Yield_Pretreatment[c, k_pret]) / total if total > 0 else 0
                    disjunct.c1.add(expr = 
                                    sum(m.F_in_conv[c, f, p, k_conv, x, j, t] for j in model.J)
                                    == mass_fraction * sum(m.F_in_conv[ct, f, p, k_conv, x, j, t] for ct in model.C for j in model.J)) 
            else:
                for c in model.C:
                    for j in model.J:
                        total = sum(m.Feedstock_Composition[ct, f] * m.Yield_Pretreatment[ct, k_pret] for ct in m.C)
                        mass_fraction = (m.Feedstock_Composition[c, f] * m.Yield_Pretreatment[c, k_pret]) / total if total > 0 else 0
                        disjunct.c1.add(expr = 
                                        m.F_in_conv[c, f, p, k_conv, x, j, t] 
                                        == mass_fraction * sum(m.F_in_conv[ct, f, p, k_conv, x, j, t] for ct in model.C))                
                
        else:
            disjunct.c0 = pe.ConstraintList() 
            for j in model.J:
                disjunct.c0.add(m.Y_conv[f, p, k_conv, x, j, t] + m.Y_pret[f, k_pret, x, t] <= 1)

    model.flow_disjunct = gdp.Disjunct(model.F, model.P, model.K_pret, model.K_conv, model.X, model.T, [0, 1], rule = rule_flow_disjunct)

    model.flow_disjunction = gdp.Disjunction(model.F, model.P, model.K_pret, model.K_conv, model.X, model.T, 
                                                    rule = lambda m, f, p, k_pret, k_conv, x, t: [
                                                        m.flow_disjunct[f, p, k_pret, k_conv, x, t, i] for i in [0, 1]
                                                    ])

    # Bioreactor Liquid
    def constraint_sl_liquid(model):
        model.sl_liquid = pe.ConstraintList()
        for t in model.T:
            for x in model.X:
                    for k_conv in model.K_conv:
                        if k_conv in model.K_conv_F:
                            model.sl_liquid.add(model.F_liquid_conv[k_conv, x, t] == 
                                                ((1 - model.Solid_Loading[k_conv]) / model.Solid_Loading[k_conv]) * sum(model.F_in_conv[c, f, p, k_conv, x, j, t] for c in model.C_Solid for j in model.J for p in model.P for f in model.F)
                                                )
                        elif k_conv in model.K_conv_CF:
                            model.sl_liquid.add(model.F_liquid_conv[k_conv, x, t] == 
                                                ((1 - model.Solid_Loading[k_conv]) / model.Solid_Loading[k_conv]) * sum(model.F_in_conv[c, f, p, k_conv, x, j, t] for c in model.C_Solid for j in model.J for f in model.F for p in model.P) 
                                                - sum(model.F_in_conv[c, f, p, k_conv, x, j, t]  for c in model.C_Liquid for j in model.J for p in model.P for f in model.F)
                                                # byrpoducts generated in pretreatment
                                                - sum(sum((model.F_in_conv[c, f, p, k_conv, x, j, t] / (model.Feedstock_Composition[c, f] * model.Yield_Pretreatment[c, k_pret]) if (model.Feedstock_Composition[c, f] * model.Yield_Pretreatment[c, k_pret]) > 0 else 0 ) for c in model.C for p in model.P for j in model.J) * (model.Byproduct_Production[bp, k_pret]) for k_pret in model.K_pret for bp in model.BP for f in model.F) 
                                                # utilities affecting mass balnace in pretreatment
                                                - sum(sum((model.F_in_conv[c, f, p, k_conv, x, j, t] / (model.Feedstock_Composition[c, f] * model.Yield_Pretreatment[c, k_pret]) if (model.Feedstock_Composition[c, f] * model.Yield_Pretreatment[c, k_pret]) > 0 else 0 ) for c in model.C for p in model.P for j in model.J) * (model.Utility_Consumption[rm, k_pret]) for k_pret in model.K_pret for rm in model.RM for f in model.F) 
                                                )
                        else:
                            model.sl_liquid.add(model.F_liquid_conv[k_conv, x, t] == 0)
                        model.sl_liquid.add(model.F_water_conv[k_conv, x, t] >= model.F_liquid_conv[k_conv, x, t])

    constraint_sl_liquid(model)

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
                                            model.ST[f, x, t - 2] + sum(model.F_in_feed[s, f, x, t] for s in model.S))
                else:
                    model.Mass_Balances.add(sum(model.F_in_prep[f, k_prep, x, t] for k_prep in model.K_prep) <= 
                                            sum(model.F_in_feed[s, f, x, t] for s in model.S))

        for p in model.P:
            model.Mass_Balances.add(sum(model.F_out_prod[m, p, x, t] for x in model.X for m in model.M) >= 0.05 * sum(model.Demand[p, m, t] for m in model.M))
            for m in model.M:
                model.Mass_Balances.add(sum(model.F_out_prod[m, p, x, t] for x in model.X) <= model.Demand[p, m, t])
                if t>1:
                    model.Mass_Balances.add(sum(model.F_out_pur[p, k_pur, x, t] for k_pur in model.K_pur) >= 
                                            (model.ST[p, x, t - 2] + sum(model.F_out_prod[m, p, x, t] for m in model.M)))
                else:
                    model.Mass_Balances.add(sum(model.F_out_pur[p, k_pur, x, t] for k_pur in model.K_pur) >= 
                                            sum(model.F_out_prod[m, p, x, t] for m in model.M))
        
        for x in model.X:

            # Storage
            if t > 1:
                for f in model.F:
                    model.Mass_Balances.add(model.ST[f, x, t] == model.ST[f, x, t - 2] + sum(model.F_in_feed[s, f, x, t] for s in model.S) - sum(model.F_in_prep[f, k_prep, x, t] for k_prep in model.K_prep))
                for p in model.P:
                    model.Mass_Balances.add(model.ST[p, x, t] == model.ST[p, x, t - 2] - sum(model.F_out_prod[m, p, x, t] for m in model.M) + sum(model.F_out_pur[p, k_pur, x, t] for k_pur in model.K_pur))
            else:
                for fp in model.FP:
                    model.Mass_Balances.add(model.ST[fp, x, t] == 0)

            for c in model.C:
                for f in model.F:
                    model.Mass_Balances.add(sum(model.F_in_conv[c, f, p, k, x, j, t] for p in model.P for k in model.K_conv for j in model.J) <= sum(model.F_out_pret[c, f, k, x, t] for k in model.K_pret))     
            

            # In Processes (Transformations)
            for f in model.F:
                for k_prep in model.K_prep:
                    model.Mass_Balances.add(model.F_out_prep[f, k_prep, x, t] == model.F_in_prep[f, k_prep, x, t])
                for k_pret in model.K_pret:
                    for c in model.C:
                        model.Mass_Balances.add(model.F_out_pret[c, f, k_pret, x, t] <= model.F_in_pret[f, k_pret, x, t]  * model.Feedstock_Composition[c, f] * model.Yield_Pretreatment[c, k_pret])
                
                for p in model.P:
                    for k_conv in model.K_conv:
                        for j in model.J:
                            if p not in model.C:
                                for k_pret in model.K_pret:
                                    model.Mass_Balances_test.add(model.F_out_conv[p, f, k_conv, x, j, t]  <= 
                                                                sum(model.F_in_conv[c, f, p, k_conv, x, j, t] *  model.Yield_Conversion[p, c, k_conv, j] * model.Yield_Hydrolysis[k_conv, k_pret, c, p] for c in model.C)
                                                                + M * (1 - model.Y_pret[f, k_pret, x, t])
                                                                )
                                model.Mass_Balances_test.add(model.F_out_conv[p, f, k_conv, x, j, t]  <= 
                                                                M * sum(model.Y_pret[f, k_pret, x, t] for k_pret in model.K_pret)
                                                            )
                            else:
                                model.Mass_Balances_test.add(model.F_out_conv[p, f, k_conv, x, j, t] == sum(model.F_in_conv[c, f, pt, k_conv, x, j, t] * model.Yield_Hydrolysis[k_conv, k_pret, c, p] for c in model.C if c in model.P for pt in model.P))
            for p in model.P:
                for k_pur in model.K_pur:
                    model.Mass_Balances.add(model.F_out_pur[p, k_pur, x, t] <= model.F_in_pur[p, k_pur, x, t] * model.Yield_Purification[p, k_pur])
    


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
                        if u in ["Water"]:
                            model.utility_consumption.add(model.F_util[u, k, x, t] == 
                                                        sum(model.F_in_conv[c, f, p, k, x, j, t] for f in model.F for c in model.C for p in model.P for j in model.J)  * model.Utility_Consumption[u, k] 
                                                        + model.F_water_conv[k, x, t])
                        else:
                            model.utility_consumption.add(model.F_util[u, k, x, t] == sum(model.F_in_conv[c, f, p, k, x, j, t] for f in model.F for c in model.C for p in model.P for j in model.J)  * model.Utility_Consumption[u, k])
                for k in model.K_pur:
                    if model.Utility_Consumption_Flow[u, k] == "IN":
                        model.utility_consumption.add(model.F_util[u, k, x, t] == sum(model.F_in_pur[p, k, x, t] for p in model.P) * model.Utility_Consumption[u, k])
                    else:
                        model.utility_consumption.add(model.F_util[u, k, x, t] == sum(model.F_out_pur[p, k, x, t] for p in model.P) * model.Utility_Consumption[u, k])
                # I cna only produce RNG or Eletricity that is consumed
                model.utility_consumption.add(sum(model.F_util[u, k, x, t] for k in model.K) >= 0)
                
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

    # Tranportation
    def rule_transportS_disjunct(disjunct, s, x, t, flag):
        m = disjunct.model()
        
        if flag:
            disjunct.c1 = pe.ConstraintList() 
            disjunct.c1.add(sum(m.Y_supply[s, f, x, t] for f in m.F) >= 1)
            disjunct.c1.add(
                m.NT_S[s, x, t] >=
                (sum(m.F_in_feed[s, f, x, t] for f in m.F) + m.Q_TR_min) / m.Q_TR_max
            )
        else:
            disjunct.c0 = pe.ConstraintList() 
            disjunct.c0.add(sum(m.Y_supply[s, f, x, t] for f in m.F) == 0)
            disjunct.c0.add(m.NT_S[s, x, t] == 0 )
            #disjunct.c0.add(sum(m.F_in_feed[s, f, x, t] for f in m.F) == 0 )

    model.transportS_disjunct = gdp.Disjunct(model.S, model.X, model.T, [0, 1], rule = rule_transportS_disjunct)
    model.transportS_disjunction = gdp.Disjunction(model.S, model.X, model.T, 
                                                    rule = lambda m, s, x, t: [
                                                        m.transportS_disjunct[s, x, t, i] for i in [0, 1]
                                                    ])
    """ 
    def rule_transportM_disjunct(disjunct, x, mk, t, flag):
        m = disjunct.model()
        
        if flag:
            disjunct.c1 = pe.ConstraintList() 
            disjunct.c1.add(sum(m.Y_market[mk, p, x, t] for p in m.P) >= 1)
            disjunct.c1.add(
                m.NT_M[x, mk, t] >=
                (sum(m.F_out_prod[mk, p, x, t] for p in m.P) + m.Q_TR_min) / m.Q_TR_max
            )
        else:
            disjunct.c0 = pe.ConstraintList() 
            disjunct.c0.add(sum(m.Y_market[mk, p, x, t] for p in m.P) == 0)
            disjunct.c0.add(m.NT_M[x, mk, t] == 0 )
            #disjunct.c0.add(sum(m.F_out_prod[mk, p, x, t] for p in m.P) == 0 )

    model.transportM_disjunct = gdp.Disjunct(model.X, model.M, model.T, [0, 1], rule = rule_transportM_disjunct)
    model.transportM_disjunction = gdp.Disjunction(model.X, model.M, model.T, 
                                                    rule = lambda m, x, mk, t: [
                                                        m.transportM_disjunct[x, mk, t, i] for i in [0, 1]
                                                    ])
    """

    model.prod_const = pe.ConstraintList()
    for x in model.X:
        for t in model.T:
            for p in model.P:
                #model.prod_const.add(sum(model.F_out_prod[m, p, x, t] for m in model.M) <= sum(model.F_out_pur[p, k_pur, x, t] for k_pur in model.K_pur))
                pass
    ####################
    # Environmental Constraints
    ####################
    

    ####################
    # Economic Constraints
    ####################
    model.LCA_calc = pe.ConstraintList()
    model.LCA_calc.add(
        model.LCA == 
        sum( 
            (
            sum(model.F_util[u, k, x, t] * model.LCI_Factor[u] for k in model.K for u in model.U) 
            + model.LCI_Factor["Transport"] * (sum(model.F_in_feed[s, f, x, t] * model.Dist_S[s, x] for f in model.F for s in model.S) + sum(model.F_out_prod[m, p, x, t] * model.Dist_M[x, m] for p in model.P for m in model.M))
            + model.LCI_Factor["Waste to landfill"] * (sum(model.F_in_prep[f, k, x, t] for f in model.F for k in model.K_prep) - sum(model.F_out_pur[p, k, x, t] for p in model.P for k in model.K_pur))
            )
            for x in model.X for t in model.T)
    )

    model.LCA_calc.add(model.LCA <= max_LCA)

    model.C_TR = pe.Var(model.X, model.T, domain = pe.NonNegativeReals)
    model.C_TR_calc = pe.ConstraintList()
    for x in model.X:
        for t in model.T:
            model.C_TR_calc.add(
                model.C_TR[x, t] == 0
                + model.TransportFactorKm * (
                    sum(model.NT_S[s, x, t] * model.Dist_S[s, x] for s in model.S)
                    + sum(model.NT_M[x, m, t] * model.Dist_M[x, m] for m in model.M)
                )
                + model.TransportFactorQty * (
                    sum(model.F_in_feed[s, f, x, t] for f in model.F for s in model.S)
                    + sum(model.F_out_prod[m, p, x, t] for p in model.P for m in model.M)
                )
            )

    model.Util_Cost = pe.Var(domain = pe.Reals)
    model.Util_Cost_calc = pe.Constraint(rule = lambda m: m.Util_Cost == sum(sum(m.F_util[u, k, x, t] for k in model.K) * m.Price_U[u] for u in model.U for x in m.X for t in m.T) )

    model.TCI_calc = pe.ConstraintList()
    for t in model.T:
        for x in model.X:
            if t > 1:
                model.TCI_calc.add(model.PC[x, t] == sum((model.Technology_Cost[k] * q**model.Technology_ScaleFactor[k] * model.CEPCI[t]/model.CEPCI[1] * (model.Y_tech[k, q, x, t] - model.Y_tech[k, q, x, t - 2])) for k in model.K for q in model.Q))
                model.TCI_calc.add(model.CAPEX[x, t] == sum((model.Technology_Cost[k] * q**model.Technology_ScaleFactor[k] * model.Technology_InstalationFactor[k] * model.CEPCI[t]/model.CEPCI[1] * (model.Y_tech[k, q, x, t] - model.Y_tech[k, q, x, t - 2])) for k in model.K for q in model.Q))
            else:
                model.TCI_calc.add(model.PC[x, 1] == sum((model.Technology_Cost[k] * q**model.Technology_ScaleFactor[k] * model.Y_tech[k, q, x, 1]) for k in model.K for q in model.Q))
                model.TCI_calc.add(model.CAPEX[x, 1] == sum((model.Technology_Cost[k] * q**model.Technology_ScaleFactor[k] * model.Technology_InstalationFactor[k] * model.Y_tech[k, q, x, 1]) for k in model.K for q in model.Q))
            
            model.TCI_calc.add(model.TCI[x, t] == 1.067556 * model.PC[x, t] +  1.5036 * model.CAPEX[x, t])

    model.Y_x = pe.Var(model.X, model.T, domain = pe.Binary)
    model.COM_calc = pe.ConstraintList()
    for t in model.T:
        for x in model.X:
            for k in model.K:
                for q in model.Q:
                    model.COM_calc.add(model.Y_x[x, t] >= model.Y_tech[k, q, x, t])
            model.COM_calc.add(
                model.COM[x,t] == 0
                + 0.1830096 * model.PC[x, t] + 0.2577600 * model.CAPEX[x, t]
                + 1.23 * sum(sum(model.F_util[u, k, x, t] for k in model.K) * model.Price_U[u] for u in model.U)
                + model.N_OP * model.Salary[x] * model.Y_x[x, t] 
                + model.DepreciatonFactor * model.PC[x, t]
            )
        
    model.EBIT_Plant = pe.Var(model.X, model.T, domain = pe.Reals)
    model.EBIT_Plant_calc = pe.Constraint(model.X, model.T,
                                    rule=lambda m, x, t: 
                                    m.EBIT_Plant[x, t] == 0
                                    # Earnings
                                    + sum(sum(m.F_out_prod[mk, p, x, t] * m.Price_P[p, mk, t] for p in model.P) for mk in model.M) 
                                    # Feedstock Cost
                                    - sum(sum(m.F_in_feed[s, f, x, t] * m.Price_F[f, s, t] for f in model.F) for s in model.S) 
                                    # Util cost
                                    # - sum(sum(m.F_util[u, k, x, t] for k in model.K) * m.Price_U[u] for u in model.U)
                                    # Cost of manufactoring
                                    - model.COM[x, t]
                                    # Transportation
                                    - model.C_TR[x, t]
                                    # Storage
                                    - m.HoldingCostFactor * sum(m.ST[p, x, t] * min(m.Price_P[p, mk, t] for mk in model.M) for p in model.P)  + sum(m.ST[f, x, t] * min(m.Price_F[f, s, t] for s in model.S) for f in model.F)
                                    )

    model.EBIT = pe.Var(model.T, domain = pe.Reals)
    model.EBIT_calc = pe.Constraint(model.T,
                                    rule=lambda m, t: 
                                    m.EBIT[t] == sum(m.EBIT_Plant[x, t] 
                                                     #
                                                     for x in m.X))

    ####################
    # Objectives
    ####################


    def OP_rule(m):      
        return sum(m.EBIT[t] 
                   #+ sum(m.DepreciatonFactor * m.PC[x, t] for x in m.X) 
                   for t in m.T)

    def NPV_rule(m):
        return sum(((1 - 0.2) * m.EBIT[t] + sum(m.DepreciatonFactor * m.PC[x, t] for x in m.X) - sum(m.TCI[x, t] for x in m.X for t in m.T)) / (1+0.03)**t for t in m.T)

    model.Objective_function = pe.Objective(rule = OP_rule, sense = pe.maximize) 

    ####################
    # Solvemfunc.print_vars(model)big
    ####################
    pe.TransformationFactory('gdp.bigm').apply_to(model, bigM=M)
    solver = po.SolverFactory("gurobi", solver_io="python")
    solver_parameters = "ResultFile=iismodel.mps" # write an ILP file to print the IIS
    solver.options['DualReductions'] = 0
    solver.options['NonConvex'] = 2
    #solver.options['Threads'] = 6
    #solver.options['NodefileStart'] = 0.5
    #solver.options['NodefileDir'] = "../Optimization/Nodes/"
    
    results = solver.solve(model, 
                        options_string=solver_parameters
                        )
    return results, model

# %%
def get_optimal(model, results):
    optimal = {}
    optimal["NPV"] = pe.value(model.Objective_function)
    optimal["Results"] = results
    for v in model.component_objects(pe.Var, active=True):
        try:
            varobject = getattr(model, str(v))
            optimal[str(v)] = {}
            for index in varobject:
                if varobject[index].value is not None:
                    if round(varobject[index].value, 4) != 0:
                        optimal[str(v)][index] = varobject[index].value
        except:
            pass
    return optimal

# %%
def solve_model():

    save = "SC"
    print("Starting!")
    results, model = create_n_solve(10**20)
    optimal = get_optimal(model, results)

    with open(f'FSC_maxNPV.pkl', mode='wb') as file:
        cloudpickle.dump(model, file)

    aux = model.LCA.value
    LCAs = [
            0.25*aux, 0.5*aux, 0.75*aux,
            0.05*aux, 0.4*aux, 0.6*aux, 0.9*aux
            ]

    all_optimal = {}
    all_optimal["maxNPV"] = optimal
    counter = 1
    print("E-Factor...")
    for max_LCA in LCAs:
        try:
            results, model = create_n_solve(max_LCA)
            optimal = get_optimal(model, results)
            with open(f'FSC_LCA{counter}.pkl', mode='wb') as file:
                cloudpickle.dump(model, file)
        except:
             pass

        
        counter +=1

        all_optimal[f"{save}/LCA{counter}"] = optimal

    with open(f'FSC_optimal.pkl', mode='wb') as file:
            pickle.dump(all_optimal, file)

    print("Done")
    return all_optimal



# %%
all_optimal = solve_model()

# %%
all_optimal["maxNPV"]

# %%
all_optimal["maxNPV"]

# %%



