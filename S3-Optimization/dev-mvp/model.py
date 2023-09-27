import pandas as pd
import numpy as np
import itertools
import pyomo.environ as pe
import pyomo.gdp as gdp
import gurobipy
import pyomo.opt as po


sets = {
    # "S"     : ["x1", "x2", ...],  # Indices for supply geogpraphies
    # "X"     : ["x1", "x2", ...],  # Indices for plant site locations
    # "M"     : ["m1", "m2", ...],  # Indices for markets
    "T"       : [1],  # Indices for periods (10 years)
    "F"       : ["Bstraw", "Wstraw"],  # Indices for feedstocks
    "C"       : [],  # Indices for components
    "C_Solid" : ["Cellulose", "Hemicellulose", "Lignin"],  # Indices for solid components
    "C_Liquid": ["Xylose", "Glucose"],  # Indices for liquid components
    
    "BP"      : ['Succinic Acid','Glycolic Acid','Formic Acid','Acetic Acid','Phenolic','Furfural','HMF'
                ],  # Indices for by-products
    
    "U"       : [
                    'Ammonia',
                    'Calcium hydroxide',
                    'Cooling water',
                    'Cyclohexane',
                    'Electricity',
                    'Ethylene Glycol',
                    'Gypsum disposal',
                    'H2SO4',
                    'HCL',
                    'Heat',
                    'Ionic Liquid mat.',
                    'Iron (II) chloride',
                    'Lime mat.',
                    'LP-Steam',
                    'MeOH',
                    'Monoethanolamine',
                    'Refrigeration -15,5°C',
                    'RNG',
                    'Steam 150 psig',
                    'Tri-n-octylamine (0,25 mol/kg) in 1-octanol',
                    'Waste to landfill',
                    'Water',
                    'Olivine',
                    'Magnesium oxide',
                    'Cooling',
                    'Enzymes'
                ],  # Indices for utilities
    
    "K"       : [],  # Indices for process options
    "K_prep"  : ['MILL3mm','MILL1.4mm','MILL10mm','MILL0.853mm','MILL0.15mm'],  # Indices for preprocessing options
    
    "K_pret"  : ['Dilute Acid','LHW','AFEX','Steam Explosion','Lime','Ionic Liquid'],  # Indices for pretreatment options
    
    "K_conv"  : ['SHF','SSF','SHCF','SSCF'],  # Indices for fermentation options
    
    "K_pur"   : ["Reactive Extraction", 'Precipitation', 'Electrodialysis', 'Direct Crystalization', 'Azeotropic Dist', 'Extraction', 'Absorption'],  # Indices for purification options
    
    "J"       : ["Strain1","Strain2","Strain3","Strain4","Strain5","Strain6 ","Strain7","Strain8","Strain9" 
                 #"8G-3", "SPT2OE", "Y2034", "Angle", "CRD51", "MD4", "TP50", "SR8/TFP3-XI76", "SR8/TFP3-XI76"
                ],  # Indices for fermentation strains
    
    "P"       : ['Ethanol', 'Lactic Acid', 'Succinic Acid'],  # Indices for final products
    
    "Q"       : [1]  # Indices for technology capacity levels
}

sets["C"] = sets["C_Solid"]  + sets["C_Liquid"] 
sets["K"] = sets["K_prep"]  + sets["K_pret"] + sets["K_conv"] + sets["K_pur"]

# Funciton for create the dataframe from the cartesian product of specific sets
def create_parametric_df(sets):
    df = pd.DataFrame(itertools.product(*sets.values()), columns=sets.keys())
    # drop duplicated rows
    df = df.drop_duplicates().reset_index(drop=True)
    # Create an empty value columns
    df.loc[:, 'Value'] = np.nan
    # Input 0 in impossible combinations (relevant for yields as you can't convert cellulose into cellulose)
    # Identify rows where a value in a column is equal to a value in another column for all columns
    equal_values_rows = df[df.apply(lambda row: len(set(row)) < len(row), axis = 1)]
    # Cahnge value for identified rows
    df.loc[equal_values_rows.index, "Value"] = 0
    return df

def create_compatability_dic():
    
    # Dictionary to store the DataFrames
    dfs = {}

    # Compatability Preprocessing with Pretreatment
    dfs['Compatibility_Pretreatment'] = create_parametric_df({'K_prep' : sets['K_prep'], 'K_pret' : sets['K_pret']})
    
    # Compatability Pretreatment, Conversion, Fermentation Strain, Component and Product   
    dfs['Conversion_Matrix_Pret_Conv'] = create_parametric_df({'K_pret': sets['K_pret'], 'K_conv': sets['K_conv']})
    dfs['Conversion_Matrix_C_Conv'] = create_parametric_df({'C': sets['C'], 'K_conv': sets['K_conv']})
    dfs['Conversion_Matrix_Conv_J'] = create_parametric_df({'K_conv': sets['K_conv'], 'J': sets['J']})
    dfs['Conversion_Matrix_J_P'] = create_parametric_df({'J': sets['J'], 'P': sets['P']})
    dfs['Conversion_Matrix_C_P'] = create_parametric_df({'C': sets['C'], 'P': sets['P']})
    dfs['Conversion_Matrix_Conv_P'] = create_parametric_df({'K_conv': sets['K_conv'], 'P': sets['P']})
    
    dfs['Compatibility_Purification'] = create_parametric_df({'P': sets['P'], 'K_pur': sets['K_pur']})

    return dfs

def create_parametric_dic(compat_dfs):
    # Dictionary to store the DataFrames
    dfs = {}

    # Generate DataFrames and store them in the dictionary
    dfs['Technology_Cost'] = create_parametric_df({'K': sets['K']})
    dfs['MinMax_Flows'] = create_parametric_df({'Flow': ['IN', 'OUT'], 'MinMax': ['Min', 'Max'], 'K': sets['K'], 'Q': sets['Q']})
    dfs['Utility_Consumption'] = create_parametric_df({'U': sets['U'], 'K': sets['K']})
    dfs['Byproduct_Production'] = create_parametric_df({'BP': sets['BP'], 'K': sets['K']})
    dfs['Feedstock_Composition'] = create_parametric_df({'C': sets['C'], 'F': sets['F']})
    
    dfs['Compatibility_Pretreatment'] = compat_dfs['Compatibility_Pretreatment']
    dfs['Yield_Pretreatment'] = create_parametric_df({'C': sets['C'], 'K_pret': sets['K_pret']})
    
    # Compatability Conversion
    df = pd.merge(compat_dfs['Conversion_Matrix_Pret_Conv'], compat_dfs['Conversion_Matrix_C_Conv'], on = "K_conv", how = "outer")
    df["Value"] = df.Value_x * df.Value_y
    df = df[df["Value"] != 0].drop(['Value_x', 'Value_y'], axis = 1)
    df = df.merge(compat_dfs['Conversion_Matrix_Conv_J'], on = "K_conv", how = "left")
    df["Value"] = df.Value_x * df.Value_y
    df = df[df["Value"] != 0].drop(['Value_x', 'Value_y'], axis = 1)
    df = df.merge(compat_dfs['Conversion_Matrix_J_P'], on = "J", how = "left")
    df["Value"] = df.Value_x * df.Value_y
    df = df[df["Value"] != 0].drop(['Value_x', 'Value_y'], axis = 1)
    df = df.merge(compat_dfs['Conversion_Matrix_C_P'], on = ["C", "P"], how = "left")
    df["Value"] = df.Value_x * df.Value_y
    df = df[df["Value"] != 0].drop(['Value_x', 'Value_y'], axis = 1)
    df = df.merge(compat_dfs['Conversion_Matrix_Conv_P'], on = ["K_conv", "P"], how = "left")
    df["Value"] = df.Value_x * df.Value_y
    df = df[df["Value"] != 0].drop(['Value_x', 'Value_y'], axis = 1)

    df = df[['C', 'K_pret', 'K_conv', 'J', 'P', 'Value']]

    dfs['Compatibility_Conversion'] = df
    
    dfs['Yield_Hydrolysis'] = df[['K_conv', 'K_pret', 'C', 'P', 'Value']].drop_duplicates().reset_index(drop=True).sort_values(['K_conv', 'K_pret', 'C', 'P'])
    dfs['Yield_Conversion'] = df[[ 'P', 'C', 'K_conv', 'J', 'Value']].drop_duplicates().reset_index(drop=True).sort_values(['J', 'P', 'C', 'K_conv'])
    
    dfs['Solid_Loading'] = create_parametric_df({'K_conv' : sets['K_conv']})

    dfs['Compatibility_Purification'] = compat_dfs['Compatibility_Purification']
    
    dfs['Yield_Purification'] = dfs['Compatibility_Purification'][dfs['Compatibility_Purification'].Value != 0].copy()
    dfs['Depreciation_Factor'] = create_parametric_df({'T' : sets['T']})
    dfs['Operating_Cost'] = dfs['Depreciation_Factor'].copy()
    dfs['Maintenance_Cost'] = dfs['Depreciation_Factor'].copy()
    dfs['Demand'] = create_parametric_df({'P': sets['P'], 'T': sets['T']})
    dfs['Supply'] = create_parametric_df({'F': sets['F'], 'T': sets['T']})
    dfs['Price_F'] = create_parametric_df({'F': sets['F'], 'T': sets['T']})
    dfs['Price_U'] = create_parametric_df({'U': sets['U'], 'T': sets['T']})
    dfs['Price_P'] = create_parametric_df({'P': sets['P'], 'T': sets['T']})

    return dfs

def create_compatability_excel(xls_path = "CompatabilityMatrices.xlsx"):
    dfs = create_compatability_dic()
    writer = pd.ExcelWriter(xls_path, engine='xlsxwriter')
    for name, df in dfs.items():
        df.to_excel(writer, sheet_name=name, index = False)
    writer.close()
    writer.handles = None
    return dfs

def read_compatability_excel(xls_path = "CompatabilityMatricesInputs.xlsx"):

    dfs = create_compatability_dic()
    excel_dfs = {}
    
    # Read excel
    xls = pd.ExcelFile(xls_path)
    for name, df in dfs.items():
        excel_df = pd.read_excel(xls, sheet_name=name)
        try:
            excel_df = excel_df.drop("Unnamed: 0", axis = 1)
        except:
            pass
        
        # concat excel df with full df
        df = pd.concat([excel_df, dfs[name]]).reset_index(drop=True)
        # keep first appearence (coming from excel df) of duplicated values in all columns except Value
        df = df.drop_duplicates(subset = df.columns.difference(['Value']), keep='first').reset_index(drop=True)
        # change all NA to 0 excpet for solid loading
        df.Value = 1 if name in ["Solid_Loading"] else df.Value.fillna(0)
        
        excel_dfs[name] = df
    
    xls.close()
    xls.handles = None
    
    return excel_dfs

def create_parametric_excel(xls_path = "ParametricTables.xlsx", compat = None):
    
    if compat is None:
        # Try to get input tables
        print("No compatabiity tables given. Attempting to retrieve from CompatabilityMatricesInputs.xlsx")
        compat_dfs = read_compatability_excel()
    elif type(compat) is str:
        compat_dfs = read_compatability_excel(compat)
    else:
        compat_dfs = compat
    
    dfs = create_parametric_dic(compat_dfs)
    writer = pd.ExcelWriter(xls_path, engine='xlsxwriter')
    for name, df in dfs.items():
        df.to_excel(writer, sheet_name=name, index = False)
    writer.close()
    writer.handles = None
    return dfs

def read_parametric_excel(xls_path, compat = None):
    
    if compat is None:
        # Try to get input tables
        print("No compatabiity tables given. Attempting to retrieve from CompatabilityMatricesInputs.xlsx")
        compat_dfs = read_compatability_excel()
    elif type(compat) is str:
        compat_dfs = read_compatability_excel(compat)
    else:
        compat_dfs = compat
        
    # Create empty dic to force parameters to be 0 when user accidently forgets certain inputs
    dfs = create_parametric_dic(compat_dfs)
    
    # Read excel
    excel_dfs = {}
    xls = pd.ExcelFile(xls_path)
    for name, df in dfs.items():
        excel_df = pd.read_excel(xls, sheet_name=name)
        try:
            excel_df = excel_df.drop("Unnamed: 0", axis = 1)
        except:
            pass
        # concat excel df with full df
        df = pd.concat([excel_df, dfs[name]]).reset_index(drop=True)
        # keep first appearence (coming from excel df) of duplicated values in all columns except Value
        df = df.drop_duplicates(subset = df.columns.difference(['Value']), keep='first').reset_index(drop=True)
        # change all NA to 0 excpet for solid loading
        df.Value = 1 if name in ["Solid_Loading"] else df.Value.fillna(0)
        
        excel_dfs[name] = df
    xls.close()
    xls.handles = None
    return excel_dfs

compat_dfs = create_compatability_excel("Inputs\CompatabilityMatrices-V4.xlsx")
dfs = create_parametric_excel("Inputs\ParametricTables-V4.xlsx", compat="Inputs\CompatabilityMatricesInputs-V4.xlsx")
dfs = read_parametric_excel("Inputs\ParametricTablesInputs-V4.xlsx", compat="Inputs\CompatabilityMatricesInputs-V4.xlsx")




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
        model.capex_calculation.add(model.CAPEX[t] == sum((model.Technology_Cost[k] * (model.Y_tech[k, q, t] - model.Y_tech[k, q, t - 1])) for k in model.K for q in model.Q))
    else:
        model.capex_calculation.add(model.CAPEX[1] == sum((model.Technology_Cost[k] * model.Y_tech[k, q, 1]) for k in model.K for q in model.Q))

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
        model.Mass_Balances.add(sum(model.F_out_pur[p, k, t] for k in model.K_pur) <= model.Demand[p, t])

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

####################
# Objecitve
####################
model.EBIT = pe.Var(model.T, domain = pe.Reals)
model.EBIT_calc = pe.Constraint(model.T, 
                                rule=lambda m, t: 
                                m.EBIT[t] == sum(sum(m.F_out_pur[p, k_pur, t] for k_pur in m.K_pur) * m.Price_P[p, t] for p in model.P) \
                                - sum(sum(m.F_in_prep[f, k, t] for k in m.K_prep) * m.Price_F[f, t] for f in model.F) 
                                - sum(sum(m.F_util[u, k, t] for k in model.K) * m.Price_U[u, t] for u in model.U) 
                                - (m.Operating_Cost[t] + m.Maintenance_Cost[t])* m.CAPEX[t]
                                )

def NPV_rule(m):      
    return sum( (m.EBIT[t] - m.CAPEX[t]) for t in m.T)
model.Objective_function = pe.Objective(rule = NPV_rule, sense = pe.maximize) 

####################
# Solve
####################

solver = po.SolverFactory("gurobi", olver_io="python")
solver_parameters = "ResultFile=iismodel.ilp" # write an ILP file to print the IIS
results = solver.solve(model, options_string=solver_parameters)
results.write()



def print_vars():
    for v in model.component_objects(pe.Var, active=True):
        try:
            varobject = getattr(model, str(v))
            print("Variable", v)
            for index in varobject:
                if varobject[index].value is not None:
                    if varobject[index].value > 0:
                        print("   ", index, varobject[index].value) 
        except:
            pass
print_vars()