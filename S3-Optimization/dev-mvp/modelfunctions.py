
import pandas as pd
import numpy as np
import itertools
import pyomo.environ as pe
import pyomo.gdp as gdp
import gurobipy
import pyomo.opt as po
import json

def create_sets(jsonfile = None):
    
    # Create default sets
    if not jsonfile:
        sets = {
            "S"     : ["DK"],  # Indices for supply geogpraphies
            "X"     : ["DK"],  # Indices for plant site locations
            "M"     : ["DK"],  # Indices for markets
            "T"       : [1],  # Indices for periods (10 years)
            "F"       : ["Bstraw", "Wstraw"],  # Indices for feedstocks
            "C"       : [],  # Indices for components
            "C_Solid" : ["Cellulose", "Hemicellulose", "Lignin"],  # Indices for solid components
            "C_Liquid": ["Xylose", "Glucose", "Lignin_Liquid"],  # Indices for liquid components
            
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
            
            "Q"       : [0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 50,  100]  # Indices for technology capacity levels
        }
        sets["C"] = sets["C_Solid"]  + sets["C_Liquid"] 
        sets["K"] = sets["K_prep"]  + sets["K_pret"] + sets["K_conv"] + sets["K_pur"]
        sets["FP"] = sets["F"]  + sets["P"] 

    # Read from file
    else:
        with open(jsonfile, 'r') as fp:
            sets = json.load(fp)

    return sets

# Funciton for create the dataframe from the cartesian product of specific sets
def create_parametric_df(sets, cols = []):
    df = pd.DataFrame(itertools.product(*sets.values()), columns=sets.keys())
    # drop duplicated rows
    df = df.drop_duplicates().reset_index(drop=True)
    # Create an empty value columns
    df.loc[:, 'Value'] = np.nan
    for col in cols:
        df.loc[:, col] = np.nan
    # Input 0 in impossible combinations (relevant for yields as you can't convert cellulose into cellulose)
    # Identify rows where a value in a column is equal to a value in another column for all columns
    equal_values_rows = df[df.apply(lambda row: len(set(row)) < len(row), axis = 1)]
    # Cahnge value for identified rows
    df.loc[equal_values_rows.index, "Value"] = 0
    return df

def create_compatability_dic(sets):
    
    # Dictionary to store the DataFrames
    dfs = {}

    # Compatability Supplier and Plant Site
    dfs['Compatibility_SupplierPlant'] = create_parametric_df({'S' : sets['S'], 'X' : sets['X']})

    # Compatability Plant Site and Market
    dfs['Compatibility_PlantMarket'] = create_parametric_df({'X' : sets['X'], 'M' : sets['M']})

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


def create_parametric_dic(sets, compat_dfs):
    # Dictionary to store the DataFrames
    dfs = {}

    # Generate DataFrames and store them in the dictionary
    dfs['Transportation_Cost_Product'] = create_parametric_df({'P': sets['P']})
    dfs['Transportation_Cost_Feedstock'] = create_parametric_df({'F': sets['F']})
    dfs['Technology_Cost'] = create_parametric_df({'K': sets['K']}, cols = ["ScaleFactor", "InstalationFactor", "BaseScale", "BaseScaleUnit", "Flow"])
    dfs['MinMax_Flows'] = create_parametric_df({'Flow': ['IN', 'OUT'], 'MinMax': ['Min', 'Max'], 'K': sets['K'], 'Q': sets['Q']})
    dfs['CEPCI'] = create_parametric_df({'T': sets['T']})

    dfs["Storage_Cost"] = create_parametric_df({'FP': sets['FP'], 'X': sets['X']})
    dfs['Utility_Consumption'] = create_parametric_df({'U': sets['U'], 'K': sets['K']})
    dfs['Byproduct_Production'] = create_parametric_df({'BP': sets['BP'], 'K': sets['K']})
    dfs['Feedstock_Composition'] = create_parametric_df({'C': sets['C'], 'F': sets['F']})
    
    # Compatability Supplier and Plant Site
    dfs['Compatibility_SupplierPlant'] = compat_dfs['Compatibility_SupplierPlant'].copy()
    # Compatability Plant Site and Market
    dfs['Compatibility_PlantMarket'] = compat_dfs['Compatibility_PlantMarket'].copy()
    
    # Distances Supplier and Plant Site
    dfs['Distances_SupplierPlant'] = compat_dfs['Compatibility_SupplierPlant'].copy()
    dfs['Distances_PlantMarket'] = compat_dfs['Compatibility_PlantMarket'].copy()

    dfs['Compatibility_Pretreatment'] = compat_dfs['Compatibility_Pretreatment'].copy()
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

    dfs['Compatibility_Conversion'] = df.copy()
    
    dfs['Yield_Hydrolysis'] = df[['K_conv', 'K_pret', 'C', 'P', 'Value']].drop_duplicates().reset_index(drop=True).sort_values(['K_conv', 'K_pret', 'C', 'P'])
    dfs['Yield_Conversion'] = df[[ 'P', 'C', 'K_conv', 'J', 'Value']].drop_duplicates().reset_index(drop=True).sort_values(['J', 'P', 'C', 'K_conv'])
    
    dfs['Solid_Loading'] = create_parametric_df({'K_conv' : sets['K_conv']})

    dfs['Compatibility_Purification'] = compat_dfs['Compatibility_Purification'].copy()
    
    dfs['Yield_Purification'] = dfs['Compatibility_Purification'][dfs['Compatibility_Purification'].Value != 0].copy()
    dfs['Depreciation_Factor'] = create_parametric_df({'T' : sets['T']})
    dfs['Operating_Cost'] = dfs['Depreciation_Factor'].copy()
    dfs['Maintenance_Cost'] = dfs['Depreciation_Factor'].copy()
    dfs['Demand'] = create_parametric_df({'P': sets['P'], 'M': sets['M'], 'T': sets['T']})
    dfs['Supply'] = create_parametric_df({'F': sets['F'], 'S': sets['S'], 'T': sets['T']})
    dfs['Price_F'] = create_parametric_df({'F': sets['F'], 'S': sets['S'], 'T': sets['T']})
    dfs['Price_U'] = create_parametric_df({'U': sets['U'], 'X': sets['X'], 'T': sets['T']})
    dfs['Price_P'] = create_parametric_df({'P': sets['P'], 'M': sets['M'], 'T': sets['T']})
    dfs['Economics'] = pd.DataFrame({'Indicator': ['SalvageValue', 'LifeCycle'], 'Value' : [0, 0]})

    return dfs

def create_compatability_excel(sets, xls_path = "Inputs\CompatabilityMatrices.xlsx"):
    dfs = create_compatability_dic(sets)
    writer = pd.ExcelWriter(xls_path, engine='xlsxwriter')
    for name, df in dfs.items():
        df.to_excel(writer, sheet_name=name, index = False)
    writer.close()
    writer.handles = None
    return dfs

def read_compatability_excel(sets, xls_path = "Inputs\CompatabilityMatricesInputs.xlsx"):

    dfs = create_compatability_dic(sets)
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
        df.Value = df.Value if name in ["Solid_Loading"] else df.Value.fillna(0)
        
        excel_dfs[name] = df
    
    xls.close()
    xls.handles = None
    
    return excel_dfs

def create_parametric_excel(sets, xls_path = "Inputs\ParametricTables.xlsx", compat = None):
    
    if compat is None:
        # Try to get input tables
        print("No compatabiity tables given. Attempting to retrieve from CompatabilityMatricesInputs.xlsx")
        compat_dfs = read_compatability_excel(sets)
    elif type(compat) is str:
        compat_dfs = read_compatability_excel(sets, compat)
    else:
        compat_dfs = compat
    
    dfs = create_parametric_dic(sets, compat_dfs)
    writer = pd.ExcelWriter(xls_path, engine='xlsxwriter')
    for name, df in dfs.items():
        df.to_excel(writer, sheet_name=name, index = False)
    writer.close()
    writer.handles = None
    return dfs

def read_parametric_excel(sets, xls_path, compat = None):
    
    if compat is None:
        # Try to get input tables
        print("No compatabiity tables given. Attempting to retrieve from CompatabilityMatricesInputs.xlsx")
        compat_dfs = read_compatability_excel(sets)
    elif type(compat) is str:
        compat_dfs = read_compatability_excel(sets, compat)
    else:
        compat_dfs = compat
        
    # Create empty dic to force parameters to be 0 when user accidently forgets certain inputs
    dfs = create_parametric_dic(sets, compat_dfs)
    
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
        if name == "Technology_Cost":
            df = df.drop_duplicates("K").reset_index(drop=True)
        # change all NA to 0 excpet for solid loading
        df.Value = df.Value if name in ["Solid_Loading"] else df.Value.fillna(0)
        
        excel_dfs[name] = df
    xls.close()
    xls.handles = None
    return excel_dfs

def print_vars(model):
    for v in model.component_objects(pe.Var, active=True):
        try:
            varobject = getattr(model, str(v))
            print("Variable", v)
            for index in varobject:
                if varobject[index].value is not None:
                    if varobject[index].value != 0:
                        print("   ", index, varobject[index].value) 
        except:
            pass