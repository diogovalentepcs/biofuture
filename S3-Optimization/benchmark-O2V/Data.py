#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 14:27:40 2021

@author: christianrappazzo
"""

import pandas as pd
import math 

file_path = 'C:\\Users\\digval\\Documents\\dtu-MasterThesis\\RappazoCode\\Superstructure_Data.xlsx'

df = pd.read_excel(file_path,header=0, index_col=0,sheet_name='Feedstock1')
Feedstocks1 = df.to_dict()   #[-]

df = pd.read_excel(file_path,header=0, index_col=0,sheet_name='Feedstock2')
Feedstocks2 = df.to_dict()   #[-]           
             
df = pd.read_excel(file_path,header=0, index_col=0,sheet_name='Bioreffinery capacity')
Bioreffinery_Capacity = df.to_dict() #[kg]

df = pd.read_excel(file_path,header=0, index_col=0,sheet_name='Preprocessing')
Preprocessing = df.to_dict() #[mm]

df = pd.read_excel(file_path,header=0, index_col=0,sheet_name='PretreatmentSize')
PretreatmentSize = df.to_dict() #[mm]

df = pd.read_excel(file_path,header=0, index_col=0,sheet_name='PretreatmentRecovery')
PretreatmentRecovery = df.to_dict()  #[-]

df = pd.read_excel(file_path,header=0, index_col=0,sheet_name='Bioconversions')
SHF_Cell_To_Glu_Conversion = df[['SHF_Cell_To_Glu']].to_dict() 
SHF_HC_To_Xyl_Conversion = df[['SHF_HC_To_Xyl']].to_dict() 
SSF_Cell_To_P_Conversion = df[['SSF_Cell_To_P']].to_dict() 
SSF_HC_To_P_Conversion = df[['SSF_HC_To_P']].to_dict() 
SHCF_Cell_To_Glu_Conversion = df[['SHCF_Cell_To_Glu']].to_dict() 
SHCF_HC_To_Xyl_Conversion = df[['SHCF_HC_To_Xyl']].to_dict() 
SSCF_Cell_To_P_Conversion = df[['SSCF_Cell_To_P']].to_dict() 
SSCF_HC_To_P_Conversion = df[['SSCF_HC_To_P']].to_dict()


df = pd.read_excel(file_path,header=0, index_col=0,sheet_name='PretreatmentByproducts')
PretreatmentByproducts = df.to_dict() #[kg of byproduct by kg of feedstock]

df = pd.read_excel(file_path,header=0, index_col=0,sheet_name='PretreatmentRawMaterials')
PretreatmentRawMaterials = df.to_dict() #[kg of raw materials by kg of feedstock]

df = pd.read_excel(file_path,header=0, index_col=0,sheet_name='FermentationStrains')
FermentationStrains = df.to_dict() 

df = pd.read_excel(file_path,header=0, index_col=0,sheet_name='Purification process')
Purification_process = df.to_dict() 

df = pd.read_excel(file_path,header=0, index_col=0,sheet_name='Manure Conversion process ')
Manure_Conversion_Process  = df.to_dict() 

df = pd.read_excel(file_path,header=0, index_col=0,sheet_name='Utility Price')
Utility_Price = df.to_dict()

df = pd.read_excel(file_path,header=0, index_col=0,sheet_name='PreprocessingRawMaterials')
PreprocessingRawMaterials = df.to_dict()

df = pd.read_excel(file_path,header=0, index_col=0,sheet_name='Lignin Conversion process')
Lignin_Conversion_Process = df.to_dict()

df = pd.read_excel(file_path,header=0, index_col=0,sheet_name='FermentationRawMaterials')
FermentationRawMaterials = df.to_dict()

df = pd.read_excel(file_path,header=0, index_col=0,sheet_name='MainProducts')
MainProducts = df.to_dict()

df = pd.read_excel(file_path,sheet_name='Efactor')
Efactor_comp = df['Components'].tolist()

df = pd.read_excel(file_path,sheet_name='Raw')
Raw = df['Raw'].tolist()
Energy = ['Cooling water  kwh','Heat kWh','Electricity kwh','LP-Steam kwh','Refrigeration  -15,5°C kwh','RNG  kWh','Steam 150 psig kg','Cooling kwh']




B1 = list(Feedstocks1.keys())
B2 = list(Feedstocks2.keys())
C = ['MILL3mm','MILL1.4mm','MILL10mm','MILL0.853mm','MILL0.15mm']   
D = ['Dilute Acid','LHW','AFEX','Steam Explosion','Lime','Ionic Liquid']  
 
E = ['Cellulose','Hemicellulose','Lignin_S','Lignin_L','Glucose','Xylose']
F = ['Succinic Acid','Glycolic Acid','Formic Acid','Acetic Acid','Phenolic','Furfural','HMF']
G = ['H2SO4 kg','Water kg','Ammonia kg','Ionic Liquid mat. kg','Lime mat. kg']
H = ['Ethanol Fermentation','Lactic Acid Fermentation','Succinic Acid Fermentation']
I = ['SHF','SSF','SHCF','SSCF']
J = list(FermentationStrains.keys())
K = list(Purification_process.keys())
L = list(Manure_Conversion_Process.keys())
Utilities = list(Utility_Price['Utility price $'].keys())
Fermentation_Utilities = list(FermentationRawMaterials['Dilute Acid'].keys())


Fermentation_SHCF_SSCF_non_compatibity = ['Lime','Ionic Liquid','LHW']

Fermentation_SSF_SSCF_non_compatibity = ['Lactic Acid Fermentation','Succinic Acid Fermentation']
                            
