#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  5 12:52:11 2021

@author: christianrappazzo
"""
import pandas as pd
import sys
from pyomo.environ import *
from pyomo.gdp import *
from Data import *

def create_n_solve():
    model = ConcreteModel()

    ##################
    #SET Declaration##
    ##################

    model.B1 = Set(initialize = B1)
    model.B2 = Set(initialize = B2)
    model.C = Set(initialize = C)
    model.D = Set(initialize = D)
    model.E = Set(initialize = E)
    model.F = Set(initialize = F)
    model.G = Set(initialize = G)
    model.H = Set(initialize = H)
    model.I = Set(initialize = I)
    model.J = Set(initialize = J)
    model.K = Set(initialize = K)
    model.L = Set(initialize = L)

    model.Utilities = Set(initialize = Utilities)
    model.Fermentation_Utilities = Set(initialize = Fermentation_Utilities)
    model.Efactor_comp = Set(initialize = Efactor_comp)


    ########################
    #PARAMETERS Declaration#
    ########################

    # Load feedstock availability into a model parameter for later access
    model.Bioreffinery_Capacity2 = Param(model.H, initialize = lambda model,h:Bioreffinery_Capacity['Capacity'][h])

    model.Bioreffinery_Capacity = Param(initialize = sum(model.Bioreffinery_Capacity2[h] for h in model.H))

    # Load output form of preprocessing option c
    model.OUT_preprocessing = Param(model.C, initialize = lambda model,c:Preprocessing[c]['Size'])

    # Load input form of preatreatment option d
    model.IN_pretreatment = Param(model.D, initialize = lambda model,d:PretreatmentSize[d]['Max Size'])

    # Load fraction of cellulose in biomass b_1
    model.alpha_cel = Param(model.B1,initialize = lambda model,b1:Feedstocks1[b1]['cellulose'] )

    # Load fraction of hemicellulose in biomass b_1
    model.alpha_HC = Param(model.B1,initialize = lambda model,b1:Feedstocks1[b1]['HC'] )

    # Load fraction of lignin in biomass b_1
    model.alpha_Lig = Param(model.B1,initialize = lambda model,b1:Feedstocks1[b1]['lignin'] )

    # Load recovery factor of component e using preatreatment option d
    model.eta_S2 = Param(model.D,model.E,initialize = lambda model,d,e:PretreatmentRecovery[d][e])

    # Load generation factor of by-product f using preatreatment option d
    model.eta_S2_BP = Param(model.D,model.F,initialize = lambda model,d,f:PretreatmentByproducts[d][f])

    def pi_S3_rule(model,d,h,i):
        x = 1
        if i == 'SSCF' or i == 'SHCF':
            if d in Fermentation_SHCF_SSCF_non_compatibity:
                x=0
        if i == 'SSF' or i == 'SSCF':
            if h in Fermentation_SSF_SSCF_non_compatibity:
                x=0       
        return x
    model.pi_S3 = Param(model.D,model.H,model.I, initialize = pi_S3_rule)

    model.SHF_Cell_To_Glu_Conversion = Param(model.D, initialize = lambda model,d:SHF_Cell_To_Glu_Conversion['SHF_Cell_To_Glu'][d])

    model.SHF_HC_To_Xyl_Conversion = Param(model.D, initialize = lambda model,d:SHF_HC_To_Xyl_Conversion['SHF_HC_To_Xyl'][d])

    def pi_S3_J_S_rule(model,h,i,j):
        x=''
        if i == 'SSF' or i == 'SHF':
            x+='_S'
        return FermentationStrains[j][h]*FermentationStrains[j][i+x]
    model.pi_S3_J_S = Param(model.H,model.I,model.J, initialize = pi_S3_J_S_rule)

    def pi_S3_J_L_rule(model,h,i,j):
        x=''
        if i == 'SSF' or i == 'SHF':
            x='_L'
        return FermentationStrains[j][h]*FermentationStrains[j][i+x]
    model.pi_S3_J_L = Param(model.H,model.I,model.J, initialize = pi_S3_J_L_rule)

    model.eta_S3_P_Glu = Param(model.D,model.J, initialize = lambda model,d,j:FermentationStrains[j]['Glucose to product'])

    model.eta_S3_P_Xyl = Param(model.D,model.J, initialize = lambda model,d,j:FermentationStrains[j]['Xylose to product'])

    model.eta_S3_P_Cel_SSF = Param(model.D,model.J, initialize = lambda model,d,j:FermentationStrains[j]['Glucose to product']*SSF_Cell_To_P_Conversion['SSF_Cell_To_P'][d])

    model.eta_S3_P_HC_SSF = Param(model.D,model.J, initialize = lambda model,d,j:FermentationStrains[j]['Xylose to product']*SSF_HC_To_P_Conversion['SSF_HC_To_P'][d])

    model.eta_S3_P_Cel_SSCF = Param(model.D,model.J, initialize = lambda model,d,j:FermentationStrains[j]['Glucose to product']*SSCF_Cell_To_P_Conversion['SSCF_Cell_To_P'][d])

    model.eta_S3_P_HC_SSCF = Param(model.D,model.J, initialize = lambda model,d,j:FermentationStrains[j]['Xylose to product']*SSCF_HC_To_P_Conversion['SSCF_HC_To_P'][d])

    model.eta_S3_P_Xyl_SSCF = Param(model.D,model.J, initialize = lambda model,d,j:FermentationStrains[j]['Xylose to product cofermentation'])

    model.SL_S3= Param(initialize = 0.2)

    model.SHCF_Cell_To_Glu_Conversion = Param(model.D, initialize = lambda model,d:SHCF_Cell_To_Glu_Conversion['SHCF_Cell_To_Glu'][d])

    model.SHCF_HC_To_Xyl_Conversion = Param(model.D, initialize = lambda model,d:SHCF_HC_To_Xyl_Conversion['SHCF_HC_To_Xyl'][d])


    def pi_S4_rule(model,h,k):
        return Purification_process[k][h]
    model.pi_S4 = Param(model.H,model.K, initialize = pi_S4_rule)

    def eta_S4_rule(model,k):
        return Purification_process[k]['Recovery']
    model.eta_S4 = Param(model.K, initialize = eta_S4_rule)

    model.Utility_Price = Param(model.Utilities,  initialize = lambda model,x:Utility_Price['Utility price $'][x])





    model.Preprocessing_Raw_Materials = Param(model.C, model.Utilities, initialize = lambda model,c,x:PreprocessingRawMaterials[c][x])

    model.Preatreatment_Raw_Materials = Param(model.D, model.Utilities, initialize = lambda model,d,x:PretreatmentRawMaterials[d][x])

    model.Purification_Raw_Materials = Param(model.K, model.Utilities, initialize = lambda model,k,x:Purification_process[k][x])

    model.Manure_Conversion_Raw_Materials = Param(model.L, model.Utilities, initialize = lambda model,l,x:Manure_Conversion_Process[l][x])

    model.Lignin_Conversion_Raw_Materials = Param(model.Utilities, initialize = lambda model,x:Lignin_Conversion_Process['Process 1'][x])

    model.Fermentation_Raw_Materials = Param(model.D,model.Fermentation_Utilities, initialize = lambda model,d,x:FermentationRawMaterials[d][x])


    #######################
    #VARIABLES Declaration#
    #######################

    #Creation of variable expressing the amount of feedstock b1 sent to the preprocessing option c  
    def F_S1_bound_rule(model,*argv):
        return (0,20*model.Bioreffinery_Capacity)
    model.F_S1 = Var(model.B1,model.C, bounds = F_S1_bound_rule)

    #Creation of variable expressing  the amount of feedstock b1 preprocessed using preprocessing option c    sent to the preatreatment option d   
    def Flow_bound_rule(model,*argv):
        return (0,50*model.Bioreffinery_Capacity)

    def RM_bound_rule(model,*argv):
        return (-50*model.Bioreffinery_Capacity,50*model.Bioreffinery_Capacity)

    model.F_S1ToS2 = Var(model.B1,model.C,model.D,bounds = Flow_bound_rule)

    #Creation of the binary variable describing the choice of the biomass of type b1
    model.Y_S1_B1 = Var(model.B1,domain = Binary)

    #Creation of the binary variable describing the preprocessing option c    
    model.Y_S1_preprocessing = Var(model.B1,model.C,domain = Binary)

    #Creation of the binary variable describing the choice of preatreatment option d
    model.Y_S2 = Var(model.B1,model.D,domain = Binary)

    #Creation of variable expressing  the amount of component e recovered from feedstock $b_1$ undergoing  preatreatment option d 
    model.F_S2 = Var(model.B1,model.D,model.E,bounds = Flow_bound_rule)

    #Creation of variable expressing  the amount of by-product f formed from feedstock $b_1$ undergoing  preatreatment option d 
    model.F_S2_BP = Var(model.B1,model.D,model.F,bounds = Flow_bound_rule)

    model.F_S2_RM = Var(model.B1,model.D,model.G,bounds = RM_bound_rule) 

    model.F_S2_Solid = Var(model.B1,model.D,bounds = Flow_bound_rule) 

    model.F_S2_Liquid = Var(model.B1,model.D,bounds = Flow_bound_rule) 

    model.F_S3_Solid = Var(model.B1,model.D,model.H,model.I,bounds = Flow_bound_rule)

    model.F_S3_Liquid = Var(model.B1,model.D,model.H,model.I,bounds = Flow_bound_rule)

    model.Y_S3_H_Liquid = Var(model.B1,model.H,domain = Binary)

    model.Y_S3_H_Solid = Var(model.B1,model.H,domain = Binary)

    model.Y_S3_H = Var(model.B1,model.H,domain = Binary)

    model.Y_S3_I = Var(model.B1,model.H,model.I,domain = Binary)

    model.F_S3_Glu = Var(model.B1,model.D,model.H,model.I,bounds = Flow_bound_rule)

    model.F_S3_Xyl = Var(model.B1,model.D,model.H,model.I,bounds = Flow_bound_rule)

    model.F_S3_Xyl_L = Var(model.B1,model.D,model.H,model.I,bounds = Flow_bound_rule)

    model.Y_S3_J_SHF1 = Var(model.B1,model.H,model.J,domain = Binary)

    model.Y_S3_J_SHF2 = Var(model.B1,model.H,model.J,domain = Binary)

    model.F_S3_Glu_F = Var(model.B1,model.D,model.H,model.I,model.J,bounds = Flow_bound_rule)

    model.F_S3_Xyl_F = Var(model.B1,model.D,model.H,model.I,model.J,bounds = Flow_bound_rule)

    model.F_S3_Xyl_L_F = Var(model.B1,model.D,model.H,model.I,model.J,bounds = Flow_bound_rule)

    model.F_S3_P_Glu = Var(model.B1,model.D,model.H,model.I,model.J,bounds = Flow_bound_rule)

    model.F_S3_P_Xyl = Var(model.B1,model.D,model.H,model.I,model.J,bounds = Flow_bound_rule)

    model.F_S3_P_Xyl_L = Var(model.B1,model.D,model.H,model.I,model.J,bounds = Flow_bound_rule)
    
    model.Y_S3_J_SSF1 = Var(model.B1,model.H,model.J,domain = Binary)

    model.Y_S3_J_SSF2 = Var(model.B1,model.H,model.J,domain = Binary)

    model.F_S3_Cel = Var(model.B1,model.D,model.H,model.I,bounds = Flow_bound_rule)

    model.F_S3_HC = Var(model.B1,model.D,model.H,model.I,bounds = Flow_bound_rule)

    model.F_S3_Cel_F = Var(model.B1,model.D,model.H,model.I,model.J, bounds = Flow_bound_rule)

    model.F_S3_HC_F = Var(model.B1,model.D,model.H,model.I,model.J, bounds = Flow_bound_rule)

    model.F_S3_P_Cel = Var(model.B1,model.D,model.H,model.I,model.J,bounds = Flow_bound_rule)

    model.F_S3_P_HC = Var(model.B1,model.D,model.H,model.I,model.J,bounds = Flow_bound_rule)

    model.F_S3_Water = Var(model.B1,model.H,model.I,bounds = Flow_bound_rule)

    model.F_S3_Tot = Var(model.B1,model.H,model.I,bounds = Flow_bound_rule)

    model.F_S3_Liquid_R = Var(model.B1,model.H,model.I,bounds = Flow_bound_rule)

    model.Y_S3_J_SHCF = Var(model.B1,model.H,model.J,domain = Binary)

    model.Y_S3_J_SSCF = Var(model.B1,model.H,model.J,domain = Binary) 

    model.F_S3_P = Var(model.H,bounds = Flow_bound_rule)

    model.F_S4 = Var(model.H,model.K,bounds = Flow_bound_rule)

    model.F_S4_P_R =Var(model.H,model.K,bounds = Flow_bound_rule)

    model.F_S4_P_Tot = Var(model.H,bounds = Flow_bound_rule)

    model.Y_S4 = Var(model.H,model.K,domain = Binary) 
            
    model.F_S4_Lignin_R = Var(bounds = Flow_bound_rule)

    model.Y_S4_Lignin = Var(domain = Binary)

    model.F_Preprocessing_RM = Var(model.Utilities, bounds = RM_bound_rule)

    model.F_Pretreatment_RM = Var(model.Utilities, bounds = RM_bound_rule)

    model.F_Purification_RM = Var(model.Utilities, bounds = RM_bound_rule)

    model.F_Fermentation_RM = Var(model.Fermentation_Utilities, bounds = RM_bound_rule)

    model.F_Manure_Conversion_RM = Var(model.Utilities, bounds = RM_bound_rule)

    model.F_Lignin_Conversion_RM = Var(model.Utilities, bounds = RM_bound_rule)

    model.Preprocessing_RM_Price = Var(bounds = RM_bound_rule)

    model.Pretreatment_RM_Price = Var(bounds = RM_bound_rule)

    model.Purification_RM_Price = Var(bounds = RM_bound_rule)

    model.Fermentation_RM_Price = Var(bounds = RM_bound_rule)

    model.Manure_Conversion_RM_Price = Var(bounds = RM_bound_rule)



    model.Y_S4_Manure = Var(domain = Binary)
    model.F_Manure = Var(model.L,bounds = F_S1_bound_rule)
    model.Y_S4_Manure_L = Var(model.L, domain = Binary)

    model.Lignin_Conversion_RM_Price = Var(bounds = RM_bound_rule)

    model.Feedstock_Price = Var(bounds = RM_bound_rule)

    model.MP_Value=Var(bounds = RM_bound_rule)

    model.E_Factor=Var(bounds = RM_bound_rule)


    ##########################
    #Preprocessing constraints
    ##########################



    #EQUATION 1
    def Feedstock_Disjuncts_rule(b,b1,flag):
        model = b.model()
        if flag ==1:
            b.c1 = ConstraintList()  
            b.c1.add(model.Y_S1_B1[b1] == 1)
            b.c1.add(sum(model.F_S1[b1,c] for c in model.C)>=0)
        if flag ==0:
            b.c0 = ConstraintList()  
            b.c0.add(model.Y_S1_B1[b1] == 0)
            b.c0.add(sum(model.F_S1[b1,c] for c in model.C)==0)

    model.Feedstock_Disjuncts = Disjunct(model.B1,[0,1], rule = Feedstock_Disjuncts_rule)
    model.Feedstock_Disjunction = Disjunction(model.B1)
    for b1 in model.B1:
        model.Feedstock_Disjunction[b1] = [model.Feedstock_Disjuncts[b1,i] for i in [0,1]]

    #EQUATION 2
    def Preprocessing_selection_rule(model,b1): 
        return sum(model.Y_S1_preprocessing[b1,c] for c in model.C) == model.Y_S1_B1[b1]
    model.Preprocessing_selection = Constraint(model.B1,rule=Preprocessing_selection_rule)

        
    #EQUATION 3
    def Preprocessing_Disjuncts_rule(b,b1,c,flag):
        model = b.model()
        if flag ==1:
            b.c1 = ConstraintList()
            b.c1.add(model.Y_S1_preprocessing[b1,c]==1)
            b.c1.add(model.F_S1[b1,c]>=0.1)    
        if flag ==0:
            b.c0 = ConstraintList()
            b.c0.add(model.Y_S1_preprocessing[b1,c]==0)
            b.c0.add(model.F_S1[b1,c]==0)    
    model.Preprocessing_Disjuncts = Disjunct(model.B1,model.C,[0,1],rule = Preprocessing_Disjuncts_rule )
    model.Preprocessing_Disjunction = Disjunction(model.B1,model.C)
    for b1 in model.B1:
        for c in model.C:
            model.Preprocessing_Disjunction[b1,c] = [model.Preprocessing_Disjuncts[b1,c,i] for i in [0,1]]
        

    #########################
    #Pretreatment constraints
    #########################


    #EQUATION 4
    def Feedstock_to_preprocessing_constraint_rule(model,b1,c):
        return sum( model.F_S1ToS2[b1,c,d] for d in model.D ) == model.F_S1[b1,c]
    model.Feedstock_utilization_constraint = Constraint(model.B1, model.C, rule=Feedstock_to_preprocessing_constraint_rule)

    #EQUATION 5
    def Pretreatment_selection_rule(model,b1): 
        return sum(model.Y_S2[b1,d] for d in model.D) == model.Y_S1_B1[b1]
    model.Pretreatment_selection = Constraint(model.B1,rule=Pretreatment_selection_rule)

    #EQUATION 6
    def Pretreatment_Disjuncts_rule(b,b1,d,flag):
        model = b.model()
        if flag ==1:
            b.c1 = ConstraintList()  
            b.c1.add(model.Y_S2[b1,d]==1)
            b.c1.add(sum(model.F_S1ToS2[b1,c,d] for c in model.C)>=0.1)
        if flag ==0:
            b.c0 = ConstraintList()  
            b.c0.add(model.Y_S2[b1,d]==0)
            b.c0.add(sum(model.F_S1ToS2[b1,c,d] for c in model.C)==0)
    model.Pretreatment_Disjuncts = Disjunct(model.B1,model.D,[0,1], rule =  Pretreatment_Disjuncts_rule)
    model.Pretreatment_Disjunction=Disjunction(model.B1,model.D)
    for b1 in model.B1:
        for d in model.D:
            model.Pretreatment_Disjunction[b1,d] = [model.Pretreatment_Disjuncts[b1,d,i] for i in [0,1]]


    #EQUATION 7
    def Compatible_S1_S2_rule(model,b1):
        OUT_S1 =  sum(model.Y_S1_preprocessing[b1,c]*model.OUT_preprocessing[c] for c in model.C )
        IN_S2 =  sum(model.Y_S2[b1,d]*model.IN_pretreatment[d] for d in model.D )
        return OUT_S1 == IN_S2
    model.compatible_S1_S2 = Constraint(model.B1,rule=Compatible_S1_S2_rule)


    #EQUATION 8
    def Pretreatment_cellulose_rule(model,b1,d):
        return model.F_S2[b1,d,'Cellulose'] == sum(model.F_S1ToS2[b1,c,d]*model.alpha_cel[b1]*model.eta_S2[d,'Cellulose'] for c in model.C)
    model.Pretreatment_cellulose = Constraint(model.B1,model.D,rule = Pretreatment_cellulose_rule)

    #EQUATION 9
    def Pretreatment_hemicellulose_rule(model,b1,d):
        return model.F_S2[b1,d,'Hemicellulose'] == sum(model.F_S1ToS2[b1,c,d]*model.alpha_HC[b1]*model.eta_S2[d,'Hemicellulose'] for c in model.C)
    model.Pretreatment_hemicellulose = Constraint(model.B1,model.D,rule = Pretreatment_hemicellulose_rule)

    #EQUATION 10
    def Pretreatment_lignin_S_rule(model,b1,d):
        return model.F_S2[b1,d,'Lignin_S'] == sum(model.F_S1ToS2[b1,c,d]*model.alpha_Lig[b1]*model.eta_S2[d,'Lignin_S'] for c in model.C)
    model.Pretreatment_lignin_S = Constraint(model.B1,model.D,rule = Pretreatment_lignin_S_rule)

    #EQUATION 11
    def Pretreatment_lignin_L_rule(model,b1,d):
        return model.F_S2[b1,d,'Lignin_L'] == sum(model.F_S1ToS2[b1,c,d]*model.alpha_Lig[b1]*model.eta_S2[d,'Lignin_L'] for c in model.C)
    model.Pretreatment_lignin_L = Constraint(model.B1,model.D,rule = Pretreatment_lignin_L_rule)

    #EQUATION 12
    def Pretreatment_glucose_rule(model,b1,d):
        return model.F_S2[b1,d,'Glucose'] == sum(model.F_S1ToS2[b1,c,d]*model.alpha_cel[b1]*model.eta_S2[d,'Glucose'] for c in model.C)
    model.Pretreatment_glucose = Constraint(model.B1,model.D,rule = Pretreatment_glucose_rule)

    #EQUATION 13
    def Pretreatment_xylose_rule(model,b1,d):
        return model.F_S2[b1,d,'Xylose'] == sum(model.F_S1ToS2[b1,c,d]*model.alpha_HC[b1]*model.eta_S2[d,'Xylose'] for c in model.C)
    model.Pretreatment_xylose = Constraint(model.B1,model.D,rule = Pretreatment_xylose_rule)

    #EQUATION 14
    def Pretreatment_byproducts_rule(model,b1,d,f):
        return model.F_S2_BP[b1,d,f] == sum(model.F_S1ToS2[b1,c,d]*model.eta_S2_BP[d,f] for c in model.C)
    model.Pretreatment_byproducts = Constraint(model.B1,model.D,model.F,rule = Pretreatment_byproducts_rule)


    #EQUATION 15
    def Pretreatment_raw_materials_rule(model,b1,d,g):
        return model.F_S2_RM[b1,d,g] == sum(model.F_S1ToS2[b1,c,d]*model.Preatreatment_Raw_Materials[d,g] for c in model.C)
    model.Pretreatment_raw_materials = Constraint(model.B1,model.D,model.G,rule = Pretreatment_raw_materials_rule)


    c
    #EQUATION 16
    def Pretreatment_solid_rule(model,b1,d):
        return model.F_S2_Solid[b1,d] == model.F_S2[b1,d,'Cellulose'] +model.F_S2[b1,d,'Hemicellulose']+model.F_S2[b1,d,'Lignin_S']
    model.Pretreatment_solid = Constraint(model.B1,model.D,rule = Pretreatment_solid_rule)

    #EQUATION 17
    def Pretreatment_liquid_rule(model,b1,d):
        return model.F_S2_Liquid[b1,d] == model.F_S2[b1,d,'Glucose'] + model.F_S2[b1,d,'Xylose']+model.F_S2[b1,d,'Lignin_L'] + sum(model.F_S2_BP[b1,d,f] for f in model.F) + sum(model.F_S2_RM[b1,d,g] for g in model.G)
    model.Pretreatment_liquid = Constraint(model.B1,model.D,rule = Pretreatment_liquid_rule)


    #########################
    #Fermentation constraints
    #########################

    #EQUATION 18
    def Fermentation_liquid_rule(model,b1,d):
        return model.F_S2_Liquid[b1,d] == sum(sum(model.F_S3_Liquid[b1,d,h,i] for i in model.I) for h in model.H)
    model.Fermentation_liquid = Constraint(model.B1,model.D,rule = Fermentation_liquid_rule)

    #EQUATION 19
    def Fermentation_solid_rule(model,b1,d):
        return model.F_S2_Solid[b1,d] == sum(sum(model.F_S3_Solid[b1,d,h,i] for i in model.I) for h in model.H)
    model.Fermentation_solid = Constraint(model.B1,model.D,rule = Fermentation_solid_rule)


    #EQUATION 20
    def D_H_I_compatibility_rule1(model,b1,d,h,i):
        return model.F_S3_Solid[b1,d,h,i] == model.pi_S3[d,h,i]*model.F_S3_Solid[b1,d,h,i]
    model.D_H_I_compatibility1 = Constraint(model.B1,model.D,model.H,model.I,rule = D_H_I_compatibility_rule1)

    #EQUATION 21
    def D_H_I_compatibility_rule2(model,b1,d,h,i):
        return model.F_S3_Liquid[b1,d,h,i] == model.pi_S3[d,h,i]*model.F_S3_Liquid[b1,d,h,i]
    model.D_H_I_compatibility2 = Constraint(model.B1,model.D,model.H,model.I,rule = D_H_I_compatibility_rule2)

    ##############################
    #SHF Fermentation constraints#
    ##############################

    #Equation 22
    def SHF_Solid_loading_rule(model,b1,h):
        return sum(model.F_S3_Solid[b1,d,h,'SHF'] for d in model.D) == model.SL_S3()*(sum(model.F_S3_Solid[b1,d,h,'SHF'] for d in model.D) +  model.F_S3_Water[b1,h,'SHF'])
    model.SHF_Solid_loading_constraint = Constraint(model.B1,model.H, rule = SHF_Solid_loading_rule)


    #Equation 23
    def SHF_Stream_rule(b,b1,d,h,flag):
        model = b.model()
        if flag ==1:
            b.c1 = ConstraintList()  
            b.c1.add(model.Y_S3_I[b1,h,'SHF']==1)
            b.c1.add(model.F_S3_Solid[b1,d,h,'SHF'] == model.F_S2_Solid[b1,d])
            b.c1.add(model.F_S3_Liquid[b1,d,h,'SHF']== model.F_S2_Liquid[b1,d])
            b.c1.add(model.F_S3_Glu[b1,d,h,'SHF'] ==  model.SHF_Cell_To_Glu_Conversion[d]*model.F_S2[b1,d,'Cellulose'])
            b.c1.add(model.F_S3_Xyl[b1,d,h,'SHF'] ==  model.SHF_HC_To_Xyl_Conversion[d]*model.F_S2[b1,d,'Hemicellulose'])
            b.c1.add(model.F_S3_Xyl_L[b1,d,h,'SHF'] == model.F_S2[b1,d,'Xylose'])
        if flag ==0:
            b.c0 = ConstraintList()  
            b.c0.add(model.Y_S3_I[b1,h,'SHF']==0)
            b.c0.add(model.F_S3_Solid[b1,d,h,'SHF']==0)
            b.c0.add(model.F_S3_Liquid[b1,d,h,'SHF']==0)
            b.c0.add(model.F_S3_Glu[b1,d,h,'SHF'] == 0)
            b.c0.add(model.F_S3_Xyl[b1,d,h,'SHF'] ==  0)
            b.c0.add(model.F_S3_Xyl_L[b1,d,h,'SHF'] == 0)

    model.SHF_Stream_Disjuncts = Disjunct(model.B1,model.D,model.H,[0,1], rule =  SHF_Stream_rule )
    model.SHF_Stream_Disjunction = Disjunction(model.B1,model.D,model.H)
    for b1 in model.B1:
        for d in model.D:
            for h in model.H:
                model.SHF_Stream_Disjunction[b1,d,h] = [model.SHF_Stream_Disjuncts[b1,d,h,i] for i in [0,1]]
        
    #EQUATION 24
    def SHF_strain_selection1_rule(model,b1,h):
        return sum(model.Y_S3_J_SHF1[b1,h,j] for j in model.J) == model.Y_S3_I[b1,h,'SHF']
    model.SHF_strain_selection1_constraint = Constraint(model.B1, model.H, rule = SHF_strain_selection1_rule)

    #EQUATION 25
    def SHF_strain_selection2_rule(model,b1,h):
        return sum(model.Y_S3_J_SHF2[b1,h,j] for j in model.J) == model.Y_S3_I[b1,h,'SHF']
    model.SHF_strain_selection2_constraint = Constraint(model.B1, model.H, rule = SHF_strain_selection2_rule)



    #EQUATION 26
    def SHF_glusose_fermentation_rule(model,b1,d,h):
        return sum(model.F_S3_Glu_F[b1,d,h,'SHF',j] for j in model.J) == model.F_S3_Glu[b1,d,h,'SHF']
    model.SHF_glusose_fermentation_constraint = Constraint(model.B1, model.D, model.H, rule = SHF_glusose_fermentation_rule)

    #EQUATION 27
    def SHF_xylose_fermentation_rule(model,b1,d,h):
        return sum(model.F_S3_Xyl_F[b1,d,h,'SHF',j] for j in model.J) == model.F_S3_Xyl[b1,d,h,'SHF']
    model.SHF_xylose_fermentation_constraint = Constraint(model.B1, model.D, model.H, rule = SHF_xylose_fermentation_rule)

    #EQUATION 28
    def SHF_xylose_L_fermentation_rule(model,b1,d,h):
        return sum(model.F_S3_Xyl_L_F[b1,d,h,'SHF',j] for j in model.J) == model.F_S3_Xyl_L[b1,d,h,'SHF']
    model.SHF_xylose_L_fermentation_constraint = Constraint(model.B1, model.D, model.H, rule = SHF_xylose_L_fermentation_rule)

    #EQUATION 29
    def SHF_Fermentation_Strain_compatibility_rule1(model,b1,d,h,j):
        return model.F_S3_Xyl_F[b1,d,h,'SHF',j] == model.pi_S3_J_S[h,'SHF',j]*model.F_S3_Xyl_F[b1,d,h,'SHF',j]
    model.SHF_Fermentation_Strain_compatibility1 = Constraint(model.B1,model.D,model.H,model.J,rule = SHF_Fermentation_Strain_compatibility_rule1)

    #EQUATION 30
    def SHF_Fermentation_Strain_compatibility_rule2(model,b1,d,h,j):
        return model.F_S3_Glu_F[b1,d,h,'SHF',j] == model.pi_S3_J_S[h,'SHF',j]*model.F_S3_Glu_F[b1,d,h,'SHF',j]
    model.SHF_Fermentation_Strain_compatibility2 = Constraint(model.B1,model.D,model.H,model.J,rule = SHF_Fermentation_Strain_compatibility_rule2)

    #EQUATION 31
    def SHF_Fermentation_Strain_compatibility_rule3(model,b1,d,h,j):
        return model.F_S3_Xyl_L_F[b1,d,h,'SHF',j] == model.pi_S3_J_L[h,'SHF',j]*model.F_S3_Xyl_L_F[b1,d,h,'SHF',j]
    model.SHF_Fermentation_Strain_compatibility3 = Constraint(model.B1,model.D,model.H,model.J,rule = SHF_Fermentation_Strain_compatibility_rule3)



    #EQUATION 32
    def SHF_Fermentation_Strain_Xyl_L_Disjuncts_rule(b,b1,d,h,j,flag):
        model = b.model()
        if flag ==1:
            b.c1 = ConstraintList()  
            b.c1.add(model.Y_S3_J_SHF1[b1,h,j]==1)
            b.c1.add(model.F_S3_P_Xyl_L[b1,d,h,'SHF',j] == model.eta_S3_P_Xyl[d,j]*model.F_S3_Xyl_L_F[b1,d,h,'SHF',j])
        if flag ==0:
            b.c0 = ConstraintList()  
            b.c0.add(model.Y_S3_J_SHF1[b1,h,j]==0)
            b.c0.add(model.F_S3_P_Xyl_L[b1,d,h,'SHF',j] == 0)

    model.SHF_Fermentation_Strain_Xyl_L_Disjuncts = Disjunct(model.B1,model.D,model.H,model.J,[0,1], rule = SHF_Fermentation_Strain_Xyl_L_Disjuncts_rule)
    model.SHF_Fermentation_Strain_Xyl_L_Disjunction = Disjunction(model.B1,model.D,model.H,model.J)
    for b1 in model.B1:
        for d in model.D:
            for h in model.H:
                for j in model.J:
                    model.SHF_Fermentation_Strain_Xyl_L_Disjunction[b1,d,h,j] = [model.SHF_Fermentation_Strain_Xyl_L_Disjuncts[b1,d,h,j,i] for i in [0,1]]
                    
    #EQUATION 33
    def SHF_Fermentation_Strain_Glu_Xyl_Disjuncts_rule(b,b1,d,h,j,flag):
        model = b.model()
        if flag ==1:
            b.c1 = ConstraintList()  
            b.c1.add(model.Y_S3_J_SHF2[b1,h,j]==1)
            b.c1.add(model.F_S3_P_Xyl[b1,d,h,'SHF',j] == model.eta_S3_P_Xyl[d,j]*model.F_S3_Xyl_F[b1,d,h,'SHF',j])
            b.c1.add(model.F_S3_P_Glu[b1,d,h,'SHF',j] == model.eta_S3_P_Glu[d,j]*model.F_S3_Glu_F[b1,d,h,'SHF',j])
        if flag ==0:
            b.c0 = ConstraintList()  
            b.c0.add(model.Y_S3_J_SHF2[b1,h,j]==0)
            b.c0.add(model.F_S3_P_Xyl[b1,d,h,'SHF',j] == 0)
            b.c0.add(model.F_S3_P_Glu[b1,d,h,'SHF',j] == 0)

    model.SHF_Fermentation_Strain_Glu_Xyl_Disjuncts = Disjunct(model.B1,model.D,model.H,model.J,[0,1], rule = SHF_Fermentation_Strain_Glu_Xyl_Disjuncts_rule)
    model.SHF_Fermentation_Strain_Glu_Xyl_Disjunction = Disjunction(model.B1,model.D,model.H,model.J)
    for b1 in model.B1:
        for d in model.D:
            for h in model.H:
                for j in model.J:
                    model.SHF_Fermentation_Strain_Glu_Xyl_Disjunction[b1,d,h,j] = [model.SHF_Fermentation_Strain_Glu_Xyl_Disjuncts[b1,d,h,j,i] for i in [0,1]]
                    
    ##############################
    #SSF Fermentation constraints
    ##############################


    #EQUATION 34
    def SSF_Solid_loading_rule(model,b1,h):
        return sum(model.F_S3_Solid[b1,d,h,'SSF'] for d in model.D) == model.SL_S3()*(sum(model.F_S3_Solid[b1,d,h,'SSF'] for d in model.D) +  model.F_S3_Water[b1,h,'SSF'])
    model.SSF_Solid_loading_constraint = Constraint(model.B1,model.H, rule = SSF_Solid_loading_rule)


    #EQUATION 35
    def SSF_Stream_rule(b,b1,d,h,flag):
        model = b.model()
        if flag ==1:
            b.c1 = ConstraintList()         
            b.c1.add(model.Y_S3_I[b1,h,'SSF']==1)
            b.c1.add(model.F_S3_Solid[b1,d,h,'SSF'] == model.F_S2_Solid[b1,d])
            b.c1.add(model.F_S3_Liquid[b1,d,h,'SSF'] == model.F_S2_Liquid[b1,d])
            b.c1.add(model.F_S3_Cel[b1,d,h,'SSF'] ==  model.F_S2[b1,d,'Cellulose'])
            b.c1.add(model.F_S3_HC[b1,d,h,'SSF'] ==  model.F_S2[b1,d,'Hemicellulose'])
            b.c1.add(model.F_S3_Xyl_L[b1,d,h,'SSF'] == model.F_S2[b1,d,'Xylose'])
        if flag ==0:
            b.c0 = ConstraintList()  
            b.c0.add(model.Y_S3_I[b1,h,'SSF']==0)
            b.c0.add(model.F_S3_Solid[b1,d,h,'SSF'] == 0)
            b.c0.add(model.F_S3_Liquid[b1,d,h,'SSF'] == 0)
            b.c0.add(model.F_S3_Cel[b1,d,h,'SSF'] == 0)
            b.c0.add(model.F_S3_HC[b1,d,h,'SSF'] ==  0)
            b.c0.add(model.F_S3_Xyl_L[b1,d,h,'SSF'] == 0)

    model.SSF_Stream_Disjuncts = Disjunct(model.B1,model.D,model.H,[0,1],rule =SSF_Stream_rule )
    model.SSF_Stream_Disjunction = Disjunction(model.B1,model.D,model.H)
    for b1 in model.B1:
        for d in model.D:
            for h in model.H:
                model.SSF_Stream_Disjunction[b1,d,h] = [model.SSF_Stream_Disjuncts[b1,d,h,i] for i in [0,1]]


    #EQUATION 36
    def SSF_strain_selection1_rule(model,b1,h):
        return sum(model.Y_S3_J_SSF1[b1,h,j] for j in model.J) == model.Y_S3_I[b1,h,'SSF']
    model.SSF_strain_selection1_constraint = Constraint(model.B1, model.H, rule = SSF_strain_selection1_rule)

    #EQUATION 37
    def SSF_xylose_L_fermentation_rule(model,b1,d,h):
        return sum(model.F_S3_Xyl_L_F[b1,d,h,'SSF',j] for j in model.J) == model.F_S3_Xyl_L[b1,d,h,'SSF']
    model.SSF_xylose_L_fermentation_constraint = Constraint(model.B1, model.D, model.H, rule = SSF_xylose_L_fermentation_rule)

    #EQUATION 38
    def SSF_Fermentation_Strain_compatibility_rule3(model,b1,d,h,j):
        return model.F_S3_Xyl_L_F[b1,d,h,'SSF',j] == model.pi_S3_J_L[h,'SSF',j]*model.F_S3_Xyl_L_F[b1,d,h,'SSF',j]
    model.SSF_Fermentation_Strain_compatibility3 = Constraint(model.B1,model.D,model.H,model.J,rule = SSF_Fermentation_Strain_compatibility_rule3)



    #EQUATION 39
    def SSF_Fermentation_Strain_Xyl_L_Disjuncts_rule(b,b1,d,h,j,flag):
        model = b.model()
        if flag ==1:
            b.c1 = ConstraintList()  
            b.c1.add(model.Y_S3_J_SSF1[b1,h,j]==1)
            b.c1.add(model.F_S3_P_Xyl_L[b1,d,h,'SSF',j] == model.eta_S3_P_Xyl[d,j]*model.F_S3_Xyl_L_F[b1,d,h,'SSF',j])
        if flag ==0:
            b.c0 = ConstraintList()  
            b.c0.add(model.Y_S3_J_SSF1[b1,h,j]==0)
            b.c0.add(model.F_S3_P_Xyl_L[b1,d,h,'SSF',j] == 0)

    model.SSF_Fermentation_Strain_Xyl_L_Disjuncts = Disjunct(model.B1,model.D,model.H,model.J,[0,1], rule = SSF_Fermentation_Strain_Xyl_L_Disjuncts_rule)
    model.SSF_Fermentation_Strain_Xyl_L_Disjunction = Disjunction(model.B1,model.D,model.H,model.J)
    for b1 in model.B1:
        for d in model.D:
            for h in model.H:
                for j in model.J:
                    model.SSF_Fermentation_Strain_Xyl_L_Disjunction[b1,d,h,j] = [model.SSF_Fermentation_Strain_Xyl_L_Disjuncts[b1,d,h,j,i] for i in [0,1]]


    #EQUATION 40
    def SSF_strain_selection2_rule(model,b1,h):
        return sum(model.Y_S3_J_SSF2[b1,h,j] for j in model.J) == model.Y_S3_I[b1,h,'SSF']
    model.SSF_strain_selection2_constraint = Constraint(model.B1, model.H, rule = SSF_strain_selection2_rule)

    #EQUATION 41
    def SSF_cellulose_fermentation_rule(model,b1,d,h):
        return sum(model.F_S3_Cel_F[b1,d,h,'SSF',j] for j in model.J) == model.F_S3_Cel[b1,d,h,'SSF']
    model.SSF_cellulose_fermentation_constraint = Constraint(model.B1, model.D, model.H, rule = SSF_cellulose_fermentation_rule)

    #EQUATION 42
    def SSF_hemicellulose_fermentation_rule(model,b1,d,h):
        return sum(model.F_S3_HC_F[b1,d,h,'SSF',j] for j in model.J) == model.F_S3_HC[b1,d,h,'SSF']
    model.SSF_hemicellulose_fermentation_constraint = Constraint(model.B1, model.D, model.H, rule = SSF_hemicellulose_fermentation_rule)

    #EQUATION 43
    def SSF_Fermentation_Strain_compatibility_rule1(model,b1,d,h,j):
        return model.F_S3_HC_F[b1,d,h,'SSF',j] == model.pi_S3_J_S[h,'SSF',j]*model.F_S3_HC_F[b1,d,h,'SSF',j]
    model.SSF_Fermentation_Strain_compatibility1 = Constraint(model.B1,model.D,model.H,model.J,rule = SSF_Fermentation_Strain_compatibility_rule1)


    #EQUATION 44
    def SSF_Fermentation_Strain_compatibility_rule2(model,b1,d,h,j):
        return model.F_S3_Cel_F[b1,d,h,'SSF',j] == model.pi_S3_J_S[h,'SSF',j]*model.F_S3_Cel_F[b1,d,h,'SSF',j]
    model.SSF_Fermentation_Strain_compatibility2 = Constraint(model.B1,model.D,model.H,model.J,rule = SSF_Fermentation_Strain_compatibility_rule2)

    #EQUATION 45
    def SSF_Fermentation_Strain_Cel_HC_Disjuncts_rule(b,b1,d,h,j,flag):
        model = b.model()
        if flag ==1:
            b.c1 = ConstraintList()  
            b.c1.add(model.Y_S3_J_SSF2[b1,h,j]==1)
            b.c1.add(model.F_S3_P_HC[b1,d,h,'SSF',j] == model.eta_S3_P_HC_SSF[d,j]*model.F_S3_HC_F[b1,d,h,'SSF',j])
            b.c1.add(model.F_S3_P_Cel[b1,d,h,'SSF',j] == model.eta_S3_P_Cel_SSF[d,j]*model.F_S3_Cel_F[b1,d,h,'SSF',j])
        if flag ==0:
            b.c0 = ConstraintList()  
            b.c0.add(model.Y_S3_J_SSF2[b1,h,j]==0)
            b.c0.add(model.F_S3_P_HC[b1,d,h,'SSF',j] == 0)
            b.c0.add(model.F_S3_P_Cel[b1,d,h,'SSF',j] == 0)

    model.SSF_Fermentation_Strain_Cel_HC_Disjuncts = Disjunct(model.B1,model.D,model.H,model.J,[0,1], rule = SSF_Fermentation_Strain_Cel_HC_Disjuncts_rule)
    model.SSF_Fermentation_Strain_Cel_HC_Disjunction = Disjunction(model.B1,model.D,model.H,model.J)
    for b1 in model.B1:
        for d in model.D:
            for h in model.H:
                for j in model.J:
                    model.SSF_Fermentation_Strain_Cel_HC_Disjunction[b1,d,h,j] = [model.SSF_Fermentation_Strain_Cel_HC_Disjuncts[b1,d,h,j,i] for i in [0,1]]

    ###############################
    #SHCF Fermentation constraints#
    ###############################

    #EQUATION 46
    def SHCF_Solid_Loading_Disjuncts_rule(b,b1,h,flag):
        model = b.model()
        if flag ==0:
            b.c0 = ConstraintList()  
            b.c0.add(sum(model.F_S3_Solid[b1,d,h,'SHCF'] for d in model.D) == model.SL_S3()*(sum(model.F_S3_Solid[b1,d,h,'SHCF'] for d in model.D) +sum(model.F_S3_Liquid[b1,d,h,'SHCF'] for d in model.D) +  model.F_S3_Water[b1,h,'SHCF']))
            b.c0.add(model.F_S3_Liquid_R[b1,h,'SHCF'] ==0)
        if flag ==1:
            b.c1 = ConstraintList()  
            b.c1.add(sum(model.F_S3_Solid[b1,d,h,'SHCF'] for d in model.D) == model.SL_S3()*(sum(model.F_S3_Solid[b1,d,h,'SHCF'] for d in model.D)+sum(model.F_S3_Liquid[b1,d,h,'SHCF'] for d in model.D) -  model.F_S3_Liquid_R[b1,h,'SHCF']))
            b.c1.add(model.F_S3_Water[b1,h,'SHCF'] == 0)

    model.SHCF_Solid_Loading_Disjuncts = Disjunct(model.B1,model.H,[0,1], rule = SHCF_Solid_Loading_Disjuncts_rule)
    model.SHCF_Solid_Loading_Disjunction = Disjunction(model.B1,model.H)
    for b1 in model.B1:     
        for h in model.H:
            model.SHCF_Solid_Loading_Disjunction[b1,h] = [model.SHCF_Solid_Loading_Disjuncts[b1,h,i] for i in [0,1]]



    #EQUATION 47
    def SHCF_Stream_rule(b,b1,d,h,flag):
        model = b.model()
        if flag ==1:
            b.c1 = ConstraintList()  
            b.c1.add(model.Y_S3_I[b1,h,'SHCF']==1)
            b.c1.add(model.F_S3_Solid[b1,d,h,'SHCF']== model.F_S2_Solid[b1,d])
            b.c1.add(model.F_S3_Liquid[b1,d,h,'SHCF']== model.F_S2_Liquid[b1,d])
            b.c1.add(model.F_S3_Glu[b1,d,h,'SHCF'] ==   model.SHCF_Cell_To_Glu_Conversion[d]*model.F_S2[b1,d,'Cellulose'])
            b.c1.add(model.F_S3_Xyl[b1,d,h,'SHCF'] ==  model.SHCF_HC_To_Xyl_Conversion[d]*model.F_S2[b1,d,'Hemicellulose'])
            b.c1.add(model.F_S3_Xyl_L[b1,d,h,'SHCF'] == model.F_S2[b1,d,'Xylose'])
        if flag ==0:
            b.c0 = ConstraintList()  
            b.c0.add(model.Y_S3_I[b1,h,'SHCF'] == 0)
            b.c0.add(model.F_S3_Solid[b1,d,h,'SHCF'] == 0)
            b.c0.add(model.F_S3_Liquid[b1,d,h,'SHCF'] == 0)
            b.c0.add(model.F_S3_Glu[b1,d,h,'SHCF'] == 0)
            b.c0.add(model.F_S3_Xyl[b1,d,h,'SHCF'] ==  0)
            b.c0.add(model.F_S3_Xyl_L[b1,d,h,'SHCF'] == 0)

    model.SHCF_Stream_Disjuncts = Disjunct(model.B1,model.D,model.H,[0,1], rule =  SHCF_Stream_rule )
    model.SHCF_Stream_Disjunction = Disjunction(model.B1,model.D,model.H)
    for b1 in model.B1:
        for d in model.D:
            for h in model.H:
                model.SHCF_Stream_Disjunction[b1,d,h] = [model.SHCF_Stream_Disjuncts[b1,d,h,i] for i in [0,1]]
        
    #EQUATION 48
    def SHCF_strain_selection_rule(model,b1,h):
        return sum(model.Y_S3_J_SHCF[b1,h,j] for j in model.J) == model.Y_S3_I[b1,h,'SHCF']
    model.SHCF_strain_selection_constraint = Constraint(model.B1, model.H, rule = SHCF_strain_selection_rule)

    #EQUATION 49
    def SHCF_glusose_fermentation_rule(model,b1,d,h):
        return sum(model.F_S3_Glu_F[b1,d,h,'SHCF',j] for j in model.J) == model.F_S3_Glu[b1,d,h,'SHCF']
    model.SHCF_glusose_fermentation_constraint = Constraint(model.B1, model.D, model.H, rule = SHCF_glusose_fermentation_rule)

    #EQUATION 50
    def SHCF_xylose_fermentation_rule(model,b1,d,h):
        return sum(model.F_S3_Xyl_F[b1,d,h,'SHCF',j] for j in model.J) == model.F_S3_Xyl[b1,d,h,'SHCF']
    model.SHCF_xylose_fermentation_constraint = Constraint(model.B1, model.D, model.H, rule = SHCF_xylose_fermentation_rule)

    #EQUATION 51
    def SHCF_xylose_L_fermentation_rule(model,b1,d,h):
        return sum(model.F_S3_Xyl_L_F[b1,d,h,'SHCF',j] for j in model.J) == model.F_S3_Xyl_L[b1,d,h,'SHCF']
    model.SHCF_xylose_L_fermentation_constraint = Constraint(model.B1, model.D, model.H, rule = SHCF_xylose_L_fermentation_rule)

    #EQUATION 52
    def SHCF_Fermentation_Strain_compatibility_rule1(model,b1,d,h,j):
        return model.F_S3_Xyl_F[b1,d,h,'SHCF',j] == model.pi_S3_J_S[h,'SHCF',j]*model.F_S3_Xyl_F[b1,d,h,'SHCF',j]
    model.SHCF_Fermentation_Strain_compatibility1 = Constraint(model.B1,model.D,model.H,model.J,rule = SHCF_Fermentation_Strain_compatibility_rule1)

    #EQUATION 53
    def SHCF_Fermentation_Strain_compatibility_rule2(model,b1,d,h,j):
        return model.F_S3_Glu_F[b1,d,h,'SHCF',j] == model.pi_S3_J_S[h,'SHCF',j]*model.F_S3_Glu_F[b1,d,h,'SHCF',j]
    model.SHCF_Fermentation_Strain_compatibility2 = Constraint(model.B1,model.D,model.H,model.J,rule = SHCF_Fermentation_Strain_compatibility_rule2)

    #EQUATION 54
    def SHCF_Fermentation_Strain_compatibility_rule3(model,b1,d,h,j):
        return model.F_S3_Xyl_L_F[b1,d,h,'SHCF',j] == model.pi_S3_J_L[h,'SHCF',j]*model.F_S3_Xyl_L_F[b1,d,h,'SHCF',j]
    model.SHCF_Fermentation_Strain_compatibility3 = Constraint(model.B1,model.D,model.H,model.J,rule = SHCF_Fermentation_Strain_compatibility_rule3)

    #Equation 55
    def SHCF_Fermentation_Strain_Glu_Xyl_Disjuncts_rule(b,b1,d,h,j,flag):
        model = b.model()
        if flag ==1:
            b.c1 = ConstraintList()  
            b.c1.add(model.Y_S3_J_SHCF[b1,h,j]==1)
            b.c1.add(model.F_S3_P_Xyl[b1,d,h,'SHCF',j] == model.eta_S3_P_Xyl[d,j]*(model.F_S3_Xyl_F[b1,d,h,'SHCF',j]+model.F_S3_Xyl_L_F[b1,d,h,'SHCF',j]))
            b.c1.add(model.F_S3_P_Glu[b1,d,h,'SHCF',j] == model.eta_S3_P_Glu[d,j]*model.F_S3_Glu_F[b1,d,h,'SHCF',j])
        if flag ==0:
            b.c0 = ConstraintList()  
            b.c0.add(model.Y_S3_J_SHCF[b1,h,j]==0)
            b.c0.add(model.F_S3_P_Xyl[b1,d,h,'SHCF',j] == 0)
            b.c0.add(model.F_S3_P_Glu[b1,d,h,'SHCF',j] == 0)

    model.SHCF_Fermentation_Strain_Glu_Xyl_Disjuncts = Disjunct(model.B1,model.D,model.H,model.J,[0,1], rule = SHCF_Fermentation_Strain_Glu_Xyl_Disjuncts_rule)
    model.SHCF_Fermentation_Strain_Glu_Xyl_Disjunction = Disjunction(model.B1,model.D,model.H,model.J)
    for b1 in model.B1:
        for d in model.D:
            for h in model.H:
                for j in model.J:
                    model.SHCF_Fermentation_Strain_Glu_Xyl_Disjunction[b1,d,h,j] = [model.SHCF_Fermentation_Strain_Glu_Xyl_Disjuncts[b1,d,h,j,i] for i in [0,1]]
        

    ###############################
    #SSCF Fermentation constraints#
    ###############################

    #EQUATION 56
    def SSCF_Solid_Loading_Disjuncts_rule(b,b1,h,flag):
        model = b.model()
        if flag ==0:
            b.c0 = ConstraintList()  
            b.c0.add(sum(model.F_S3_Solid[b1,d,h,'SSCF'] for d in model.D) == model.SL_S3()*(sum(model.F_S3_Solid[b1,d,h,'SSCF'] for d in model.D) + sum(model.F_S3_Liquid[b1,d,h,'SSCF'] for d in model.D) +  model.F_S3_Water[b1,h,'SSCF']))
            b.c0.add(model.F_S3_Liquid_R[b1,h,'SSCF'] ==0)
        if flag ==1:
            b.c1 = ConstraintList()  
            b.c1.add(sum(model.F_S3_Solid[b1,d,h,'SSCF'] for d in model.D) == model.SL_S3()*(sum(model.F_S3_Solid[b1,d,h,'SSCF'] for d in model.D)+ sum(model.F_S3_Liquid[b1,d,h,'SSCF'] for d in model.D) -  model.F_S3_Liquid_R[b1,h,'SSCF']))
            b.c1.add(model.F_S3_Water[b1,h,'SSCF'] == 0)

    model.SSCF_Solid_Loading_Disjuncts = Disjunct(model.B1,model.H,[0,1], rule = SSCF_Solid_Loading_Disjuncts_rule)
    model.SSCF_Solid_Loading_Disjunction = Disjunction(model.B1,model.H)
    for b1 in model.B1:     
        for h in model.H:
            model.SSCF_Solid_Loading_Disjunction[b1,h] = [model.SSCF_Solid_Loading_Disjuncts[b1,h,i] for i in [0,1]]

    #EQUATION 57
    def SSCF_Stream_rule(b,b1,d,h,flag):
        model = b.model()
        if flag ==1:
            b.c1 = ConstraintList()  
            b.c1.add(model.Y_S3_I[b1,h,'SSCF']==1)
            b.c1.add(model.F_S3_Solid[b1,d,h,'SSCF'] == model.F_S2_Solid[b1,d])
            b.c1.add(model.F_S3_Liquid[b1,d,h,'SSCF'] == model.F_S2_Liquid[b1,d])
            b.c1.add(model.F_S3_Cel[b1,d,h,'SSCF'] ==  model.F_S2[b1,d,'Cellulose'])
            b.c1.add(model.F_S3_HC[b1,d,h,'SSCF'] ==  model.F_S2[b1,d,'Hemicellulose'])
            b.c1.add(model.F_S3_Xyl_L[b1,d,h,'SSCF'] == model.F_S2[b1,d,'Xylose'])
        if flag ==0:
            b.c0 = ConstraintList()  
            b.c0.add(model.Y_S3_I[b1,h,'SSCF']==0)
            b.c0.add(model.F_S3_Solid[b1,d,h,'SSCF'] == 0)
            b.c0.add(model.F_S3_Liquid[b1,d,h,'SSCF'] == 0)
            b.c0.add(model.F_S3_Cel[b1,d,h,'SSCF'] == 0)
            b.c0.add(model.F_S3_HC[b1,d,h,'SSCF'] ==  0)
            b.c0.add(model.F_S3_Xyl_L[b1,d,h,'SSCF'] == 0)

    model.SSCF_Stream_Disjuncts = Disjunct(model.B1,model.D,model.H,[0,1], rule =  SSCF_Stream_rule )
    model.SSCF_Stream_Disjunction = Disjunction(model.B1,model.D,model.H)
    for b1 in model.B1:
        for d in model.D:
            for h in model.H:
                model.SSCF_Stream_Disjunction[b1,d,h] = [model.SSCF_Stream_Disjuncts[b1,d,h,i] for i in [0,1]]

    #EQUATION 58
    def SSCF_strain_selection_rule(model,b1,h):
        return sum(model.Y_S3_J_SSCF[b1,h,j] for j in model.J) == model.Y_S3_I[b1,h,'SSCF']
    model.SSCF_strain_selection_constraint = Constraint(model.B1, model.H, rule = SSCF_strain_selection_rule)

    #EQUATION 59
    def SSCF_cellulose_fermentation_rule(model,b1,d,h):
        return sum(model.F_S3_Cel_F[b1,d,h,'SSCF',j] for j in model.J) == model.F_S3_Cel[b1,d,h,'SSCF']
    model.SSCF_cellulose_fermentation_constraint = Constraint(model.B1, model.D, model.H, rule = SSCF_cellulose_fermentation_rule)

    #EQUATION 60
    def SSCF_hemicellulose_fermentation_rule(model,b1,d,h):
        return sum(model.F_S3_HC_F[b1,d,h,'SSCF',j] for j in model.J) == model.F_S3_HC[b1,d,h,'SSCF']
    model.SSCF_hemicellulose_fermentation_constraint = Constraint(model.B1, model.D, model.H, rule = SSCF_hemicellulose_fermentation_rule)

    #EQUATION 61
    def SSCF_xylose_L_fermentation_rule(model,b1,d,h):
        return sum(model.F_S3_Xyl_L_F[b1,d,h,'SSCF',j] for j in model.J) == model.F_S3_Xyl_L[b1,d,h,'SSCF']
    model.SSCF_xylose_L_fermentation_constraint = Constraint(model.B1, model.D, model.H, rule = SSCF_xylose_L_fermentation_rule)

    #EQUATION 62
    def SSCF_Fermentation_Strain_compatibility_rule1(model,b1,d,h,j):
        return model.F_S3_HC_F[b1,d,h,'SSCF',j] == model.pi_S3_J_S[h,'SSCF',j]*model.F_S3_HC_F[b1,d,h,'SSCF',j]
    model.SSCF_Fermentation_Strain_compatibility1 = Constraint(model.B1,model.D,model.H,model.J,rule = SSCF_Fermentation_Strain_compatibility_rule1)


    #EQUATION 63
    def SSCF_Fermentation_Strain_compatibility_rule2(model,b1,d,h,j):
        return model.F_S3_Cel_F[b1,d,h,'SSCF',j] == model.pi_S3_J_S[h,'SSCF',j]*model.F_S3_Cel_F[b1,d,h,'SSCF',j]
    model.SSCF_Fermentation_Strain_compatibility2 = Constraint(model.B1,model.D,model.H,model.J,rule = SSCF_Fermentation_Strain_compatibility_rule2)

    #EQUATION 64
    def SSCF_Fermentation_Strain_compatibility_rule3(model,b1,d,h,j):
        return model.F_S3_Xyl_L_F[b1,d,h,'SSCF',j] == model.pi_S3_J_L[h,'SSCF',j]*model.F_S3_Xyl_L_F[b1,d,h,'SSCF',j]
    model.SSCF_Fermentation_Strain_compatibility3 = Constraint(model.B1,model.D,model.H,model.J,rule = SSCF_Fermentation_Strain_compatibility_rule3)

    #EQUATION 65 
    def SSCF_Fermentation_Strain_Cel_HC_Disjuncts_rule(b,b1,d,h,j,flag):
        model = b.model()
        if flag ==1:
            b.c1 = ConstraintList()  
            b.c1.add(model.Y_S3_J_SSCF[b1,h,j]==1)
            b.c1.add(model.F_S3_P_HC[b1,d,h,'SSCF',j] == model.eta_S3_P_HC_SSCF[d,j]*model.F_S3_HC_F[b1,d,h,'SSCF',j])
            b.c1.add(model.F_S3_P_Cel[b1,d,h,'SSCF',j] == model.eta_S3_P_Cel_SSCF[d,j]*model.F_S3_Cel_F[b1,d,h,'SSCF',j])        
            b.c1.add(model.F_S3_P_Xyl_L[b1,d,h,'SSCF',j] == model.eta_S3_P_Xyl_SSCF[d,j]*model.F_S3_Xyl_L_F[b1,d,h,'SSCF',j])
        if flag ==0:
            b.c0 = ConstraintList()  
            b.c0.add(model.Y_S3_J_SSCF[b1,h,j]==0)
            b.c0.add(model.F_S3_P_HC[b1,d,h,'SSCF',j] == 0)
            b.c0.add(model.F_S3_P_Cel[b1,d,h,'SSCF',j] == 0)
            b.c0.add(model.F_S3_P_Xyl_L[b1,d,h,'SSCF',j] == 0)

    model.SSCF_Fermentation_Strain_Cel_HC_Disjuncts = Disjunct(model.B1,model.D,model.H,model.J,[0,1], rule = SSCF_Fermentation_Strain_Cel_HC_Disjuncts_rule)
    model.SSCF_Fermentation_Strain_Cel_HC_Disjunction = Disjunction(model.B1,model.D,model.H,model.J)    
    for b1 in model.B1:
        for d in model.D:
            for h in model.H:
                for j in model.J:
                    model.SSCF_Fermentation_Strain_Cel_HC_Disjunction[b1,d,h,j] = [model.SSCF_Fermentation_Strain_Cel_HC_Disjuncts[b1,d,h,j,i] for i in [0,1]]

    ##########################
    #Fermentation constraints#
    ##########################

    #EQUATION 65 - Bis (not numbered in the report, between eq 65 and 66)
    def Fermentation_Product_rule(model,h):
        return model.F_S3_P[h] == sum(sum(sum(model.F_S3_P_Xyl_L[b1,d,h,'SHF',j] + model.F_S3_P_Xyl[b1,d,h,'SHF',j] + model.F_S3_P_Glu[b1,d,h,'SHF',j]  for j in model.J) for d in model.D) for b1 in model.B1)  + sum(sum(sum(model.F_S3_P_Xyl_L[b1,d,h,'SSF',j] + model.F_S3_P_HC[b1,d,h,'SSF',j] + model.F_S3_P_Cel[b1,d,h,'SSF',j]  for j in model.J) for d in model.D) for b1 in model.B1) + sum(sum(sum(model.F_S3_P_Xyl_L[b1,d,h,'SSCF',j] + model.F_S3_P_HC[b1,d,h,'SSCF',j] + model.F_S3_P_Cel[b1,d,h,'SSCF',j]  for j in model.J) for d in model.D) for b1 in model.B1) +sum(sum(sum(model.F_S3_P_Xyl[b1,d,h,'SHCF',j] + model.F_S3_P_Glu[b1,d,h,'SHCF',j]  for j in model.J) for d in model.D) for b1 in model.B1) 
    model.Fermentation_Product = Constraint(model.H,rule =  Fermentation_Product_rule)

    ###########################################
    #Concentration and purification constraints#
    ###########################################


    #EQUATION 66
    def Purification_distribution_rule(model,h):
        return model.F_S3_P[h] == sum(model.F_S4[h,k] for k in model.K)
    model.Purification_distibution = Constraint(model.H, rule = Purification_distribution_rule)

    #EQUATION 67
    def Purification_selection_rule(model,h):
        return sum(model.Y_S4[h,k] for k in model.K) <= sum(sum(model.Y_S3_I[b1,h,i] for b1 in model.B1) for i in model.I)
    model.Purification_selection = Constraint(model.H, rule = Purification_selection_rule)


    #EQUATION 68
    def Purification_selection_rule2(model,h):
        return sum(model.Y_S4[h,k] for k in model.K) <= 1
    model.Purification_selection2 = Constraint(model.H, rule = Purification_selection_rule2)


    #EQUATION 69
    def Purification_Compatibility_rule(model,h,k):
        return model.F_S4[h,k] == model.pi_S4[h,k]*model.F_S4[h,k]
    model.Purification_Compatibility = Constraint(model.H,model.K, rule = Purification_Compatibility_rule)

    #EQUATION 70
    def Purification_Disjuncts_rule(b,h,k,flag):
        model = b.model()
        if flag ==1:
            b.c1 = ConstraintList()  
            b.c1.add(model.Y_S4[h,k]==1)
            b.c1.add(model.F_S4_P_R[h,k] == model.eta_S4[k]*model.F_S4[h,k])
        if flag ==0:
            b.c0 = ConstraintList()  
            b.c0.add(model.Y_S4[h,k]==0)
            b.c0.add(model.F_S4_P_R[h,k] ==0)
    model.Purification_Disjuncts = Disjunct(model.H,model.K,[0,1], rule =  Purification_Disjuncts_rule)
    model.Purification_Disjunction = Disjunction(model.H,model.K)
    for h in model.H:
        for k in model.K:
            model.Purification_Disjunction[h,k] = [model.Purification_Disjuncts[h,k,i] for i in [0,1]]
            
    #EQUATION 71
    def Product_Purification_rule(model,h):
        return sum(model.F_S4_P_R[h,k] for k in model.K) == model.F_S4_P_Tot[h]
    model.Product_Purification = Constraint(model.H,rule= Product_Purification_rule)

    ####################
    #Lignin conversion #
    ####################

    #EQUATION 72
    def Lignin_Conversion_Disjuncts_rule(b,flag):
        model = b.model()
        if flag ==1:
            b.c1 = ConstraintList()  
            b.c1.add(model.Y_S4_Lignin==1)
            b.c1.add(model.F_S4_Lignin_R == sum(sum(model.F_S2[b1,d,'Lignin_S']  for d in model.D) for b1 in model.B1))
        if flag ==0:
            b.c0 = ConstraintList()  
            b.c0.add(model.Y_S4_Lignin == 0)
            b.c0.add(model.F_S4_Lignin_R == 0 )
    model.Lignin_Conversion_Disjuncts = Disjunct([0,1], rule =  Lignin_Conversion_Disjuncts_rule)
    model.Lignin_Conversion_Disjunction = Disjunction(expr=[model.Lignin_Conversion_Disjuncts[i] for i in [0,1]])

    #####################
    # Manure conversion #
    #####################

    #EQUATION 73
    def Manure_Conversion_selection_rule(model):
        return model.Y_S4_Manure == sum(model.Y_S4_Manure_L[l] for l in model.L)
    model.Manure_Conversion_selection = Constraint(rule = Manure_Conversion_selection_rule)

    #EQUATION 74
    def Manure_Conversion_Disjuncts_rule(b,l,flag):
        model = b.model()
        if flag ==1:
            b.c1 = ConstraintList()  
            b.c1.add(model.Y_S4_Manure_L[l]==1 )
            b.c1.add(model.F_Manure[l]  >= 0.1)
        if flag ==0:
            b.c0 = ConstraintList()  
            b.c0.add(model.Y_S4_Manure_L[l] == 0)
            b.c0.add(model.F_Manure[l] == 0 )
    model.Manure_Conversion_Disjuncts = Disjunct(model.L,[0,1], rule =  Manure_Conversion_Disjuncts_rule)
    model.Manure_Conversion_Disjunction = Disjunction(model.L)
    for l in model.L:
        model.Manure_Conversion_Disjunction[l] = [model.Manure_Conversion_Disjuncts[l,i] for i in [0,1]]
        
        
    ######################
    #Economic Constraints#
    ######################

    def Preprocessing_RM_rule(model,x):
        return model.F_Preprocessing_RM[x] == sum(sum(model.F_S1[b1,c] for b1 in model.B1)*model.Preprocessing_Raw_Materials[c,x] for c in model.C)
    model.Preprocessing_RM = Constraint(model.Utilities, rule = Preprocessing_RM_rule)

    def Pretreatment_RM_rule(model,x):
        return model.F_Pretreatment_RM[x] == sum(sum(sum(model.F_S1ToS2[b1,c,d] for b1 in model.B1) for c in model.C)*model.Preatreatment_Raw_Materials[d,x] for d in model.D)
    model.Pretreatment_RM = Constraint(model.Utilities, rule = Pretreatment_RM_rule)

    def Fermentation_RM_rule(model,x):
        if x=='Enzyme kg':
            return  model.F_Fermentation_RM[x] == sum(sum(model.F_S2[b1,d,'Cellulose'] for b1 in model.B1)*model.Fermentation_Raw_Materials[d,x] for d in model.D)
        if x=='Water kg':
            return model.F_Fermentation_RM[x] == sum(sum(sum(model.F_S3_Water[b1,h,i] for i in model.I) for h in model.H) for b1 in model.B1)
    model.Fermentation_RM = Constraint(model.Fermentation_Utilities, rule = Fermentation_RM_rule)

    def Purification_RM_rule(model,x):
        return model.F_Purification_RM[x] == sum((sum(model.F_S4_P_R[h,k] for h in model.H))*model.Purification_Raw_Materials[k,x] for k in model.K)
    model.Purification_RM = Constraint(model.Utilities, rule = Purification_RM_rule)

    def Manure_Conversion_RM_rule(model,x):
        return model.F_Manure_Conversion_RM[x] == sum(model.F_Manure[l]*model.Manure_Conversion_Raw_Materials[l,x] for l in model.L)
    model.Manure_Conversion_RM = Constraint(model.Utilities, rule = Manure_Conversion_RM_rule)



    def Lignin_Conversion_RM_rule(model,x):
        a=0
        if x == 'Waste to landfill kg':
            a=sum(sum(model.F_S2[b1,d,'Lignin_S']  for d in model.D) for b1 in model.B1) - model.F_S4_Lignin_R 
        return model.F_Lignin_Conversion_RM[x] == model.F_S4_Lignin_R*model.Lignin_Conversion_Raw_Materials[x]  + a
    model.Lignin_Conversion_RM = Constraint(model.Utilities, rule = Lignin_Conversion_RM_rule)






    def Preprocessing_RM_Price_rule(model):
        return model.Preprocessing_RM_Price == sum(model.F_Preprocessing_RM[x]*model.Utility_Price[x] for x in model.Utilities)
    model.Preprocessing_RM_Cost = Constraint(rule = Preprocessing_RM_Price_rule)
        
    def Pretreatment_RM_Price_rule(model):
        return model.Pretreatment_RM_Price == sum(model.F_Pretreatment_RM[x]*model.Utility_Price[x] for x in model.Utilities)
    model.Pretreatment_RM_Cost = Constraint(rule = Pretreatment_RM_Price_rule)

    def Purification_RM_Price_rule(model):
        return model.Purification_RM_Price == sum(model.F_Purification_RM[x]*model.Utility_Price[x] for x in model.Utilities)
    model.Purification_RM_Cost = Constraint(rule =Purification_RM_Price_rule)

    def Fermentation_RM_Price_rule(model):
        return model.Fermentation_RM_Price == sum(model.F_Fermentation_RM[x]*FermentationRawMaterials['Price $'][x] for x in model.Fermentation_Utilities)
    model.Fermentation_RM_Cost = Constraint(rule = Fermentation_RM_Price_rule)

    def Manure_Conversion_RM_Price_rule(model):
        return model.Manure_Conversion_RM_Price == sum(model.F_Manure_Conversion_RM[x]*model.Utility_Price[x] for x in model.Utilities)
    model.Manure_Conversion_RM_Cost= Constraint(rule  = Manure_Conversion_RM_Price_rule)

    def Lignin_Conversion_RM_Price_rule(model):
        return model.Lignin_Conversion_RM_Price == sum(model.F_Lignin_Conversion_RM[x]*model.Utility_Price[x] for x in model.Utilities)
    model.Lignin_Conversion_RM_Cost = Constraint(rule = Lignin_Conversion_RM_Price_rule)

    def Feedstock_Price_rule(model):
        return model.Feedstock_Price == sum(sum(model.F_S1[b1,c] for c in model.C)*Feedstocks1[b1]['Price $'] for b1 in model.B1) + sum(sum(model.F_Manure[l] for l in model.L)*Feedstocks2[b2]['Price $'] for b2 in model.B2)
    model.Feedstock_Cost = Constraint(rule = Feedstock_Price_rule)

    def MP_Value_rule(model):
        return model.MP_Value == model.F_S4_P_Tot['Ethanol Fermentation']*MainProducts['Ethanol']['Price'] +model.F_S4_P_Tot['Succinic Acid Fermentation']*MainProducts['Succinic Acid']['Price'] + model.F_S4_P_Tot['Lactic Acid Fermentation']*MainProducts['Lactic Acid']['Price']
    model.MP_Value_C = Constraint(rule = MP_Value_rule)


    def EFACTOR_rule(model):
        return model.E_Factor*sum(model.F_S4_P_Tot[h] for h in model.H) == (sum(model.F_Lignin_Conversion_RM[x] + model.F_Manure_Conversion_RM[x] + model.F_Purification_RM[x] + model.F_Pretreatment_RM[x] + model.F_Preprocessing_RM[x] for x in model.Efactor_comp)+ sum(sum(model.F_S1[b1,c] for c in model.C) for b1 in model.B1) + sum(sum(model.F_Manure[l] for l in model.L) for b2 in model.B2)+sum(sum(sum(model.F_S2_BP[b1,d,f] for b1 in model.B1) for d in model.D) for f in model.F))
    model.EFACTOR_rule_C = Constraint(rule = EFACTOR_rule)

    #################################################
    ##Constraints used to impose production targets##
    #################################################

    def Biorefinery_Capacity_rule(model,h):
        return model.F_S4_P_Tot[h] == model.Bioreffinery_Capacity2[h]
    model.Biorefinery_Capacity= Constraint(model.H,rule = Biorefinery_Capacity_rule)

    def Biorefinery_Capacity_rule2(model,h):
        return model.F_Manure_Conversion_RM['RNG  kWh'] == -1000
    model.Biorefinery_Capacity2= Constraint(model.H,rule = Biorefinery_Capacity_rule2)


    ########################
    #Objective function#####
    ########################

    def obj_rule(model):      
        return model.MP_Value - model.Feedstock_Price - model.Preprocessing_RM_Price - model.Pretreatment_RM_Price - model.Purification_RM_Price - model.Fermentation_RM_Price -model.Manure_Conversion_RM_Price - model.Lignin_Conversion_RM_Price
    model.Objective_function = Objective(rule = obj_rule, sense = maximize) 


    """
    def obj_rule(model):      
        return model.E_Factor
    model.Objective_function = Objective(rule = obj_rule, sense = minimize) 

    """


    #Solving the model and displaying the solution
    TransformationFactory('gdp.bigm').apply_to(model)
    opt = SolverFactory("gurobi", solver_io="python") 
    opt.options['NonConvex'] = 2
    res = opt.solve(model)
    return model, res

print("Iterations Started!")
#model, res = create_n_solve()
#print(res.solver.wallclock_time)
#res.write()


def print_vars(model):
    for v in model.component_objects(Var, active=True):
        try:
            varobject = getattr(model, str(v))
            print("Variable", v)
            for index in varobject:
                if varobject[index].value is not None:
                    if varobject[index].value != 0:
                        print("   ", index, varobject[index].value) 
        except:
            pass
#print_vars(model)
""" time_list = []
for i in range(10):
    model, res = create_n_solve()
    print("Iteration", i + 1, "|", res.solver.wallclock_time)
    time_list.append(res.solver.wallclock_time)
    ds = pd.Series(time_list, name = "JulienRunTime")
    ds.to_pickle("runtimes.pickle") """

ds = pd.read_pickle("runtimes.pickle")
df.to_excel("runtimes-2.xlsx")