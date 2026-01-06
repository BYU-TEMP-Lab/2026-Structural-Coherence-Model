import numpy as np
from numpy import array
from matplotlib import pyplot as plt
import math
from scipy.optimize import curve_fit
import xlrd
import xlwt
import pandas as pd
import tkinter as tk
from tkinter import IntVar
import csv

def Kinetic_Theory_Emp(df,MSTDB_df,compound_input,mol_fracs,Temp_Range, V_m=0, C_p=0, alpha=0, expon=0):
    
    #Define smaller dataframes of each variable contained in the spreadsheet
    Compound_df=df['Compound']
    M_i_df=df['M_i (g/mol)']
    T_i_m_df=df['T_i_m (K)']
    V_i_m_df=df['V_i_m (m^3/mol)']
    alpha_i_m_df=df['alpha_i_m (K^-1)']
    C_i_0_df=df['C_i_0 (m/s)']
    C_i_p_df=df['C_i_p (J/mol/K)']
    r_c_df=df['r_c (m)']
    r_a_df=df['r_a (m)']
    rho_i_m_df=df['rho_m (g/m^3)']

    #initialize arrays of indices of compounds in the dataframe, and boolean values of whether they are in the dataframe
    indices=np.zeros(len(compound_input))
    indices_exist=np.full(len(compound_input),False)


    #find indices in dataframe and store index values
    for i in range(len(compound_input)):
        for j in range(len(Compound_df)):
            if compound_input[i]==Compound_df[j]:
                indices[i]=j
                indices_exist[i]=True

    #check to make sure all the input compounds actually exist. Quit program if not
    for i in range(len(compound_input)):
        problem=False
        if indices_exist[i]==False:
            print('Warning: ' + str(compound_input[i]) + ' is not included in the compound data spreadsheet' )
            problem=True
        if problem==True:
            print('Force quit.')
            quit()

    #Generate arrays of all parameters that already exist in the dataframe, but in the order of the input compound array
    Compound=[]
    for i in range(len(compound_input)):
        Compound.append('')
    M_i=np.zeros(len(compound_input))
    T_i_m=np.zeros(len(compound_input))
    V_i_m=np.zeros(len(compound_input))
    alpha_i_m=np.zeros(len(compound_input))
    C_i_0=np.zeros(len(compound_input))
    C_i_p=np.zeros(len(compound_input))
    r_c=np.zeros(len(compound_input))
    r_a=np.zeros(len(compound_input))
    rho_i_m=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        Compound[i]=Compound_df[indices[i]]
        M_i[i]=M_i_df[indices[i]]*0.001         # Convert to kg/mol
        T_i_m[i]=T_i_m_df[indices[i]]
        V_i_m[i]=V_i_m_df[indices[i]]
        alpha_i_m[i]=alpha_i_m_df[indices[i]]
        C_i_0[i]=C_i_0_df[indices[i]]
        C_i_p[i]=C_i_p_df[indices[i]]
        r_c[i]=r_c_df[indices[i]]
        r_a[i]=r_a_df[indices[i]]
        rho_i_m[i]=rho_i_m_df[indices[i]]*0.001 # Convert to kg/m^3

    #Calculate Gruneisen parameter
    gamma_i_m=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        gamma_i_m[i]=((M_i[i]*alpha_i_m[i]*C_i_0[i]**2)/C_i_p[i])

    #initialize a temperature range
    T_melt = Temp_Range[0]
    T=np.linspace(T_melt,Temp_Range[1],num=100)

    
    #calculate constant volume heat capacities at melting temp
    C_i_v=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        C_i_v[i]=C_i_p[i]/(1+alpha_i_m[i]*gamma_i_m[i]*T_i_m[i])


    #Calculate number of ions per compound
    n_i=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        for c in Compound[i]:
            if c.isupper():
                n_i[i] = n_i[i] + 1
            elif c.isnumeric():
                c = int(c)
                n_i[i] = n_i[i] + (c-1)

    print("KT_emp, number of ions: ",n_i)
            
 
    #Calculate compound thermal conductivites at melting points
    lambda_i_m=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        lambda_i_m[i] = 4.33 * ( (C_i_0[i] * C_i_v[i]) / (3*n_i[i]*V_i_m[i]) ) * (r_a[i]+r_c[i])


    #Calculate volume fractions
    phi_i_m=np.zeros(len(compound_input))
    denom=0
    for i in range(len(compound_input)):
        phi_i_m[i]=V_i_m[i]*mol_fracs[i]
        denom=denom+V_i_m[i]*mol_fracs[i]
    phi_i_m=phi_i_m/denom  

    #Calculate weight fractions
    kappa_i_m=np.zeros(len(compound_input))
    denom=0
    for i in range(len(compound_input)):
        kappa_i_m[i]=rho_i_m[i]*mol_fracs[i]
        denom=denom+rho_i_m[i]*mol_fracs[i] #denom=denom+rho_i_m[i]*mol_fracs[i]
    kappa_i_m=kappa_i_m/denom 

    #Calculate taus for compounds
    Tau_i_m=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        Tau_i_m[i] = ( C_i_v[i]*C_i_0[i] ) / n_i[i]


    #Calculate mixture speed of sound       <<<< Should try this as a function of temperature
    C_0_mix_m=0
    for i in range(len(compound_input)):
        C_0_mix_m=C_0_mix_m+phi_i_m[i]**2/(kappa_i_m[i]*C_i_0[i]**2)
    C_0_mix_m=1/np.sqrt(C_0_mix_m)

    #Calculate mixture specific heat
    #Calculate average number of ions
    C_p_mix=0
    n_mix=0
    if C_p==0:
        for i in range(len(compound_input)):
            C_p_mix=C_p_mix+mol_fracs[i]*C_i_p[i]
            n_mix=n_mix+mol_fracs[i]*n_i[i]
    else: 
        C_p_mix=C_p
        for i in range(len(compound_input)):
            n_mix=n_mix+mol_fracs[i]*n_i[i]
    

    #Calculate mixture thermal expansion
    alpha_mix_m=0
    if alpha==0:
        for i in range(len(compound_input)):
            alpha_mix_m=alpha_mix_m+phi_i_m[i]*alpha_i_m[i]
    else:alpha_mix_m=alpha
    
    

    #Calculate mixture molecular weight
    M_mix=0
    for i in range(len(compound_input)):
        M_mix=M_mix+mol_fracs[i]*M_i[i]
        

    #Calculate mixture gamma
    gamma_mix_m=M_mix*alpha_mix_m*C_0_mix_m**2*(1/C_p_mix)
    

    #Calculate Tau for the mixture
    Tau_mix_m=(C_p_mix/(1+alpha_mix_m*gamma_mix_m*T_melt)) *C_0_mix_m *(1/n_mix)
    
    
    #Calculate ideal thermal conductivity
    sum_term=0
    for i in range(len(compound_input)):
        sum_term=sum_term+mol_fracs[i]*(lambda_i_m[i]/Tau_i_m[i])**(-3/2)
    lambda_ideal_mix_m=( Tau_mix_m**(-3/2) *sum_term ) **(-2/3)

    #Calculate average molar volume
    V_mix_m=0
    if V_m==0:
        for i in range(len(compound_input)):
        
            V_mix_m=V_mix_m+mol_fracs[i]*M_i[i]*(1/rho_i_m[i])
    else: V_mix_m=V_m



    sum_term=0
    for i in range(len(compound_input)):
        sum_term=sum_term+mol_fracs[i]*((M_i[i]/M_mix)-1)**2
    delta_mix_m=0.4872*(lambda_ideal_mix_m/(1.380649e-23*C_0_mix_m*(n_mix*6.0221408e23*(1/V_mix_m))**(2/3)))*sum_term
    

    lambda_mix=lambda_ideal_mix_m*(1-delta_mix_m)

    lambda_mix_T=np.zeros(len(T))
    for j in range(len(T)):
        lambda_mix_T[j]=lambda_mix*(1-alpha_mix_m*(gamma_mix_m+(1/3))*(T[j]-T_melt))    #<<--T_melt must be the melting point of the liquid!

    return(T,lambda_mix_T)

def Kinetic_Theory(df,MSTDB_df,compound_input,mol_fracs,Temp_Range, V_m=0, C_p=0, alpha=0, expon=0):
    
    #Define smaller dataframes of each variable contained in the spreadsheet
    Compound_df=df['Compound']
    M_i_df=df['M_i (g/mol)']
    T_i_m_df=df['T_i_m (K)']
    V_i_m_df=df['V_i_m (m^3/mol)']
    alpha_i_m_df=df['alpha_i_m (K^-1)']
    C_i_0_df=df['C_i_0 (m/s)']
    C_i_p_df=df['C_i_p (J/mol/K)']
    r_c_df=df['r_c (m)']
    r_a_df=df['r_a (m)']
    rho_i_m_df=df['rho_m (g/m^3)']

    #initialize arrays of indices of compounds in the dataframe, and boolean values of whether they are in the dataframe
    indices=np.zeros(len(compound_input))
    indices_exist=np.full(len(compound_input),False)


    #find indices in dataframe and store index values
    for i in range(len(compound_input)):
        for j in range(len(Compound_df)):
            if compound_input[i]==Compound_df[j]:
                indices[i]=j
                indices_exist[i]=True

    #check to make sure all the input compounds actually exist. Quit program if not
    for i in range(len(compound_input)):
        problem=False
        if indices_exist[i]==False:
            print('Warning: ' + str(compound_input[i]) + ' is not included in the compound data spreadsheet' )
            problem=True
        if problem==True:
            print('Force quit.')
            quit()

    #Generate arrays of all parameters that already exist in the dataframe, but in the order of the input compound array
    Compound=[]
    for i in range(len(compound_input)):
        Compound.append('')
    M_i=np.zeros(len(compound_input))
    T_i_m=np.zeros(len(compound_input))
    V_i_m=np.zeros(len(compound_input))
    alpha_i_m=np.zeros(len(compound_input))
    C_i_0=np.zeros(len(compound_input))
    C_i_p=np.zeros(len(compound_input))
    r_c=np.zeros(len(compound_input))
    r_a=np.zeros(len(compound_input))
    rho_i_m=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        Compound[i]=Compound_df[indices[i]]
        M_i[i]=M_i_df[indices[i]]*0.001         # Convert to kg/mol
        T_i_m[i]=T_i_m_df[indices[i]]
        V_i_m[i]=V_i_m_df[indices[i]]
        alpha_i_m[i]=alpha_i_m_df[indices[i]]
        C_i_0[i]=C_i_0_df[indices[i]]
        C_i_p[i]=C_i_p_df[indices[i]]
        r_c[i]=r_c_df[indices[i]]
        r_a[i]=r_a_df[indices[i]]
        rho_i_m[i]=rho_i_m_df[indices[i]]*0.001 # Convert to kg/m^3

    #Calculate Gruneisen parameter
    gamma_i_m=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        gamma_i_m[i]=((M_i[i]*alpha_i_m[i]*C_i_0[i]**2)/C_i_p[i])

    #initialize a temperature range
    T_melt = Temp_Range[0]
    T=np.linspace(T_melt,Temp_Range[1],num=100)


    #calculate constant volume heat capacities
    C_i_v=np.zeros((len(compound_input),len(T)))
    for i in range(len(compound_input)):
        for j in range(len(T)):
            C_i_v[i][j]=C_i_p[i]/(1+alpha_i_m[i]*gamma_i_m[i]*T[j])

    
    #Calculate number of ions per compound
    n_i=np.zeros(len(compound_input))
    n_i_c=np.zeros(len(compound_input))
    n_i_a=np.zeros(len(compound_input))
    comps = 0

    #Calculate number of ions per compound
    n_i=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        for c in Compound[i]:
            if c.isupper():
                n_i[i] = n_i[i] + 1
            elif c.isnumeric():
                c = int(c)
                n_i[i] = n_i[i] + (c-1)
            if comps == 0:
                n_i_c[i] = n_i[i]
                if "NO3" in Compound[i]:
                    n_i_c[i] = n_i_c[i] +1
                    n_i_a[i] = 3
                comps = 1
            else:
                if "NO3" in Compound[i]:
                    continue
                n_i_a[i] = n_i[i] - n_i_c[i]
        comps = 0

    print("KT, number of ions: ",n_i)
    print("KT, number of cations: ",n_i_c)
    print("KT, number of anions: ",n_i_a)

    #Calculate melting point thermal conductivities based on compound minimum thermal conductivites
    lambda_i_m=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        lambda_i_m[i]= (1 + n_i_c[i]/n_i_a[i]) * 1.380649e-23 * (6.0221408e23 * n_i[i]/M_i[i] * rho_i_m[i])**(2/3) * C_i_0[i]

    #Calculate volume fractions
    phi_i_m=np.zeros(len(compound_input))
    denom=0
    for i in range(len(compound_input)):
        phi_i_m[i]=V_i_m[i]*mol_fracs[i]
        denom=denom+V_i_m[i]*mol_fracs[i]
    phi_i_m=phi_i_m/denom  

    #Calculate weight fractions
    kappa_i_m=np.zeros(len(compound_input))
    denom=0
    for i in range(len(compound_input)):
        kappa_i_m[i]=rho_i_m[i]*mol_fracs[i]
        denom=denom+rho_i_m[i]*mol_fracs[i]
    kappa_i_m=kappa_i_m/denom 

    #Calculate taus for compounds
    Tau_i_m=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        Tau_i_m[i] = ( C_i_v[i]*C_i_0[i] ) / n_i[i]


    #Calculate mixture speed of sound
    C_0_mix_m=0
    for i in range(len(compound_input)):
        C_0_mix_m=C_0_mix_m+phi_i_m[i]**2/(kappa_i_m[i]*C_i_0[i]**2)
    C_0_mix_m=1/np.sqrt(C_0_mix_m)

    #Calculate mixture specific heat
    #Calculate average number of ions
    C_p_mix=0
    n_mix=0
    if C_p==0:
        for i in range(len(compound_input)):
            C_p_mix=C_p_mix+mol_fracs[i]*C_i_p[i]
            n_mix=n_mix+mol_fracs[i]*n_i[i]
    else: 
        C_p_mix=C_p
        for i in range(len(compound_input)):
            n_mix=n_mix+mol_fracs[i]*n_i[i]
    

    #Calculate mixture thermal expansion
    alpha_mix_m=0
    if alpha==0:
        for i in range(len(compound_input)):
            alpha_mix_m=alpha_mix_m+phi_i_m[i]*alpha_i_m[i]
    else:alpha_mix_m=alpha
    
    

    #Calculate mixture molecular weight
    M_mix=0
    for i in range(len(compound_input)):
        M_mix=M_mix+mol_fracs[i]*M_i[i]
        

    #Calculate mixture gamma
    gamma_mix_m=M_mix*alpha_mix_m*C_0_mix_m**2*(1/C_p_mix)
    

    #Calculate Tau for the mixture
    Tau_mix_m=(C_p_mix/(1+alpha_mix_m*gamma_mix_m*T_melt)) *C_0_mix_m *(1/n_mix)
    
    
    #Calculate ideal thermal conductivity
    sum_term=0
    for i in range(len(compound_input)):
        sum_term=sum_term+mol_fracs[i]*(lambda_i_m[i]/Tau_i_m[i])**(-3/2)
    lambda_ideal_mix_m=( Tau_mix_m**(-3/2) *sum_term ) **(-2/3)

    #Calculate average molar volume
    V_mix_m=0
    if V_m==0:
        for i in range(len(compound_input)):
        
            V_mix_m=V_mix_m+mol_fracs[i]*M_i[i]*(1/rho_i_m[i])
        V_mix_m=V_mix_m
    else: V_mix_m=V_m

    sum_term=0
    for i in range(len(compound_input)):
        sum_term=sum_term+mol_fracs[i]*((M_i[i]/M_mix)-1)**2
    delta_mix_m=0.4872*(lambda_ideal_mix_m/(1.380649e-23*C_0_mix_m*(n_mix*6.0221408e23*(1/V_mix_m))**(2/3)))*sum_term
    

    lambda_mix=lambda_ideal_mix_m*(1-delta_mix_m)

    lambda_mix_T=np.zeros(len(T))
    for j in range(len(T)):
        lambda_mix_T[j]=lambda_mix*(1-alpha_mix_m*(gamma_mix_m+(1/3))*(T[j]-T_melt))    #<<--T_melt must be the melting point of the liquid!

    return(T,lambda_mix_T)


def Kinetic_Theory_RhoT(df,MSTDB_df,compound_input,mol_fracs,Temp_Range, V_m=0, C_p=0, alpha=0, expon=0):
    
    #Define smaller dataframes of each variable contained in the spreadsheet
    Compound_df=df['Compound']
    M_i_df=df['M_i (g/mol)']
    T_i_m_df=df['T_i_m (K)']
    V_i_m_df=df['V_i_m (m^3/mol)']
    alpha_i_m_df=df['alpha_i_m (K^-1)']
    C_i_0_df=df['C_i_0 (m/s)']
    C_i_p_df=df['C_i_p (J/mol/K)']
    r_c_df=df['r_c (m)']
    r_a_df=df['r_a (m)']
    rho_i_m_df=df['rho_m (g/m^3)']

    #initialize arrays of indices of compounds in the dataframe, and boolean values of whether they are in the dataframe
    indices=np.zeros(len(compound_input))
    indices_exist=np.full(len(compound_input),False)


    #find indices in dataframe and store index values
    for i in range(len(compound_input)):
        for j in range(len(Compound_df)):
            if compound_input[i]==Compound_df[j]:
                indices[i]=j
                indices_exist[i]=True

    #check to make sure all the input compounds actually exist. Quit program if not
    for i in range(len(compound_input)):
        problem=False
        if indices_exist[i]==False:
            print('Warning: ' + str(compound_input[i]) + ' is not included in the compound data spreadsheet' )
            problem=True
        if problem==True:
            print('Force quit.')
            quit()

    #Generate arrays of all parameters that already exist in the dataframe, but in the order of the input compound array
    Compound=[]
    for i in range(len(compound_input)):
        Compound.append('')
    M_i=np.zeros(len(compound_input))
    T_i_m=np.zeros(len(compound_input))
    V_i_m=np.zeros(len(compound_input))
    alpha_i_m=np.zeros(len(compound_input))
    C_i_0=np.zeros(len(compound_input))
    C_i_p=np.zeros(len(compound_input))
    r_c=np.zeros(len(compound_input))
    r_a=np.zeros(len(compound_input))
    rho_i_m=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        Compound[i]=Compound_df[indices[i]]
        M_i[i]=M_i_df[indices[i]]*0.001         # Convert to kg/mol
        T_i_m[i]=T_i_m_df[indices[i]]
        V_i_m[i]=V_i_m_df[indices[i]]
        alpha_i_m[i]=alpha_i_m_df[indices[i]]
        C_i_0[i]=C_i_0_df[indices[i]]
        C_i_p[i]=C_i_p_df[indices[i]]
        r_c[i]=r_c_df[indices[i]]
        r_a[i]=r_a_df[indices[i]]
        rho_i_m[i]=rho_i_m_df[indices[i]]*0.001 # Convert to kg/m^3

    #Calculate number of ions per compound
    n_i=np.zeros(len(compound_input))
    n_i_c=np.zeros(len(compound_input))
    n_i_a=np.zeros(len(compound_input))
    comps = 0
    for i in range(len(compound_input)):
        n_i[i]=(sum(1 for c in compound_input[i] if c.isupper())+ sum(1 for c in compound_input[i] if c.isnumeric()))
        for c in compound_input[i]:
            if c.isupper() and comps == 0:
                n_i_c[i] = n_i_c[i] + 1
                comps = 1
            if c.isnumeric() and comps == 1:
                n_i_c[i] = n_i_c[i]*float(c)
            if c.isupper() and comps == 1:
                n_i_a[i] = n_i_a[i] + 1
                comps = 2
            if c.isnumeric() and comps == 2:
                n_i_a[i] = n_i_a[i]*float(c)
        comps = 0

    #Calculate compound thermal conductivites at melting points
    lambda_i_m=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        lambda_i_m[i]= (1 + n_i_c[i]/n_i_a[i]) * 1.380649e-23 * (6.0221408e23 * n_i[i]/M_i[i] * rho_i_m[i])**(2/3) * C_i_0[i]

    #Calculate Gruneisen parameter
    gamma_i_m=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        gamma_i_m[i]=((M_i[i]*alpha_i_m[i]*C_i_0[i]**2)/C_i_p[i])
        # Molecular weight, heat capacity, thermal expansion, speed of sound (all at melting temp)

    #initialize a temperature range
    T_melt = Temp_Range[0]
    T=np.linspace(T_melt,Temp_Range[1],num=100)


    # Find the temperature dependent density from data
    rho_liquid = np.zeros((len(Compound),len(T)))
    for i in range(len(Compound)):
        density_A, density_B = density(df,compound_input[i])
        for j in range(len(T)):
            rho_liquid[i][j] = density_A - density_B*T[j] 
    
    
    #Calculate number of ions per compound
    n_i=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        n_i[i]=(sum(1 for c in Compound[i] if c.isupper())+ sum(1 for c in Compound[i] if c.isnumeric()))

 
    #Calculate compound thermal conductivites at melting points
    lambda_i_m=np.zeros(len(T))
    for i in range(len(compound_input)):
        for j in range(len(T)):
            l = (1 + n_i_c[i]/n_i_a[i]) * 1.380649e-23 * (6.0221408e23 * n_i[i]/M_i[i] * rho_i_m[i])**(2/3) * C_i_0[i] * (rho_liquid[i][j]/rho_i_m[i])**(gamma_i_m[i] + 1/3)
            lambda_i_m[j] = l

    return(T,lambda_i_m)

def Kinetic_Theory_Unary(df,MSTDB_df,compound_input,mol_fracs,Temp_Range, V_m=0, C_p=0, alpha=0, expon=0):
    
    #Define smaller dataframes of each variable contained in the spreadsheet
    Compound_df=df['Compound']
    M_i_df=df['M_i (g/mol)']
    T_i_m_df=df['T_i_m (K)']
    V_i_m_df=df['V_i_m (m^3/mol)']
    alpha_i_m_df=df['alpha_i_m (K^-1)']
    C_i_0_df=df['C_i_0 (m/s)']
    C_i_p_df=df['C_i_p (J/mol/K)']
    r_c_df=df['r_c (m)']
    r_a_df=df['r_a (m)']
    rho_i_m_df=df['rho_m (g/m^3)']

    #initialize arrays of indices of compounds in the dataframe, and boolean values of whether they are in the dataframe
    indices=np.zeros(len(compound_input))
    indices_exist=np.full(len(compound_input),False)


    #find indices in dataframe and store index values
    for i in range(len(compound_input)):
        for j in range(len(Compound_df)):
            if compound_input[i]==Compound_df[j]:
                indices[i]=j
                indices_exist[i]=True

    #check to make sure all the input compounds actually exist. Quit program if not
    for i in range(len(compound_input)):
        problem=False
        if indices_exist[i]==False:
            print('Warning: ' + str(compound_input[i]) + ' is not included in the compound data spreadsheet' )
            problem=True
        if problem==True:
            print('Force quit.')
            quit()

    #Generate arrays of all parameters that already exist in the dataframe, but in the order of the input compound array
    Compound=[]
    for i in range(len(compound_input)):
        Compound.append('')
    M_i=np.zeros(len(compound_input))
    T_i_m=np.zeros(len(compound_input))
    V_i_m=np.zeros(len(compound_input))
    alpha_i_m=np.zeros(len(compound_input))
    C_i_0=np.zeros(len(compound_input))
    C_i_p=np.zeros(len(compound_input))
    r_c=np.zeros(len(compound_input))
    r_a=np.zeros(len(compound_input))
    rho_i_m=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        Compound[i]=Compound_df[indices[i]]
        M_i[i]=M_i_df[indices[i]]
        T_i_m[i]=T_i_m_df[indices[i]]
        V_i_m[i]=V_i_m_df[indices[i]]
        alpha_i_m[i]=alpha_i_m_df[indices[i]]
        C_i_0[i]=C_i_0_df[indices[i]]
        C_i_p[i]=C_i_p_df[indices[i]]
        r_c[i]=r_c_df[indices[i]]
        r_a[i]=r_a_df[indices[i]]
        rho_i_m[i]=rho_i_m_df[indices[i]]

    #Calculate Gruneisen parameter
    gamma_i_m=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        gamma_i_m[i]=((M_i[i]*alpha_i_m[i]*C_i_0[i]**2)/C_i_p[i])*(1/1000)
        # Molecular weight, heat capacity, thermal expansion, speed of sound (all at melting temp)

    #initialize a temperature range
    T_melt = Temp_Range[0]
    T=np.linspace(T_melt,Temp_Range[1],num=100)

    
    #Calculate number of ions per compound
    n_i=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        n_i[i]=(sum(1 for c in Compound[i] if c.isupper())+ sum(1 for c in Compound[i] if c.isnumeric()))

 
    #Calculate compound thermal conductivites at melting points
    lambda_i_m=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        lambda_i_m[i]=4.33*1.5*   ( (  C_i_0[i]*   (C_i_p[i]/(1+alpha_i_m[i]*gamma_i_m[i]*T_i_m[i]))    )   /   (3*n_i[i]*V_i_m[i])    )    *     (r_a[i]+r_c[i])


    lambda_mix_T=np.zeros(len(T))
    for j in range(len(T)):
        lambda_mix_T[j]=lambda_i_m[i]*(1-alpha_i_m[i]*(gamma_i_m[i]+(1/3))*(T[j]-T_melt))

    return(T,lambda_mix_T)

def Phonon_Gas_Model(df,MSTDB_df,compound_input,mol_fracs,Temp_Range, V_m=0, C_p=0, alpha=0, expon=0):
    
    #Define smaller dataframes of each variable contained in the spreadsheet
    Compound_df=df['Compound']
    M_i_df=df['M_i (g/mol)']
    T_i_m_df=df['T_i_m (K)']
    V_i_m_df=df['V_i_m (m^3/mol)']
    alpha_i_m_df=df['alpha_i_m (K^-1)']
    C_i_0_df=df['C_i_0 (m/s)']
    C_i_p_df=df['C_i_p (J/mol/K)']
    C_i_p_sp_df=df['C_i_p_sp (J/g/K)']
    r_c_df=df['r_c (m)']
    r_a_df=df['r_a (m)']
    rho_i_m_df=df['rho_m (g/m^3)']
    C_0_T_A_df=df['SoS(T)_A (A+B*T)']
    C_0_T_B_df=df['SoS(T)_B']
    rho_T_A_df=df['A (Density (g/cm3):   A - BT(K))']
    rho_T_B_df=df['B (Density (g/cm3):   A - BT(K))'] 

    #initialize arrays of indices of compounds in the dataframe, and boolean values of whether they are in the dataframe
    indices=np.zeros(len(compound_input))
    indices_exist=np.full(len(compound_input),False)


    #find indices in dataframe and store index values
    for i in range(len(compound_input)):
        for j in range(len(Compound_df)):
            if compound_input[i]==Compound_df[j]:
                indices[i]=j
                indices_exist[i]=True

    #check to make sure all the input compounds actually exist. Quit program if not
    for i in range(len(compound_input)):
        problem=False
        if indices_exist[i]==False:
            print('Warning: ' + str(compound_input[i]) + ' is not included in the compound data spreadsheet' )
            problem=True
        if problem==True:
            print('Force quit.')
            quit()

    #Generate arrays of all parameters that already exist in the dataframe, but in the order of the input compound array
    Compound=[]
    for i in range(len(compound_input)):
        Compound.append('')
    M_i=np.zeros(len(compound_input))
    T_i_m=np.zeros(len(compound_input))
    V_i_m=np.zeros(len(compound_input))
    alpha_i_m=np.zeros(len(compound_input))
    C_i_0=np.zeros(len(compound_input))
    C_i_p=np.zeros(len(compound_input))
    C_i_p_sp=np.zeros(len(compound_input))
    r_c=np.zeros(len(compound_input))
    r_a=np.zeros(len(compound_input))
    rho_i_m=np.zeros(len(compound_input))
    C_0_T_A=np.zeros(len(compound_input))
    C_0_T_B=np.zeros(len(compound_input))
    rho_T_A=np.zeros(len(compound_input))
    rho_T_B=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        Compound[i]=Compound_df[indices[i]]
        M_i[i]=M_i_df[indices[i]]*0.001         # Convert to kg/mol
        T_i_m[i]=T_i_m_df[indices[i]]
        V_i_m[i]=V_i_m_df[indices[i]]
        alpha_i_m[i]=alpha_i_m_df[indices[i]]
        C_i_0[i]=C_i_0_df[indices[i]]
        C_i_p[i]=C_i_p_df[indices[i]]
        C_i_p_sp[i]=C_i_p_sp_df[indices[i]]/0.001   # Convert to J/kg/K
        r_c[i]=r_c_df[indices[i]]
        r_a[i]=r_a_df[indices[i]]
        rho_i_m[i]=rho_i_m_df[indices[i]]*0.001 # Convert to kg/m^3
        C_0_T_A[i]=C_0_T_A_df[indices[i]]
        C_0_T_B[i]=C_0_T_B_df[indices[i]]
        rho_T_A[i]=rho_T_A_df[indices[i]]
        rho_T_B[i]=rho_T_B_df[indices[i]]


    #Calculate Gruneisen parameter
    gamma_i_m=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        gamma_i_m[i]=((M_i[i]*alpha_i_m[i]*C_i_0[i]**2)/C_i_p[i])#*(1/1000)
        # Molecular weight, heat capacity, thermal expansion, speed of sound (all at melting temp)


    #initialize a temperature range
    T_melt = Temp_Range[0]
    T=np.linspace(T_melt,Temp_Range[1],num=100)


    #Calculate volume fractions
    phi_i_m=np.zeros(len(compound_input))
    denom=0
    for i in range(len(compound_input)):
        phi_i_m[i]=V_i_m[i]*mol_fracs[i]
        denom=denom+V_i_m[i]*mol_fracs[i]
    phi_i_m=phi_i_m/denom  


    # Find the temperature dependent density from data
    rho_liquid = np.zeros((len(Compound),len(T)))
    for i in range(len(Compound)):     
        density_A, density_B = density(df,compound_input[i])
        for j in range(len(T)):
            rho_liquid[i][j] = density_A + density_B*T[j] 


    #Calculate weight fractions
    kappa_i_m=np.zeros((len(compound_input),len(T)))
    denom=0
    for i in range(len(compound_input)):
        for j in range(len(T)):
            kappa_i_m[i][j]=rho_liquid[i][j]*mol_fracs[i]
            denom=denom+rho_liquid[i][j]*mol_fracs[i]
    kappa_i_m=kappa_i_m/denom 


    #Calculate number of ions per compound
    n_i=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        for c in Compound[i]:
            if c.isupper():
                n_i[i] = n_i[i] + 1
            elif c.isnumeric():
                c = int(c)
                n_i[i] = n_i[i] + (c-1)
    
    print("PGM, number of ions: ",n_i)


    #Calculate average number of ions
    C_p_mix=0
    C_p_sp_mix=0
    n_mix=0
    if C_p==0:
        for i in range(len(compound_input)):
            C_p_mix=C_p_mix+mol_fracs[i]*C_i_p[i]
            C_p_sp_mix=C_p_sp_mix+mol_fracs[i]*C_i_p_sp[i]
            n_mix=n_mix+mol_fracs[i]*n_i[i]
    else: 
        C_p_mix=C_p
        C_p_sp_mix=C_p
        for i in range(len(compound_input)):
            n_mix=n_mix+mol_fracs[i]*n_i[i]
 

    #Calculate mixture molecular weight
    M_mix=0
    for i in range(len(compound_input)):
        M_mix=M_mix+mol_fracs[i]*M_i[i]


    #Calculate average molar volume
    V_mix_m=np.zeros(len(T))
    if V_m==0:
        for i in range(len(compound_input)):
            for j in range(len(T)):
                V_mix_m[j]=V_mix_m[j]+mol_fracs[i]*M_i[i]*(1/rho_liquid[i][j])
    else: V_mix_m=V_m


    #Calculate mixture speed of sound
    C_0_mix_T = np.zeros(len(T))
    i = 0

    # If measurement data for speed of sound with temp does not exist
    for j in range(len(T)):
        for i in range(len(compound_input)):
            if C_0_T_A[i] == 0:
                C_0_mix_T[j]=1/np.sqrt(C_0_mix_T[j]+phi_i_m[i]**2/(kappa_i_m[i][j]*C_i_0[i]**2))
            else:
                SoS = C_0_T_B[i]*T[j]+C_0_T_A[i]
                C_0_mix_T[j]=1/np.sqrt(C_0_mix_T[j]+phi_i_m[i]**2/(kappa_i_m[i][j]*SoS**2))


    #Calculate compound thermal conductivity
    lambda_i_m=np.zeros(len(T))
    lambda_i_mb=np.zeros(len(T))
    lambda_i_mg=np.zeros(len(T))
    for i in range(len(compound_input)):
        for j in range(len(T)):
            if all(i == 0 for i in C_0_T_A):
                SoS = C_0_T_B[i]*T[j]+C_0_T_A[i]    # Calculate sound velocity with temp
            else:
                SoS = C_i_0
            #lambda_i_m[j] = 1/3 * C_p_mix * SoS*1.4 * rho_liquid[i][j] * 1/M_mix * (1/(V_mix_m/(n_mix*6.022*10**23)))**(-1/3)
            #lambda_i_m[j] = 1/3 * C_p_mix * SoS * rho_liquid[i][j] * 1/M_mix * (1/(V_mix_m/(n_mix*6.022*10**23)))**(-1/3)
            #lambda_i_m[j] = 1/3 * C_p_mix * 1/V_mix_m * SoS * (6.022*10**23/V_mix_m)**(-1/3)
            #lambda_i_m[j] = 1/3 * C_p_mix * 1/V_mix_m[j] * SoS * (r_a[i]+r_c[i])
            lambda_i_mb[j] = 1/3 * C_p_mix * 1/V_mix_m[j] * C_i_0[i] * (r_a[i]+r_c[i])   # Not consistent with Zhao
            lambda_i_m[j] = 1/3 * C_p_sp_mix * rho_liquid[i][j] * C_i_0[i] * (r_a[i]+r_c[i])    # Verified with Zhao's results, uses Zhao's data and specific heat capacity
            lambda_i_mg[j] = 1/3 * C_p_mix * 1/M_mix * rho_liquid[i][j] * C_i_0[i] * (r_a[i]+r_c[i])  # Verified with Zhao's results, uses Zhao's data and but MSTDB heat capacity

    # lambda_mix_T=np.zeros(len(T))
    # for j in range(len(T)):
    #     lambda_mix_T[j]=lambda_i_m[i]*(1-alpha_i_m[i]*(gamma_i_m[i]+(1/3))*(T[j]-T_melt))
    nan_check = np.isnan(lambda_i_m)
    contains_nan = nan_check.any()
    if contains_nan:
        return(T,lambda_i_mb)
    else:
        return(T,lambda_i_m)

def Bridgman(df,MSTDB_df,compound_input,mol_fracs,Temp_Range, V_m=0, C_p=0, alpha=0, expon=0):

    #Define smaller dataframes of each variable contained in the spreadsheet
    Compound_df=df['Compound']
    M_i_df=df['M_i (g/mol)']
    T_i_m_df=df['T_i_m (K)']
    V_i_m_df=df['V_i_m (m^3/mol)']
    alpha_i_m_df=df['alpha_i_m (K^-1)']
    C_i_0_df=df['C_i_0 (m/s)']
    C_i_p_df=df['C_i_p (J/mol/K)']
    r_c_df=df['r_c (m)']
    r_a_df=df['r_a (m)']
    rho_i_m_df=df['rho_m (g/m^3)']
    C_0_T_A_df=df['SoS(T)_A (A+B*T)']
    C_0_T_B_df=df['SoS(T)_B']

    #initialize arrays of indices of compounds in the dataframe, and boolean values of whether they are in the dataframe
    indices=np.zeros(len(compound_input))
    indices_exist=np.full(len(compound_input),False)

    #find indices in dataframe and store index values
    for i in range(len(compound_input)):
        for j in range(len(Compound_df)):
            if compound_input[i]==Compound_df[j]:
                indices[i]=j
                indices_exist[i]=True

    #check to make sure all the input compounds actually exist. Quit program if not
    for i in range(len(compound_input)):
        problem=False
        if indices_exist[i]==False:
            print('Warning: ' + str(compound_input[i]) + ' is not included in the compound data spreadsheet' )
            problem=True
        if problem==True:
            print('Force quit.')
            quit()

    #Generate arrays of all parameters that already exist in the dataframe, but in the order of the input compound array
    Compound=[]
    for i in range(len(compound_input)):
        Compound.append('')
    M_i=np.zeros(len(compound_input))
    T_i_m=np.zeros(len(compound_input))
    V_i_m=np.zeros(len(compound_input))
    alpha_i_m=np.zeros(len(compound_input))
    C_i_0=np.zeros(len(compound_input))
    C_i_p=np.zeros(len(compound_input))
    r_c=np.zeros(len(compound_input))
    r_a=np.zeros(len(compound_input))
    rho_i_m=np.zeros(len(compound_input))
    C_0_T_A=np.zeros(len(compound_input))
    C_0_T_B=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        Compound[i]=Compound_df[indices[i]]
        M_i[i]=M_i_df[indices[i]]*0.001         # Convert to kg/mol
        T_i_m[i]=T_i_m_df[indices[i]]
        V_i_m[i]=V_i_m_df[indices[i]]
        alpha_i_m[i]=alpha_i_m_df[indices[i]]
        C_i_0[i]=C_i_0_df[indices[i]]
        C_i_p[i]=C_i_p_df[indices[i]]
        r_c[i]=r_c_df[indices[i]]
        r_a[i]=r_a_df[indices[i]]
        rho_i_m[i]=rho_i_m_df[indices[i]]*0.001 # Convert to kg/m^3
        C_0_T_A[i]=C_0_T_A_df[indices[i]]
        C_0_T_B[i]=C_0_T_B_df[indices[i]]

    #initialize a temperature range
    T_melt = Temp_Range[0]
    T=np.linspace(T_melt,Temp_Range[1],num=100)

    #Calculate Gruneisen parameter
    gamma_i_m=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        gamma_i_m[i]=((M_i[i]*alpha_i_m[i]*C_i_0[i]**2)/C_i_p[i])

    #Calculate constant volume heat capacities
    C_i_v=np.zeros((len(compound_input),len(T)))
    for i in range(len(compound_input)):
        for j in range(len(T)):
            C_i_v[i][j]=C_i_p[i]/(1+alpha_i_m[i]*gamma_i_m[i]*T[j])       # <-- This isn't supposed to use the melting point parameters

    #Calculate volume fractions
    phi_i_m=np.zeros(len(compound_input))
    denom=0
    for i in range(len(compound_input)):
        phi_i_m[i]=V_i_m[i]*mol_fracs[i]
        denom=denom+V_i_m[i]*mol_fracs[i]
    phi_i_m=phi_i_m/denom  

    #Calculate weight fractions
    kappa_i_m=np.zeros(len(compound_input))
    denom=0
    for i in range(len(compound_input)):
        kappa_i_m[i]=rho_i_m[i]*mol_fracs[i]
        denom=denom+rho_i_m[i]*mol_fracs[i]
    kappa_i_m=kappa_i_m/denom 

    #Calculate mixture speed of sound at melting temp
    C_0_mix_m = 0
    C_0_mix_T = np.zeros(len(T))
    i = 0
    # If measurement data for speed of sound with temp exists
    for i in range(len(compound_input)):
        C_0_mix_m=C_0_mix_m+phi_i_m[i]**2/(kappa_i_m[i]*C_i_0[i]**2)
    C_0_mix_m=1/np.sqrt(C_0_mix_m)
    # If measurement data for speed of sound with temp does not exist
    for j in range(len(T)):
        for i in range(len(compound_input)):
            if C_0_T_A[i] == 0:
                C_0_mix_T[j]=C_0_mix_T[j]+phi_i_m[i]**2/(kappa_i_m[i]*C_i_0[i]**2)
            else:
                SoS = C_0_T_B[i]*T[j]+C_0_T_A[i]
                C_0_mix_T[j]=C_0_mix_T[j]+phi_i_m[i]**2/(kappa_i_m[i]*SoS**2)
        C_0_mix_T[j]=1/np.sqrt(C_0_mix_T[j])

    # Find the temperature dependent density from data
    rho_liquid = np.zeros((len(Compound),len(T)))
    for i in range(len(Compound)):
        density_A, density_B = density(df,compound_input[i])
        for j in range(len(T)):
            rho_liquid[i][j] = density_A - density_B*T[j] 

    #Calculate number of ions per compound
    n_i=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        n_i[i]=(sum(1 for c in Compound[i] if c.isupper())+ sum(1 for c in Compound[i] if c.isnumeric()))

    #Calculate average number of ions
    C_p_mix=0
    C_v_mix=np.zeros(len(T))    #Const vol heat capacity for mixture <-- NOT SURE IF THIS IS OKAY
    n_mix=0
    rho_mix=np.zeros(len(T))
    M_i_mix=0
    if C_p==0:
        for i in range(len(compound_input)):
            C_p_mix=C_p_mix+mol_fracs[i]*C_i_p[i]
            M_i_mix=M_i_mix+mol_fracs[i]*M_i[i]
            for j in range(len(T)):
                C_v_mix[j]=C_v_mix[j]+mol_fracs[i]*C_i_v[i][j]  # Linear mixing
                rho_mix[j]=rho_mix[j]+mol_fracs[i]*rho_liquid[i][j]
            n_mix=n_mix+mol_fracs[i]*n_i[i]
    else: 
        C_p_mix=C_p
        for i in range(len(compound_input)):
            M_i_mix=M_i_mix+mol_fracs[i]*M_i[i]
            for j in range(len(T)):
                C_v_mix[j]=C_v_mix[j]+mol_fracs[i]*C_i_v[i][j]  # Linear mixing
            n_mix=n_mix+mol_fracs[i]*n_i[i]

    #Calculate average molar volume
    V_mix_m=0
    if V_m==0:
        for i in range(len(compound_input)):       
            V_mix_m=V_mix_m+mol_fracs[i]*M_i[i]*(1/rho_i_m[i])
    else: V_mix_m=V_m

    #Calculate thermal conductivity of mix
    lambda_mix_T=np.zeros(len(T))

    if expon == 0:
        for j in range(len(T)):
            lambda_mix_T[j] = 3*1.380649e-23*C_0_mix_m*((n_mix*6.0221408e23)/V_mix_m)**(2/3)*(0.931/np.sqrt(C_p_mix/C_v_mix[j]))     
    else:
        if all(i == 0 for i in C_0_T_A):
            for j in range(len(T)):
                #lambda_mix_T[j] = 3*1.380649e-23*C_0_mix_m/(M_i_mix/(6.0221408e23*rho_mix[j]))**(2/3)*(0.931/np.sqrt(C_p_mix/C_v_mix[j]))       # Liquid density temp dependence
                lambda_mix_T[j] = 3*1.380649e-23*C_0_mix_m*((n_mix*6.0221408e23)/V_mix_m)**(2/3)*(0.931/np.sqrt(C_p_mix/C_v_mix[j]))           
        else:
            for j in range(len(T)):
                #lambda_mix_T[j] = 3*1.380649e-23*C_0_mix_T[j]/(M_i_mix/(6.0221408e23*rho_mix[j]))**(2/3)*(0.931/np.sqrt(C_p_mix/C_v_mix[j]))    # Liquid density temp dependence
                lambda_mix_T[j] = 3*1.380649e-23*C_0_mix_T[j]*((n_mix*6.0221408e23)/V_mix_m)**(2/3)*(0.931/np.sqrt(C_p_mix/C_v_mix[j]))        

    # ^^^ The heat capacity at constant volume (C_v_mix) is where density is changing, so that's how it's temp-dependent. 

    return(T,lambda_mix_T)

def ChapmanEnskog(df,MSTDB_df,compound_input,mol_fracs,Temp_Range, V_m=0, C_p=0, alpha=0, expon=0):
    
    #Define smaller dataframes of each variable contained in the spreadsheet
    Compound_df=df['Compound']
    M_i_df=df['M_i (g/mol)']
    T_i_m_df=df['T_i_m (K)']
    V_i_m_df=df['V_i_m (m^3/mol)']
    alpha_i_m_df=df['alpha_i_m (K^-1)']
    C_i_0_df=df['C_i_0 (m/s)']
    C_i_p_df=df['C_i_p (J/mol/K)']
    r_c_df=df['r_c (m)']
    r_a_df=df['r_a (m)']
    rho_i_m_df=df['rho_m (g/m^3)']

    # viscosityA_df=MSTDB_df['Viscosity (mN*s/m2): A*exp(B/(R*T(K)))']
    # viscosityB_df=MSTDB_df['Viscosity (mN*s/m2): A*exp(B/(R*T(K)))_B']

    #initialize arrays of indices of compounds in the dataframe, and boolean values of whether they are in the dataframe
    indices=np.zeros(len(compound_input))
    indices_exist=np.full(len(compound_input),False)


    #find indices in dataframe and store index values
    for i in range(len(compound_input)):
        for j in range(len(Compound_df)):
            if compound_input[i]==Compound_df[j]:
                indices[i]=j
                indices_exist[i]=True

    #check to make sure all the input compounds actually exist. Quit program if not
    for i in range(len(compound_input)):
        problem=False
        if indices_exist[i]==False:
            print('Warning: ' + str(compound_input[i]) + ' is not included in the compound data spreadsheet' )
            problem=True
        if problem==True:
            print('Force quit.')
            quit()

    #Generate arrays of all parameters that already exist in the dataframe, but in the order of the input compound array
    Compound=[]
    for i in range(len(compound_input)):
        Compound.append('')
    M_i=np.zeros(len(compound_input))
    T_i_m=np.zeros(len(compound_input))
    V_i_m=np.zeros(len(compound_input))
    alpha_i_m=np.zeros(len(compound_input))
    C_i_0=np.zeros(len(compound_input))
    C_i_p=np.zeros(len(compound_input))
    r_c=np.zeros(len(compound_input))
    r_a=np.zeros(len(compound_input))
    rho_i_m=np.zeros(len(compound_input))
    viscosityA=np.zeros(len(compound_input))
    viscosityB=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        Compound[i]=Compound_df[indices[i]]
        M_i[i]=M_i_df[indices[i]]   / 1000
        T_i_m[i]=T_i_m_df[indices[i]]
        V_i_m[i]=V_i_m_df[indices[i]]
        alpha_i_m[i]=alpha_i_m_df[indices[i]]
        C_i_0[i]=C_i_0_df[indices[i]]
        C_i_p[i]=C_i_p_df[indices[i]]
        r_c[i]=r_c_df[indices[i]]
        r_a[i]=r_a_df[indices[i]]
        rho_i_m[i]=rho_i_m_df[indices[i]]   / 1000

        n =(type(MSTDB_df.iloc[15,0]))

        comp = Compound[i]

        print(type(Compound[i]))

        

        viscosityA[i]=get_cell_value('Molten_Salt_Thermophysical_Properties.csv', Compound[i], 'Viscosity (mN*s/m2): A*exp(B/(R*T(K)))')
        viscosityB[i]=get_cell_value('Molten_Salt_Thermophysical_Properties.csv', Compound[i], 'Viscosity (mN*s/m2): A*exp(B/(R*T(K)))_B')


    #initialize a temperature range
    T_melt = Temp_Range[0]
    T=np.linspace(T_melt,Temp_Range[1],num=100)

    # Find the temperature dependent density from data
    rho_liquid = np.zeros((len(T),len(Compound)))
    for i in range(len(Compound)):
        density_A, density_B = density(df,compound_input[i])
        for j in range(len(T)):
            rho_liquid[j][i] = density_A - density_B*T[j] 

    
    # Mass of molecule
    m = np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        m[i] = M_i[i]/6.0221408e23

    # Diameter of molecule
    d = np.zeros(len(compound_input))
    b = np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        d[i] = 2 * (r_c[i]+r_a[i])  

    D = np.zeros((len(T),len(compound_input)))  # Effective scattering diameter
    lambda_0 = np.zeros((len(T),len(compound_input)))
    k_0 = np.zeros((len(T),len(compound_input)))
    #Calculate dilute gas thermal conductivities
    for i in range(len(compound_input)):
        for j in range(len(T)):
            mu = viscosityA[i]*np.exp(viscosityB[i]/(8.314*T[j]))/1000
            lambda_0[j][i] = 15/4 * (1.380649e-23)/(m[i]) * mu
            D[j][i] = np.sqrt((5 * np.sqrt(np.pi*m[i]*1.380649e-23*T[j]) / (16*mu*np.pi)))
            k_0[j][i] = lambda_0[j][i]
        
    #Calculate thermal conductivities
    b = np.zeros((len(T),len(compound_input)))
    g = np.zeros((len(T),len(compound_input)))
    PI = np.zeros((len(T),len(compound_input)))
    lambda_compound = np.zeros((len(T),len(compound_input)))
    for i in range(len(compound_input)):
        for j in range(len(T)):
            b[j][i] = 2/3*np.pi*D[j][i]**3
            g[j][i] = (1 - 11*rho_liquid[j][i]*b[i]/9) / (1 - 2*rho_liquid[j][i]*b[j][i])           
            PI[j][i] = 4/9 * np.sqrt(np.pi*m[i]*1.380649e-23*T[j]) * rho_liquid[j][i]**2 * D[j][i]**4 * g[j][i]
            lambda_compound[j][i] = k_0[j][i] / g[j][i] * (1 + 2/5*np.pi*rho_liquid[j][i]*D[j][i]**3*g[j][i])**2 + (3*1.380649e-23) / (2*m[i]) * PI[j][i]

    # Linear mixing of compound thermal conductivities
    lambda_mix_T=np.zeros(len(T))
    for i in range(len(compound_input)):
        for j in range(len(T)):
            lambda_mix_T[j]=lambda_compound[j][i]*mol_fracs[i]

    return(T,lambda_mix_T)

def density(df,sheet_name):

    #Define smaller dataframes of each variable contained in the spreadsheet
    Compound_df=df['Compound']

    #Define smaller dataframes of each variable contained in the spreadsheet
    rho_T_A_df=df['A (Density (g/cm3):   A - BT(K))']
    rho_T_B_df=df['B (Density (g/cm3):   A - BT(K))'] 

    #find indices in dataframe and store index values
    for j in range(len(Compound_df)):
        if sheet_name==Compound_df[j]:
            index=j

    density_A=rho_T_A_df[index]*100**3/1000
    density_B=rho_T_B_df[index]*100**3/1000

    if pd.isna(density_A):
        # Read the Excel file
        df_density = pd.read_excel('Density.xlsx', sheet_name=sheet_name)
                
        # Extract the column and calculate the average
        density_column = df_density['A (Density (g/cm3):   A - BT(K))']
        density_A = density_column.mean()*100**3/1000

        # Extract the column and calculate the average
        density_column = df_density['B (Density (g/cm3):   A - BT(K))']
        density_B = density_column.mean()*100**3/1000
    else:
        pass
    
    return density_A, density_B


def get_cell_value(csv_file, target_row_header, target_column_header):
    with open(csv_file, 'r') as file:
        reader = csv.DictReader(file)
        
        for row in reader:
            if row['Formula'] == target_row_header:
                return row[target_column_header]
    
    # If the target row or column header is not found, return None or handle it as needed
    return None

def functionlibrary():
    functions = {
    'Kinetic Theory' : Kinetic_Theory,
    'Kinetic Theory, Empirical Fit' : Kinetic_Theory_Emp,
    'Kinetic Theory Unary Rho(T)' : Kinetic_Theory_RhoT,
    'Kinetic Theory Unary Salts' : Kinetic_Theory_Unary,
    'Phonon Gas Model' : Phonon_Gas_Model,
    'Bridgman' : Bridgman,
    'Chapman-Enskog' : ChapmanEnskog,
    }
    return functions
    