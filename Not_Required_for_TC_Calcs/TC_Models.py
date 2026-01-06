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
from scipy.constants import Boltzmann as k_B
from scipy.constants import Avogadro as Avog
import inspect

def Kinetic_Theory_Emp(df,MSTDB_df,SCL_PDF_df,compound_input,mol_fracs,Temp_Range, V_m=0, density_mix=0, C_p_mix=0, alpha=0, expon=0):
    
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

    # print("KT_emp, number of ions: ",n_i)
          
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
        kappa_i_m[i]=M_i[i]*mol_fracs[i]   
        denom=denom+M_i[i]*mol_fracs[i]
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
    if C_p_mix==0:
        for i in range(len(compound_input)):
            C_p_mix=C_p_mix+mol_fracs[i]*C_i_p[i]
            n_mix=n_mix+mol_fracs[i]*n_i[i]
    else: 
        C_p_mix=C_p_mix
        for i in range(len(compound_input)):
            n_mix=n_mix+mol_fracs[i]*n_i[i]


    #Calculate average molar volume
    V_mix_m=0
    if V_m==0:
        for i in range(len(compound_input)):    
            V_mix_m=V_mix_m+mol_fracs[i]*M_i[i]*(1/rho_i_m[i])
    else: V_mix_m=V_m  


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



    sum_term=0
    for i in range(len(compound_input)):
        sum_term=sum_term+mol_fracs[i]*((M_i[i]/M_mix)-1)**2
    #delta_mix_m=0.4872*(lambda_ideal_mix_m/(1.380649e-23*C_0_mix_m*(n_mix*6.0221408e23*(1/V_mix_m))**(3/2)))*sum_term
    delta_mix_m=0.4872*(lambda_ideal_mix_m/(1.380649e-23*C_0_mix_m*(n_mix*6.0221408e23*(1/V_mix_m))**(2/3)))*sum_term
    

    lambda_mix=lambda_ideal_mix_m*(1-delta_mix_m)

    lambda_mix_T=np.zeros(len(T))
    for j in range(len(T)):
        lambda_mix_T[j]=lambda_mix*(1-alpha_mix_m*(gamma_mix_m+(1/3))*(T[j]-T_melt))    #<<--T_melt must be the melting point of the liquid!

    print("KT_emp: ",lambda_mix_T[0])
    return(T,lambda_mix_T)

def Kinetic_Theory_Emp_Mix(df,MSTDB_df,SCL_PDF_df,compound_input,mol_fracs,Temp_Range, V_m=0, density_mix=0, C_p_mix=0, alpha=0, expon=0):
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

    print('')
    print('################ Kinetic Theory Empirical Fit, Mix Data ############################')
    specific_heat_mix = [0,0]
    density_mix = [0,0]
    sound_velocity_mix = [0,0]    
    molar_volume_mix = [0,0] 
    
    # Obtain measurement data for function
    func_name = inspect.currentframe().f_code.co_name
    r_ac_mix, density_mix, sound_velocity_mix, specific_heat_mix = prop_lookup(func_name,Compound,mol_fracs)


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

    # print("KT_emp, number of ions: ",n_i)
          
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
        kappa_i_m[i]=M_i[i]*mol_fracs[i]        # <<< Shouldn't this be molar mass?
        denom=denom+M_i[i]*mol_fracs[i]
    kappa_i_m=kappa_i_m/denom 

    #Calculate taus for compounds
    Tau_i_m=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        Tau_i_m[i] = ( C_i_v[i]*C_i_0[i] ) / n_i[i]


    #Calculate mixture speed of sound
    C_0_mix_m=0
    C_0_mix_est = 0
    C_0_mix_data = 0 
    for i in range(len(compound_input)):
        C_0_mix_est += phi_i_m[i]**2/(kappa_i_m[i]*C_i_0[i]**2)
    C_0_mix_data = sound_velocity_mix[1]*T_melt+sound_velocity_mix[0]
    C_0_mix_est=1/np.sqrt(C_0_mix_est)
    if sound_velocity_mix == [0,0]:
        C_0_mix_m = C_0_mix_est
    else:
        C_0_mix_m = C_0_mix_data

        C_0_mix_est_avg = np.average(C_0_mix_est)
        C_0_mix_data_avg = np.average(C_0_mix_data)
        error = 100 * (C_0_mix_data_avg-C_0_mix_est_avg)/C_0_mix_data_avg
        print("Sound Velocity - Estimated: ", C_0_mix_est_avg)
        print("Sound Velocity - Mix Data: ", C_0_mix_data_avg)
        print("Sound Velocity % Difference: ", error, " %")


    #Calculate mixture molecular weight
    M_mix=0
    for i in range(len(compound_input)):
        M_mix += mol_fracs[i]*M_i[i]

    #Calculate average molar volume
    V_mix_m=0
    V_mix_est = 0
    V_mix_data = 0
    V_mix_data_Vm = 0
    for i in range(len(compound_input)):    
        V_mix_est += mol_fracs[i]*M_i[i]*(1/rho_i_m[i])
    V_mix_data = M_mix/(density_mix[0] + density_mix[1]*T_melt)
    V_mix_data_Vm = molar_volume_mix[1]*T_melt+molar_volume_mix[0]
    if molar_volume_mix == [0,0]:
        if density_mix == [0,0]:
            V_mix_m = V_mix_est
        else:
            V_mix_m = V_mix_data
    else:
        V_mix_m = V_mix_data_Vm

        V_mix_est_avg = np.average(V_mix_est)
        V_mix_data_avg = np.average(V_mix_data)
        error_densityVm = 100 * (V_mix_data_avg-V_mix_est_avg)/V_mix_data_avg
        error_Vm = 100 * (V_mix_data_Vm-V_mix_est_avg)/V_mix_data_Vm
        #print("")
        #print("Molar volume - From estimated density: ", V_mix_est_avg)
        #print("Molar volume - From mix data density: ", V_mix_data_avg)
        #print("Molar volume - Mix data Molar Volume: ", V_mix_data_Vm)
        #print("Molar volume % Difference, Density Data/Estimation: ", error_densityVm, " %")
        #print("Molar volume % Difference, Molar Volume Data/Estimation: ", error_Vm, " %")

    #Calculate mixture specific heat
    #Calculate average number of ions
    C_p_mix=0
    C_p_mix_est = 0
    C_p_mix_data = 0
    n_mix=0
    for i in range(len(compound_input)):
        C_p_mix_est += C_p_mix+mol_fracs[i]*C_i_p[i]
        n_mix = n_mix + mol_fracs[i]*n_i[i]
    C_p_mix_data = specific_heat_mix
    if specific_heat_mix[0] == 0:
        C_p_mix = C_p_mix_est
    else:
        if specific_heat_mix[2] == 'm':
            C_p_mix = C_p_mix_data[1]*T_melt+C_p_mix_data[0]
        elif specific_heat_mix[2] == 'g':   # Converts units from J/g/K to J/mol/K
            C_p_mix = (C_p_mix_data[0]*T_melt+C_p_mix_data[1])*(M_mix)

        C_p_mix_est_avg = np.average(C_p_mix_est)
        C_p_mix_data_avg = np.average(C_p_mix)
        error = 100 * (C_p_mix_data_avg-C_p_mix_est_avg)/C_p_mix_data_avg
        #print("")
        print("Specific Heat - Estimated: ", C_p_mix_est_avg)
        print("Specific Heat - Mix Data: ", C_p_mix_data_avg)
        print("Specific Heat % Difference: ", error, " %")


    #Calculate mixture thermal expansion
    alpha_mix_m=0
    if alpha==0:
        for i in range(len(compound_input)):
            alpha_mix_m=alpha_mix_m+phi_i_m[i]*alpha_i_m[i]
    else:alpha_mix_m=alpha
      

    #Calculate mixture gamma
    gamma_mix_m=M_mix*alpha_mix_m*C_0_mix_m**2*(1/C_p_mix)
    

    #Calculate Tau for the mixture
    Tau_mix_m=(C_p_mix/(1+alpha_mix_m*gamma_mix_m*T_melt)) *C_0_mix_m *(1/n_mix)
    
    
    #Calculate ideal thermal conductivity
    sum_term=0
    for i in range(len(compound_input)):
        sum_term=sum_term+mol_fracs[i]*(lambda_i_m[i]/Tau_i_m[i])**(-3/2)
    lambda_ideal_mix_m=( Tau_mix_m**(-3/2) *sum_term ) **(-2/3)



    sum_term=0
    for i in range(len(compound_input)):
        sum_term=sum_term+mol_fracs[i]*((M_i[i]/M_mix)-1)**2
    #delta_mix_m=0.4872*(lambda_ideal_mix_m/(1.380649e-23*C_0_mix_m*(n_mix*6.0221408e23*(1/V_mix_m))**(3/2)))*sum_term
    delta_mix_m=0.4872*(lambda_ideal_mix_m/(1.380649e-23*C_0_mix_m*(n_mix*6.0221408e23*(1/V_mix_m))**(2/3)))*sum_term
    

    lambda_mix=lambda_ideal_mix_m*(1-delta_mix_m)

    lambda_mix_T=np.zeros(len(T))
    for j in range(len(T)):
        lambda_mix_T[j]=lambda_mix*(1-alpha_mix_m*(gamma_mix_m+(1/3))*(T[j]-T_melt))    #<<--T_melt must be the melting point of the liquid!

    return(T,lambda_mix_T)    

def Kinetic_Theory(df,MSTDB_df,SCL_PDF_df,compound_input,mol_fracs,Temp_Range, V_m=0, density_mix=0, C_p_mix=0, alpha=0, expon=0):
    
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


    # #calculate constant volume heat capacities
    # C_i_v=np.zeros((len(compound_input),len(T)))
    # for i in range(len(compound_input)):
    #     for j in range(len(T)):
    #         C_i_v[i][j]=C_i_p[i]/(1+alpha_i_m[i]*gamma_i_m[i]*T[j])

    #calculate constant volume heat capacities at melting temp
    C_i_v=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        C_i_v[i]=C_i_p[i]/(1+alpha_i_m[i]*gamma_i_m[i]*T_i_m[i])

    
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

    # print("KT, number of ions: ",n_i)
    # print("KT, number of cations: ",n_i_c)
    # print("KT, number of anions: ",n_i_a)

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
        kappa_i_m[i]=M_i[i]*mol_fracs[i]
        denom=denom+M_i[i]*mol_fracs[i]
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
    if C_p_mix==0:
        for i in range(len(compound_input)):
            C_p_mix=C_p_mix+mol_fracs[i]*C_i_p[i]
            n_mix=n_mix+mol_fracs[i]*n_i[i]
    else: 
        C_p_mix=C_p_mix
        for i in range(len(compound_input)):
            n_mix=n_mix+mol_fracs[i]*n_i[i]
    

    #Calculate average molar volume
    V_mix_m=0
    if V_m==0:
        for i in range(len(compound_input)):
            V_mix_m=V_mix_m+mol_fracs[i]*M_i[i]*(1/rho_i_m[i])
    else: V_mix_m=V_m


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


    sum_term=0
    for i in range(len(compound_input)):
        sum_term=sum_term+mol_fracs[i]*((M_i[i]/M_mix)-1)**2
    #delta_mix_m=0.4872*(lambda_ideal_mix_m/(1.380649e-23*C_0_mix_m*(n_mix*6.0221408e23*(1/V_mix_m))**(3/2)))*sum_term
    delta_mix_m = 0.4872*(lambda_ideal_mix_m/(1.380649e-23*C_0_mix_m*(n_mix*6.0221408e23*(1/V_mix_m))**(2/3)))*sum_term

    # if lambda_ideal_mix_m >= (1.380649e-23*C_0_mix_m*(n_mix*6.0221408e23*(1/V_mix_m))**(2/3)):
    #     delta_mix_m = 0.4872*(lambda_ideal_mix_m/(1.380649e-23*C_0_mix_m*(n_mix*6.0221408e23*(1/V_mix_m))**(2/3)))*sum_term
    # else:
    #     delta_mix_m = 0
    

    lambda_mix=lambda_ideal_mix_m*(1-delta_mix_m)

    lambda_mix_T=np.zeros(len(T))
    for j in range(len(T)):
        lambda_mix_T[j]=lambda_mix*(1-alpha_mix_m*(gamma_mix_m+(1/3))*(T[j]-T_melt))    #<<--T_melt must be the melting point of the liquid!

    print("KT_CP: ",lambda_mix_T[0])
    return(T,lambda_mix_T)

def Kinetic_Theory_Mix(df,MSTDB_df,SCL_PDF_df,compound_input,mol_fracs,Temp_Range, V_m=0, density_mix=0, C_p_mix=0, alpha=0, expon=0):
    
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

    print('')
    print('################ Kinetic Theory CP, Mix Data ############################')
    specific_heat_mix = C_p_mix
    density_mix = [0,0]
    sound_velocity_mix = [0,0]    
    molar_volume_mix = [0,0] 
    
    # Obtain measurement data for function
    func_name = inspect.currentframe().f_code.co_name
    r_ac_mix, density_mix, sound_velocity_mix, specific_heat_mix = prop_lookup(func_name,Compound,mol_fracs)

    #Calculate Gruneisen parameter
    gamma_i_m=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        gamma_i_m[i]=((M_i[i]*alpha_i_m[i]*C_i_0[i]**2)/C_i_p[i])

    #initialize a temperature range
    T_melt = Temp_Range[0]
    T=np.linspace(T_melt,Temp_Range[1],num=100)


    # #calculate constant volume heat capacities
    # C_i_v=np.zeros((len(compound_input),len(T)))
    # for i in range(len(compound_input)):
    #     for j in range(len(T)):
    #         C_i_v[i][j]=C_i_p[i]/(1+alpha_i_m[i]*gamma_i_m[i]*T[j])

    #calculate constant volume heat capacities at melting temp
    C_i_v=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        C_i_v[i]=C_i_p[i]/(1+alpha_i_m[i]*gamma_i_m[i]*T_i_m[i])

    
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

    # print("KT, number of ions: ",n_i)
    # print("KT, number of cations: ",n_i_c)
    # print("KT, number of anions: ",n_i_a)

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
        kappa_i_m[i]=M_i[i]*mol_fracs[i]
        denom=denom+M_i[i]*mol_fracs[i]
    kappa_i_m=kappa_i_m/denom 

    #Calculate taus for compounds
    Tau_i_m=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        Tau_i_m[i] = ( C_i_v[i]*C_i_0[i] ) / n_i[i]


    #Calculate mixture speed of sound
    C_0_mix_m=0
    C_0_mix_est = 0
    C_0_mix_data = 0 
    for i in range(len(compound_input)):
        C_0_mix_est += phi_i_m[i]**2/(kappa_i_m[i]*C_i_0[i]**2)
    C_0_mix_data = sound_velocity_mix[1]*T_melt+sound_velocity_mix[0]
    C_0_mix_est=1/np.sqrt(C_0_mix_est)
    if sound_velocity_mix == [0,0]:
        C_0_mix_m = C_0_mix_est
    else:
        C_0_mix_m = C_0_mix_data

        C_0_mix_est_avg = np.average(C_0_mix_est)
        C_0_mix_data_avg = np.average(C_0_mix_data)
        error = 100 * (C_0_mix_data_avg-C_0_mix_est_avg)/C_0_mix_data_avg
        print("Sound Velocity - Estimated: ", C_0_mix_est_avg)
        print("Sound Velocity - Mix Data: ", C_0_mix_data_avg)
        print("Sound Velocity % Difference: ", error, " %")
    

    #Calculate mixture molecular weight
    M_mix=0
    for i in range(len(compound_input)):
        M_mix=M_mix+mol_fracs[i]*M_i[i]

    #Calculate average molar volume
    V_mix_m=0
    V_mix_est = 0
    V_mix_data = 0
    V_mix_data_Vm = 0
    for i in range(len(compound_input)):    
        V_mix_est += mol_fracs[i]*M_i[i]*(1/rho_i_m[i])
    V_mix_data = M_mix/(density_mix[0] + density_mix[1]*T_melt)
    V_mix_data_Vm = molar_volume_mix[1]*T_melt+molar_volume_mix[0]
    if molar_volume_mix == [0,0]:
        if density_mix == [0,0]:
            V_mix_m = V_mix_est
        else:
            V_mix_m = V_mix_data
    else:
        V_mix_m = V_mix_data_Vm

        V_mix_est_avg = np.average(V_mix_est)
        V_mix_data_avg = np.average(V_mix_data)
        error_densityVm = 100 * (V_mix_data_avg-V_mix_est_avg)/V_mix_data_avg
        error_Vm = 100 * (V_mix_data_Vm-V_mix_est_avg)/V_mix_data_Vm
        #print("")
        #print("Molar volume - From estimated density: ", V_mix_est_avg)
        #print("Molar volume - From mix data density: ", V_mix_data_avg)
        #print("Molar volume - Mix data Molar Volume: ", V_mix_data_Vm)
        #print("Molar volume % Difference, Density Data/Estimation: ", error_densityVm, " %")
        #print("Molar volume % Difference, Molar Volume Data/Estimation: ", error_Vm, " %")

    #Calculate mixture specific heat
    #Calculate average number of ions
    C_p_mix=0
    C_p_mix_est = 0
    C_p_mix_data = 0
    n_mix=0
    for i in range(len(compound_input)):
        C_p_mix_est += C_p_mix+mol_fracs[i]*C_i_p[i]
        n_mix = n_mix + mol_fracs[i]*n_i[i]
    C_p_mix_data = specific_heat_mix
    if specific_heat_mix[0] == 0:
        C_p_mix = C_p_mix_est
    else:
        if specific_heat_mix[2] == 'm':
            C_p_mix = C_p_mix_data[1]*T_melt+C_p_mix_data[0]
        elif specific_heat_mix[2] == 'g':   # Converts units from J/g/K to J/mol/K
            C_p_mix = (C_p_mix_data[0]*T_melt+C_p_mix_data[1])*(M_mix)

        C_p_mix_est_avg = np.average(C_p_mix_est)
        C_p_mix_data_avg = np.average(C_p_mix)
        error = 100 * (C_p_mix_data_avg-C_p_mix_est_avg)/C_p_mix_data_avg
        #print("")
        print("Specific Heat - Estimated: ", C_p_mix_est_avg)
        print("Specific Heat - Mix Data: ", C_p_mix_data_avg)
        print("Specific Heat % Difference: ", error, " %")

    #Calculate mixture thermal expansion
    alpha_mix_m=0
    if alpha==0:
        for i in range(len(compound_input)):
            alpha_mix_m=alpha_mix_m+phi_i_m[i]*alpha_i_m[i]
    else:alpha_mix_m=alpha

        

    #Calculate mixture gamma
    gamma_mix_m=M_mix*alpha_mix_m*C_0_mix_m**2*(1/C_p_mix)
    

    #Calculate Tau for the mixture
    Tau_mix_m=(C_p_mix/(1+alpha_mix_m*gamma_mix_m*T_melt)) *C_0_mix_m *(1/n_mix)
    
    
    #Calculate ideal thermal conductivity
    sum_term=0
    for i in range(len(compound_input)):
        sum_term=sum_term+mol_fracs[i]*(lambda_i_m[i]/Tau_i_m[i])**(-3/2)
    lambda_ideal_mix_m=( Tau_mix_m**(-3/2) *sum_term ) **(-2/3)


    sum_term=0
    for i in range(len(compound_input)):
        sum_term=sum_term+mol_fracs[i]*((M_i[i]/M_mix)-1)**2
    #delta_mix_m=0.4872*(lambda_ideal_mix_m/(1.380649e-23*C_0_mix_m*(n_mix*6.0221408e23*(1/V_mix_m))**(3/2)))*sum_term
    delta_mix_m=0.4872*(lambda_ideal_mix_m/(1.380649e-23*C_0_mix_m*(n_mix*6.0221408e23*(1/V_mix_m))**(2/3)))*sum_term
    

    lambda_mix=lambda_ideal_mix_m*(1-delta_mix_m)

    lambda_mix_T=np.zeros(len(T))
    for j in range(len(T)):
        lambda_mix_T[j]=lambda_mix*(1-alpha_mix_m*(gamma_mix_m+(1/3))*(T[j]-T_melt))    #<<--T_melt must be the melting point of the liquid!

    return(T,lambda_mix_T)

def Kinetic_Theory24(df,MSTDB_df,SCL_PDF_df,compound_input,mol_fracs,Temp_Range, V_m=0, density_mix=0, C_p_mix=0, alpha=0, expon=0):
    print('')
    print('################ Kinetic Theory CP, Mix Data ############################')

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

    # print("KT, number of ions: ",n_i)
    # print("KT, number of cations: ",n_i_c)
    # print("KT, number of anions: ",n_i_a)
        
    # Calculate compound psi term
    psi_i=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        psi_i[i] = 1 + n_i_c[i]/n_i_a[i]

    # Calculate compound number density
    n_dens_i = np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        n_dens_i[i] = Avog *n_i[i]/V_i_m[i]

    #Calculate compound thermal conductivity(T) based on compound minimum thermal conductivites
    lambda_i = np.zeros((len(Compound),len(T)))
    for i in range(len(compound_input)):
        for j in range(len(T)):
            lambda_i[i][j] = (1 + n_i_c[i]/n_i_a[i]) * k_B * n_dens_i[i]**(2/3) * C_i_0[i] * (1 - alpha_i_m[i]*(gamma_i_m[i] + 1/3)*(T[j] - T_i_m[i]))

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
        kappa_i_m[i]=M_i[i]*mol_fracs[i]
        denom=denom+M_i[i]*mol_fracs[i]
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
    if C_p_mix==0:
        for i in range(len(compound_input)):
            C_p_mix=C_p_mix+mol_fracs[i]*C_i_p[i]
            n_mix=n_mix+mol_fracs[i]*n_i[i]
    else: 
        C_p_mix=C_p_mix
        for i in range(len(compound_input)):
            n_mix=n_mix+mol_fracs[i]*n_i[i]
 

    #Calculate average molar volume
    V_mix_m=0
    if V_m==0:
        for i in range(len(compound_input)):
            V_mix_m=V_mix_m+mol_fracs[i]*M_i[i]*(1/rho_i_m[i])
    else: V_mix_m=V_m

        
    #Calculate mixture number density
    n_dens_mix = Avog * n_mix/V_mix_m


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
    

    #Calculate kinetic thermal conductivity
    lambda_k = np.zeros(len(T))
    for j in range(len(T)):
        # if Compound == ['NaCl','UCl3'] and mol_fracs == [0.63,0.37]:
        #     lambda_k[j] = 1.240955586 * k_B * n_dens_mix**(2/3) * C_0_mix_m * (1 - alpha_mix_m*(gamma_mix_m + 1/3)*(T[j] - T_melt))
        # elif Compound == ['NaCl','UCl3'] and mol_fracs == [0.658,0.342]:
        #     lambda_k[j] = 1.258369925 * k_B * n_dens_mix**(2/3) * C_0_mix_m * (1 - alpha_mix_m*(gamma_mix_m + 1/3)*(T[j] - T_melt))
        # else:
        lambda_k[j] = (1 + np.sum(n_i_c)/np.sum(n_i_a)) * k_B * n_dens_mix**(2/3) * C_0_mix_m * (1 - alpha_mix_m*(gamma_mix_m + 1/3)*(T[j] - T_melt))
    
    
    #Calculate ideal thermal conductivity
    lambda_ideal = np.zeros(len(T))
    for i in range(len(compound_input)):
        for j in range(len(T)):
            lambda_ideal[j] += lambda_i[i][j]*mol_fracs[i]


    # Calculate delta of mix
    delta = np.zeros(len(T))
    for i in range(len(compound_input)):
        for j in range(len(T)):
            # if Compound == ['NaCl','UCl3'] and mol_fracs == [0.63,0.37]:
            #     delta[j] += (0.193646716)/M_mix * lambda_i[i][j]/lambda_ideal[j] * mol_fracs[i] * (1 - M_i[i]/M_mix)**2
            # elif Compound == ['NaCl','UCl3'] and mol_fracs == [0.658,0.342]:
            #     delta[j] += (0.171475397)/M_mix * lambda_i[i][j]/lambda_ideal[j] * mol_fracs[i] * (1 - M_i[i]/M_mix)**2
            # else:
            delta[j] += lambda_i[i][j]/lambda_ideal[j] * mol_fracs[i] * (1 - M_i[i]/M_mix)**2
    

    lambda_mix_T = lambda_k * (1 - delta)


    print("k: ",lambda_mix_T[0])
    # print("lambda_BC: ",lambda_mix_T[0]/(1/3*C_p_mix*C_0_mix_m))
    print("cp: ",C_p_mix/V_mix_m)
    print("vs: ",C_0_mix_m)
    print("lambda: ", 1 + np.sum(n_i_c)/np.sum(n_i_a))
    return(T,lambda_mix_T)

def Kinetic_Theory24_Mix(df,MSTDB_df,SCL_PDF_df,compound_input,mol_fracs,Temp_Range, V_m=0, density_mix=0, C_p_mix=0, alpha=0, expon=0):
    
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

    print('')
    print('################ Kinetic Theory 2024, Mix Data ############################')
    specific_heat_mix = C_p_mix
    density_mix = [0,0]
    sound_velocity_mix = [0,0]    
    molar_volume_mix = [0,0] 
    
    # Obtain measurement data for function
    func_name = inspect.currentframe().f_code.co_name
    r_ac_mix, density_mix, sound_velocity_mix, specific_heat_mix = prop_lookup(func_name,Compound,mol_fracs)

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

    # print("KT, number of ions: ",n_i)
    # print("KT, number of cations: ",n_i_c)
    # print("KT, number of anions: ",n_i_a)
        
    # Calculate compound psi term
    psi_i=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        psi_i[i] = 1 + n_i_c[i]/n_i_a[i]

    # Calculate compound number density
    n_dens_i = np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        n_dens_i[i] = Avog *n_i[i]/V_i_m[i]

    #Calculate compound thermal conductivity(T) based on compound minimum thermal conductivites
    lambda_i = np.zeros((len(Compound),len(T)))
    for i in range(len(compound_input)):
        for j in range(len(T)):
            lambda_i[i][j] = (1 + n_i_c[i]/n_i_a[i]) * k_B * n_dens_i[i]**(2/3) * C_i_0[i] * (1 - alpha_i_m[i]*(gamma_i_m[i] + 1/3)*(T[j] - T_i_m[i]))

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
        kappa_i_m[i]=M_i[i]*mol_fracs[i]
        denom=denom+M_i[i]*mol_fracs[i]
    kappa_i_m=kappa_i_m/denom 

    #Calculate taus for compounds
    Tau_i_m=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        Tau_i_m[i] = ( C_i_v[i]*C_i_0[i] ) / n_i[i]


    #Calculate mixture speed of sound
    C_0_mix_m=0
    C_0_mix_est = 0
    C_0_mix_data = 0 
    for i in range(len(compound_input)):
        C_0_mix_est += phi_i_m[i]**2/(kappa_i_m[i]*C_i_0[i]**2)
    C_0_mix_data = sound_velocity_mix[1]*T_melt+sound_velocity_mix[0]
    C_0_mix_est=1/np.sqrt(C_0_mix_est)
    if sound_velocity_mix == [0,0]:
        C_0_mix_m = C_0_mix_est
    else:
        C_0_mix_m = C_0_mix_data

        C_0_mix_est_avg = np.average(C_0_mix_est)
        C_0_mix_data_avg = np.average(C_0_mix_data)
        error = 100 * (C_0_mix_data_avg-C_0_mix_est_avg)/C_0_mix_data_avg
        print("Sound Velocity - Estimated: ", C_0_mix_est_avg)
        print("Sound Velocity - Mix Data: ", C_0_mix_data_avg)
        print("Sound Velocity % Difference: ", error, " %")

    #Calculate mixture molecular weight
    M_mix=0
    for i in range(len(compound_input)):
        M_mix=M_mix+mol_fracs[i]*M_i[i]


    #Calculate average molar volume
    V_mix_m=0
    V_mix_est = 0
    V_mix_data = 0
    V_mix_data_Vm = 0
    for i in range(len(compound_input)):    
        V_mix_est += mol_fracs[i]*M_i[i]*(1/rho_i_m[i])
    V_mix_data = M_mix/(density_mix[0] + density_mix[1]*T_melt)
    V_mix_data_Vm = molar_volume_mix[1]*T_melt+molar_volume_mix[0]
    if molar_volume_mix == [0,0]:
        if density_mix == [0,0]:
            V_mix_m = V_mix_est
        else:
            V_mix_m = V_mix_data
    else:
        V_mix_m = V_mix_data_Vm

        V_mix_est_avg = np.average(V_mix_est)
        V_mix_data_avg = np.average(V_mix_data)
        error_densityVm = 100 * (V_mix_data_avg-V_mix_est_avg)/V_mix_data_avg
        error_Vm = 100 * (V_mix_data_Vm-V_mix_est_avg)/V_mix_data_Vm
        #print("")
        #print("Molar volume - From estimated density: ", V_mix_est_avg)
        #print("Molar volume - From mix data density: ", V_mix_data_avg)
        #print("Molar volume - Mix data Molar Volume: ", V_mix_data_Vm)
        #print("Molar volume % Difference, Density Data/Estimation: ", error_densityVm, " %")
        #print("Molar volume % Difference, Molar Volume Data/Estimation: ", error_Vm, " %")
    
    #Calculate mixture specific heat
    #Calculate average number of ions
    C_p_mix=0
    C_p_mix_est = 0
    C_p_mix_data = 0
    n_mix=0
    for i in range(len(compound_input)):
        C_p_mix_est += C_p_mix+mol_fracs[i]*C_i_p[i]
        n_mix = n_mix + mol_fracs[i]*n_i[i]
    C_p_mix_data = specific_heat_mix
    if specific_heat_mix[0] == 0:
        C_p_mix = C_p_mix_est
    else:
        if specific_heat_mix[2] == 'm':
            C_p_mix = C_p_mix_data[1]*T_melt+C_p_mix_data[0]
        elif specific_heat_mix[2] == 'g':   # Converts units from J/g/K to J/mol/K
            C_p_mix = (C_p_mix_data[0]*T_melt+C_p_mix_data[1])*(M_mix)

        C_p_mix_est_avg = np.average(C_p_mix_est)
        C_p_mix_data_avg = np.average(C_p_mix)
        error = 100 * (C_p_mix_data_avg-C_p_mix_est_avg)/C_p_mix_data_avg
        #print("")
        print("Specific Heat - Estimated: ", C_p_mix_est_avg)
        print("Specific Heat - Mix Data: ", C_p_mix_data_avg)
        print("Specific Heat % Difference: ", error, " %")
        
    #Calculate mixture number density
    n_dens_mix = Avog * n_mix/V_mix_m


    #Calculate mixture thermal expansion
    alpha_mix_m=0
    if alpha==0:
        for i in range(len(compound_input)):
            alpha_mix_m=alpha_mix_m+phi_i_m[i]*alpha_i_m[i]
    else:alpha_mix_m=alpha
        

    #Calculate mixture gamma
    gamma_mix_m=M_mix*alpha_mix_m*C_0_mix_m**2*(1/C_p_mix)
    

    #Calculate kinetic thermal conductivity
    lambda_k = np.zeros(len(T))
    for j in range(len(T)):
        lambda_k[j] = (1 + np.sum(n_i_c)/np.sum(n_i_a)) * k_B * n_dens_mix**(2/3) * C_0_mix_m * (1 - alpha_mix_m*(gamma_mix_m + 1/3)*(T[j] - T_melt))
    
    
    #Calculate ideal thermal conductivity
    lambda_ideal = np.zeros(len(T))
    for i in range(len(compound_input)):
        for j in range(len(T)):
            lambda_ideal[j] += lambda_i[i][j]*mol_fracs[i]


    # Calculate delta of mix
    delta = np.zeros(len(T))
    for i in range(len(compound_input)):
        for j in range(len(T)):
            delta[j] += lambda_i[i][j]/lambda_ideal[j] * mol_fracs[i] * (1 - M_i[i]/M_mix)**2
    

    lambda_mix_T = lambda_k * (1 - delta)


    print("k: ",lambda_mix_T[0])
    # print("lambda_BC: ",lambda_mix_T[0]/(1/3*C_p_mix*C_0_mix_m))
    print("cp: ",C_p_mix/V_mix_m)
    print("vs: ",C_0_mix_m)
    print("lambda: ", 1 + np.sum(n_i_c)/np.sum(n_i_a))
    return(T,lambda_mix_T)

def Zhao_PGM(df,MSTDB_df,SCL_PDF_df,compound_input,mol_fracs,Temp_Range, V_m=0, density_mix=0, sound_velocity_mix=0, C_p_mix=0, alpha=0, expon=0):
    print("")
    print("Zhao_PGM ######################################")
    print(compound_input)
    #Define smaller dataframes of each variable contained in the spreadsheet
    Compound_df=df['Compound']
    M_i_df=df['M_i (g/mol)']
    T_i_m_df=df['T_i_m (K)']
    V_i_m_df=df['V_i_m (m^3/mol)']
    C_i_p_df=df['C_i_p (J/mol/K)']
    C_i_p_sp_df=df['C_i_p_sp (J/g/K)']
    r_c_df=df['r_c (m)']
    r_a_df=df['r_a (m)']
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
    C_i_p=np.zeros(len(compound_input))
    C_i_p_sp=np.zeros(len(compound_input))
    r_c=np.zeros(len(compound_input))
    r_a=np.zeros(len(compound_input))
    C_0_T_A=np.zeros(len(compound_input))
    C_0_T_B=np.zeros(len(compound_input))
    rho_T_A=np.zeros(len(compound_input))
    rho_T_B=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        Compound[i]=Compound_df[indices[i]]
        M_i[i]=M_i_df[indices[i]]*0.001         # Convert to kg/mol
        T_i_m[i]=T_i_m_df[indices[i]]
        V_i_m[i]=V_i_m_df[indices[i]]
        C_i_p[i]=C_i_p_df[indices[i]]
        C_i_p_sp[i]=C_i_p_sp_df[indices[i]]/0.001   # Convert to J/kg/K
        r_c[i]=r_c_df[indices[i]]
        r_a[i]=r_a_df[indices[i]]
        C_0_T_A[i]=C_0_T_A_df[indices[i]]
        C_0_T_B[i]=C_0_T_B_df[indices[i]]
        rho_T_A[i]=rho_T_A_df[indices[i]]/0.001
        rho_T_B[i]=rho_T_B_df[indices[i]]/0.001

    # print('')
    # print('####### Phonon Gas Model ############################')

    specific_heat_mix = C_p_mix


    #initialize a temperature range
    T_melt = Temp_Range[0]
    T=np.linspace(T_melt,Temp_Range[1],num=100)

    #Calculate number of ions per compound
    n_i=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        for c in Compound[i]:
            if c.isupper():
                n_i[i] = n_i[i] + 1
            elif c.isnumeric():
                c = int(c)
                n_i[i] = n_i[i] + (c-1)
    
    #print("PGM, number of ions: ",n_i)


    # if Compound == ['NaCl','UCl3'] and mol_fracs == [0.63,0.37]:
    #     density_mix = [4.22*100**3/1000,-0.00113*100**3/1000]   # Desyatnik, 1975
    # elif Compound == ['LiF','BeF2']:
    #     density_mix = [2518.15,-0.424]  # Vidrio, 2022, FLiBe 33.59mol%BeF2
    #     sound_velocity_mix = [4272.309646,-1.101929884] # 66-34% Cantor, 1968
    #     specific_heat_mix = 2414.7*0.033  # 67-33% Sohal, 2010  
    #     #density_mix = [2413,-0.488]  # Janz, 1974, FLiBe 33mol%BeF2
    #     #2110 * (1.885 + 2.762*mol_fracs[1] + mol_fracs[1]**2) / (1.773 + 2.663*mol_fracs[0] )
    # elif Compound == ['LiF','NaF','KF']:
    #     density_mix = [2729.3,-0.73]  # 46.5-11.5-42%, Vriesema [1979], Ingersoll et al. [2007], and Williams et al. [2006]
    #     sound_velocity_mix = [3241.15,-1.20]    # 46.5-11.5-42% Robertson, 2022
    #     specific_heat_mix = 1882.8*0.0413  # 46.5-11.5-42% Sohal, 2010
    # else:
    density_mix = [0,0]
    sound_velocity_mix = [0,0]    
    specific_heat_mix = [0,0]  

    #Find the temperature dependent density from data
    rho_i = np.zeros((len(Compound),len(T)))
    for i in range(len(Compound)):
        if pd.isna(rho_T_A[i]):     
            rho_T_A[i], rho_T_B[i] = density(df,compound_input[i])
        else:
            pass
        for j in range(len(T)):
            rho_i[i][j] = rho_T_A[i] + rho_T_B[i]*T[j] 

    
    # Calculate the temp-dependent sound velocity of compounds
    C_0_i = np.zeros((len(Compound),len(T)))
    for i in range(len(compound_input)):
        for j in range(len(T)):
            C_0_i[i][j] = C_0_T_B[i]*T[j]+C_0_T_A[i]


    # Calculate molar volume of compounds from temp-dependent density
    V_i = np.zeros((len(Compound),len(T)))
    for i in range(len(compound_input)):
        for j in range(len(T)):
            V_i[i][j] = M_i[i] / rho_i[i][j]

    
    # Calculate volume fractions
    phi_i=np.zeros((len(Compound),len(T)))
    for j in range(len(T)):
        denom_phi = 0
        for i in range(len(compound_input)):
            phi_i[i][j] = V_i[i][j]*mol_fracs[i]
            denom_phi += V_i[i][j]*mol_fracs[i]
        for i in range(len(compound_input)):
            phi_i[i][j]=phi_i[i][j]/denom_phi

    # Calculate mass (fractions
    kappa_i=np.zeros(len(Compound))
    denom_kappa = 0
    for i in range(len(compound_input)):
        kappa_i[i] = M_i[i]*mol_fracs[i]
        denom_kappa += M_i[i]*mol_fracs[i]
    kappa_i=kappa_i/denom_kappa           


    #Calculate mixture molecular weight
    M_mix=0
    for i in range(len(compound_input)):
        M_mix=M_mix+mol_fracs[i]*M_i[i]


    # Calculate the average temp-dependent molar volume of mixture
    V_mix = np.zeros(len(T))      
    V_mix_est = np.zeros(len(T)) 
    V_mix_data = np.zeros(len(T)) 
    for i in range(len(compound_input)):
        for j in range(len(T)):
            V_mix_est[j] += mol_fracs[i]*V_i[i][j]
            if density_mix == [0,0]:
                pass
            else:
                V_mix_data[j] = M_mix/(density_mix[0] + density_mix[1]*T[j])
    if density_mix == [0,0]:
        V_mix = V_mix_est
    else:
        V_mix = V_mix_data

        V_mix_est_avg = np.average(V_mix_est)
        V_mix_data_avg = np.average(V_mix_data)
        error = 100 * (V_mix_data_avg-V_mix_est_avg)/V_mix_data_avg
        #print("Molar volume from estimated mix density: ", V_mix_est_avg)
        #print("Molar volume from data mix density: ", V_mix_data_avg)
        #print("Molar volume % Difference: ", error, " %")


    #Calculate the temp-dependent density of mixture
    rho_mix = np.zeros(len(T))
    rho_mix_est = np.zeros(len(T)) 
    rho_mix_data = np.zeros(len(T)) 
    for i in range(len(compound_input)):
        for j in range(len(T)):
            rho_mix_est[j] += rho_i[i][j]*mol_fracs[i]   #M_mix/V_mix[j]
            rho_mix_data[j] = density_mix[0] + density_mix[1]*T[j]   # If mixture temp-dependent density data exists
    if density_mix == [0,0]:
        rho_mix = rho_mix_est
    else:
        rho_mix = rho_mix_data

        rho_mix_est_avg = np.average(rho_mix_est)
        rho_mix_data_avg = np.average(rho_mix_data)
        error = 100 * (rho_mix_data_avg-rho_mix_est_avg)/rho_mix_data_avg
        #print("Estimated mix density: ", rho_mix_est_avg)
        #print("Data mix density: ", rho_mix_data_avg)
        print("Density % Difference: ", error, " %")


    # Calculate the temp-dependent sound velocity of mixture
    C_0_mix = np.zeros(len(T))
    C_0_mix_est = np.zeros(len(T)) 
    C_0_mix_data = np.zeros(len(T)) 
    for i in range(len(compound_input)):
        for j in range(len(T)):
            C_0_mix_est[j] += phi_i[i][j]**2 / (kappa_i[i] * C_0_i[i][j]**2)
            C_0_mix_data[j] = sound_velocity_mix[1]*T[j]+sound_velocity_mix[0]    # If mixture temp-dependent sound velocity data exists
    for j in range(len(T)):
        C_0_mix_est[j] = 1 / np.sqrt(C_0_mix_est[j])
    if sound_velocity_mix == [0,0]:
        C_0_mix = C_0_mix_est
    else:
        C_0_mix = C_0_mix_data

        C_0_mix_est_avg = np.average(C_0_mix_est)
        C_0_mix_data_avg = np.average(C_0_mix_data)
        error = 100 * (C_0_mix_data_avg-C_0_mix_est_avg)/C_0_mix_data_avg
        #print("Estimated mix sound velocity: ", C_0_mix_est_avg)
        #print("Data mix sound velocity: ", C_0_mix_data_avg)
        print("sound velocity % Difference: ", error, " %")


    # Calculate mixture specific heat & # of ions
    C_p_mix=0
    C_p_mix_est = 0
    C_p_mix_data = 0
    n_mix=0
    for i in range(len(compound_input)):
        C_p_mix_est += C_p_mix+mol_fracs[i]*C_i_p[i]
        C_p_mix_data=specific_heat_mix
        n_mix=n_mix+mol_fracs[i]*n_i[i]
    if specific_heat_mix[0] == 0:
        C_p_mix = C_p_mix_est
    else:
        C_p_mix = C_p_mix_data

        C_p_mix_est_avg = np.average(C_p_mix_est)
        C_p_mix_data_avg = np.average(C_p_mix)
        error = 100 * (C_p_mix_data_avg-C_p_mix_est_avg)/C_p_mix_data_avg
        #print("Estimated mix specific heat: ", C_p_mix_est_avg)
        #print("Data mix specific heat: ", C_p_mix_data_avg)
        print("specific heat % Difference: ", error, " %")



    # Calculate mixture radial distance
    r_ac_mix = 0
    r_ac_mix1 = 0
    for i in range(len(compound_input)):
        r_ac_mix=r_ac_mix1+mol_fracs[i]*(r_a[i]+r_c[i])

    
    print("Calculated mixture specific heat (melting point): ",C_p_mix/V_mix[0])
    print("Calculated mixtur sound velocity (melting point): ",C_0_mix[0])
    print("Calculated average mean free path: ",r_ac_mix)

    #Calculate compound thermal conductivity
    lambda_i_m=np.zeros(len(T))
    lambda_i_mg=np.zeros(len(T))
    lambda_i_mb=np.zeros(len(T))
    sound_w_time = 0
    for j in range(len(T)):
        lambda_i_m[j] = 1/3 * C_i_p_sp[i] * rho_mix[j] * C_0_mix[j] * r_ac_mix    # Verified with Zhao's results, uses Zhao's data and specific heat capacity
        #lambda_i_mg[j] = 1/3 * C_p_mix * 1/M_mix * rho_mix[j] * C_0_mix[j] * r_ac_mix  # Verified with Zhao's results, uses Zhao's data and but MSTDB heat capacity
        lambda_i_mb[j] = 1/3 * C_p_mix * 1/V_mix[j] * C_0_mix[j] * r_ac_mix  # Uses molar volume at melting point only 
    

    # if sound_w_time == 1:
    #     print("Calculated with temp-dependent sound velocity data.")
    # else:
    #     print("No temp-dependent sound velocity data available. Calculated with melting temp sound velocity only.")
    

    nan_check = np.isnan(lambda_i_m)
    contains_nan = nan_check.any()
    if contains_nan or lambda_i_m[0] == 0:
        print("k: ",lambda_i_mb[0])
        # print("lambda_BC: ",lambda_i_mb[0]/(1/3*C_p_mix*C_0_mix))
        print("cp: ",C_p_mix)
        print("vs: ",C_0_mix[0])
        return(T,lambda_i_mb)
    else:
        print("k: ",lambda_i_m[0])
        # print("lambda_BC: ",lambda_i_m[0]/(1/3*C_p_mix*C_0_mix))
        print("cp: ",C_p_mix)
        print("vs: ",C_0_mix[0])
        return(T,lambda_i_m)

def GECM(df,MSTDB_df,SCL_PDF_df,compound_input,mol_fracs,Temp_Range, V_m=0, density_mix=0, sound_velocity_mix=0, C_p_mix=0, alpha=0, expon=0):
    print("")
    print("GECM ######################################")
    print(compound_input)

    #Define smaller dataframes of each variable contained in the spreadsheet
    Compound_df=df['Compound']
    M_i_df=df['M_i (g/mol)']
    T_i_m_df=df['T_i_m (K)']
    V_i_m_df=df['V_i_m (m^3/mol)']
    C_i_p_df=df['C_i_p (J/mol/K)']
    C_i_p_sp_df=df['C_i_p_sp (J/g/K)']
    r_c_df=df['r_c (m)']
    r_a_df=df['r_a (m)']
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
    C_i_p=np.zeros(len(compound_input))
    C_i_p_sp=np.zeros(len(compound_input))
    r_c=np.zeros(len(compound_input))
    r_a=np.zeros(len(compound_input))
    C_0_T_A=np.zeros(len(compound_input))
    C_0_T_B=np.zeros(len(compound_input))
    rho_T_A=np.zeros(len(compound_input))
    rho_T_B=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        Compound[i]=Compound_df[indices[i]]
        M_i[i]=M_i_df[indices[i]]*0.001         # Convert to kg/mol
        T_i_m[i]=T_i_m_df[indices[i]]
        V_i_m[i]=V_i_m_df[indices[i]]
        C_i_p[i]=C_i_p_df[indices[i]]
        C_i_p_sp[i]=C_i_p_sp_df[indices[i]]/0.001   # Convert to J/kg/K
        r_c[i]=r_c_df[indices[i]]
        r_a[i]=r_a_df[indices[i]]
        C_0_T_A[i]=C_0_T_A_df[indices[i]]
        C_0_T_B[i]=C_0_T_B_df[indices[i]]
        rho_T_A[i]=rho_T_A_df[indices[i]]/0.001
        rho_T_B[i]=rho_T_B_df[indices[i]]/0.001

    # print('')
    # print('####### Phonon Gas Model ############################')

    specific_heat_mix = C_p_mix


    #initialize a temperature range
    T_melt = Temp_Range[0]
    T=np.linspace(T_melt,Temp_Range[1],num=100)

    #Calculate number of ions per compound
    n_i=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        for c in Compound[i]:
            if c.isupper():
                n_i[i] = n_i[i] + 1
            elif c.isnumeric():
                c = int(c)
                n_i[i] = n_i[i] + (c-1)
    
    #print("PGM, number of ions: ",n_i)
    density_mix = [0,0]
    sound_velocity_mix = [0,0]
    specific_heat_mix = [0,0]

    # Obtain measurement data for function
    func_name = inspect.currentframe().f_code.co_name
    r_ac_mix, density_mix, sound_velocity_mix, specific_heat_mix = prop_lookup(func_name,Compound,mol_fracs)
    print("SCL: ",r_ac_mix)
    
    #Find the temperature dependent density from data
    rho_i = np.zeros((len(Compound),len(T)))
    for i in range(len(Compound)):
        if pd.isna(rho_T_A[i]):     
            rho_T_A[i], rho_T_B[i] = density(df,compound_input[i])
        else:
            pass
        for j in range(len(T)):
            rho_i[i][j] = rho_T_A[i] + rho_T_B[i]*T[j] 

    
    # Calculate the temp-dependent sound velocity of compounds
    C_0_i = np.zeros((len(Compound),len(T)))
    for i in range(len(compound_input)):
        for j in range(len(T)):
            C_0_i[i][j] = C_0_T_B[i]*T[j]+C_0_T_A[i]


    # Calculate molar volume of compounds from temp-dependent density
    V_i = np.zeros((len(Compound),len(T)))
    for i in range(len(compound_input)):
        for j in range(len(T)):
            V_i[i][j] = M_i[i] / rho_i[i][j]

    
    # Calculate volume fractions
    phi_i=np.zeros((len(Compound),len(T)))
    for j in range(len(T)):
        denom_phi = 0
        for i in range(len(compound_input)):
            phi_i[i][j] = V_i[i][j]*mol_fracs[i]
            denom_phi += V_i[i][j]*mol_fracs[i]
        for i in range(len(compound_input)):
            phi_i[i][j]=phi_i[i][j]/denom_phi

    # Calculate mass (fractions
    kappa_i=np.zeros(len(Compound))
    denom_kappa = 0
    for i in range(len(compound_input)):
        kappa_i[i] = M_i[i]*mol_fracs[i]
        denom_kappa += M_i[i]*mol_fracs[i]
    kappa_i=kappa_i/denom_kappa           


    #Calculate mixture molecular weight
    M_mix=0
    for i in range(len(compound_input)):
        M_mix=M_mix+mol_fracs[i]*M_i[i]


    # Calculate the average temp-dependent molar volume of mixture
    V_mix = np.zeros(len(T))      
    V_mix_est = np.zeros(len(T)) 
    V_mix_data = np.zeros(len(T)) 
    for i in range(len(compound_input)):
        for j in range(len(T)):
            V_mix_est[j] += mol_fracs[i]*V_i[i][j]
            if density_mix == [0,0]:
                pass
            else:
                V_mix_data[j] = M_mix/(density_mix[0] + density_mix[1]*T[j])
    if density_mix == [0,0]:
        V_mix = V_mix_est
    else:
        V_mix = V_mix_data

        V_mix_est_avg = np.average(V_mix_est)
        V_mix_data_avg = np.average(V_mix_data)
        error = 100 * (V_mix_data_avg-V_mix_est_avg)/V_mix_data_avg
        #print("Molar volume from estimated mix density: ", V_mix_est_avg)
        #print("Molar volume from data mix density: ", V_mix_data_avg)
        #print("Molar volume % Difference: ", error, " %")


    #Calculate the temp-dependent density of mixture
    rho_mix = np.zeros(len(T))
    rho_mix_est = np.zeros(len(T)) 
    rho_mix_data = np.zeros(len(T)) 
    for i in range(len(compound_input)):
        for j in range(len(T)):
            rho_mix_est[j] += rho_i[i][j]*mol_fracs[i]   #M_mix/V_mix[j]
            rho_mix_data[j] = density_mix[0] + density_mix[1]*T[j]   # If mixture temp-dependent density data exists
    if density_mix == [0,0]:
        rho_mix = rho_mix_est
    else:
        rho_mix = rho_mix_data

        rho_mix_est_avg = np.average(rho_mix_est)
        rho_mix_data_avg = np.average(rho_mix_data)
        error = 100 * (rho_mix_data_avg-rho_mix_est_avg)/rho_mix_data_avg
        #print("Estimated mix density: ", rho_mix_est_avg)
        #print("Data mix density: ", rho_mix_data_avg)
        print("Density % Difference: ", error, " %")


    # Calculate the temp-dependent sound velocity of mixture
    C_0_mix = np.zeros(len(T))
    C_0_mix_est = np.zeros(len(T)) 
    C_0_mix_data = np.zeros(len(T)) 
    for i in range(len(compound_input)):
        for j in range(len(T)):
            C_0_mix_est[j] += phi_i[i][j]**2 / (kappa_i[i] * C_0_i[i][j]**2)
            C_0_mix_data[j] = sound_velocity_mix[1]*T[j]+sound_velocity_mix[0]    # If mixture temp-dependent sound velocity data exists
    for j in range(len(T)):
        C_0_mix_est[j] = 1 / np.sqrt(C_0_mix_est[j])
    if sound_velocity_mix == [0,0]:
        C_0_mix = C_0_mix_est
    else:
        C_0_mix = C_0_mix_data

        C_0_mix_est_avg = np.average(C_0_mix_est)
        C_0_mix_data_avg = np.average(C_0_mix_data)
        error = 100 * (C_0_mix_data_avg-C_0_mix_est_avg)/C_0_mix_data_avg
        #print("Estimated mix sound velocity: ", C_0_mix_est_avg)
        #print("Data mix sound velocity: ", C_0_mix_data_avg)
        print("sound velocity % Difference: ", error, " %")


    # Calculate mixture specific heat & # of ions
    C_p_mix=0
    C_p_mix_est = 0
    C_p_mix_data = 0
    n_mix=0
    for i in range(len(compound_input)):
        C_p_mix_est += C_p_mix+mol_fracs[i]*C_i_p[i]
        C_p_mix_data=specific_heat_mix
        n_mix=n_mix+mol_fracs[i]*n_i[i]
    if specific_heat_mix[0] == 0:
        C_p_mix = C_p_mix_est
    else:
        C_p_mix = C_p_mix_data

        C_p_mix_est_avg = np.average(C_p_mix_est)
        C_p_mix_data_avg = np.average(C_p_mix_data)
        error = 100 * (C_p_mix_data_avg-C_p_mix_est_avg)/C_p_mix_data_avg
        #print("Estimated mix specific heat: ", C_p_mix_est_avg)
        #print("Data mix specific heat: ", C_p_mix_data_avg)
        print("specific heat % Difference: ", error, " %")



    # Calculate mixture radial distance
    r_ac_mix_r_sum = 0
    # Calculate mixture radial distance
    r_ac_mix1 = 0
    for i in range(len(compound_input)):
        r_ac_mix_r_sum=r_ac_mix1+mol_fracs[i]*(r_a[i]+r_c[i]) 
    #print("PDF CatAn peak MFP: ", r_ac_mix)
    print("Avg sum radii: ", r_ac_mix_r_sum)
    #print("MFP % Difference: ", error, " %")


    # Find matching composition in SCL_PDF_df
    matching_row = find_SCL_composition(SCL_PDF_df, compound_input, mol_fracs)

    if matching_row is not None:
        #print("Matching composition found in SCL_PDF_df")
        # Extract data from the matching row
        scl_data = matching_row.iloc[0]
        
        # Extract pair data
        pairs = []
        for i in range(1, 6):  # Assuming up to 5 pairs
            pair_prefix = f'Pair{i}'
            if f'{pair_prefix}_x-SCL' in scl_data:
                pair = {
                    'name': scl_data[f'{pair_prefix}'],
                    'x_scl': scl_data[f'{pair_prefix}_x-SCL'],
                    'peak_x': scl_data[f'{pair_prefix}_Peak-x'],
                    'peak_y': scl_data[f'{pair_prefix}_Peak-y'],
                    'min_x': scl_data[f'{pair_prefix}_Min-x'],
                    'min_y': scl_data[f'{pair_prefix}_Min-y']
                }
                pairs.append(pair)

        # Use the extracted data in your calculations
        # For example, you can use the first pair's x_scl as r_ac_mix:
        r_ac_mix = pairs[0]['x_scl'] if pairs else 0
    else:
        print("No matching composition found in SCL_PDF_df")
        r_ac_mix = r_ac_mix_r_sum

    
    # print("Calculated mixture specific heat (melting point): ",C_p_mix/V_mix[0])
    # print("Calculated mixtur sound velocity (melting point): ",C_0_mix[0])
    # print("Calculated average mean free path: ",r_ac_mix)

    #Calculate compound thermal conductivity
    lambda_i_m=np.zeros(len(T))
    lambda_i_mg=np.zeros(len(T))
    lambda_i_mb=np.zeros(len(T))
    sound_w_time = 0
    for j in range(len(T)):
        lambda_i_m[j] = 1/3 * C_i_p_sp[i] * rho_mix[j] * C_0_mix[j] * r_ac_mix    # Verified with Zhao's results, uses Zhao's data and specific heat capacity
        #lambda_i_mg[j] = 1/3 * C_p_mix * 1/M_mix * rho_mix[j] * C_0_mix[j] * r_ac_mix  # Verified with Zhao's results, uses Zhao's data and but MSTDB heat capacity
        lambda_i_mb[j] = 1/3 * C_p_mix * 1/V_mix[j] * C_0_mix[j] * r_ac_mix  # Uses molar volume at melting point only 
    

    # if sound_w_time == 1:
    #     print("Calculated with temp-dependent sound velocity data.")
    # else:
    #     print("No temp-dependent sound velocity data available. Calculated with melting temp sound velocity only.")
    

    nan_check = np.isnan(lambda_i_m)
    contains_nan = nan_check.any()
    if contains_nan or lambda_i_m[0] == 0:
        print("k: ",lambda_i_mb[0])
        # print("PGM-PDF_Avg (1/3*Cp*vs): ",lambda_i_mb[0]/r_ac_mix)
        # print("PGM-PDF_Avg (1/3*Cp*MFP): ",lambda_i_mb[0]/C_0_mix[0])
        # print("PGM-PDF_Avg (1/3*vs*MFP): ",lambda_i_mb[0]/(C_p_mix/V_mix[0]))
        # print("lambda_BC: ",lambda_i_mb[0]/(1/3*C_p_mix*C_0_mix[0]))
        # print("lambda_SCL: ",r_ac_mix)
        print("cp: ",C_p_mix/V_mix[0])
        print("vs: ",C_0_mix[0])
        return(T,lambda_i_mb)
    else:
        print("k: ",lambda_i_m[0])
        # print("PGM-PDF_Avg (1/3*Cp*vs): ",lambda_i_mb[0]/r_ac_mix)
        # print("PGM-PDF_Avg (1/3*Cp*MFP): ",lambda_i_mb[0]/C_0_mix[0])
        # print("PGM-PDF_Avg (1/3*vs*MFP): ",lambda_i_mb[0]/C_i_p_sp[0])
        # print("lambda_BC: ",lambda_i_m[0]/(1/3*C_p_mix*C_0_mix[0]))
        # print("lambda_SCL: ",r_ac_mix)
        print("cp: ",C_i_p_sp[0]* rho_mix[0])
        print("vs: ",C_0_mix[0])
        print("lambda: ", r_ac_mix)
        return(T,lambda_i_m)

        return(T,lambda_i_m)


def GECM_Mix(df,MSTDB_df,SCL_PDF_df,compound_input,mol_fracs,Temp_Range, V_m=0, density_mix=0, sound_velocity_mix=0, C_p_mix=0, alpha=0, expon=0):
    print(" ")
    print("GECM_Mix ######################################")
    print(compound_input)
    #Define smaller dataframes of each variable contained in the spreadsheet
    Compound_df=df['Compound']
    M_i_df=df['M_i (g/mol)']
    T_i_m_df=df['T_i_m (K)']
    V_i_m_df=df['V_i_m (m^3/mol)']
    C_i_p_df=df['C_i_p (J/mol/K)']
    C_i_p_sp_df=df['C_i_p_sp (J/g/K)']
    r_c_df=df['r_c (m)']
    r_a_df=df['r_a (m)']
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
    C_i_p=np.zeros(len(compound_input))
    C_i_p_sp=np.zeros(len(compound_input))
    r_c=np.zeros(len(compound_input))
    r_a=np.zeros(len(compound_input))
    C_0_T_A=np.zeros(len(compound_input))
    C_0_T_B=np.zeros(len(compound_input))
    rho_T_A=np.zeros(len(compound_input))
    rho_T_B=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        Compound[i]=Compound_df[indices[i]]
        M_i[i]=M_i_df[indices[i]]*0.001         # Convert to kg/mol
        T_i_m[i]=T_i_m_df[indices[i]]
        V_i_m[i]=V_i_m_df[indices[i]]
        C_i_p[i]=C_i_p_df[indices[i]]
        C_i_p_sp[i]=C_i_p_sp_df[indices[i]]/0.001   # Convert to J/kg/K
        r_c[i]=r_c_df[indices[i]]
        r_a[i]=r_a_df[indices[i]]
        C_0_T_A[i]=C_0_T_A_df[indices[i]]
        C_0_T_B[i]=C_0_T_B_df[indices[i]]
        rho_T_A[i]=rho_T_A_df[indices[i]]/0.001
        rho_T_B[i]=rho_T_B_df[indices[i]]/0.001

    # print('')
    # print('####### Phonon Gas Model ############################')

    specific_heat_mix = C_p_mix


    #initialize a temperature range
    T_melt = Temp_Range[0]
    T=np.linspace(T_melt,Temp_Range[1],num=100)

    #Calculate number of ions per compound
    n_i=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        for c in Compound[i]:
            if c.isupper():
                n_i[i] = n_i[i] + 1
            elif c.isnumeric():
                c = int(c)
                n_i[i] = n_i[i] + (c-1)
    
    #print("PGM, number of ions: ",n_i)
    density_mix = [0,0]
    sound_velocity_mix = [0,0]
    specific_heat_mix = [0,0]

    # Obtain measurement data for function
    func_name = inspect.currentframe().f_code.co_name
    r_ac_mix, density_mix, sound_velocity_mix, specific_heat_mix = prop_lookup(func_name,Compound,mol_fracs)
    print("SCL: ",r_ac_mix)

    #Find the temperature dependent density from data
    rho_i = np.zeros((len(Compound),len(T)))
    for i in range(len(Compound)):
        if pd.isna(rho_T_A[i]):     
            rho_T_A[i], rho_T_B[i] = density(df,compound_input[i])
        else:
            pass
        for j in range(len(T)):
            rho_i[i][j] = rho_T_A[i] + rho_T_B[i]*T[j] 

    
    # Calculate the temp-dependent sound velocity of compounds
    C_0_i = np.zeros((len(Compound),len(T)))
    for i in range(len(compound_input)):
        for j in range(len(T)):
            C_0_i[i][j] = C_0_T_B[i]*T[j]+C_0_T_A[i]


    # Calculate molar volume of compounds from temp-dependent density
    V_i = np.zeros((len(Compound),len(T)))
    for i in range(len(compound_input)):
        for j in range(len(T)):
            V_i[i][j] = M_i[i] / rho_i[i][j]

    
    # Calculate volume fractions
    phi_i=np.zeros((len(Compound),len(T)))
    for j in range(len(T)):
        denom_phi = 0
        for i in range(len(compound_input)):
            phi_i[i][j] = V_i[i][j]*mol_fracs[i]
            denom_phi += V_i[i][j]*mol_fracs[i]
        for i in range(len(compound_input)):
            phi_i[i][j]=phi_i[i][j]/denom_phi

    # Calculate mass (fractions
    kappa_i=np.zeros(len(Compound))
    denom_kappa = 0
    for i in range(len(compound_input)):
        kappa_i[i] = M_i[i]*mol_fracs[i]
        denom_kappa += M_i[i]*mol_fracs[i]
    kappa_i=kappa_i/denom_kappa           


    #Calculate mixture molecular weight
    M_mix=0
    for i in range(len(compound_input)):
        M_mix=M_mix+mol_fracs[i]*M_i[i]


    # Calculate the average temp-dependent molar volume of mixture
    V_mix = np.zeros(len(T))      
    V_mix_est = np.zeros(len(T)) 
    V_mix_data = np.zeros(len(T)) 
    for i in range(len(compound_input)):
        for j in range(len(T)):
            V_mix_est[j] += mol_fracs[i]*V_i[i][j]
            if density_mix == [0,0]:
                pass
            else:
                V_mix_data[j] = M_mix/(density_mix[0] + density_mix[1]*T[j])
    if density_mix == [0,0]:
        V_mix = V_mix_est
    else:
        V_mix = V_mix_data

        V_mix_est_avg = np.average(V_mix_est)
        V_mix_data_avg = np.average(V_mix_data)
        error = 100 * (V_mix_data_avg-V_mix_est_avg)/V_mix_data_avg
        #print("Molar volume from estimated mix density: ", V_mix_est_avg)
        #print("Molar volume from data mix density: ", V_mix_data_avg)
        #print("Molar volume % Difference: ", error, " %")


    #Calculate the temp-dependent density of mixture
    rho_mix = np.zeros(len(T))
    rho_mix_est = np.zeros(len(T)) 
    rho_mix_data = np.zeros(len(T)) 
    for i in range(len(compound_input)):
        for j in range(len(T)):
            rho_mix_est[j] += rho_i[i][j]*mol_fracs[i]   #M_mix/V_mix[j]
            rho_mix_data[j] = density_mix[0] + density_mix[1]*T[j]   # If mixture temp-dependent density data exists
    if density_mix == [0,0]:
        rho_mix = rho_mix_est
    else:
        rho_mix = rho_mix_data

        rho_mix_est_avg = np.average(rho_mix_est)
        rho_mix_data_avg = np.average(rho_mix_data)
        error = 100 * (rho_mix_data_avg-rho_mix_est_avg)/rho_mix_data_avg
        #print("Estimated mix density: ", rho_mix_est_avg)
        #print("Data mix density: ", rho_mix_data_avg)
        print("Density % Difference: ", error, " %")


    # Calculate the temp-dependent sound velocity of mixture
    C_0_mix = np.zeros(len(T))
    C_0_mix_est = np.zeros(len(T)) 
    C_0_mix_data = np.zeros(len(T)) 
    for i in range(len(compound_input)):
        for j in range(len(T)):
            C_0_mix_est[j] += phi_i[i][j]**2 / (kappa_i[i] * C_0_i[i][j]**2)
            C_0_mix_data[j] = sound_velocity_mix[1]*T[j]+sound_velocity_mix[0]    # If mixture temp-dependent sound velocity data exists
    for j in range(len(T)):
        C_0_mix_est[j] = 1 / np.sqrt(C_0_mix_est[j])
    if sound_velocity_mix == [0,0]:
        C_0_mix = C_0_mix_est
    else:
        C_0_mix = C_0_mix_data

        C_0_mix_est_avg = np.average(C_0_mix_est)
        C_0_mix_data_avg = np.average(C_0_mix_data)
        error = 100 * (C_0_mix_data_avg-C_0_mix_est_avg)/C_0_mix_data_avg
        #print("Estimated mix sound velocity: ", C_0_mix_est_avg)
        #print("Data mix sound velocity: ", C_0_mix_data_avg)
        print("sound velocity % Difference: ", error, " %")


    # Calculate mixture specific heat & # of ions
    C_p_mix=0
    C_p_mix_est = 0
    C_p_mix_data = 0
    n_mix=0
    for i in range(len(compound_input)):
        C_p_mix_est += C_p_mix+mol_fracs[i]*C_i_p[i]
        C_p_mix_data=specific_heat_mix
        n_mix=n_mix+mol_fracs[i]*n_i[i]
    if specific_heat_mix[0] == 0:
        C_p_mix = C_p_mix_est
    else:
        if specific_heat_mix[2] == 'm':
            C_p_mix = C_p_mix_data[1]*T_melt+C_p_mix_data[0]
        elif specific_heat_mix[2] == 'g':   # Converts units from J/g/K to J/mol/K
            C_p_mix = (C_p_mix_data[0]*T_melt+C_p_mix_data[1])*(M_mix)

        C_p_mix_est_avg = np.average(C_p_mix_est)
        C_p_mix_data_avg = np.average(C_p_mix)
        error = 100 * (C_p_mix_data_avg-C_p_mix_est_avg)/C_p_mix_data_avg
        #print("Estimated mix specific heat: ", C_p_mix_est_avg)
        #print("Data mix specific heat: ", C_p_mix_data_avg)
        print("specific heat % Difference: ", error, " %")



    # Calculate mixture radial distance
    if r_ac_mix == 0:
        # Calculate mixture radial distance
        r_ac_mix1 = 0
        for i in range(len(compound_input)):
            r_ac_mix=r_ac_mix1+mol_fracs[i]*(r_a[i]+r_c[i])
    else:
        r_ac_mix = r_ac_mix
        r_ac_mix_calc = 0
        # Calculate mixture radial distance
        r_ac_mix1 = 0
        for i in range(len(compound_input)):
            r_ac_mix_calc=r_ac_mix1+mol_fracs[i]*(r_a[i]+r_c[i]) 
        #print("PDF CatAn peak MFP: ", r_ac_mix)
        print("Avg sum radii: ", r_ac_mix_calc)
        #print("MFP % Difference: ", error, " %")

    
    # print("Calculated mixture specific heat (melting point): ",C_p_mix/V_mix[0])
    # print("Calculated mixtur sound velocity (melting point): ",C_0_mix[0])
    # print("Calculated average mean free path: ",r_ac_mix)

    #Calculate compound thermal conductivity
    lambda_i_m=np.zeros(len(T))
    lambda_i_mg=np.zeros(len(T))
    lambda_i_mb=np.zeros(len(T))
    sound_w_time = 0
    for j in range(len(T)):
        lambda_i_m[j] = 1/3 * C_i_p_sp[i] * rho_mix[j] * C_0_mix[j] * r_ac_mix    # Verified with Zhao's results, uses Zhao's data and specific heat capacity
        #lambda_i_mg[j] = 1/3 * C_p_mix * 1/M_mix * rho_mix[j] * C_0_mix[j] * r_ac_mix  # Verified with Zhao's results, uses Zhao's data and but MSTDB heat capacity
        lambda_i_mb[j] = 1/3 * C_p_mix * 1/V_mix[j] * C_0_mix[j] * r_ac_mix  # Uses molar volume at melting point only 
    

    # if sound_w_time == 1:
    #     print("Calculated with temp-dependent sound velocity data.")
    # else:
    #     print("No temp-dependent sound velocity data available. Calculated with melting temp sound velocity only.")
    

    nan_check = np.isnan(lambda_i_m)
    contains_nan = nan_check.any()
    if contains_nan or lambda_i_m[0] == 0:
        print("k: ",lambda_i_mb[0])
        # print("PGM-PDF_Avg (1/3*Cp*vs): ",lambda_i_mb[0]/r_ac_mix)
        # print("PGM-PDF_Avg (1/3*Cp*MFP): ",lambda_i_mb[0]/C_0_mix[0])
        # print("PGM-PDF_Avg (1/3*vs*MFP): ",lambda_i_mb[0]/(C_p_mix/V_mix[0]))
        # print("lambda_BC: ",lambda_i_mb[0]/(1/3*C_p_mix*C_0_mix[0]))
        # print("lambda_SCL: ",r_ac_mix)
        print("cp: ",C_p_mix/V_mix[0])
        print("vs: ",C_0_mix[0])
        return(T,lambda_i_mb)
    else:
        print("k: ",lambda_i_m[0])
        # print("PGM-PDF_Avg (1/3*Cp*vs): ",lambda_i_mb[0]/r_ac_mix)
        # print("PGM-PDF_Avg (1/3*Cp*MFP): ",lambda_i_mb[0]/C_0_mix[0])
        # print("PGM-PDF_Avg (1/3*vs*MFP): ",lambda_i_mb[0]/C_i_p_sp[0])
        # print("lambda_BC: ",lambda_i_m[0]/(1/3*C_p_mix*C_0_mix[0]))
        # print("lambda_SCL: ",r_ac_mix)
        print("cp: ",C_i_p_sp[0]* rho_mix[0])
        print("vs: ",C_0_mix[0])
        print("lambda: ", r_ac_mix)
        return(T,lambda_i_m)

def GECM_VelCorr(df,MSTDB_df,SCL_PDF_df,compound_input,mol_fracs,Temp_Range, V_m=0, density_mix=0, sound_velocity_mix=0, C_p_mix=0, alpha=0, expon=0):
    print(" ")
    print("GECM_VelCorr ######################################")
    print(compound_input)
    #Define smaller dataframes of each variable contained in the spreadsheet
    Compound_df=df['Compound']
    mr_i_df=df['m_r (g/mol)']
    Mc_i_df=df['Mc (g/mol)']
    M_i_df=df['M_i (g/mol)']
    T_i_m_df=df['T_i_m (K)']
    V_i_m_df=df['V_i_m (m^3/mol)']
    C_i_p_df=df['C_i_p (J/mol/K)']
    C_i_p_sp_df=df['C_i_p_sp (J/g/K)']
    r_c_df=df['r_c (m)']
    r_a_df=df['r_a (m)']
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
    mr_i=np.zeros(len(compound_input))
    Mc_i=np.zeros(len(compound_input))
    T_i_m=np.zeros(len(compound_input))
    V_i_m=np.zeros(len(compound_input))
    C_i_p=np.zeros(len(compound_input))
    C_i_p_sp=np.zeros(len(compound_input))
    r_c=np.zeros(len(compound_input))
    r_a=np.zeros(len(compound_input))
    C_0_T_A=np.zeros(len(compound_input))
    C_0_T_B=np.zeros(len(compound_input))
    rho_T_A=np.zeros(len(compound_input))
    rho_T_B=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        Compound[i]=Compound_df[indices[i]]
        M_i[i]=M_i_df[indices[i]]*0.001         # Convert to kg/mol
        mr_i[i]=mr_i_df[indices[i]]*0.001         # Convert to kg/mol
        Mc_i[i]=Mc_i_df[indices[i]]*0.001         # Convert to kg/mol
        T_i_m[i]=T_i_m_df[indices[i]]
        V_i_m[i]=V_i_m_df[indices[i]]
        C_i_p[i]=C_i_p_df[indices[i]]
        C_i_p_sp[i]=C_i_p_sp_df[indices[i]]/0.001   # Convert to J/kg/K
        r_c[i]=r_c_df[indices[i]]
        r_a[i]=r_a_df[indices[i]]
        C_0_T_A[i]=C_0_T_A_df[indices[i]]
        C_0_T_B[i]=C_0_T_B_df[indices[i]]
        rho_T_A[i]=rho_T_A_df[indices[i]]/0.001
        rho_T_B[i]=rho_T_B_df[indices[i]]/0.001

    # print('')
    # print('####### Phonon Gas Model ############################')

    specific_heat_mix = C_p_mix


    #initialize a temperature range
    T_melt = Temp_Range[0]
    T=np.linspace(T_melt,Temp_Range[1],num=100)

    #Calculate number of ions per compound
    n_i=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        for c in Compound[i]:
            if c.isupper():
                n_i[i] = n_i[i] + 1
            elif c.isnumeric():
                c = int(c)
                n_i[i] = n_i[i] + (c-1)
    
    #print("PGM, number of ions: ",n_i)


    # if Compound == ['NaCl','UCl3'] and mol_fracs == [0.63,0.37]:
    #     density_mix = [4.22*100**3/1000,-0.00113*100**3/1000]   # Desyatnik, 1975
    # elif Compound == ['LiF','BeF2']:
    #     density_mix = [2518.15,-0.424]  # Vidrio, 2022, FLiBe 33.59mol%BeF2
    #     sound_velocity_mix = [4272.309646,-1.101929884] # 66-34% Cantor, 1968
    #     specific_heat_mix = 2414.7*0.033  # 67-33% Sohal, 2010  
    #     #density_mix = [2413,-0.488]  # Janz, 1974, FLiBe 33mol%BeF2
    #     #2110 * (1.885 + 2.762*mol_fracs[1] + mol_fracs[1]**2) / (1.773 + 2.663*mol_fracs[0] )
    # elif Compound == ['LiF','NaF','KF']:
    #     density_mix = [2729.3,-0.73]  # 46.5-11.5-42%, Vriesema [1979], Ingersoll et al. [2007], and Williams et al. [2006]
    #     sound_velocity_mix = [3241.15,-1.20]    # 46.5-11.5-42% Robertson, 2022
    #     specific_heat_mix = 1882.8*0.0413  # 46.5-11.5-42% Sohal, 2010
    # else:
    density_mix = [0,0]
    sound_velocity_mix = [0,0]    
    specific_heat_mix = [0,0]  


    #Find the temperature dependent density from data
    rho_i = np.zeros((len(Compound),len(T)))
    for i in range(len(Compound)):
        if pd.isna(rho_T_A[i]):     
            rho_T_A[i], rho_T_B[i] = density(df,compound_input[i])
        else:
            pass
        for j in range(len(T)):
            rho_i[i][j] = rho_T_A[i] + rho_T_B[i]*T[j] 

    
    # Calculate the temp-dependent sound velocity of compounds
    C_0_correction_i = np.zeros(len(Compound))
    C_0_i = np.zeros((len(Compound),len(T)))
    for i in range(len(compound_input)):
        C_0_correction_i[i] = np.sqrt(Mc_i[i]/(2*mr_i[i]))  # Correction for sound velocity with fast velocity of metal ions, Hosokawa ,2009   
        for j in range(len(T)):
            C_0_i[i][j] = (C_0_T_B[i]*T[j]+C_0_T_A[i])/C_0_correction_i[i] 


    # Calculate molar volume of compounds from temp-dependent density
    V_i = np.zeros((len(Compound),len(T)))
    for i in range(len(compound_input)):
        for j in range(len(T)):
            V_i[i][j] = M_i[i] / rho_i[i][j]

    
    # Calculate volume fractions
    phi_i=np.zeros((len(Compound),len(T)))
    for j in range(len(T)):
        denom_phi = 0
        for i in range(len(compound_input)):
            phi_i[i][j] = V_i[i][j]*mol_fracs[i]
            denom_phi += V_i[i][j]*mol_fracs[i]
        for i in range(len(compound_input)):
            phi_i[i][j]=phi_i[i][j]/denom_phi

    # Calculate mass (fractions
    kappa_i=np.zeros(len(Compound))
    denom_kappa = 0
    for i in range(len(compound_input)):
        kappa_i[i] = M_i[i]*mol_fracs[i]
        denom_kappa += M_i[i]*mol_fracs[i]
    kappa_i=kappa_i/denom_kappa           


    #Calculate mixture molecular weight
    M_mix=0
    for i in range(len(compound_input)):
        M_mix=M_mix+mol_fracs[i]*M_i[i]


    # Calculate the average temp-dependent molar volume of mixture
    V_mix = np.zeros(len(T))      
    V_mix_est = np.zeros(len(T)) 
    V_mix_data = np.zeros(len(T)) 
    for i in range(len(compound_input)):
        for j in range(len(T)):
            V_mix_est[j] += mol_fracs[i]*V_i[i][j]
            if density_mix == [0,0]:
                pass
            else:
                V_mix_data[j] = M_mix/(density_mix[0] + density_mix[1]*T[j])
    if density_mix == [0,0]:
        V_mix = V_mix_est
    else:
        V_mix = V_mix_data

        V_mix_est_avg = np.average(V_mix_est)
        V_mix_data_avg = np.average(V_mix_data)
        error = 100 * (V_mix_data_avg-V_mix_est_avg)/V_mix_data_avg
        #print("Molar volume from estimated mix density: ", V_mix_est_avg)
        #print("Molar volume from data mix density: ", V_mix_data_avg)
        #print("Molar volume % Difference: ", error, " %")


    #Calculate the temp-dependent density of mixture
    rho_mix = np.zeros(len(T))
    rho_mix_est = np.zeros(len(T)) 
    rho_mix_data = np.zeros(len(T)) 
    for i in range(len(compound_input)):
        for j in range(len(T)):
            rho_mix_est[j] += rho_i[i][j]*mol_fracs[i]   #M_mix/V_mix[j]
            rho_mix_data[j] = density_mix[0] + density_mix[1]*T[j]   # If mixture temp-dependent density data exists
    if density_mix == [0,0]:
        rho_mix = rho_mix_est
    else:
        rho_mix = rho_mix_data

        rho_mix_est_avg = np.average(rho_mix_est)
        rho_mix_data_avg = np.average(rho_mix_data)
        error = 100 * (rho_mix_data_avg-rho_mix_est_avg)/rho_mix_data_avg
        #print("Estimated mix density: ", rho_mix_est_avg)
        #print("Data mix density: ", rho_mix_data_avg)
        print("Density % Difference: ", error, " %")


    # Calculate the temp-dependent sound velocity of mixture
    C_0_mix = np.zeros(len(T))
    C_0_mix_est = np.zeros(len(T)) 
    C_0_mix_data = np.zeros(len(T)) 
    for i in range(len(compound_input)):
        for j in range(len(T)):
            C_0_mix_est[j] += phi_i[i][j]**2 / (kappa_i[i] * C_0_i[i][j]**2)
            C_0_mix_data[j] = sound_velocity_mix[1]*T[j]+sound_velocity_mix[0]    # If mixture temp-dependent sound velocity data exists
    for j in range(len(T)):
        C_0_mix_est[j] = 1 / np.sqrt(C_0_mix_est[j])
    if sound_velocity_mix == [0,0]:
        C_0_mix = C_0_mix_est
    else:
        C_0_mix = C_0_mix_data

        C_0_mix_est_avg = np.average(C_0_mix_est)
        C_0_mix_data_avg = np.average(C_0_mix_data)
        error = 100 * (C_0_mix_data_avg-C_0_mix_est_avg)/C_0_mix_data_avg
        #print("Estimated mix sound velocity: ", C_0_mix_est_avg)
        #print("Data mix sound velocity: ", C_0_mix_data_avg)
        print("sound velocity % Difference: ", error, " %")


    # Calculate mixture specific heat & # of ions
    C_p_mix=0
    C_p_mix_est = 0
    C_p_mix_data = 0
    n_mix=0
    for i in range(len(compound_input)):
        C_p_mix_est += C_p_mix+mol_fracs[i]*C_i_p[i]
        C_p_mix_data=specific_heat_mix
        n_mix=n_mix+mol_fracs[i]*n_i[i]
    if specific_heat_mix[0] == 0:
        C_p_mix = C_p_mix_est
    else:
        C_p_mix = C_p_mix_data

        C_p_mix_est_avg = np.average(C_p_mix_est)
        C_p_mix_data_avg = np.average(C_p_mix_data)
        error = 100 * (C_p_mix_data_avg-C_p_mix_est_avg)/C_p_mix_data_avg
        #print("Estimated mix specific heat: ", C_p_mix_est_avg)
        #print("Data mix specific heat: ", C_p_mix_data_avg)
        print("specific heat % Difference: ", error, " %")



    # Calculate mixture radial distance
    r_ac_mix = 0
    r_ac_mix1 = 0
    for i in range(len(compound_input)):
        r_ac_mix=r_ac_mix1+mol_fracs[i]*(r_a[i]+r_c[i])


    #print("Calculated average mean free path: ",r_ac_mix)

    #Calculate compound thermal conductivity
    lambda_i_m=np.zeros(len(T))
    lambda_i_mg=np.zeros(len(T))
    lambda_i_mb=np.zeros(len(T))
    sound_w_time = 0
    for j in range(len(T)):
        lambda_i_m[j] = 1/3 * C_i_p_sp[i] * rho_mix[j] * C_0_mix[j] * r_ac_mix    # Verified with Zhao's results, uses Zhao's data and specific heat capacity
        #lambda_i_mg[j] = 1/3 * C_p_mix * 1/M_mix * rho_mix[j] * C_0_mix[j] * r_ac_mix  # Verified with Zhao's results, uses Zhao's data and but MSTDB heat capacity
        lambda_i_mb[j] = 1/3 * C_p_mix * 1/V_mix[j] * C_0_mix[j] * r_ac_mix  # Uses molar volume at melting point only 
    

    # if sound_w_time == 1:
    #     print("Calculated with temp-dependent sound velocity data.")
    # else:
    #     print("No temp-dependent sound velocity data available. Calculated with melting temp sound velocity only.")
    

    nan_check = np.isnan(lambda_i_m)
    contains_nan = nan_check.any()
    if contains_nan or lambda_i_m[0] == 0:
        print("PGM-VelCorr: ",lambda_i_mb[0])
        return(T,lambda_i_mb)
    else:
        print("PGM-VelCorr: ",lambda_i_m[0])
        return(T,lambda_i_m)

def Zhao_PGM_Mix(df,MSTDB_df,SCL_PDF_df,compound_input,mol_fracs,Temp_Range, V_m=0, density_mix=0, sound_velocity_mix=0, C_p_mix=0, alpha=0, expon=0):
    print(" ")
    print("Zhao_PGM_Mix ######################################")
    print(compound_input)
    #Define smaller dataframes of each variable contained in the spreadsheet
    Compound_df=df['Compound']
    M_i_df=df['M_i (g/mol)']
    T_i_m_df=df['T_i_m (K)']
    V_i_m_df=df['V_i_m (m^3/mol)']
    C_i_p_df=df['C_i_p (J/mol/K)']
    C_i_p_sp_df=df['C_i_p_sp (J/g/K)']
    r_c_df=df['r_c (m)']
    r_a_df=df['r_a (m)']
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
    C_i_p=np.zeros(len(compound_input))
    C_i_p_sp=np.zeros(len(compound_input))
    r_c=np.zeros(len(compound_input))
    r_a=np.zeros(len(compound_input))
    C_0_T_A=np.zeros(len(compound_input))
    C_0_T_B=np.zeros(len(compound_input))
    rho_T_A=np.zeros(len(compound_input))
    rho_T_B=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        Compound[i]=Compound_df[indices[i]]
        M_i[i]=M_i_df[indices[i]]*0.001         # Convert to kg/mol
        T_i_m[i]=T_i_m_df[indices[i]]
        V_i_m[i]=V_i_m_df[indices[i]]
        C_i_p[i]=C_i_p_df[indices[i]]
        C_i_p_sp[i]=C_i_p_sp_df[indices[i]]/0.001   # Convert to J/kg/K
        r_c[i]=r_c_df[indices[i]]
        r_a[i]=r_a_df[indices[i]]
        C_0_T_A[i]=C_0_T_A_df[indices[i]]
        C_0_T_B[i]=C_0_T_B_df[indices[i]]
        rho_T_A[i]=rho_T_A_df[indices[i]]/0.001
        rho_T_B[i]=rho_T_B_df[indices[i]]/0.001

    print('')
    print('################ Phonon Gas Model, Mix Data ############################')

    specific_heat_mix = C_p_mix
    density_mix = [0,0]
    sound_velocity_mix = [0,0]    
    molar_volume_mix = [0,0] 
    
    # Obtain measurement data for function
    func_name = inspect.currentframe().f_code.co_name
    r_ac_mix, density_mix, sound_velocity_mix, specific_heat_mix = prop_lookup(func_name,Compound,mol_fracs)


    #initialize a temperature range
    T_melt = Temp_Range[0]
    T=np.linspace(T_melt,Temp_Range[1],num=100)

    #Calculate number of ions per compound
    n_i=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        for c in Compound[i]:
            if c.isupper():
                n_i[i] = n_i[i] + 1
            elif c.isnumeric():
                c = int(c)
                n_i[i] = n_i[i] + (c-1)
    
    #print("PGM, number of ions: ",n_i)


    #Find the temperature dependent density from data
    rho_i = np.zeros((len(Compound),len(T)))
    for i in range(len(Compound)):
        if pd.isna(rho_T_A[i]):     
            rho_T_A[i], rho_T_B[i] = density(df,compound_input[i])
        else:
            pass
        for j in range(len(T)):
            rho_i[i][j] = rho_T_A[i] + rho_T_B[i]*T[j] 

    
    # Calculate the temp-dependent sound velocity of compounds
    C_0_i = np.zeros((len(Compound),len(T)))
    for i in range(len(compound_input)):
        for j in range(len(T)):
            C_0_i[i][j] = C_0_T_B[i]*T[j]+C_0_T_A[i]


    # Calculate molar volume of compounds from temp-dependent density
    V_i = np.zeros((len(Compound),len(T)))
    for i in range(len(compound_input)):
        for j in range(len(T)):
            V_i[i][j] = M_i[i] / rho_i[i][j]

    
    # Calculate volume fractions
    phi_i=np.zeros((len(Compound),len(T)))
    for j in range(len(T)):
        denom_phi = 0
        for i in range(len(compound_input)):
            phi_i[i][j] = V_i[i][j]*mol_fracs[i]
            denom_phi += V_i[i][j]*mol_fracs[i]
        for i in range(len(compound_input)):
            phi_i[i][j]=phi_i[i][j]/denom_phi

    # Calculate mass (fractions
    kappa_i=np.zeros(len(Compound))
    denom_kappa = 0
    for i in range(len(compound_input)):
        kappa_i[i] = M_i[i]*mol_fracs[i]
        denom_kappa += M_i[i]*mol_fracs[i]
    kappa_i=kappa_i/denom_kappa           


    #Calculate mixture molecular weight
    M_mix=0
    for i in range(len(compound_input)):
        M_mix += mol_fracs[i]*M_i[i]


    # Calculate the average temp-dependent molar volume of mixture
    V_mix = np.zeros(len(T))      
    V_mix_est = np.zeros(len(T)) 
    V_mix_data = np.zeros(len(T)) 
    V_mix_data_Vm = np.zeros(len(T)) 
    for i in range(len(compound_input)):
        for j in range(len(T)):
            V_mix_est[j] += mol_fracs[i]*V_i[i][j]
            V_mix_data[j] = M_mix/(density_mix[0] + density_mix[1]*T[j])
            V_mix_data_Vm[j] = molar_volume_mix[0] + molar_volume_mix[1]*T[j]
    if molar_volume_mix == [0,0]:
        if density_mix == [0,0]:
            V_mix = V_mix_est
        else:
            V_mix = V_mix_data
    else:
        V_mix = V_mix_data_Vm


        V_slope_est,V_int_est = np.polyfit(T,V_mix_est,1)
        V_data = molar_volume_mix[1]*T_melt+molar_volume_mix[0]
        V_est = V_slope_est*T_melt+V_int_est
        
        error_C0_slope = 100 * (molar_volume_mix[1]-V_slope_est)/molar_volume_mix[1]
        error_C0_melt = 100 * (V_data-V_est)/V_data
        # error_C0_slope = 100 * (sound_velocity_mix[1]-V_slope_est)/sound_velocity_mix[1]
        # error_C0_int = 100 * (sound_velocity_mix[0]-V_int_est)/sound_velocity_mix[0]
        #print("Molar volume Melt Temp- Estimated:", V_est)
        #print("Molar volume Melt Temp- Data:", V_data)
        #print("Molar volume Slope % Difference: ", error_C0_slope, " %")
        #print("Molar volume Melt Temp % Difference: ", error_C0_melt, " %")
        # print("Sound Velocity Intercept % Difference: ", error_C0_int, " %")
        # error_V_slope = 100 * (molar_volume_mix[1]-V_slope_est)/molar_volume_mix[1]
        # error_V_int = 100 * (molar_volume_mix[0]-V_int_est)/molar_volume_mix[0]
        # print("Molar Volume Slope % Difference: ", error_V_slope, " %")
        # print("Molar Volume Intercept % Difference: ", error_V_int, " %")

        # V_mix_est_avg = np.average(V_mix_est)
        # V_mix_data_avg = np.average(V_mix_data)
        # V_mix_data_Vm_avg = np.average(V_mix_data_Vm)
        # error_densityVm = 100 * (V_mix_data_avg-V_mix_est_avg)/V_mix_data_avg
        # error_Vm = 100 * (V_mix_data_Vm_avg-V_mix_est_avg)/V_mix_data_Vm_avg
        # #print("")
        # print("Molar Volume - From estimated density: ", V_mix_est_avg)
        # print("Molar volume - From mix data density: ", V_mix_data_avg)
        # print("Molar volume - Mix data Molar Volume: ", V_mix_data_Vm_avg)
        # print("Molar Volume % Difference, Density Data/Estimation: ", error_densityVm, " %")
        # print("Molar Volume % Difference, Molar Volume Data/Estimation: ", error_Vm, " %")


    #Calculate the temp-dependent density of mixture
    rho_mix = np.zeros(len(T))
    rho_mix_est = np.zeros(len(T)) 
    rho_mix_data = np.zeros(len(T)) 
    for i in range(len(compound_input)):
        for j in range(len(T)):
            rho_mix_est[j] += rho_i[i][j]*mol_fracs[i]   #M_mix/V_mix[j]
            rho_mix_data[j] = density_mix[0] + density_mix[1]*T[j]   # If mixture temp-dependent density data exists
    if density_mix == [0,0]:
        rho_mix = rho_mix_est
    else:
        rho_mix = rho_mix_data

        rho_slope_est,rho_int_est = np.polyfit(T,rho_mix_est,1)
        error_rho_slope = 100 * (density_mix[1]-rho_slope_est)/density_mix[1]
        error_rho_int = 100 * (density_mix[0]-rho_int_est)/density_mix[0]
        print("Density Slope % Difference: ", error_rho_slope, " %")
        print("Density Intercept % Difference: ", error_rho_int, " %")
        # rho_mix_est_avg = np.average(rho_mix_est)
        # rho_mix_data_avg = np.average(rho_mix_data)
        # error = 100 * (rho_mix_data_avg-rho_mix_est_avg)/rho_mix_data_avg
        # print("Estimated mix density: ", rho_mix_est_avg)
        # print("Data mix density: ", rho_mix_data_avg)
        # print("Density % Difference: ", error, " %")


    # Calculate the temp-dependent sound velocity of mixture
    C_0_mix = np.zeros(len(T))
    C_0_mix_est = np.zeros(len(T)) 
    C_0_mix_data = np.zeros(len(T)) 
    for i in range(len(compound_input)):
        for j in range(len(T)):
            C_0_mix_est[j] += phi_i[i][j]**2 / (kappa_i[i] * C_0_i[i][j]**2)
            C_0_mix_data[j] = sound_velocity_mix[1]*T[j]+sound_velocity_mix[0]    # If mixture temp-dependent sound velocity data exists
    for j in range(len(T)):
        C_0_mix_est[j] = 1 / np.sqrt(C_0_mix_est[j])
    if sound_velocity_mix == [0,0]:
        C_0_mix = C_0_mix_est
    else:
        C_0_mix = C_0_mix_data

        C0_slope_est,C0_int_est = np.polyfit(T,C_0_mix_est,1)
        vsound_data = sound_velocity_mix[1]*T_melt+sound_velocity_mix[0]
        vsound_est = C0_slope_est*T_melt+C0_int_est
        
        error_C0_slope = 100 * (sound_velocity_mix[1]-C0_slope_est)/sound_velocity_mix[1]
        error_C0_melt = 100 * (vsound_data-vsound_est)/vsound_data
        # error_C0_slope = 100 * (sound_velocity_mix[1]-C0_slope_est)/sound_velocity_mix[1]
        # error_C0_int = 100 * (sound_velocity_mix[0]-C0_int_est)/sound_velocity_mix[0]
        print("Sound Velocity Slope % Difference: ", error_C0_slope, " %")
        print("Sound Velocity Melt Temp % Difference: ", error_C0_melt, " %")
        # print("Sound Velocity Intercept % Difference: ", error_C0_int, " %")

        # C_0_mix_est_avg = np.average(C_0_mix_est)
        # C_0_mix_data_avg = np.average(C_0_mix_data)
        # error = 100 * (C_0_mix_data_avg-C_0_mix_est_avg)/C_0_mix_data_avg
        # print("Estimated mix sound velocity: ", C_0_mix_est_avg)
        # print("Data mix sound velocity: ", C_0_mix_data_avg)
        # print("sound velocity % Difference: ", error, " %")


    # Calculate mixture specific heat & # of ions
    C_p_mix=0
    C_p_mix_est = 0
    C_p_mix_data = 0
    n_mix=0
    for i in range(len(compound_input)):
        C_p_mix_est += C_p_mix+mol_fracs[i]*C_i_p[i]
        C_p_mix_data=specific_heat_mix
        n_mix=n_mix+mol_fracs[i]*n_i[i]
    if specific_heat_mix[0] == 0:
        C_p_mix = C_p_mix_est
    else:
        if specific_heat_mix[2] == 'm':
            C_p_mix = C_p_mix_data[1]*T_melt+C_p_mix_data[0]
        elif specific_heat_mix[2] == 'g':   # Converts units from J/g/K to J/mol/K
            C_p_mix = (C_p_mix_data[0]*T_melt+C_p_mix_data[1])*(M_mix)

        C_p_mix_est_avg = np.average(C_p_mix_est)
        C_p_mix_data_avg = np.average(C_p_mix)
        error = 100 * (C_p_mix_data_avg-C_p_mix_est_avg)/C_p_mix_data_avg
        print("Specific Heat - Estimated: ", C_p_mix_est_avg)
        print("Specific Heat - Mix Data: ", C_p_mix_data_avg)
        print("Specific Heat % Difference: ", error, " %")



    # Calculate mixture radial distance
    r_ac_mix = 0
    r_ac_mix1 = 0
    for i in range(len(compound_input)):
        r_ac_mix=r_ac_mix1+mol_fracs[i]*(r_a[i]+r_c[i])

    #print("Calculated average mean free path: ",r_ac_mix)


    #Calculate compound thermal conductivity
    lambda_i_m=np.zeros(len(T))
    lambda_i_mg=np.zeros(len(T))
    lambda_i_mb=np.zeros(len(T))
    sound_w_time = 0
    for j in range(len(T)):
        lambda_i_m[j] = 1/3 * C_i_p_sp[i] * rho_mix[j] * C_0_mix[j] * r_ac_mix    # Verified with Zhao's results, uses Zhao's data and specific heat capacity
        #lambda_i_mg[j] = 1/3 * C_p_mix * 1/M_mix * rho_mix[j] * C_0_mix[j] * r_ac_mix  # Verified with Zhao's results, uses Zhao's data and but MSTDB heat capacity
        lambda_i_mb[j] = 1/3 * C_p_mix * 1/V_mix[j] * C_0_mix[j] * r_ac_mix  # Uses molar volume at melting point only 
    

    # if sound_w_time == 1:
    #     print("Calculated with temp-dependent sound velocity data.")
    # else:
    #     print("No temp-dependent sound velocity data available. Calculated with melting temp sound velocity only.")

    nan_check = np.isnan(lambda_i_m)
    contains_nan = nan_check.any()
    if contains_nan or lambda_i_m[0] == 0:
        print("k: ",lambda_i_mb[0])
        print("cp: ",C_p_mix)
        print("vs: ",C_0_mix[0])
        # print("lambda_BC: ",lambda_i_mb[0]/(1/3*C_p_mix*C_0_mix[0]))
        return(T,lambda_i_mb)
    else:
        print("k: ",lambda_i_m[0])
        print("cp: ",C_p_mix)
        print("vs: ",C_0_mix[0])
        # print("lambda_BC: ",lambda_i_m[0]/(1/3*C_p_mix*C_0_mix[0]))
        return(T,lambda_i_m)
    

def Bridgman(df,MSTDB_df,SCL_PDF_df,compound_input,mol_fracs,Temp_Range, V_m=0, C_p=0, alpha=0, expon=0):

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

def ChapmanEnskog(df,MSTDB_df,SCL_PDF_df,compound_input,mol_fracs,Temp_Range, V_m=0, C_p=0, alpha=0, expon=0):
    
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


def Ideal(df,MSTDB_df,SCL_PDF_df,compound_input,mol_fracs,Temp_Range, V_m=0, density_mix=0, C_p_mix=0, alpha=0, expon=0):
    
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
        
    # Calculate compound psi term
    psi_i=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        psi_i[i] = 1 + n_i_c[i]/n_i_a[i]

    # Calculate compound number density
    n_dens_i = np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        n_dens_i[i] = Avog *n_i[i]/V_i_m[i]

    #Calculate compound thermal conductivity(T) based on compound minimum thermal conductivites
    lambda_i = np.zeros((len(Compound),len(T)))
    for i in range(len(compound_input)):
        for j in range(len(T)):
            lambda_i[i][j] = (1 + n_i_c[i]/n_i_a[i]) * k_B * n_dens_i[i]**(2/3) * C_i_0[i] * (1 - alpha_i_m[i]*(gamma_i_m[i] + 1/3)*(T[j] - T_i_m[i]))
 
    #Calculate ideal thermal conductivity
    lambda_mix_T = np.zeros(len(T))
    for i in range(len(compound_input)):
        for j in range(len(T)):
            lambda_mix_T[j] += lambda_i[i][j]*mol_fracs[i]


    print("Ideal Mixing: ",lambda_mix_T[0])
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

def find_SCL_composition(SCL_PDF_df, compound_input, mol_fracs):
    composition_str = '-'.join([f"{mol:.4f}{comp}" for mol, comp in zip(mol_fracs, compound_input)])
    matching_row = SCL_PDF_df[SCL_PDF_df['Composition'] == composition_str]
    return matching_row if not matching_row.empty else None

def get_cell_value(csv_file, target_row_header, target_column_header):
    with open(csv_file, 'r') as file:
        reader = csv.DictReader(file)
        
        for row in reader:
            if row['Formula'] == target_row_header:
                return row[target_column_header]
    
    # If the target row or column header is not found, return None or handle it as needed
    return None

def prop_lookup(method,Compound, mol_fracs):
    r_ac_mix = 0
    density_mix = [0,0]
    sound_velocity_mix = [0,0]
    specific_heat_mix = [0,0]
    if 'GECM' in method:
        if Compound == ['LiCl']:
            r_ac_mix = 4.118098594123714E-10    # 100% Walz, 878 K, 2019
        elif Compound == ['NaCl']:
            r_ac_mix = 4.599211434277691E-10    # 100% Walz, 1074.15 K, 2019 # (other) 100% Andersson, 1250 K, 2022
        elif Compound == ['KCl']:
            r_ac_mix = 5.061494421764213E-10    # 100% Walz, 1043 K, 2019
        elif Compound == ['LiF']:
            r_ac_mix = 3.315990988491789E-10    # 100% Walz, 1121 K, 2019
        elif Compound == ['NaF']:
            r_ac_mix = 3.7969068224247864E-10   # 100% Walz, 1266.15 K, 2019
        elif Compound == ['KF']:
            r_ac_mix = 4.254924144410876E-10    # 100% Walz, 1131.15 K, 2019
        elif Compound == ['RbF']:
            r_ac_mix = 3.143356787E-10          # 100% Walz, 1068.15 K, 2019
        elif Compound == ['CsF']:
            r_ac_mix = 3.392004301E-10          # 100% Walz, 955.15 K, 2019
        elif Compound == ['MgCl2']:
            r_ac_mix = 4.7852054075804915E-10   # 100% McGreevy, 998 K, 1987
        elif Compound == ['CaCl2']:
            r_ac_mix = 5.033965984936765E-10    # 100% McGreevy, 1093 K, 1987
        elif Compound == ['SrCl2']:
            r_ac_mix = 5.10679E-10              # 100% McGreevy, 1198 K, 1987
        elif Compound == ['NaNO3']:
            r_ac_mix = 2.921000675E-10          # 100% 
        elif sorted(Compound) == sorted(['NaCl','UCl3']): # and mol_fracs == [0.63,0.37]:
            r_ac_mix = 2.89258E-10 #3.00685145E-10  # 64-36% Andersson, 2022
        elif sorted(Compound) == sorted(['LiF','BeF2']) and mol_fracs == [0.5,0.5]:
            r_ac_mix = 1.9572e-10  # 50-50% Sun, 2024
        elif sorted(Compound) == sorted(['LiF','BeF2']) and mol_fracs == [0.66,0.34]:
            r_ac_mix = 2.080697E-10 # 2.090906E-10 #1.7449105E-10  first peak   # #1.660250E-10     # 50-50% Sun, 2024
        elif sorted(Compound) == sorted(['LiF','NaF']):
            r_ac_mix = 2.28638E-10  # 60-40% Grizzi, 2024
        elif sorted(Compound) == sorted(['LiF','NaF','KF']):
            r_ac_mix = 2.32906E-10 # 2.618998138E-10  # 46.5-11.5-42% Frandsen, 2020
        elif sorted(Compound) == sorted(['NaF','KF','MgF2']):
            r_ac_mix = 2.56459E-10  # 34.5-59-6.5%, Rudenko, 2024
        elif sorted(Compound) == sorted(['MgCl2','NaCl','KCl']):
            r_ac_mix = 2.74842E-10  # 20.47-41.3-38.23%, Jiang, 2024
        # elif sorted(Compound) == sorted(['LiF','KF','UF4']):
        #     r_ac_mix = 2.733382581251662e-10  # 0.2727LiF-0.1818NaF-0.091UF4 Grizzi, 2024
        elif sorted(Compound) == sorted(['LiF','NaF','UF4']):
            r_ac_mix = 2.19714e-10  # 54.54LiF-336.36NaF-9.1UF4 Grizzi, 2024
        elif sorted(Compound) == sorted(['NaCl','KCl','ZnCl']):
            r_ac_mix = 2.65503e-10  # 0.22NaCl-0.393KCl-0.387ZnCl2",'Xi, 2024; 1073K
        else:
            r_ac_mix = 0

    # Properties of mixtures (if known)
            # density_mix = [A,B]           <-- kg/m^3; A: Intercept, B: Slope
            # sound_velocity_mix = [A,B]    <-- m/s;    A: Intercept, B: Slope
            # specific_heat_mix = [A,B,u]   <-- J/molK; A: Intercept, B: Slope, u: Units ('g'=J/g-K, 'm'=J/mol-K)
    if 'Mix' in method:
        if sorted(Compound) == sorted(['NaCl','UCl3']) and mol_fracs == [0.63,0.37]:
            density_mix = [4.22*100**3/1000,-0.00113*100**3/1000]   # Desyatnik, 1975
        elif sorted(Compound) == sorted(['LiF','BeF2']):
            density_mix = [2518.15,-0.424]  # Vidrio, 2022, FLiBe 33.59mol%BeF2
            sound_velocity_mix = [4272.309646,-1.101929884] # 66-34% Cantor, 1968
            specific_heat_mix = [2414.7*0.033,0,'m']  # 67-33% Sohal, 2010  
            #density_mix = [2413,-0.488]  # Janz, 1974, FLiBe 33mol%BeF2
            #2110 * (1.885 + 2.762*mol_fracs[1] + mol_fracs[1]**2) / (1.773 + 2.663*mol_fracs[0] )
        elif sorted(Compound) == sorted(['LiF','NaF']):
            sound_velocity_mix = [3054,-0.698]    # 63-37% Minchenko, 1985
        elif sorted(Compound) == sorted(['LiF','NaF','KF']):
            density_mix = [2729.3,-0.73]  # 46.5-11.5-42%, Vriesema [1979], Ingersoll et al. [2007], and Williams et al. [2006]
            sound_velocity_mix = [3241.15,-1.20]    # 46.5-11.5-42% Robertson, 2022
            specific_heat_mix = [1882.8*0.0413,0,'m']  # 46.5-11.5-42% Sohal, 2010
            molar_volume_mix = [1.34991e-5,7.55e-9] # Kubikova,2013 
        elif sorted(Compound) == sorted(['NaF','KF','MgF2']):
            density_mix = [2518.951923,-0.451923077]  # 34.5-59-6.5%, Rudenko, 2024
            specific_heat_mix = [71.29715445,0,'m']  # 34.5-59-6.5%, Rudenko, 2024
        elif sorted(Compound) == sorted(['MgCl2','NaCl','KCl']):
            density_mix = [1958.8438,-0.56355]  # 45.98-15.11-38.91%, Wang, 2021
            specific_heat_mix = [1.30138,-0.0005,'g']  # 45.98-15.11-38.91%, Wang, 2021
        elif Compound == ['LiCl']:
            sound_velocity_mix = [3500, 0] 
        elif Compound == ['NaCl']:
            sound_velocity_mix = [2720, 0]
        elif Compound == ['KCl']:
            sound_velocity_mix = [2360, 0]
        elif Compound == ['LiF']:
            sound_velocity_mix = [3500, 0]
        elif Compound == ['NaF']:
            sound_velocity_mix = [2720, 0]
        elif Compound == ['KF']:
            sound_velocity_mix = [2360, 0]
        elif Compound == ['MgCl2']: 
            sound_velocity_mix = [2101, 0] # Gheribi, Theoretical #[3400, 0]
        elif Compound == ['CaCl2']:
            sound_velocity_mix = [2200, 0]
    else:
        density_mix = [0,0]
        sound_velocity_mix = [0,0]
        specific_heat_mix = [0,0]

    return r_ac_mix, density_mix, sound_velocity_mix, specific_heat_mix

def functionlibrary():
    functions = {
    'Gheribi-KTCPλ' : Kinetic_Theory,
    'Gheribi-KTCPλ, Mix Data' : Kinetic_Theory_Mix,
    'Gheribi-KT-Emp' : Kinetic_Theory_Emp,
    'Gheribi-KT-Emp-Mix' : Kinetic_Theory_Emp_Mix,
    'Gheribi-KT24' : Kinetic_Theory24,
    'Gheribi-KT24, Mix Data' : Kinetic_Theory24_Mix,
    'Zhao-PGM' : Zhao_PGM,
    'Zhao-PGM, Mix Data' : Zhao_PGM_Mix,
    'Present Model' : GECM,
    'Present Model, Mix Data' : GECM_Mix,
    'Present Model-VelCorr' : GECM_VelCorr,
    'Bridgman' : Bridgman,
    'Chapman-Enskog' : ChapmanEnskog,
    'Ideal' : Ideal,
    }
    return functions

# def functionlibraryUnabbreviated():
#     functions = {
#     'Kinetic Theory Cahill-Poll λ' : Kinetic_Theory,
#     'Kinetic Theory Cahill-Poll λ, Mixture Data ' : Kinetic_Theory_Mix,
#     'Kinetic Theory Emp-Fit' : Kinetic_Theory_Emp,
#     'Kinetic Theory Emp-Fit, Mixture Data' : Kinetic_Theory_Emp_Mix,
#     'Kinetic Theory 2024' : Kinetic_Theory24,
#     'Kinetic Theory 2024, Mixture Data' : Kinetic_Theory24_Mix,
#     'Phonon Gas Model' : Zhao_PGM,
#     'Phonon Gas Model, PDF-MaxPeak' : Zhao_PGM_PDF_Max,
#     'Phonon Gas Model, PDF-Avg1stPeaks' : Zhao_PGM_PDF_Avg1stPeaks,
#     'Phonon Gas Model, PDF-MaxPeak, Mix Data' : Zhao_PGM_PDF_Max_Mix,
#     'Phonon Gas Model, PDF-Avg1stPeaks, Mix Data' : Zhao_PGM_PDF__Avg1stPeaks_Mix,
#     'Phonon Gas Model, Velocity Correction' : GECM_VelCorr,
#     'Phonon Gas Model, Mix Data' : Zhao_PGM_Mix,
#     'Phonon Gas Model, Melting Point Data' : Zhao_PGM_MP,
#     'Bridgman' : Bridgman,
#     'Chapman-Enskog' : ChapmanEnskog,
#     'Ideal Mixing' : Ideal,
#     }
#     return functions
    