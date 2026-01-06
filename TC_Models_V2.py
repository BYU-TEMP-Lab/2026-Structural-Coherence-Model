import os
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
        print("gmass_Tm: ", (1 - M_i[i]/M_mix)**2, " of compound ", compound_input[i])
        print("k_kin/k_i: ", lambda_i[i][j]/lambda_ideal[j])

    print("delta_Tm: ", delta[0])
    

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
        print("gmass_Tm: ", (1 - M_i[i]/M_mix)**2, " of compound ", compound_input[i])
        print("k_kin/k_i: ", lambda_i[i][j]/lambda_ideal[j])

    print("delta_Tm: ", delta[0])

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
        r_ac_mix += mol_fracs[i]*(r_a[i]+r_c[i])

    
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
        print("cp: ",C_i_p_sp[0])
        print("vs: ",C_0_mix[0])
        return(T,lambda_i_mb)
    else:
        print("k: ",lambda_i_m[0])
        # print("lambda_BC: ",lambda_i_m[0]/(1/3*C_p_mix*C_0_mix))
        print("cp: ",C_p_mix/V_mix[0])
        print("vs: ",C_0_mix[0])
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
        r_ac_mix += mol_fracs[i]*(r_a[i]+r_c[i])

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
        print("cp: ",C_i_p_sp[0])
        print("vs: ",C_0_mix[0])
        # print("lambda_BC: ",lambda_i_mb[0]/(1/3*C_p_mix*C_0_mix))
        return(T,lambda_i_mb)
    else:
        print("k: ",lambda_i_m[0])
        print("cp: ",C_p_mix/V_mix[0])
        print("vs: ",C_0_mix[0])
        # print("lambda_BC: ",lambda_i_m[0]/(1/3*C_p_mix*C_0_mix))
        return(T,lambda_i_m)
   
def GECM(df,MSTDB_df,SCL_PDF_df,compound_input,mol_fracs,Temp_Range, V_m=0, density_mix=0, sound_velocity_mix=0, C_p_mix=0, alpha=0, expon=0):
    print("")
    print("GECM ######################################")
    print(compound_input)
    
    # Find matching SCL data row
    scl_row = find_matching_scl_row(SCL_PDF_df, compound_input, mol_fracs)
    scl_data = scl_row.to_dict() if scl_row is not None else {}

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

    # Get SCL properties for all pairs
    scl_props = get_scl_properties(compound_input, mol_fracs, scl_data)
    pair_scl = scl_props['pair_scl']
    pair_peak_y = scl_props['pair_peak_y']
    pair_min_y = scl_props['pair_min_y']
    
    # Calculate average SCL if available, otherwise use prop_lookup
    if np.any(pair_scl > 0):
        r_ac_mix_scl = np.mean(pair_scl)
        print(f"Using SCL results, average SCL: {r_ac_mix_scl*1e10:.3f} A")
        r_ac_mix = r_ac_mix_scl
    else:
        print("No SCL data available, falling back to prop_lookup")
        r_ac_mix_prop, _, _, _ = prop_lookup(func_name, compound_input, mol_fracs)
        if r_ac_mix_prop == 0:
            print("    No SCL data found in prop_lookup either. Using avg sum radii.")
            # Calculate mixture radial distance
            r_sum = 0
            for i in range(len(compound_input)):
                r_sum += mol_fracs[i]*(r_a[i]+r_c[i])
                print(f"        {compound_input[i]}: r_a + r_c = {(r_a[i]+r_c[i])*1e10:.3f} A")
            r_ac_mix = r_sum
            print(f"    Using average sum of radii: {r_ac_mix*1e10:.3f} A")
        else:
            print(f"    Prop_lookup SCL: {r_ac_mix_prop*1e10:.3f} A")
            r_ac_mix = r_ac_mix_prop

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
    r_ac_mix = 0
    r_ac_mix1 = 0
    for i in range(len(compound_input)):
        r_ac_mix += mol_fracs[i]*(r_a[i]+r_c[i])

    
    print("Calculated mixture specific heat (melting point): ",C_p_mix/V_mix[0])
    print("Calculated mixtur sound velocity (melting point): ",C_0_mix[0])
    print("Calculated average mean free path: ",r_ac_mix)

    # Calculate temperature-dependent slopes for specific heat and sound velocity
    # First, ensure C_p_mix is an array with the same length as T
    if not hasattr(C_p_mix, '__len__') or len(C_p_mix) == 1:
        # If C_p_mix is a scalar or single value, create an array with that value
        C_p_mix_array = np.full_like(T, C_p_mix[0] if hasattr(C_p_mix, '__len__') else C_p_mix)
    else:
        C_p_mix_array = C_p_mix
    
    # Calculate specific heat per volume for all temperatures
    C_p_vol = np.zeros_like(T)
    for j in range(len(T)):
        C_p_vol[j] = C_p_mix_array[j] / V_mix[j]  # J/m³/K
    
    # Use linear regression to find the slope of C_p_vol vs T
    if len(T) > 1:
        specific_heat_slope = np.polyfit(T - T[0], C_p_vol, 1)[0]  # dC_p/dT in J/m³/K²
        sound_velocity_slope = np.polyfit(T - T[0], C_0_mix, 1)[0]  # dC_0/dT in m/s/K
    else:
        specific_heat_slope = 0
        sound_velocity_slope = 0
    
    # Calculate thermal conductivity for both methods
    lambda_i_m_pos = np.zeros(len(T))
    lambda_i_m_sop = np.zeros(len(T))
    
    # k as product of sums
    for j in range(len(T)):
        if V_mix[j] > 0:
            lambda_i_m_pos[j] = 1/3 * C_p_vol[j] * C_0_mix[j] * r_ac_mix  # Uses molar volume at melting point only 
    

    # k as sum of products
    
    # Calculate partial specific heat for each compound
    C_p_part = [np.zeros(len(T)) for _ in range(len(Compound))]
    f_K = np.zeros(len(Compound))
    for i in range(len(Compound)):
        # Calculate bond strength factor from PDF peaks
        f_K[i] = (pair_peak_y[i] - pair_min_y[i]) / pair_peak_y[i]
        for j in range(len(T)):
            C_p_part[i][j] += mol_fracs[i]*C_p_vol[j]*f_K[i]
            C_p_part_denom += mol_fracs[i]*f_K[i]
    C_p_part = C_p_part / C_p_part_denom

    # Calculate partial sound velocity
    C_0_part = [np.zeros(len(T)) for _ in range(len(Compound))]
    for i in range(len(Compound)):
        for j in range(len(T)):
            C_0_part[i][j] += mol_fracs[i]*C_0_mix[j]
            C_0_part_denom += mol_fracs[i]
    C_0_part = C_0_part / C_0_part_denom

    # Calculate thermal conductivity
    for i in range(len(Compound)):
        for j in range(len(T)):
            lambda_i_m_sop[j] += C_p_part[i][j]*C_0_part[i][j]*pair_scl[i] 
    lambda_i_m_sop = 1/3*lambda_i_m_sop
    
    cp_vol_at_melt = C_p_vol[0] if len(C_p_vol) > 0 else 0 
    print(f"k_pos: {lambda_i_m_pos[0]:.3f} W/m/K")
    print(f"k_sop: {lambda_i_m_sop[0]:.3f} W/m/K")
    print(f"cp: {cp_vol_at_melt:.2f} J/m³/K")
    print(f"vs: {C_0_mix[0]:.2f} m/s")
    
    return (T, {
        'thermal_conductivity': lambda_i_m_pos,
        'specific_heat_m': cp_vol_at_melt,  # J/m³/K at melting point
        'specific_heat_prime': specific_heat_slope,  # dC_p/dT in J/m³/K²
        'sound_velocity_m': C_0_mix[0],  # m/s at melting point
        'sound_velocity_prime': sound_velocity_slope  # dC_0/dT in m/s/K
    })


def GECM_Mix(df, MSTDB_df, SCL_PDF_df, compound_input, mol_fracs, Temp_Range, V_m=0, density_mix=0, sound_velocity_mix=0, C_p_mix=0, alpha=0, expon=0):
    print(" ")
    print("GECM_Mix ######################################")
    print(compound_input)
    
    # Find matching SCL data row
    scl_row = find_matching_scl_row(SCL_PDF_df, compound_input, mol_fracs)
    scl_data = scl_row.to_dict() if scl_row is not None else {}
    
    # Define smaller dataframes of each variable contained in the spreadsheet
    Compound_df = df['Compound']
    M_i_df = df['M_i (g/mol)']
    T_i_m_df = df['T_i_m (K)']
    V_i_m_df = df['V_i_m (m^3/mol)']
    C_i_p_df = df['C_i_p (J/mol/K)']
    C_i_p_sp_df = df['C_i_p_sp (J/g/K)']
    r_c_df = df['r_c (m)']
    r_a_df = df['r_a (m)']
    C_0_T_A_df = df['SoS(T)_A (A+B*T)']
    C_0_T_B_df = df['SoS(T)_B']
    rho_T_A_df = df['A (Density (g/cm3):   A - BT(K))']
    rho_T_B_df = df['B (Density (g/cm3):   A - BT(K))']

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

    # Get SCL properties for all pairs
    scl_props = get_scl_properties(compound_input, mol_fracs, scl_data)
    pair_scl = scl_props['pair_scl']
    pair_peak_y = scl_props['pair_peak_y']
    pair_min_y = scl_props['pair_min_y']
    
    # Calculate average SCL if available, otherwise use prop_lookup
    if np.any(pair_scl > 0):
        r_ac_mix_scl = np.mean(pair_scl)
        print(f"Using SCL results, average SCL: {r_ac_mix_scl*1e10:.3f} A")
        r_ac_mix = r_ac_mix_scl
    else:
        print("No SCL data available, falling back to prop_lookup")
        r_ac_mix_prop, _, _, _ = prop_lookup(func_name, compound_input, mol_fracs)
        if r_ac_mix_prop == 0:
            print("    No SCL data found in prop_lookup either. Using avg sum radii.")
            # Calculate mixture radial distance
            r_sum = 0
            for i in range(len(compound_input)):
                r_sum += mol_fracs[i]*(r_a[i]+r_c[i])
                print(f"        {compound_input[i]}: r_a + r_c = {(r_a[i]+r_c[i])*1e10:.3f} A")
            r_ac_mix = r_sum
            print(f"    Using average sum of radii: {r_ac_mix*1e10:.3f} A")
        else:
            print(f"    Prop_lookup SCL: {r_ac_mix_prop*1e10:.3f} A")
            r_ac_mix = r_ac_mix_prop

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
        C_p_mix_data_avg = np.average(C_p_mix_data)
        error = 100 * (C_p_mix_data_avg-C_p_mix_est_avg)/C_p_mix_data_avg
        #print("Estimated mix specific heat: ", C_p_mix_est_avg)
        #print("Data mix specific heat: ", C_p_mix_data_avg)
        print("specific heat % Difference: ", error, " %")
    # print("Calculated mixture specific heat (melting point): ",C_p_mix/V_mix[0])
    # print("Calculated mixtur sound velocity (melting point): ",C_0_mix[0])
    # print("Calculated average mean free path: ",r_ac_mix)

    # Calculate temperature-dependent slopes for specific heat and sound velocity
    # Use linear regression to find the slope of C_p_mix vs T
    if len(T) > 1:
        # Calculate specific heat per volume for all temperatures
        C_p_vol = np.zeros(len(T))
        for j in range(len(T)):
            C_p_vol[j] = C_p_mix / V_mix[j]  # Convert to J/m³/K
        
        # Calculate specific heat slope (dC_p/dT) in J/m³/K²
        specific_heat_slope = np.polyfit(T - T[0], C_p_vol, 1)[0]
        
        # Calculate sound velocity slope (dC_0/dT) in m/s/K
        sound_velocity_slope = np.polyfit(T - T[0], C_0_mix, 1)[0]
    else:
        specific_heat_slope = 0
        sound_velocity_slope = 0

    # Calculate compound thermal conductivity
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
        # print("lambda_BC: ",lambda_i_mb[0]/(1/3*C_p_mix*C_0_mix))
        # print("lambda_SCL: ",r_ac_mix)
        print("cp: ",C_p_mix/V_mix[0])
        print("vs: ",C_0_mix[0])
        return (T, {
            'thermal_conductivity': lambda_i_mb,
            'specific_heat_m': C_p_mix/V_mix[0],
            'specific_heat_prime': specific_heat_slope,  # dC_p/dT in J/m³/K²
            'sound_velocity_m': C_0_mix[0],
            'sound_velocity_prime': sound_velocity_slope  # dC_0/dT in m/s/K
        })
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
        return (T, {
            'thermal_conductivity': lambda_i_m,
            'specific_heat_m': C_i_p_sp[0] * rho_mix[0],
            'specific_heat_prime': specific_heat_slope,  # dC_p/dT in J/m³/K²
            'sound_velocity_m': C_0_mix[0],
            'sound_velocity_prime': sound_velocity_slope  # dC_0/dT in m/s/K
        })

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

def find_matching_scl_row(scl_pdf_df, compound_input, mol_fracs, tol=0.05):
    """Find the SCL data row that matches the given composition.
    
    Args:
        scl_pdf_df: DataFrame containing SCL data
        compound_input: List of compound names (e.g., ['NaCl', 'KCl'])
        mol_fracs: List of mole fractions corresponding to compound_input
        tol: Tolerance for mole fraction comparison (default: 0.05 or 5%)
        
    Returns:
        pandas.Series: The matching SCL data row, or None if no match found
    """
    if scl_pdf_df is None or scl_pdf_df.empty:
        print("No SCL data available in the provided DataFrame")
        return None
        
    # Normalize the input mole fractions
    mol_fracs = np.array(mol_fracs, dtype=float)
    mol_fracs = mol_fracs / np.sum(mol_fracs)  # Ensure they sum to 1
    
    # Convert compound list to a sorted tuple of (compound, fraction) for comparison
    input_comp = sorted(zip(compound_input, mol_fracs), key=lambda x: x[0])
    
    # Try to find a matching composition in the SCL data
    for _, row in scl_pdf_df.iterrows():
        if pd.isna(row['Composition']) or not row['Composition']:
            continue
            
        try:
            # Parse the composition string (e.g., '0.5NaCl-0.5KCl')
            comp_parts = row['Composition'].split('-')
            scl_compounds = []
            scl_fracs = []
            
            for part in comp_parts:
                # Find the index where the compound name starts
                idx = 0
                while idx < len(part) and (part[idx].isdigit() or part[idx] in '.-'):
                    idx += 1
                
                if idx == 0:  # No fraction, assume 1.0
                    scl_frac = 1.0 / len(comp_parts)  # Distribute evenly if no fractions given
                    scl_comp = part
                else:
                    scl_frac = float(part[:idx])
                    scl_comp = part[idx:]
                
                scl_compounds.append(scl_comp)
                scl_fracs.append(scl_frac)
            
            # Normalize the SCL fractions
            scl_fracs = np.array(scl_fracs, dtype=float)
            scl_fracs = scl_fracs / np.sum(scl_fracs)
            
            # Create sorted list of (compound, fraction) for comparison
            scl_comp = sorted(zip(scl_compounds, scl_fracs), key=lambda x: x[0])
            
            # Check if the compounds match exactly
            if len(input_comp) != len(scl_comp):
                continue
                
            match = True
            for (in_comp, in_frac), (scl_c, scl_f) in zip(input_comp, scl_comp):
                if in_comp != scl_c or abs(in_frac - scl_f) > tol:
                    match = False
                    break
            
            if match:
                print(f"Found matching SCL data for composition: {row['Composition']}")
                return row
                
        except Exception as e:
            print(f"Error processing SCL row {row.name}: {e}")
            continue
    
    print(f"No matching SCL data found for composition: {'-'.join([f'{f}{c}' for c, f in zip(compound_input, mol_fracs)])}")
    return None

def split_compound_to_pair(compound):
    """Split a compound into its element components and create a pair string.
    
    Args:
        compound: Compound name (e.g., 'NaCl', 'UCl3')
        
    Returns:
        str: Element pair string (e.g., 'Na-Cl', 'U-Cl')
    """
    # Handle common compound patterns
    if not isinstance(compound, str):
        return None
        
    # Remove any numbers from the compound name
    clean_compound = ''.join([c for c in compound if not c.isdigit()])
    
    # Handle common compound patterns
    if clean_compound in ['LiF', 'NaF', 'KF', 'RbF', 'CsF', 'BeF2', 'MgF2', 'CaF2', 'SrF2', 'BaF2', 'AlF3']:
        return f"{clean_compound[0]}-{clean_compound[1:]}"
    elif clean_compound in ['NaCl', 'KCl', 'RbCl', 'CsCl', 'NaBr', 'KBr', 'NaI', 'KI']:
        return f"{clean_compound[:2]}-{clean_compound[2:]}"
    elif clean_compound in ['UCl3', 'PuCl3', 'LaCl3', 'CeCl3', 'PrCl3', 'NdCl3', 'GdCl3']:
        return f"{clean_compound[0]}-{clean_compound[1:3]}"
    elif clean_compound in ['ZrF4', 'ThF4', 'UF4']:
        return f"{clean_compound[:2]}-{clean_compound[2]}"
    
    # Default case - try to split at the first lowercase letter after an uppercase letter
    for i in range(1, len(clean_compound)):
        if clean_compound[i].islower() and i+1 < len(clean_compound) and clean_compound[i+1].isupper():
            return f"{clean_compound[:i+1]}-{clean_compound[i+1:]}"
    
    # If no pattern matches, return the original compound
    return clean_compound

def get_scl_properties(compound_input, mol_fracs, scl_data):
    """
    Extract SCL properties for the given composition in the same order as compound_input.
    
    Args:
        compound_input: List of compound names (e.g., ['LiF', 'NaF', 'KF'])
        mol_fracs: List of mole fractions corresponding to compound_input
        scl_data: Dictionary containing SCL data from find_matching_scl_row()
        
    Returns:
        dict: Dictionary containing arrays of properties for each pair in compound_input order
    """
    # Initialize output arrays
    num_compounds = len(compound_input)
    pair_scl = np.zeros(num_compounds)
    pair_peak_y = np.zeros(num_compounds)
    pair_min_y = np.zeros(num_compounds)
    
    # If no SCL data, return zeros
    if not scl_data:
        return {
            'pair_scl': pair_scl,
            'pair_peak_y': pair_peak_y,
            'pair_min_y': pair_min_y
        }
    
    # Generate all possible pairs from compound_input
    pairs = []
    for i in range(num_compounds):
        for j in range(i, num_compounds):
            comp1 = compound_input[i].split('_')[0]  # Remove any _suffix
            comp2 = compound_input[j].split('_')[0]  # Remove any _suffix
            pair = f"{comp1}-{comp2}" if i != j else comp1
            pairs.append((i, j, pair))
    
    # Extract SCL data for each pair
    for i, j, pair in pairs:
        # Look for matching pair in SCL data
        for pair_idx in range(1, 7):  # Check for up to 6 pairs
            pair_col = f'Pair {pair_idx} Label'
            if pair_col in scl_data and pd.notna(scl_data[pair_col]) and scl_data[pair_col] == pair:
                scl_col = f'Pair {pair_idx} SCL_i (A)'
                peak_col = f'Pair {pair_idx} Peak Y'
                min_col = f'Pair {pair_idx} Min Y'
                
                if scl_col in scl_data and pd.notna(scl_data[scl_col]):
                    pair_scl[i] = float(scl_data[scl_col]) * 1e-10  # Convert A to m
                    pair_peak_y[i] = float(scl_data[peak_col]) if peak_col in scl_data and pd.notna(scl_data[peak_col]) else 0.0
                    pair_min_y[i] = float(scl_data[min_col]) if min_col in scl_data and pd.notna(scl_data[min_col]) else 0.0
                
                # For diagonal pairs (i==j), we're done
                if i == j:
                    break
                    
                # For off-diagonal pairs, also set the symmetric element
                if j < num_compounds:  # Safety check
                    pair_scl[j] = pair_scl[i]
                    pair_peak_y[j] = pair_peak_y[i]
                    pair_min_y[j] = pair_min_y[i]
                break
    
    return {
        'pair_scl': pair_scl,
        'pair_peak_y': pair_peak_y,
        'pair_min_y': pair_min_y
    }

def get_weighted_scl(compound_input, mol_fracs, scl_data):
    """
    Calculate weighted SCL for the given composition.
    
    Args:
        compound_input: List of compound names (e.g., ['LiF', 'NaF', 'KF'])
        mol_fracs: List of mole fractions corresponding to compound_input
        scl_data: List of SCL data entries from load_scl_data()
    """
    print(f"\nCalculating weighted SCL for composition: {compound_input} with mole fractions: {mol_fracs}")
    print(f"Number of SCL data entries: {len(scl_data)}")
    
    # Normalize input mole fractions to sum to 1
    total = sum(mol_fracs)
    if abs(total - 1.0) > 1e-6:
        print(f"Normalizing mole fractions from {mol_fracs} to {[f/total for f in mol_fracs]}")
        mol_fracs = [f/total for f in mol_fracs]
    
    # Find matching composition in SCL data (allowing for different order and small fraction differences)
    matched_entry = None
    print("\nSearching for matching composition in SCL data...")
    for i, entry in enumerate(scl_data):
        print(f"\nChecking SCL entry {i+1}:")
        print(f"  Components: {entry['components']}")
        print(f"  Fractions: {entry['fractions']}")
        print(f"  Available pairs: {list(entry['pairs'].keys())}")
        
        if len(entry['components']) == len(compound_input) and \
           set(entry['components']) == set(compound_input):
            print("  Components match! Checking fractions...")
            if compositions_match(entry['components'], entry['fractions'], 
                                compound_input, mol_fracs, tol=0.01):
                print("  Fractions match within tolerance! Using this entry.")
                matched_entry = entry
                break
            else:
                print(f"  Fractions don't match: {entry['fractions']} vs {mol_fracs}")
        else:
            print("  Components don't match")
    
    if not matched_entry:
        print(f"\nERROR: No SCL data found for composition: {compound_input} with fractions {mol_fracs}")
        print("Available SCL entries:")
        for i, entry in enumerate(scl_data):
            print(f"  {i+1}. {entry['components']} with fractions {entry['fractions']}")
        return None, None
    
    # If no pair data is available, return the average SCL
    if not matched_entry['pairs']:
        print(f"\nWARNING: No pair data available for composition: {compound_input}")
        print(f"Using average SCL: {matched_entry['avg_scl']:.3f} Å")
        return matched_entry['avg_scl'] * 1e-10, None  # Convert from A to m
    
    # Calculate weighted average of SCL_i values
    total_weight = 0
    weighted_sum = 0
    pair_scls = {}
    
    # Get the order of compounds as they appear in the SCL data
    # This ensures we match the pair ordering from the original data
    comp_order = matched_entry['components']
    frac_order = matched_entry['fractions']
    
    # For each pair in the composition
    print("\nCalculating weighted SCL for pairs:")
    for i in range(len(comp_order)):
        for j in range(i, len(comp_order)):
            # Create pair string using the same logic as in load_scl_data
            comp1 = comp_order[i].split('_')[0]  # Remove any _suffix
            comp2 = comp_order[j].split('_')[0]  # Remove any _suffix
            
            # Split compounds into element pairs using the helper function
            pair1 = split_compound_to_pair(comp1)
            pair2 = split_compound_to_pair(comp2)
            
            # Create pair string (e.g., 'Na-Cl')
            pair = f"{pair1}-{pair2}" if pair1 != pair2 else pair1
            
            # Look up SCL_i for this pair
            if pair in matched_entry['pairs']:
                scl_i = matched_entry['pairs'][pair]
                # Weight by mole fractions from the SCL data
                weight = frac_order[i] * frac_order[j] if i != j else frac_order[i] * frac_order[j]
                weighted_sum += scl_i * weight
                total_weight += weight
                pair_scls[pair] = scl_i
                print(f"  {comp_order[i]}-{comp_order[j]} → {pair}: SCL={scl_i:.3f} Å, weight={weight:.4f}, contribution={scl_i * weight:.4f}")
            else:
                print(f"  WARNING: No SCL data found for pair: {pair} (from {comp_order[i]}-{comp_order[j]})")
    
    if total_weight > 0:
        avg_scl = (weighted_sum / total_weight) * 1e-10  # Convert from A to m
        print(f"\nWeighted SCL calculation:")
        print(f"  Total weighted sum: {weighted_sum:.6f}")
        print(f"  Total weight: {total_weight:.6f}")
        print(f"  Weighted average: {weighted_sum/total_weight:.3f} Å")
        print(f"  Converted to meters: {avg_scl:.3e} m")
        print(f"\nUsing weighted average SCL of {weighted_sum/total_weight:.3f} Å ({avg_scl:.3e} m) for composition: {compound_input}")
        print(f"Mole fractions: {mol_fracs}")
        print(f"Individual pair SCLs: {pair_scls}")
        return avg_scl, pair_scls
    
    # Fall back to average SCL if no pairs were found
    print(f"\nWARNING: No valid pairs found for SCL calculation. Using average SCL as fallback.")
    print(f"Average SCL: {matched_entry['avg_scl']:.3f} Å")
    return matched_entry['avg_scl'] * 1e-10, None  # Convert from A to m
    return matched_entry['avg_scl'] * 1e-10, None  # Convert from A to m

def parse_composition(comp_str):
    """Parse composition string like '0.465LiF-0.115NaF-0.42KF' into components and fractions."""
    components = []
    fractions = []
    
    # Split by '-' and process each component
    for part in comp_str.split('-'):
        # Find the index where the compound name starts
        comp_start = 0
        for i, c in enumerate(part):
            if c.isalpha():
                comp_start = i
                break
        
        # Extract fraction and compound
        frac = float(part[:comp_start])
        comp = part[comp_start:]
        
        components.append(comp)
        fractions.append(frac)
    
    # Normalize fractions to sum to 1 (in case of rounding errors)
    total = sum(fractions)
    if abs(total - 1.0) > 1e-6:  # Only normalize if there's a significant difference
        fractions = [f/total for f in fractions]
    
    return components, fractions

def compositions_match(comp1, fracs1, comp2, fracs2, tol=0.01):
    """Check if two compositions match (allowing for different order and small fraction differences)."""
    if set(comp1) != set(comp2):
        return False
    
    # Create a mapping of compound to fraction for both compositions
    comp_frac1 = {c: f for c, f in zip(comp1, fracs1)}
    comp_frac2 = {c: f for c, f in zip(comp2, fracs2)}
    
    # Check all fractions are within tolerance
    for comp in comp_frac1:
        if abs(comp_frac1[comp] - comp_frac2[comp]) > tol:
            return False
    
    return True

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
            sound_velocity_mix = [3400, 0]
        elif Compound == ['CaCl2']:
            sound_velocity_mix = [2200, 0]
    else:
        density_mix = [0,0]
        sound_velocity_mix = [0,0]
        specific_heat_mix = [0,0]

    return r_ac_mix, density_mix, sound_velocity_mix, specific_heat_mix

def functionlibrary():
    functions = {
    'Gheribi-KT24' : Kinetic_Theory24,
    'Gheribi-KT24, Mix Data' : Kinetic_Theory24_Mix,
    'Zhao-PGM' : Zhao_PGM,
    'Zhao-PGM, Mix Data' : Zhao_PGM_Mix,
    'Present Model' : GECM,
    'Present Model, Mix Data' : GECM_Mix,
    'Ideal' : Ideal,
    }
    return functions

    