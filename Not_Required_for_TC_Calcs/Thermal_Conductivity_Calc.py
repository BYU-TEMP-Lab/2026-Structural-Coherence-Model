import numpy as np
from numpy import array
from matplotlib import pyplot as plt
import math
from scipy.optimize import curve_fit
import xlrd
import xlwt
import pandas as pd

#Import Excel spreadsheet as a data frame
df_read = pd.read_excel('C:/Users/4un/Code/TC_Comparison_Code/TC_Compound_Data.xlsx')


def calc_TC(df,compound_input,mol_fracs,T_melt, V_m=0, C_p=0, alpha=0):

    #Define smaller dataframes of each variable contained in the spreadsheet
    Compound_df=df_read['Compound']
    M_i_df=df_read['M_i (g/mol)']
    T_i_m_df=df_read['T_i_m (K)']
    V_i_m_df=df_read['V_i_m (m^3/mol)']
    alpha_i_df=df_read['alpha_i (K^-1)']
    C_i_0_df=df_read['C_i_0 (m/s)']
    C_i_p_df=df_read['C_i_p (J/mol/K)']
    r_c_df=df_read['r_c (m)']
    r_a_df=df_read['r_a (m)']
    rho_i_m_df=df_read['rho_m (g/m^3)']

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
    alpha_i=np.zeros(len(compound_input))
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
        alpha_i[i]=alpha_i_df[indices[i]]
        C_i_0[i]=C_i_0_df[indices[i]]
        C_i_p[i]=C_i_p_df[indices[i]]
        r_c[i]=r_c_df[indices[i]]
        r_a[i]=r_a_df[indices[i]]
        rho_i_m[i]=rho_i_m_df[indices[i]]

    #Calculate Gruneisen parameter
    gamma_i_m=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        gamma_i_m[i]=((M_i[i]*alpha_i[i]*C_i_0[i]**2)/C_i_p[i])*(1/1000)

    #initialize a temperature range
    T=np.linspace(T_melt,1500,num=100)

    # #calculate constant volume heat capacities
    # C_i_v=np.zeros((len(compound_input),len(T)))
    # for i in range(len(compound_input)):
    #     for j in range(len(T)):
    #         C_i_v[i][j]=C_i_p[i]/(1+alpha_i[i]*gamma_i_m[i]*T[j])

    
    
    
    #Calculate number of ions per compound
    n_i=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        n_i[i]=(sum(1 for c in Compound[i] if c.isupper())+ sum(1 for c in Compound[i] if c.isnumeric()))

 
    #Calculate compound thermal conductivites at melting points
    lambda_i_m=np.zeros(len(compound_input))
    for i in range(len(compound_input)):
        lambda_i_m[i]=4.33*   ( (  C_i_0[i]*   (C_i_p[i]/(1+alpha_i[i]*gamma_i_m[i]*T_i_m[i]))    )   /   (3*n_i[i]*V_i_m[i])    )    *     (r_a[i]+r_c[i])


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
        Tau_i_m[i]=(C_i_p[i]/(1+alpha_i[i]*gamma_i_m[i]*T_i_m[i]))*C_i_0[i]*(1/n_i[i])


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
            alpha_mix_m=alpha_mix_m+phi_i_m[i]*alpha_i[i]
    else:alpha_mix_m=alpha
    
    

    #Calculate mixture molecular weight
    M_mix=0
    for i in range(len(compound_input)):
        M_mix=M_mix+mol_fracs[i]*M_i[i]
        

    #Calculate mixture gamma
    gamma_mix_m=M_mix*alpha_mix_m*C_0_mix_m**2*(1/C_p_mix)*(1/1000)
    

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
        V_mix_m=V_mix_m*1E-6
    else: V_mix_m=V_m
   

    print("V of the mix is: " + str(V_mix_m))

    sum_term=0
    for i in range(len(compound_input)):
        sum_term=sum_term+mol_fracs[i]*((M_i[i]/M_mix)-1)**2
    delta_mix_m=0.4872*(lambda_ideal_mix_m/(1.380649e-23*C_0_mix_m*(n_mix*6.0221408e23*(1/V_mix_m))**(2/3)))*sum_term
    
    
    plt.xlim([800,1500])
    plt.ylim([0,1.5])
    

    lambda_mix=lambda_ideal_mix_m*(1-delta_mix_m)
    print(lambda_mix)

    lambda_mix_T=np.zeros(len(T))
    for j in range(len(T)):
        lambda_mix_T[j]=lambda_mix*(1-alpha_mix_m*(gamma_mix_m+(1/3))*(T[j]-T_melt))



    return(T,lambda_mix_T)


#NaCl-UCl3 mp:790

T_kin, lambda_kin=calc_TC(df_read,['LiF','BeF2'],[0.73,0.27],520+273)
T_kin2, lambda_kin2=calc_TC(df_read,['LiF','BeF2','UF4'],[0.725,0.265,0.01],520+273)
T_kin3, lambda_kin3=calc_TC(df_read,['LiF','BeF2','UF4'],[0.720,0.260,0.02],520+273)

plt.xlabel('Temperature (K)')
plt.ylabel('Thermal Conductivity (W/m K)')
plt.plot(T_kin,lambda_kin,label='0.73LiF–0.27BeF2 (kinetic theory)',color='blue')
plt.plot(T_kin2,lambda_kin2,label='(0.73LiF–0.27BeF2) + 1% UF4 (kinetic theory)',color='green')
plt.plot(T_kin3,lambda_kin3,label='(0.73LiF–0.27BeF2) + 2% UF4 (kinetic theory)',color='purple')
plt.scatter([884,933,979,1032,1078],[0.974,1.03,1.028,1.058,1.08],color='blue',label='0.73LiF–0.27BeF2 (Bobrova)')
plt.scatter([882,933,982,1031,1078],[0.89,0.946,0.945,0.999,1.05],color='green',label='(0.73LiF–0.27BeF2) + 1% UF4 (Bobrova)')
plt.scatter([832,883,931,979,1032,1073],[0.79,0.847,0.898,0.945,0.947,0.998],color='purple',label='(0.73LiF–0.27BeF2) + 2% UF4 (Bobrova)')

plt.legend()

plt.show()

#Bobrova: https://link.springer.com/article/10.1134/S0036029523020039

T_kin, lambda_kin=calc_TC(df_read,['LiF','NaF','ZrF4'],[0.55,0.22,0.23],845)
plt.plot(T_kin,lambda_kin,label='0.55LiF–0.22NaF-0.23ZrF4 (kinetic theory)',color='orange')
plt.scatter([917,980,997],[0.59,0.598,0.63],color='orange',label='0.55LiF–0.22NaF-0.23ZrF4 (Khokhlov)')

T_kin, lambda_kin=calc_TC(df_read,['LiF','ThF4'],[0.77,0.23],841)
plt.plot(T_kin,lambda_kin,label='0.77LiF–0.23ThF4 (kinetic theory)',color='pink')
plt.scatter([1081,1083],[0.83,0.83],color='pink',label='0.77LiF–0.23ThF4 (Khokhlov)')

T_kin, lambda_kin=calc_TC(df_read,['LiF','NaF','KF'],[0.465,0.115,0.42],727)
plt.plot(T_kin,lambda_kin,label='0.465LiF–0.115NaF-0.42KF (kinetic theory)',color='cyan')
plt.scatter([790,825,887,954,1032,1079],[0.82,0.84,0.87,0.91,0.93,0.96],color='cyan',label='0.465LiF–0.115NaF-0.42KF (Khokhlov)')

plt.legend()

plt.show()

#Khokhlov: https://www.sciencedirect.com/science/article/abs/pii/S002231151100002X?via%3Dihub

T_kin, lambda_kin=calc_TC(df_read,['LiF','BeF2','UF4','ThF4'],[0.675,0.20,0.005,0.12],475+273)
plt.plot(T_kin,lambda_kin,label='0.675LiF–0.20BeF2-0.005UF4-0.12ThF4 (kinetic theory)',color='black')
plt.scatter([811,817,855,889,945,1003,1079,1142],[1.21,1.18,1.23,1.23,1.22,1.14,1.10,1.05],color='black',label='0.675LiF–0.20BeF2-0.005UF4-0.12ThF4 (Cooke)')

#https://www.semanticscholar.org/paper/Transport-Properties-of-Molten-salt-Reactor-Fuel-of-Afonichkin-Khokhlov/b9ee1b2fcba1c701f0dae7fe2f0229a60dcba07e


T_kin, lambda_kin=calc_TC(df_read,['LiF','BeF2','ZrF4','UF4'],[0.712,0.23,0.05,0.008],440+273)
plt.plot(T_kin,lambda_kin,label='0.712LiF–0.23BeF2-0.05ZrF4-0.008UF4 (kinetic theory)',color='grey')
plt.scatter([787,793,891,924,965,973,1050,1043,1049,1057,1125],[1.19,1.15,1.32,1.34,1.28,1.35,1.23,1.22,1.19,1.25,1.14],color='grey',label='0.712LiF–0.23BeF2-0.05ZrF4-0.008UF4 (Cooke)')


T_kin, lambda_kin=calc_TC(df_read,['LiF','BeF2'],[0.66,0.34],475+273)
plt.plot(T_kin,lambda_kin,label='0.66LiF–0.34BeF2(kinetic theory)',color='brown')
plt.scatter([823,826,933,1026,1137],[0.91,0.95,1.03,1.02,0.97],color='brown',label='0.66LiF–0.34BeF2 (Cooke)')

#Cooke: https://moltensalt.org/references/static/downloads/pdf/ORNL-4449.pdf




#Grimes: https://moltensalt.org/references/static/downloads/pdf/FFR_chap12.pdf

#https://pdf.sciencedirectassets.com/271445/1-s2.0-S0306454917X0008X/1-s2.0-S0306454917301391/am.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEEYaCXVzLWVhc3QtMSJHMEUCIHU7eqATtjtL175IsvrqvoZeuLFHCyUiu7YvgTK2BvBhAiEAlZKrkYJ0%2FgzF2aiFOd1%2BlkMXhX%2BlKzXFr5cLzZzDEe0qvAUIrv%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FARAFGgwwNTkwMDM1NDY4NjUiDO%2FH%2BsDRvNkW2J3VTCqQBXh2244XholpjqaJyCh43ZPJnIMeWfUiGC5d0oY1PTIWQsO%2FjBgNwLG9aNsdX0PH45%2B4gkVzFqSF9X20IWjIBN8H4Hptp0%2Fpks77rQt3VRVu3lYAUHrqZsIEDi9D8jRiTCGczQ7wavtYwuxwbtQGoTT1WQElvyArFe0OcRONwzEX5vc86xwVA47Vc3l2KwYyNvoQ4zjZqGAtSpq2kQxphJUzHfS%2Bl3wN1CzTw36mglQCVZh4jPd0xaLtcZ3VYG2dMK8v4QdOvVxshT43F5IGhmWGf3YNdxceYwiRDDJYsn9ERQNqxx0KU1dBn6Gctz2ZKa6LDtvHWDv%2FYJMikR0Sl9jQbelyPau%2FuGjPfh3xv2QGoz8Ea6%2FnGaIu4YNl50ydDMq9Kedq5RM7b1GU45nnXDxOKpsXDjwJdRpdcw8ja6G1hXcH3iVeaGE0bwFsR%2Blu8yNk0o40Oowe12sCrvhceeL4rOQezKT8AWxm4ZIx%2FHccM8y8wjFpeuztxrDMq4nl%2F3%2FIptj0AGT53kyjBSOi7zMccw8yj7UvYYln%2Fc2PE66sc9u5O5GFq7bXtE%2BDo6XH83qunh6VDssJESJxT6nEOkfZk8j9wF%2Fm9LLjg8TcAh7m0zJ8t1JYDOWS7QBWT7TluadBwtBoFrhEEC8i%2FY%2Bb2W7zBQDts8mRgICtP9HYNhM1xCyJWJHQIECfgeE0%2Fyu2Vz7b8nIXp55zylaawWzslpiG77Kc3yrGG7G0CCgwGjcqxaI5ol1%2F6X7MBmZGd36ZfSL0T3ZTTS3VApBTG0fgQ1ihiJ9p3svQoNdzQ%2BFJE9PLfBIvv0FtZjbWftSFO7yIqjt5GngoU9pCpfb5sKOoh6x3P91kaTadDhgVKc0v6EqgMLTtyKsGOrEB8DWdghHsdaWDSJTgY7SCEzxuSzoN3KrOZTJ8iUyh1vjCchho71d8akooAtH7QC3MPqHLbIEsTx6XshzeRVoUbpWB7FQl3Xl1UlU7CF4vSQYlKhCWItfyOa%2FB6zd%2B4JgSKRYzbzzU07SEGvwnhrw0KM8%2Bch63KvIyXxnFpTRPHFR%2ByHyg%2F4L1K7j112HQJUGH7hyJydfLtPjq4gOGkr5i12CBMlYvS%2B9WMySODmzWOs13&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20231207T213825Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTYTZ73QENJ%2F20231207%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=bbd9fd34fb7ebd871cc17054c32fe3ef79f2e4dafb9412f5edea8f1afcedae6f&hash=15cb5e9e1a6b7d4cefb21336ef5af0ad880f5a44fdf068d026601a6c7d6e4fe3&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S0306454917301391&tid=pdf-89061768-0923-452c-a673-bbfbee283a77&sid=b5a9bac23b5ee44f429891c1aa725b6d5658gxrqa&type=client

plt.legend()

plt.show()

T_kin, lambda_kin=calc_TC(df_read,['NaF','ZrF4'],[0.81,0.19],1023)
plt.plot(T_kin,lambda_kin,label='0.81NaF–0.19ZrF4 (kinetic theory)',color='red')

T_kin, lambda_kin=calc_TC(df_read,['NaF','ZrF4'],[0.422,0.588],842)
plt.plot(T_kin,lambda_kin,label='0.422NaF–0.588ZrF4 (kinetic theory)',color='blue')

T_kin, lambda_kin=calc_TC(df_read,['LiF'],[1],1120)
plt.plot(T_kin,lambda_kin,label='LiF (kinetic theory)',color='orange')

T_kin, lambda_kin=calc_TC(df_read,['ThF4'],[1],1400)
plt.plot(T_kin,lambda_kin,label='ThF4 (kinetic theory)',color='green')

plt.legend()

plt.show()


T_kin, lambda_kin=calc_TC(df_read,['LiF','BeF2','ThF4'],[0.717,0.160,0.123],730)
plt.plot(T_kin,lambda_kin,label='0.717LiF–0.16BeF2-0.123ThF4 (kinetic theory)',color='red')

plt.legend()

plt.show()