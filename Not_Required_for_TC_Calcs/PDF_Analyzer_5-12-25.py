import pandas as pd
import numpy as np
import re
import os
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from datetime import datetime
from scipy.special import expit  # Sigmoid function for smooth transitions
from scipy.signal import find_peaks
from mendeleev import element
import csv
from scipy.integrate import quad

def read_and_process_csv(file_path):
    # Read the CSV file
    df = pd.read_csv(file_path, header=None)
    
    ion_pairs = df.iloc[0, ::2]  # Get ion pair labels from the first row, every 3rd cell
    data_dict = {ion: {} for ion in ion_pairs}
    
    for i, ion in enumerate(ion_pairs):
        base_idx = i * 2
        x_label = df.iloc[1, base_idx]
        y_label = df.iloc[1, base_idx + 1]
        x_values = df.iloc[2:, base_idx].dropna().astype(float).values
        y_values = df.iloc[2:, base_idx + 1].dropna().astype(float).values
        
        # Extend x_values and y_values
        x_extended = np.linspace(0, x_values[0], num=int(x_values[0]) + 1)
        y_extended = np.zeros_like(x_extended)  # Zero y-values for the extended x-range
        
        # Combine the extended x and y with the original data
        x_combined = np.concatenate((x_extended, x_values))
        y_combined = np.concatenate((y_extended, y_values))
        
        data_dict[ion]['x'] = x_combined
        data_dict[ion]['y'] = y_combined
    
    return data_dict

def parse_composition(composition_str):
    # Example: "0.63NaCl-0.37UCl3"
    components = composition_str.split('-')
    composition = {}
    for component in components:
        match = re.match(r"([0-9.]+)([A-Za-z0-9]+)", component)
        if match:
            fraction, salt = match.groups()
            composition[salt] = float(fraction)
    # Split the input string into fractions and compounds
    compound_split = re.findall(r'([0-9.]+)([A-Za-z0-9]+)', composition_str)

    # Initialize dictionary for ion counts
    ion_counts = {}

    # Iterate over each compound
    for fraction, compound in compound_split:
        elements = re.findall(r'([A-Z][a-z]*)([0-9]*)', compound)
        ion_count = {}
        
        # Count each ion in the compound
        for element, count in elements:
            count = int(count) if count else 1  # Default to 1 if no count is given
            if element in ion_count:
                ion_count[element] += count
            else:
                ion_count[element] = count
        
        # Add to the ion_counts dictionary
        ion_counts[compound] = ion_count

    return composition, ion_counts

def calculate_weights(composition, ion_counts, splines):
    # Initialize a dictionary to hold the total molar concentration of each element
    element_concentration = {}

    # Calculate the total contribution of each element
    for salt, fraction in composition.items():
        for element, count in ion_counts[salt].items():
            if element in element_concentration:
                element_concentration[element] += fraction * count
            else:
                element_concentration[element] = fraction * count

    # Calculate the sum of all element concentrations
    total_concentration = sum(element_concentration.values())

    # Normalize to get the relative molar concentration
    rel_conc = {element: concentration / total_concentration 
                            for element, concentration in element_concentration.items()}

    weights = {}

    pairs_with_splines = []
    for ion1, conc1 in rel_conc.items():
        for ion2, conc2 in rel_conc.items():
            ion_pair = f"{ion1}-{ion2}"
            rev_ion_pair = f"{ion2}-{ion1}"
            if ion_pair in splines or rev_ion_pair in splines:
                pairs_with_splines.append(1)
            else:
                pairs_with_splines.append(0)
            if ion_pair in weights or rev_ion_pair in weights:
                if ion_pair == rev_ion_pair:
                    pass  
                else:
                    ion_pair = rev_ion_pair
            else:
                weights[ion_pair] =+ (conc1 * conc2)

    sumweights = sum(weights.values())
    weights = {key: value/sumweights for key, value in weights.items()}
    
    return weights, pairs_with_splines

def interpolate_data(data_dict,num_points):
    splines = {}
    interpolated_splines = {}

    x_min = min(min(data["x"]) for data in data_dict.values())
    x_max = min(max(data["x"]) for data in data_dict.values())

    x_new = np.linspace(0, x_max, num_points)

    for ion, data in data_dict.items():
        x = data['x']
        y = data['y']

        # Create the spline function
        spline = interp1d(x, y, kind='linear', bounds_error=False, fill_value=0)
        splines[ion] = spline  # Store original spline

        # Generate a fixed number of x points and compute interpolated values
        interpolated_splines[ion] = interp1d(x_new, spline(x_new), kind='linear', bounds_error=False, fill_value=0)

    return interpolated_splines, x_new


def average_x_position(x, y, x_min, x_max):
    # Ensure the inputs are numpy arrays for easier slicing
    x = np.array(x)
    y = np.array(y)
    
    # Mask the x and y values to only consider the range between x_min and x_max
    mask = (x >= x_min) & (x <= x_max)
    x_filtered = x[mask]
    y_filtered = y[mask]
    
    # Calculate the total area under the curve (integral)
    area = np.trapz(y_filtered, x_filtered)
    
    # Calculate the weighted average x-position
    weighted_x = np.trapz(x_filtered * y_filtered, x_filtered) / area
    
    return weighted_x



def identify_cat_an(ion_pair):   # Detect if cation-anion pair
    element1, element2 = ion_pair.split('-')
    el1 = element(element1)
    el2 = element(element2)
    el1_states = el1.oxistates
    el2_states = el2.oxistates
    blank = "none"
    
   
    if el1_states[0] > 0 and el2_states[0] < 0:    #is_cation_anion:
        return 1, element1, element2
    elif el1_states[0] < 0 and el2_states[0] > 0:  # is_anion_cation:
        return 1, element2, element1
    elif el1_states[0] > 0 and el2_states[0] > 0:  # is_cation_cation:
        return 2, element1, element2
    else:                                          # is_anion_anion:
        return 3, element1, element2

def Save_fig(title):
    # Variables for the file name and folder
    folder_path = 'J:\Research\PDF_Analysis_kModel\SCL_Figures'
    file_name = f'{title}.png'  # Combine string and variable

    # Full path to save the file
    save_path = os.path.join(folder_path, file_name)

    # Save the plot
    plt.savefig(save_path)
    return

def cdf(x_val,spline, total_integral):              # Define the Cumulative Distribution function
    integral, _ = quad(spline, spline.x[0], x_val)
    return integral / total_integral


def Nu_integral(splines, weights, pairs_with_splines, composition_str,source,weighttype,MFP_calc,num_points,x_range):

    nu = 0.5

    for pair, spline in splines.items():
    
        P_pair = np.zeros(num_points)
        beta_pair = np.zeros(num_points)
        beta_integral_pair = np.zeros(num_points)
        sum_y_arrays = np.zeros(num_points)
        mix_c = 0
        
        b_RT = 0
        b_NI = 0
        b_SRO = 0


        if pair in weights:
            pass
        else:
            elements = pair.split('-')
            pair = '-'.join(reversed(elements))
        cat_an_weight = weights.get(pair)

        # Identify if cation-anion pair
        catan, cation, anion = identify_cat_an(pair)
        if catan == 1:
            y_spline = spline(x_range)
            First_peak_index_ca = np.argmax(y_spline)
            r_c_ca = x_range[First_peak_index_ca]
            w_c_ca = y_spline[First_peak_index_ca]
            c_i_ca = cat_an_weight * w_c_ca * r_c_ca**-1
            print("Cat-An Pair: ",pair)
        else:
            continue  # Only stay in loop if it's a cat-an pair


        # Weight and plot splines
        for pair_iter, spline_iter in splines.items():
            if pair_iter in weights:
                pass
            else:
                elements = pair_iter.split('-')
                pair_iter = '-'.join(reversed(elements))
            concentration = weights.get(pair_iter)
            y_spline_iter = spline_iter(x_range)
            avg_conc = np.linspace(concentration,concentration,num_points)
            line1, = plt.plot(x_range, y_spline_iter, label=pair_iter)
            previous_color = line1.get_color()
            plt.plot(x_range, avg_conc, linestyle=':', color=previous_color)

            catan_iter, cation_iter, anion_iter = identify_cat_an(pair_iter)
            if catan_iter == 3:
                pass
            else:
                if cation == cation_iter or cation == anion_iter:
                    sum_y_arrays += y_spline_iter    # all possible arrays including principal cation besides anion-anion

            if catan_iter == 1:
                First_peak_index = np.argmax(y_spline_iter)
                r_c = x_range[First_peak_index]
                w_c = y_spline_iter[First_peak_index]
                c_i = concentration * w_c * r_c**-1
                mix_c += c_i
                

        if mix_c == 0:
            b_RT = 0
        else:
            b_RT = 1 - c_i_ca/mix_c

        try:
            cat_an_pair = cation + '-' + anion
            cat_an_spline = splines[cat_an_pair]
        except KeyError:
            cat_an_pair = anion + '-' + cation
            cat_an_spline = splines[cat_an_pair]
        cat_an_y_array = cat_an_spline(x_range)

        try:
            cat_cat_pair = cation + '-' + cation    # The corresponding cat-cat spline is essential to know transition pointsn (actually now we just use cat-an spline crossing concentration)
            cat_cat_spline = splines[cat_cat_pair]
            cat_cat_y_array = cat_cat_spline(x_range)
            cat_cat_weight = weights.get(cat_cat_pair)
        except KeyError:
            print("No cation-cation spline data available. Estimating based on cation-anion spline (inverse when below concentration).")
            weight = weights.get(cat_cat_pair,100)
            cat_cat_y_array = -cat_an_y_array*weight
            for i in range(len(x_range)):
                if cat_an_y_array[i] > weight:
                    cat_cat_y_array[i] = 0
            cat_cat_weight = weights.get(cat_cat_pair)

        first_peak = True
        for i in range(len(x_range)):

            if cat_an_y_array[i] > cat_an_weight or first_peak == True:
                b_SRO_denom = np.max(cat_an_y_array[0:i+1])
                if b_SRO_denom == 0:
                    b_SRO = 0
                else:
                    b_SRO = 1 - cat_an_y_array[i]/b_SRO_denom
                b_NI_denom = (sum_y_arrays[i])
                if b_NI_denom == 0:
                    b_NI = 0
                else:
                    b_NI = 1 - cat_an_y_array[i]/b_NI_denom 
                
            elif cat_an_y_array[i] < cat_an_weight:
                b_SRO_denom = np.max(cat_an_y_array[0:i+1])
                if b_SRO_denom == 0:
                    b_SRO = 0
                else:
                    b_SRO = 1 - cat_cat_y_array[i]/b_SRO_denom
                b_NI_denom = (sum_y_arrays[i])
                if b_NI_denom == 0:
                    b_NI = 0
                else:
                    b_NI = 1 - cat_cat_y_array[i]/b_NI_denom
                

            if w_c_ca == cat_an_y_array[i]:
                first_peak = False

            beta_pair[i] = ((b_RT/(1-b_RT)) + b_NI/(1-b_NI) + b_SRO/(1-b_SRO))/3


        for i in range(len(x_range)):
            beta_x = x_range[:i+1] 
            beta_y = beta_pair[:i+1]
            beta_integral_pair[i] = np.trapz(beta_y,beta_x)

            P_pair[i] = np.exp(-beta_integral_pair[i])

        # Expectation value of Nu
        P_exp = np.trapz(P_pair*x_range,P_pair)/np.trapz(x_range,P_pair)

        # Find the x value corresponding to the closest y to y_avg
        index = np.argmin(np.abs(P_pair - P_exp))
        x_SCL_pair = x_range[index]

        #x_SCL_all += x_SCL_pair*cat_an_weight

        plt.plot(x_range, beta_y, linestyle='-.', label=fr"$\beta(r)$ - {cat_an_pair}")

        plt.plot(x_range, P_pair, linestyle='-.',  label=f"P(r) - {cat_an_pair}")

        plt.axvline(x=x_SCL_pair, color='g', linestyle='--', label=f"$\lambda_{{P(r)}} - {cat_an_pair} = {round(x_SCL_pair,2)}$")

        if MFP_calc > 0:
            plt.axvline(x=MFP_calc, color='k', linestyle='--', label=f"$\lambda_{{BC}} - {cat_an_pair} = {round(MFP_calc,2)}$")
        #plt.plot(x_combined, y_combined, label='Combined RDF',color='cornflowerblue')
            
        print(f"lambda_{{Nu(r)}} - {cat_an_pair} = {round(x_SCL_pair,5)}")
    
    # Get the current date
    current_date = datetime.now().strftime("%Y%m%d")
    save_title = f'PartPDF_CatAn_{weighttype}_{composition_str}'

    plt.xlabel('Radius (Å)')
    plt.ylabel('RDF')
    plt.title(f'Partial PDF, {weighttype} : {composition_str} ({source})', fontsize=10)
    plt.legend()
    Save_fig(save_title)
    plt.close()
    print("")
    print(f"lambda_{{Nu(r)}} = {round(x_SCL_pair,5)}")
    #print(f"P(r),min = {Nu_pair[999]:e}")
    print("")
    return


def Read_RDF(file_name,composition_str,source,MFP_calc):

    folder_name = "RDF_Plots\RDF_CSV"

    # Get the current working directory
    current_dir = os.getcwd()

    # Construct the path to the RDF_Plots folder
    rdf_plots_folder = os.path.join(current_dir, folder_name)
    file_path = os.path.join(rdf_plots_folder, file_name)
    print("###  ",composition_str, "   #############")
    
    data_dict = read_and_process_csv(file_path)

    num_points = 100
    splines,x_range = interpolate_data(data_dict,num_points)

    peaks = 'off'

    composition, ion_counts = parse_composition(composition_str)
    
    weights, pairs_with_splines = calculate_weights(composition, ion_counts,splines)

    # Weight and plot splines
    for pair, spline in splines.items():
        if pair in weights:
            pass
        else:
            elements = pair.split('-')
            pair = '-'.join(reversed(elements))
        weight = weights.get(pair, 100.0)  # Get the weight or use 1.0 as default
        avg_conc = np.linspace(weight,weight,num_points)
        y_weighted = weight * spline(x_range)
        spline.y = y_weighted

    Nu_integral(splines, weights, pairs_with_splines, composition_str,source,'Conc',MFP_calc,num_points,x_range)

    
Read_RDF('PDF_LiCl.csv',"1.0LiCl",'Walz, 2019; 878K',3.893446)   # Walz, 2019; 878K
Read_RDF('PDF_NaCl.csv',"1.0NaCl",'Walz, 2019; 1074.15K',4.352426)   # PDF 'Walz, 2019; 1074.15K, #'Andersson, 2022; 1250K')   # PDF Andersson, 1250 K, 2022
Read_RDF('PDF_KCl.csv',"1.0KCl",'Walz, 2019; 1043K',4.453593)   # Walz, 2019; 1043K
Read_RDF('PDF_LiF.csv',"1.0LiF",'Walz, 2019; 1121K',3.305228)   # Walz, 2019; 1121K
Read_RDF('PDF_NaF.csv',"1.0NaF",'Walz, 2019; 1266.15K',4.005403)   # Walz, 2019; 1266.15K
Read_RDF('PDF_KF.csv',"1.0KF",'Walz, 2019; 1131.15K',4.225170)   # Walz, 2019; 1131.15K
Read_RDF('PDF_MgCl2.csv',"1.0MgCl2",'McGreevy, 1987; 998K',4.049956)   # McGreevy, 1987; 998K
Read_RDF('PDF_CaCl2.csv',"1.0CaCl2",'McGreevy, 1987; 1093K',3.372692)   # McGreevy, 1987; 1093K
Read_RDF('PDF_SrCl2.csv',"1.0SrCl2",'McGreevy, 1987; 1198K',3.388477)   # McGreevy, 1987; 998K
Read_RDF('PDF_NaCl-UCl3.csv',"0.64NaCl-0.36UCl3",'Andersson, 2022; 1250K',2.455254)   # PDF Andersson, 2022
Read_RDF('PDF_NaCl-UCl3.csv',"0.64NaCl-0.36UCl3",'ANL-Andersson, 2022; 1250K',2.717874)   # PDF Andersson, 2022
Read_RDF('PDF_FLiBe_Grizzi.csv',"0.5LiF-0.5BeF2",'Sun, 2024; 900K',1.922959)   # PDF Sun, 900 K, 2024
Read_RDF('PDF_FLiBe_Langford.csv',"0.66LiF-0.34BeF2",'Langford, 2022; 973K',1.922959)   # PDF Langford, 2022 (Cylindrical)
Read_RDF('PDF_FLiBe_Fayfar.csv',"0.66LiF-0.34BeF2",'Fayfar, 2023; 973K',1.922959)   # PDF Fayfar
Read_RDF('PDF_FLiNa.csv',"0.6LiF-0.4NaF",'Grizzi, 2024; ``973``K',2.378058)   # PDF Grizzi, 900 K, 2024
Read_RDF('PDF_FLiNaK.csv',"0.465LiF-0.115NaF-0.42KF",'Frandsen, 2020; 873K',2.421498)   # PDF Frandsen, 940 K, 2020
Read_RDF('PDF_FMgNaK.csv',"0.345NaF-0.59KF-0.065MgF2",'Solano, 2021; 1073K',0)   # PDF Frandsen, 940 K, 2020
Read_RDF('PDF_LiF-NaF-UF4.csv',"0.5454LiF-0.3636NaF-0.091UF4",'Grizzi, 2024; 1473K',0)   # PDF Grizzi, 2024
Read_RDF('PDF_NaCl-KCl-ZnCl2_1073.csv',"0.22NaCl-0.393KCl-0.387ZnCl2",'Xi, 2024; 1073K',0)   # PDF Xi, 2024
