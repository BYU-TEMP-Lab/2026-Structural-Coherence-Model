import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import interp1d
from scipy.interpolate import CubicSpline
from scipy.interpolate import UnivariateSpline
from scipy.signal import find_peaks
from scipy.integrate import quad
import re
import periodictable as pt
from datetime import datetime
import matplotlib.lines as mlines
#from CoherenceLength import *
import csv
from scipy.signal import find_peaks
from mendeleev import element
from scipy.optimize import root_scalar


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

    for ion1, conc1 in rel_conc.items():
        for ion2, conc2 in rel_conc.items():
            ion_pair = f"{ion1}-{ion2}"
            rev_ion_pair = f"{ion2}-{ion1}"
            if ion_pair in splines or rev_ion_pair in splines:
                if ion_pair in weights or rev_ion_pair in weights:
                    if ion_pair == rev_ion_pair:
                        pass  
                    else:
                        ion_pair = rev_ion_pair
                else:
                    weights[ion_pair] =+ (conc1 * conc2)

    sumweights = sum(weights.values())
    weights = {key: value/sumweights for key, value in weights.items()}
    
    return weights


def interpolate_data(data_dict):
    splines = {}
    for ion, data in data_dict.items():
        x = data['x']
        y = data['y']
        spline = interp1d(x, y,kind='linear', bounds_error=False, fill_value=0)
        #spline = Akima1DInterpolator(x, y, axis=0)
        #spline = CubicSpline(x,y,bc_type='natural')
        splines[ion] = spline

        # Add a y=0 line up to the first x-position
    return splines


def identify_cat_an(ion_pair):   # Detect if cation-anion pair
    element1, element2 = ion_pair.split('-')
    el1 = element(element1)
    el2 = element(element2)
    el1_states = el1.oxistates
    el2_states = el2.oxistates
    
   
    if el1_states[0] > 0 and el2_states[0] < 0:    #is_cation_anion:
        return 1
    elif el1_states[0] < 0 and el2_states[0] > 0:  # is_anion_cation:
        return 1
    else:
        return 0
    

def cdf(x_val,spline, total_integral):              # Define the Cumulative Distribution function
    integral, _ = quad(spline, spline.x[0], x_val)
    return integral / total_integral


def average_x_position(spline, a, b):
    """
    Computes the average x-position of a function represented by a spline over the interval [a, b].
    
    Parameters:
    spline : callable
        A spline function from scipy.interpolate.UnivariateSpline or similar.
    a : float
        Lower bound of the integration interval.
    b : float
        Upper bound of the integration interval.
    
    Returns:
    float
        The average x-position.
    """
    numerator = quad(lambda x: x * spline(x), a, b)[0]
    denominator = quad(spline, a, b)[0]
    
    if denominator == 0:
        raise ValueError("Integral of the function is zero; cannot compute average x-position.")
    
    return numerator / denominator


def plot_cat_an_splines_weighted(splines, weights, composition_str,source,weighttype,MFP_calc):
    x_min = min(min(spline.x) for spline in splines.values())
    x_max = max(max(spline.x) for spline in splines.values())
    y_max = max(max(spline.y) for spline in splines.values())
    x_combined = np.linspace(x_min, x_max, 500)
    y_combined = np.zeros_like(x_combined)

    peaks_list = []
    peaks_list_csv = []

    salt_name = 0

    catan_sum_weights = 0

    sum_MFP = 0
    splinecount = 0

    for pair, spline in splines.items():
        # Identify if cation-anion pair
        catan = identify_cat_an(pair)  

        if catan == 0:
            continue     
        
        x_new = np.linspace(min(spline.x), max(spline.x), 500)

        if pair in weights:
            pass
        else:
            elements = pair.split('-')
            pair = '-'.join(reversed(elements))
        weight = weights.get(pair, 100.0)  # Get the weight or use 1.0 as default
        if catan == 1:
            catan_sum_weights += weight

    for pair, spline in splines.items():
        # Identify if cation-anion pair
        catan = identify_cat_an(pair)

        if catan == 0:
            continue    

        splinecount += 1
        x_new = np.linspace(min(spline.x), max(spline.x), 500)
        y_weighted = spline(x_new)

        if pair in weights:
            pass
        else:
            elements = pair.split('-')
            pair = '-'.join(reversed(elements))
        weight = weights.get(pair, 100.0)  # Get the weight or use 1.0 as default
        avg_conc = np.linspace(weight,weight,500)

        if pair in weights:
            pass
        else:
            elements = pair.split('-')
            pair = '-'.join(reversed(elements))
        weight = weights.get(pair, 100.0)  # Get the weight or use 1.0 as default
        y_weighted = weight * spline(x_new)# * 1/prefactor**2

        # Find all peaks
        peak_indices, _ = find_peaks(y_weighted)

        # Only get the correct peaks
        peaks_to_delete = []
        for loop_idx,peak_idx in enumerate(peak_indices):
            # Get rid of peaks below concentration line
            if y_weighted[peak_idx] < weight:
                peaks_to_delete.append(loop_idx)
                continue

            # Find lower bound for peak
            for j in range(peak_idx - 1, -1, -1):
                if y_weighted[j] < weight:
                    thresh_x_low = x_new[j]
                    idx_low = j
                    break
            else:
                thresh_x_low = None
            
            # Find upper bound for peak
            for j in range(peak_idx + 1, len(x_new)):
                if y_weighted[j] < weight:
                    thresh_x_up = x_new[j]
                    idx_high = j
                    break
            else:
                if y_weighted[-1]>0:
                    idx_high = j
                else:
                    thresh_x_up = None

            # Get rid of peaks below concentration line
            if thresh_x_low is None or thresh_x_up is None:
                peaks_to_delete.append(loop_idx)
                continue


            peaks_in_range_idx = [p for p in peak_indices if idx_low <= p <= idx_high]
            peaks_in_range_y = y_weighted[peaks_in_range_idx]

            if y_weighted[peak_idx] == max(peaks_in_range_y):
                pass
            else:
                peaks_to_delete.append(loop_idx)
                continue
        peak_indices = np.delete(peak_indices, peaks_to_delete)


        # Check if the spline ends in a peak (find_peaks only gets peaks with values on either side)
        if len(peak_indices) >= 2 and y_weighted[peak_indices[-1]] < y_weighted[-1]:
            peak_indices[-1] = len(y_weighted) -1

        # Calculate cumulative distribution function
        x_cdf = x_new
        # Total integral of the PDF for normalization
        total_integral = quad(spline, min(spline.x), max(spline.x))[0]
        cdf_values = [cdf(val,spline, total_integral) for val in x_cdf] 
        plt.plot(x_cdf, cdf_values, label=f'CDF,{pair}', linestyle='--')


        """Integrates (spline(x) - y_target) over the given x range."""       
        integral_about_conc, _ = quad(lambda x: spline(x) - weight, x_min, x_max)

       # integral_about_0, _ = quad(spline, x_limits[0], x_limits[-1])   

        print(integral_about_conc)


        below_conc = 0
        for i,idx in enumerate(peak_indices):
            x_at_peak = x_new[idx]
            y_at_peak = y_weighted[idx]
            peaks_list.append([pair, [x_at_peak, y_at_peak]])

            plt.plot(x_at_peak, y_at_peak,'ro')
            plt.text(x_at_peak, y_at_peak, f'{x_at_peak:.2f},{y_at_peak:.2f}')

            peak_prob = 0
            if x_at_peak > weight and below_conc == 0:
                below_conc = 1
                if len(peak_indices) > i+1:
                    y_conc = y_at_peak - weight
                    yn_conc = y_weighted[idx + 1] - weight
                    peak_prob = (y_conc - yn_conc)/y_weighted[0]
                    sum_MFP += peak_prob*x_at_peak
                else:
                    y_conc = y_at_peak - weight
                    peak_prob = (y_conc)/y_weighted[0]
                    sum_MFP += peak_prob*x_at_peak

            # Start a new line in the CSV file
            if i == len(peak_indices)-1 and splinecount == len(splines):
                with open('pair_peaks_list_catan.csv', mode='a', newline='') as file:
                    writer = csv.writer(file)
                    if catan == 1:
                        writer.writerow([composition_str, pair, i, x_at_peak, y_at_peak, weight, weight/catan_sum_weights, peak_prob, sum_MFP])
                    else:
                        writer.writerow([composition_str, pair, i, x_at_peak, y_at_peak, weight, '', '', sum_MFP])
            else:
                with open('pair_peaks_list_catan.csv', mode='a', newline='') as file:
                    writer = csv.writer(file)
                    if catan == 1:
                        writer.writerow([composition_str, pair, i, x_at_peak, y_at_peak, weight, weight/catan_sum_weights, peak_prob, ''])
                    else:
                        writer.writerow([composition_str, pair, i, x_at_peak, y_at_peak, weight, '', '',''])       
        
        # Output new CSV of splines
        """ # Define the folder name
        folder_name = "Splines"

        # Create the folder if it doesn't exist
        os.makedirs(folder_name, exist_ok=True)                

        filename = f"output_{pair}_{composition_str}.csv"    
        with open(filename, mode='w', newline='') as file:
            writer = csv.writer(file)
            
            # Write the header
            writer.writerow(["x", "y"])
            
            # Write the data
            for x, y in zip(x_new, y_weighted):
                writer.writerow([x, y]) """

        line1, = plt.plot(x_new, y_weighted, label=pair)
        previous_color = line1.get_color()
        plt.plot(x_new, avg_conc, linestyle=':', color=previous_color)
        salt_name += 1

    if MFP_calc > 0:
        plt.axvline(x=MFP_calc, color='r', linestyle='-', label="bc-MFP")

    #plt.plot(x_combined, y_combined, label='Combined RDF',color='cornflowerblue')
    
    # Get the current date
    current_date = datetime.now().strftime("%Y%m%d")
    save_title = f'PartPDF_CatAn_{weighttype}_{composition_str}'

    plt.xlabel('Radius (Å)')
    plt.ylabel('RDF')
    plt.title(f'Partial PDF, {weighttype} : {composition_str} ({source})', fontsize=10)
    plt.legend()
    Save_fig(save_title)
    plt.close()
    return peaks_list


def PopDF(df,salt,column,new_value):
    # Check if the salt composition exists
    if salt in df['Salt Composition'].values:
        # Update the value for the specified column
        df.loc[df['Salt Composition'] == salt, column] = new_value
    else:
        # Create a new entry with NaN for other columns
        new_row = {'Salt Composition': salt, **{attr: (new_value if attr == column else None) for attr in attributes}}
        df = df.append(new_row, ignore_index=True)
    return df
    

def Save_fig(title):
    # Variables for the file name and folder
    folder_path = 'J:\Research\PDF_Analysis_kModel\Figures'
    file_name = f'{title}.png'  # Combine string and variable

    # Full path to save the file
    save_path = os.path.join(folder_path, file_name)

    # Save the plot
    plt.savefig(save_path)
    return


def MFPcalc(peaklist,glist):

    prod = 1
    n = len(9)
    for i in range(0,n):
        for j in range(1,i):
            prod = prod * (1- (glist[j]*glist[j+1])/glist[j+1])
    
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
    splines = interpolate_data(data_dict)

    peaks = 'off'

    composition, ion_counts = parse_composition(composition_str)
    weights = calculate_weights(composition, ion_counts,splines)
    peaks_list = plot_cat_an_splines_weighted(splines, weights, composition_str,source,'Conc',MFP_calc)

    print(peaks_list)

    print("")



attributes = ['Peak Sharpness','Peak Intensity','Decay Intensity','+- Pos_1','+- Pos_2','+- Pos_3','++ Pos_1','++ Pos_2','++ Pos_3','-- Pos_1','-- Pos_2','-- Pos_3']    

# Save the peaks_list to a CSV file
csv_filename = f'pair_peaks_list_catan.csv'
with open(csv_filename, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Salt','Pair', 'Peak Number', 'Distance (A)', 'Weighted g(r)', 'Concentration', 'Cat-An Weighted g(r)', 'MFP Probability', 'Avg. MFP'])  # Header row


Read_RDF('PDF_LiCl.csv',"1.0LiCl",'Walz, 2019; 878K',3.893446)   # Walz, 2019; 878K
Read_RDF('PDF_NaCl.csv',"1.0NaCl",'Walz, 2019; 1074.15K',4.352426)   # PDF 'Walz, 2019; 1074.15K, #'Andersson, 2022; 1250K')   # PDF Andersson, 1250 K, 2022
Read_RDF('PDF_KCl.csv',"1.0KCl",'Walz, 2019; 1043K',4.453593)   # Walz, 2019; 1043K
Read_RDF('PDF_LiF.csv',"1.0LiF",'Walz, 2019; 1121K',3.305228)   # Walz, 2019; 1121K
Read_RDF('PDF_NaF.csv',"1.0NaF",'Walz, 2019; 1266.15K',4.005403)   # Walz, 2019; 1266.15K
Read_RDF('PDF_KF.csv',"1.0KF",'Walz, 2019; 1131.15K',4.225170)   # Walz, 2019; 1131.15K
Read_RDF('PDF_RbF.csv',"1.0RbF",'Walz, 2019; 1068.15K',0)   # Walz, 2019; 1068.15K
Read_RDF('PDF_CsF.csv',"1.0CsF",'Walz, 2019; 955.15K',0)   # Walz, 2019; 955.15K
Read_RDF('PDF_MgCl2.csv',"1.0MgCl2",'McGreevy, 1987; 998K',4.049956)   # McGreevy, 1987; 998K
Read_RDF('PDF_CaCl2.csv',"1.0CaCl2",'McGreevy, 1987; 1093K',3.372692)   # McGreevy, 1987; 1093K
Read_RDF('PDF_NaNO3.csv',"1.0NaNO3",'Tahara, 2017; 623K',0)   # Tahara, 2017; 623K
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

#Read_RDF('PDF_NaCl-UCl3.csv',"0.64NaCl-0.36UCl3",'Andersson, 2022; 1250K',2.455254)   # PDF Andersson, 2022



plt.show()

#  if el1_states[0] > 0 and el2_states[0] < 0:    #is_cation_anion:
#         return 1, element1, element2
#     elif el1_states[0] < 0 and el2_states[0] > 0:  # is_anion_cation:
#         return 1, element2, element1
#     elif el1_states[0] > 0 and el2_states[0] > 0:  # is_cation_cation:
#         return 2, element1, element1
#     else:                                          # is_anion_anion:
#         return 3, element1, element2