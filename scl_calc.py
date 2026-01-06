import os
import pandas as pd
import numpy as np
from scipy.integrate import trapezoid
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import re
from mendeleev import element
from scipy.signal import find_peaks
from datetime import datetime
from scipy.ndimage import gaussian_filter1d
import csv

def standardize_ion_pair(ion_pair):
    # Convert to string if it's not already (handles numpy.float64 and other types)
    if not isinstance(ion_pair, str):
        ion_pair = str(ion_pair)
    
    # Split the ion pair into two parts
    parts = ion_pair.split('-')
    if len(parts) != 2:
        raise ValueError(f"Invalid ion pair format: {ion_pair}. Expected format: 'A-B'")
    
    # Extract elements from each part (handling cases like 'LiCl' -> 'Li' and 'Cl')
    def extract_elements(compound):
        # This regex matches element symbols (1-2 letters) followed by optional numbers
        elements = re.findall(r'([A-Z][a-z]?\d*)', compound)
        return [re.sub(r'\d', '', el) for el in elements]  # Remove any numbers
    
    # Get elements from each part
    elements1 = extract_elements(parts[0])
    elements2 = extract_elements(parts[1])
    
    if not elements1 or not elements2:
        raise ValueError(f"Could not extract elements from ion pair: {ion_pair}")
    
    # For now, just take the first element from each part
    # You might want to enhance this based on your specific needs
    element1 = elements1[0]
    element2 = elements2[0]
    
    # Get element properties
    el1 = element(element1)
    el2 = element(element2)
    el1_states = el1.oxistates or [0]  # Default to 0 if no oxidation states
    el2_states = el2.oxistates or [0]  # Default to 0 if no oxidation states

    # Determine if the elements are cations or anions
    if el1_states[0] > 0 and el2_states[0] < 0:
        return f"{element1}-{element2}"  # cation-anion
    elif el1_states[0] < 0 and el2_states[0] > 0:
        return f"{element2}-{element1}"  # anion-cation
    elif el1_states[0] > 0 and el2_states[0] > 0:
        return '-'.join(sorted([element1, element2]))  # cation-cation, sorted alphabetically
    else:
        # If both are cations or both are anions, sort alphabetically
        return '-'.join(sorted([element1, element2]))
    
# Function to interpolate splines over a common x-range
def interval_cut(splines, filtered_crossing_points, x_range, num_pnts_interval):
    for ion_pair in splines.keys():
        splines[ion_pair]['intervals'] = []  # Add a new key for intervals

    for index, (start, end) in enumerate(zip(filtered_crossing_points[:-1], filtered_crossing_points[1:])):
        x_interval = np.linspace(x_range[start], x_range[end], num_pnts_interval)
        for ion_pair, spline in splines.items():
            y_interpolated = spline['spline'](x_interval)
            splines[ion_pair]['intervals'].append({'interval_number': index, 'x': x_interval, 'y': y_interpolated})

    return splines

def extract_and_average(x, y, x_0, width):
    """
    Extract values from y that span a width centered at x_0 and compute their average.

    Parameters:
    - x: numpy array of x-values
    - y: numpy array of y-values
    - x_0: center point
    - width: total width of the range to consider

    Returns:
    - average of the extracted y-values
    """
    # Calculate the range
    x_min = x_0 - width / 2
    x_max = x_0 + width / 2

    # Create a mask for the range
    mask = (x >= x_min) & (x <= x_max)

    # Extract the y-values within the range
    y_extracted = y[mask]

    # Compute the average
    average = np.mean(y_extracted) if len(y_extracted) > 0 else None

    return average

def format_composition_with_subscripts(composition_str):
    """
    Format composition string with proper subscript notation for chemical formulas.
    Converts numbers after element symbols to subscripts.
    Example: "0.5NaCl-0.5UCl3" -> "0.5NaCl-0.5UCl₃"
    """
    import re

    def format_compound(compound):
        # Pattern to match element symbol followed by numbers
        # This matches patterns like: NaCl2, UCl3, MgF2, etc.
        pattern = r'([A-Z][a-z]?)(\d*)'
        def replace_with_subscript(match):
            element = match.group(1)
            number = match.group(2)
            if number:
                # Convert number to subscript using Unicode subscript characters
                subscript_number = ''.join([f'₀₁₂₃₄₅₆₇₈₉'[int(digit)] for digit in number])
                return f"{element}{subscript_number}"
            return element

        return re.sub(pattern, replace_with_subscript, compound)

    # Split by '-' and format each compound
    compounds = composition_str.split('-')
    formatted_compounds = [format_compound(comp) for comp in compounds]

    return '-'.join(formatted_compounds)
    
def select_elements_at_intervals(array, interval):    
    selected_elements = []
    current_index = 0
    while current_index < len(array):
        selected_elements.append(array[current_index])
        # Calculate the next index by finding the closest element to the current index + interval
        next_index = np.argmin(np.abs(array - (array[current_index] + interval)))
        if next_index <= current_index:
            break
        current_index = next_index
    return selected_elements

class IonPairPDF:
    def __init__(self, ion_pair, x_values, y_values, weights):
        self.ion_pair = standardize_ion_pair(ion_pair)  # Standardize the ion pair name
        self.x = np.array(x_values)
        self.y = np.array(y_values)
        self.weights = weights  # Store the weights
        self.spline = self.create_spline()
        self.type = self.type_ion()
        self.peak,self.minima = self.find_first_peak_and_minimum(self.x, self.y, any_type=False)

    def create_spline(self):
        return interp1d(self.x, self.y, kind='linear', bounds_error=False, fill_value=0)
    
    def type_ion(self):
        element1, element2 = self.ion_pair.split('-')
        el1 = element(element1)
        el2 = element(element2)
        el1_states = el1.oxistates
        el2_states = el2.oxistates
        if el1_states[0] > 0 and el2_states[0] < 0:    #is_cation_anion:
            return "ca"
        elif el1_states[0] < 0 and el2_states[0] > 0:  # is_anion_cation:
            return "ac"
        elif (el1_states[0] > 0 and el2_states[0] > 0) and el1 == el2 :  # is_similar_cation_cation:
            return "cc_sim"
        elif (el1_states[0] > 0 and el2_states[0] > 0) and el1!= el2 :  # is_different_cation_cation:
            return "cc_diff"
        elif (el1_states[0] < 0 and el2_states[0] < 0) and el1 == el2 :  # is_similar_anion_anion:
            return "aa_sim"
        elif (el1_states[0] < 0 and el2_states[0] < 0) and el1!= el2 :  # is_different_anion_anion:
            return "aa_diff"
        
    def find_first_peak_and_minimum(self, x, y, any_type=False):
        if self.type == "ca" or any_type == True:
            # Get the weight for this ion pair
            weight = self.weights.get(self.ion_pair, 1.0)  # Default to 1.0 if not found
            
            # Weight the y values
            y_weighted = y * weight

            # Find all significant peaks with some minimum prominence and width
            # prominence: minimum height difference between peak and its lowest contour line
            # width: minimum width of the peak in samples
            # distance: minimum distance between peaks in samples (about 0.5 Angstrom)
            x_spacing = np.mean(np.diff(x))  # average x spacing
            min_peak_distance = int(0.5 / x_spacing)  # convert 0.5 Angstrom to samples
            
            peaks, properties = find_peaks(
                y_weighted,
                prominence=0.1 * np.max(y_weighted),  # at least 10% of max height
                width=2,  # minimum width of 2 samples
                distance=min_peak_distance  # minimum distance between peaks
            )

            # If we found peaks, take the first one
            if len(peaks) > 0:
                first_peak_index = peaks[0]
                first_peak_x = x[first_peak_index]
                first_peak_y = y_weighted[first_peak_index]
            else:
                first_peak_x, first_peak_y = None, None

            # Find the first minimum after the first peak
            if first_peak_x is not None:
                after_peak = y_weighted[first_peak_index:]
                minima, _ = find_peaks(
                    -after_peak,
                    prominence=0.1 * np.max(y_weighted),
                    width=2,
                    distance=min_peak_distance
                )
                if len(minima) > 0:
                    first_min_index = minima[0] + first_peak_index
                    first_min_x = x[first_min_index]
                    first_min_y = y_weighted[first_min_index]
                else:
                    first_min_x, first_min_y = None, None
            else:
                first_min_x, first_min_y = None, None
        else:
            # For non-cation-anion pairs, return None
            first_peak_x, first_peak_y = None, None
            first_min_x, first_min_y = None, None

        return (first_peak_x, first_peak_y), (first_min_x, first_min_y)


class MoltenSaltPDF:
    def __init__(self, comp, pdf_file, source, temp, gamma_bc):
        self.comp = comp
        self.source = source
        self.temp = temp
        self.gamma_bc = gamma_bc
        self.composition, self.ion_counts, self.comp = self.parse_composition(comp)        
        self.weights = self.calculate_weights()
        self.ion_pairs = self.load_pdf_data(pdf_file)
        self.interpolated_splines, self.weighted_splines, self.x_new = self.interpolate_data()
        self.plot_data = []  # Store plot data here
        self.ion_pair_results = {}  # Store ion pair analysis results

    def load_pdf_data(self, pdf_file):
        folder_name = "RDF_CSV"
        # Use the directory where the script is located
        script_dir = os.path.dirname(os.path.abspath(__file__))
        rdf_plots_folder = os.path.join(script_dir, folder_name)
        file_path = os.path.join(rdf_plots_folder, pdf_file)
        
        df = pd.read_csv(file_path, header=None)
        ion_pairs = {}

        for i in range(0, df.shape[1], 2):
            # Get the ion pair name and handle potential NaN or non-string values
            raw_ion_pair = df.iloc[0, i]
            if pd.isna(raw_ion_pair):
                continue  # Skip if the value is NaN
            ion_pair = standardize_ion_pair(raw_ion_pair)  # Standardize the ion pair name
            # Get and clean the data
            x_values = df.iloc[2:, i].dropna().astype(float).values
            y_values = df.iloc[2:, i+1].dropna().astype(float).values
            
            # Sort values by x
            sort_idx = np.argsort(x_values)
            x_sorted = x_values[sort_idx]
            y_sorted = y_values[sort_idx]
            
            # Remove duplicate x-values (keep first occurrence)
            unique_mask = np.concatenate(([True], np.diff(x_sorted) > 0))
            x_unique = x_sorted[unique_mask]
            y_unique = y_sorted[unique_mask]
            
            # Check for non-monotonic points that weren't duplicates
            if len(x_unique) < len(x_sorted):
                print(f"Warning: Removed {len(x_sorted) - len(x_unique)} non-monotonic points for {ion_pair}")
            
            # Extend x_values and y_values
            x_extended = np.linspace(0, x_unique[0], num=int(x_unique[0]) + 1)
            y_extended = np.zeros_like(x_extended)
            
            x_combined = np.concatenate((x_extended, x_unique))
            y_combined = np.concatenate((y_extended, y_unique))

            ion_pairs[ion_pair] = IonPairPDF(ion_pair, x_combined, y_combined, self.weights)

        return ion_pairs

    def parse_composition(self, composition_str):
        components = composition_str.split('-')
        composition = {}
        for component in components:
            match = re.match(r"([0-9.]+)([A-Za-z0-9]+)", component)
            if match:
                fraction, salt = match.groups()
                composition[salt] = float(fraction)
        
        compound_split = re.findall(r'([0-9.]+)([A-Za-z0-9]+)', composition_str)
        ion_counts = {}
        for fraction, compound in compound_split:
            elements = re.findall(r'([A-Z][a-z]*)([0-9]*)', compound)
            ion_count = {}
            for element, count in elements:
                count = int(count) if count else 1
                if element in ion_count:
                    ion_count[element] += count
                else:
                    ion_count[element] = count
            ion_counts[compound] = ion_count

        # Function to get cation atomic number for sorting
        def get_cation_atomic_number(salt):
            try:
                from mendeleev import element
                import re

                # Extract cation from salt (e.g., 'NaCl' -> 'Na', 'LiF' -> 'Li')
                # Look for element symbol followed by optional number
                elements = re.findall(r'([A-Z][a-z]?)(\d*)', salt)
                if elements:
                    cation = elements[0][0]  # First element is typically the cation
                    el = element(cation)
                    return el.atomic_number
            except:
                pass
            return 999  # Default high number for unknown elements

        # Sort salts by cation atomic number instead of alphabetically
        sorted_salts = sorted(composition.keys(), key=get_cation_atomic_number)
        sorted_composition_str = '-'.join([f"{composition[salt]}{salt}" for salt in sorted_salts])

        return composition, ion_counts, sorted_composition_str

    def calculate_weights(self):
        element_concentration = {}
        for salt, fraction in self.composition.items():
            for element, count in self.ion_counts[salt].items():
                if element in element_concentration:
                    element_concentration[element] += fraction * count
                else:
                    element_concentration[element] = fraction * count

        total_concentration = sum(element_concentration.values())
        rel_conc = {element: concentration / total_concentration 
                                for element, concentration in element_concentration.items()}

        weights = {}
        for ion1, conc1 in rel_conc.items():
            for ion2, conc2 in rel_conc.items():
                ion_pair = standardize_ion_pair(f"{ion1}-{ion2}")  # Standardize the ion pair name
                if ion_pair in weights:
                    weights[ion_pair] = conc1 * conc2
                else:
                    weights[ion_pair] = conc1 * conc2

        sumweights = sum(weights.values())
        weights = {key: value/sumweights for key, value in weights.items()}
        
        return weights

    def interpolate_data(self, num_pnts_interval=1000):
        interpolated_splines = {}
        weighted_splines = {}
        x_min = min(min(pdf.x) for pdf in self.ion_pairs.values())
        x_max = min(max(pdf.x) for pdf in self.ion_pairs.values())
        x_new = np.linspace(0, x_max, num_pnts_interval)

        for ion_pair, pdf in self.ion_pairs.items():
            interpolated_spline = interp1d(x_new, pdf.spline(x_new), kind='linear', bounds_error=False, fill_value=0)
            interpolated_splines[ion_pair] = interpolated_spline
            
            if ion_pair in self.weights:
                weight = self.weights[ion_pair]
                weighted_splines[ion_pair] = lambda x, s=interpolated_spline, w=weight: s(x) * w
            else:
                weighted_splines[ion_pair] = interpolated_spline

        return interpolated_splines, weighted_splines, x_new

    def analyze_pdf(self):
        print("")
        print(f"###  {self.comp}   #############")

        # Prepare data for CSV
        csv_data = [self.comp]  # Start with the composition

        # Check if the salt is a unary divalent salt
        unary_divalent_salts = ["MgCl2", "CaCl2"]
        spline_noise = any(salt in self.comp for salt in unary_divalent_salts)

        weighted_SCLs = 0
        weights_pair = 0

        for ion_pair_ca_i, pdf_data_ca_i in self.ion_pairs.items():
            if pdf_data_ca_i.type == "ca":
                print(f"Cation-anion pair: {ion_pair_ca_i}")

                # Access peak and minima
                peak_i = pdf_data_ca_i.peak
                minima_i = pdf_data_ca_i.minima
                print(f"Peak X: {peak_i[0]}, Minima X: {minima_i[0]}")

                # Use the interpolated x-range for all calculations
                x_range_i = self.x_new
                
                cation_i, anion_i = ion_pair_ca_i.split('-')
                y_combined = np.zeros_like(x_range_i)
                y_weighted_ca_total = np.zeros_like(x_range_i)

                conc_ca_i = self.weights[ion_pair_ca_i]
                y_weighted_ca_i = self.weighted_splines[ion_pair_ca_i](x_range_i)

                c_ca_i = conc_ca_i# * peak_i[1] * peak_i[0]**-1
                sum_c_ca_j = 0

                ion_pair_cc_i = '-'.join(sorted([cation_i, cation_i]))
                no_cc_spline = False
                try:
                    y_weighted_cc_i = self.weighted_splines[ion_pair_cc_i](x_range_i)
                except KeyError:
                    y_weighted_cc_i = np.linspace(self.weights[ion_pair_cc_i], self.weights[ion_pair_cc_i], len(x_range_i))
                    no_cc_spline = True
                    #print("No cation-cation spline data available. Approximating based on concentration.")

                y_ideal = np.zeros_like(x_range_i)
                sum_conc_ca = 0
                sum_conc_cc = 0
                conc_cc_j = 0
                sum_c_cc_j = 0
                y_weighted_cc_total = np.zeros_like(x_range_i)
                # Calculate peak properties, combined spline, etc.
                for ion_pair_j, pdf_data_j in self.ion_pairs.items():
                    ion1, ion2 = ion_pair_j.split('-')
                    y_weighted_j = self.weighted_splines[ion_pair_j](x_range_i)

                    if pdf_data_j.type == "ca":
                        # Access peak and minima
                        peak_j = pdf_data_j.peak
                        minima_j = pdf_data_j.minima
                        conc_ca_j = self.weights[ion_pair_j]
                        sum_c_ca_j += conc_ca_j# * peak_j[1] * peak_j[0]**-1
                        sum_conc_ca += conc_ca_j
                        if cation_i in [ion1, ion2]:
                            y_weighted_ca_total += y_weighted_j

                    if cation_i in [ion1, ion2]:
                        y_combined += y_weighted_j

                    if pdf_data_j.type == "cc_sim" or pdf_data_j.type == "cc_diff":
                        conc_cc_j = self.weights[ion_pair_j]
                        sum_c_cc_j += conc_cc_j
                        sum_conc_cc += conc_cc_j
                        if cation_i in [ion1, ion2]:
                            y_weighted_cc_total += y_weighted_j
                    
                # Create a spline for y_combined using the interpolated x-range
                y_combined_spline = interp1d(x_range_i, y_combined, kind='linear', bounds_error=False, fill_value=0)
                
                # Use the interpolated x-range for transfer points
                transfer_points = select_elements_at_intervals(self.x_new, peak_i[0])
                transfer_points_indices = [np.argmin(np.abs(self.x_new - point)) for point in transfer_points]
                transfer_points = self.x_new[transfer_points_indices]  # Ensure exact x-values from the grid

                # Calculate cation-cation peak position and check if it exists
                cc_transfer = False
                DF = 0
                if ion_pair_cc_i in self.ion_pairs:
                    pdf_cc_i = self.ion_pairs[ion_pair_cc_i]
                    cc_peak_i, _ = pdf_cc_i.find_first_peak_and_minimum(pdf_cc_i.x, pdf_cc_i.y, any_type=True)
                    print(f"Cation-cation pair: {ion_pair_cc_i}")
                    if cc_peak_i[0] is not None:
                        print(f"First cc peak: x = {cc_peak_i[0]:.2f}, y = {cc_peak_i[1]:.2f}")
                        if cc_peak_i[0] > peak_i[0] and cc_peak_i[0] < 2*transfer_points[1]:
                            # DF = ((peak_i[0] / cc_peak_i[0]) ** 2 - (1 / 2) ** 2) / (1 - (1 / 2) ** 2)
                            DF = ((peak_i[0] / cc_peak_i[0]) - (1 / 2)) / (1 - (1 / 2))
                            cc_transfer = True
                    else:
                        print("No peak found for cation-cation pair")
                else:
                    print(f"No data found for cation-cation pair: {ion_pair_cc_i}")

                b_KF = np.zeros(len(transfer_points))
                b_NI = np.zeros(len(transfer_points))
                b_PC = np.zeros(len(transfer_points))
                b_SRO = np.zeros(len(transfer_points))
                #b_CS = np.zeros(len(transfer_points))   # Charge screening from anion density

                # Bond strength/diffusion factor: How well ideal recipient can receive energy
                KF = (abs(peak_i[1]-minima_i[1])/peak_i[1]) 
                # Non-ideal recipient factor: Probability of energy transfer to ideal recipient
                NI = (y_weighted_ca_i[transfer_points_indices[1]]/y_combined[transfer_points_indices[1]]) # First transfer
                # Phonon transfer factor: Probability of succeeding ideal recipient in structure
                PH = conc_ca_i/sum_conc_ca 

                # Relative transfer participation of pair
                RTE = c_ca_i/sum_c_ca_j

                next_not_last = True
                for i in range(0,len(transfer_points)):
                    r_i = transfer_points_indices[i]
                    if i == len(transfer_points)-1:
                        next_not_last = False
                    else:
                        r_i_next = transfer_points_indices[i+1]

                    if i == 0:           # Assume all energy travels minimum distance
                        continue

                    # if i % 2 == 1:      # Anion-to-cation transfer (ca peak)
                    #     if y_combined[r_i_next] == 0 and next_not_last:
                    #         PH = 1
                    #     elif next_not_last:
                    #         PH = y_weighted_cc_i[r_i_next] / y_weighted_cc_total[r_i_next]
                    #     else:
                    #         PH = 0
                    # else:               # Cation-to-anion transfer (cc peak)
                    #     if y_combined[r_i_next] == 0 and next_not_last:
                    #         PH = 1
                    #     elif next_not_last:
                    #         PH = y_weighted_ca_i[r_i_next] / y_weighted_ca_total[r_i_next]
                    #     else:
                    #         PH = 0
                    # PHi = PH + (1 - PH)*0.5  # Effective phonon transfer probability considering bond angle from 180-90

                    # b_KF[i] += 1 - KF
                    # b_NI[i] += 1 - NI
                    # b_PC[i] += 1 - PHi

                    if i % 2 == 1:      # Anion-to-cation transfer
                        if y_combined[r_i] == 0:
                            NI = 1
                        else:
                            NI = (y_weighted_ca_i[r_i]/y_combined[r_i])
                    else:               # Cation-to-anion transfer
                        if y_combined[r_i] == 0:
                            NI = 1
                        else:
                            NI = (y_weighted_cc_i[r_i]/y_combined[r_i])
                    if cc_transfer:
                        # PH_f = PH + (1 - PH)*KF**12*DF**(1/12)
                        PH_f = PH #+ (1 - PH)*KF**6*DF # <--- Good
                        # PH_f = PH + DF*KF - PH*DF*KF
                    else:
                        PH_f = PH

                    b_KF[i] += 1 - KF
                    b_NI[i] += 1 - NI
                    b_PC[i] += 1 - PH_f

                    # if i % 2 == 1:      # Anion-to-cation transfer
                    #     if y_combined[r_i] == 0:
                    #         b_NI[i] += 0
                    #     else:
                    #         b_NI[i] += 1 - (y_weighted_ca_i[r_i]/y_combined[r_i])
                    #     if y_combined[r_i_next] == 0 and next_not_last:
                    #         b_PC[i] += 0
                    #     elif next_not_last:
                    #         b_PC[i] += 1 - y_weighted_cc_i[r_i_next] / y_combined[r_i_next]
                    #     else:
                    #         b_PC[i] += 1
                    # else:               # Cation-to-anion transfer
                    #     if y_combined[r_i] == 0:
                    #         b_NI[i] += 0
                    #     else:
                    #         b_NI[i] += 1 - (y_weighted_cc_i[r_i]/y_combined[r_i])
                    #     if y_combined[r_i_next] == 0 and next_not_last:
                    #         b_PC[i] += 0
                    #     elif next_not_last:
                    #         b_PC[i] += 1 - y_weighted_ca_i[r_i_next] / y_weighted_ca_total[r_i_next]
                    #     else:
                    #         b_PC[i] += 1

                    # Ensure b_factors are within [0, 1]
                    b_KF[i] = min(max(b_KF[i], 0), 1)
                    b_NI[i] = min(max(b_NI[i], 0), 1)
                    b_PC[i] = min(max(b_PC[i], 0), 1)
                    b_SRO[i] = min(max(b_SRO[i], 0), 1)

                beta_i = np.zeros(len(transfer_points))
                beta_i_integral = 0#np.zeros(len(transfer_points))

                b_array = [b_KF,b_NI,b_PC,b_SRO]

                S_i = np.zeros(len(transfer_points))   # Initialize cumulative survival function
                
                for k in range(len(transfer_points_indices)):
                    b_i_sum = 0
                    for b_i in b_array:
                        if b_i[k] == 1:
                            b_i_sum = (float('inf'))
                        else:
                            b_i_sum += b_i[k] / (1 - b_i[k])
                    beta_i[k] += (b_i_sum)

                #print(f'beta_i: {beta_i}')

                # Calculate cumulative survival function S_i
                for k in range(0,len(transfer_points_indices)-1):
                    beta_i_integral += beta_i[k]*(transfer_points[k+1] - transfer_points[k])
                    S_i[k] = np.exp(-beta_i_integral)
                #S_i = S_i*b_RTO[0]
                #print(f'S_i: {S_i}')

                S_i_y = np.zeros(len(x_range_i))   # Initialize S_i_y with zeros
                # Iterate over x_range_i and assign the corresponding S_i value
                for i in range(len(x_range_i)):
                    for k in range(0, len(transfer_points_indices)-1):
                        if transfer_points_indices[k+1] > i >= transfer_points_indices[k]:
                            S_i_y[i] = S_i[k]
                            break  # Exit the loop once the correct S_i is found
                x_SCL_pair = np.trapz(S_i_y, x=x_range_i)
                print(f'x_SCL_pair: {x_SCL_pair}')  
                
                weights_pair += RTE
                weighted_SCLs += x_SCL_pair * RTE

                # Store results for this ion pair
                self.ion_pair_results[ion_pair_ca_i] = {
                    'type': pdf_data_ca_i.type,
                    'peak_x': peak_i[0],
                    'peak_y': peak_i[1],
                    'minima_x': minima_i[0],
                    'minima_y': minima_i[1],
                    'transfer_factors': [b_KF, b_NI, b_PC, b_SRO],
                    'scl': x_SCL_pair,
                    'weight': RTE,
                    'ion_pair_ca_i': ion_pair_ca_i,
                    'ion_pair_cc_i': ion_pair_cc_i,
                    'cc_transfer': cc_transfer,
                    'cc_peak_x': cc_peak_i[0] if 'cc_peak_i' in locals() and cc_peak_i[0] is not None else None,
                    'cc_peak_y': cc_peak_i[1] if 'cc_peak_i' in locals() and cc_peak_i[0] is not None else None,
                    'DF': DF if 'DF' in locals() else None
                }

                # Store plot data
                self.plot_data.append({
                    # 'x_range': x_range_i,
                    # 'y_ideal': y_ideal,
                    # 'y_combined': y_combined,
                    # 'beta_y': beta_y,
                    # 'beta_i': beta_i,
                    # 'beta_i_integral': beta_i_integral,
                    # 'S_i': S_i,
                    'x_range': x_range_i,
                    'x_SCL_pair': x_SCL_pair,
                    'S_i': S_i_y,
                    'ion_pair_ca_i': ion_pair_ca_i,
                    'MFP_bc': self.gamma_bc,
                })

                # print(f"lambda_P(r) - {ion_pair_ca_i} = {round(x_SCL_pair, 5)}")
        avg_SCL = weighted_SCLs / weights_pair
        self.plot_data.append({'avg_SCL': round(avg_SCL, 5)})
        print(f"Average SCL: {round(avg_SCL, 5)}")
        print(f"lambda_BC = {round(self.gamma_bc, 5)}")

        # Define CSV filename and check if it exists
        csv_filename = 'scl_results.csv'
        file_exists = os.path.isfile(csv_filename)
        
        # Define the base headers that are always present
        base_headers = ['Composition', 'Source', 'Temperature (K)', 'Average SCL (A)']
        
        # Prepare data row with base information
        base_data = [self.comp, self.source, self.temp, round(avg_SCL, 5)]
        
        # Get the current ion pairs and sort them for consistent ordering
        current_pairs = sorted(self.ion_pair_results.items())
        
        # Create a dictionary to map pair labels to their data
        pair_data = {}
        for pair_num, (ion_pair, result) in enumerate(current_pairs, 1):
            if pair_num > 6:  # Limit to 6 pairs as per requirement
                break
                
            pair_data[pair_num] = {
                'label': ion_pair,
                'data': [
                    round(result['scl'], 5) if result['scl'] is not None else '',
                    round(result['peak_x'], 5) if result['peak_x'] is not None else '',
                    round(result['peak_y'], 5) if result['peak_y'] is not None else '',
                    round(result['minima_x'], 5) if result['minima_x'] is not None else '',
                    round(result['minima_y'], 5) if result['minima_y'] is not None else '',
                    round(result.get('cc_peak_x', ''), 5) if result.get('cc_peak_x') is not None else '',
                    round(result.get('cc_peak_y', ''), 5) if result.get('cc_peak_y') is not None else ''
                ]
            }
        
        # If the file doesn't exist, create it with headers
        if not file_exists:
            headers = base_headers.copy()
            
            # Generate headers for all pairs (up to 6 as per requirement)
            for pair_num in range(1, 7):
                # Add pair label header
                headers.append(f'Pair {pair_num} Label')
                
                # Add pair-specific headers
                pair_headers = [
                    f'Pair {pair_num} SCL_i (A)',
                    f'Pair {pair_num} Peak X (A)',
                    f'Pair {pair_num} Peak Y',
                    f'Pair {pair_num} Min X (A)',
                    f'Pair {pair_num} Min Y',
                    f'Pair {pair_num} CC Peak X (A)',
                    f'Pair {pair_num} CC Peak Y'
                ]
                headers.extend(pair_headers)
            
            # Write headers to the new file
            with open(csv_filename, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(headers)
        else:
            # Read existing data to check for duplicate entries
            with open(csv_filename, 'r', newline='') as f:
                reader = csv.reader(f)
                headers = next(reader, [])
                existing_data = list(reader)
            
            # Filter out any existing entry with the same composition AND source
            existing_data = [row for row in existing_data if len(row) > 1 and not (row[0] == self.comp and row[1] == self.source)]
            
            # Reconstruct the file with filtered data
            with open(csv_filename, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(headers)
                writer.writerows(existing_data)
        
        # Prepare the new data row according to headers
        new_data_row = base_data.copy()
        
        # Fill in the data according to the headers
        i = len(base_headers)  # Start after base headers
        pair_num = 1
        
        while i < len(headers):
            # If this is a pair label header
            if f'Pair {pair_num} Label' in headers[i]:
                # Add the pair label if we have data for this pair number
                if pair_num in pair_data:
                    new_data_row.append(pair_data[pair_num]['label'])
                else:
                    new_data_row.append('')
                
                # Add the pair data if we have it
                if pair_num in pair_data:
                    new_data_row.extend(pair_data[pair_num]['data'])
                else:
                    # Add empty values for all pair data columns
                    new_data_row.extend([''] * 7)  # 7 data columns per pair
                
                i += 8  # Move to next pair (1 label + 7 data columns)
                pair_num += 1
            else:
                # If we encounter an unexpected header format, just add an empty value
                new_data_row.append('')
                i += 1
        
        # Append the new data row to the CSV
        with open(csv_filename, 'a', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(new_data_row)

    def save_plot_data(self, folder='plot_data', overwrite=True):
        """Save the plot data for this salt to a CSV file.
        
        Args:
            folder (str): Directory to save the plot data files
            overwrite (bool): If True, overwrite existing files with the same name
                             If False, skip saving if file exists
        """
        os.makedirs(folder, exist_ok=True)
        # Create a safe filename that includes both composition and source
        safe_source = ''.join(c if c.isalnum() else '_' for c in self.source.split(',')[0].strip())
        filename = os.path.join(folder, f'{self.comp.replace("-", "_")}_{safe_source}_plot_data.csv')
        
        # Skip if file exists and we're not overwriting
        if not overwrite and os.path.exists(filename):
            print(f"Skipping {filename} - file already exists (use overwrite=True to replace)")
            return
        
        # Get the x-range used in the PDF analysis
        if not hasattr(self, 'x_new') or len(self.x_new) == 0:
            print(f"No PDF data available for {self.comp}")
            return
            
        # Create a DataFrame with the interpolated x-values
        df = pd.DataFrame({'r (A)': self.x_new})
        
        # Add weighted PDF g(r) for each ion pair using the interpolated splines
        for ion_pair, spline_func in self.weighted_splines.items():
            # Evaluate the weighted spline at the interpolated x-values
            weighted_g_r = spline_func(self.x_new)
            # Add to DataFrame with 'weighted_' prefix
            df[f'g(r)_weighted_{ion_pair}'] = weighted_g_r
            
            # Also save the unweighted PDF if available
            if ion_pair in self.interpolated_splines:
                unweighted_g_r = self.interpolated_splines[ion_pair](self.x_new)
                df[f'g(r)_unweighted_{ion_pair}'] = unweighted_g_r
        
        # Add S_i values, ensuring they're on the same x-grid
        for i, data in enumerate(self.plot_data):
            if 'S_i' in data:
                # Get the ion pair name or use a default
                ion_pair = data.get('ion_pair_ca_i', f'pair_{i}')
                series_name = f'S_i_{ion_pair}'
                
                # If S_i is already on the same grid, use it directly
                if 'x_range' in data and len(data['x_range']) == len(self.x_new) and np.allclose(data['x_range'], self.x_new):
                    df[series_name] = data['S_i']
                # Otherwise, interpolate to the common grid
                elif 'x_range' in data and 'S_i' in data and len(data['S_i']) > 0:
                    # Create interpolation function for this S_i series
                    interp_func = interp1d(
                        data['x_range'], 
                        data['S_i'], 
                        kind='linear', 
                        bounds_error=False, 
                        fill_value=np.nan
                    )
                    # Add interpolated values to DataFrame
                    df[series_name] = interp_func(self.x_new)
        
        # Add metadata
        if any('avg_SCL' in data for data in self.plot_data):
            df['Average_SCL'] = next((data['avg_SCL'] for data in self.plot_data if 'avg_SCL' in data), np.nan)
        if hasattr(self, 'gamma_bc'):
            df['lambda_BC'] = self.gamma_bc
        
        # Save to CSV
        df.to_csv(filename, index=False)
        print(f"Plot data saved to {filename}")
        return filename

    def plot_pdf(self, show_plot=True, save_plot=True, output_dir='scl_plots'):
        """Plot the PDF and SCL analysis.
        
        Args:
            show_plot (bool): Whether to display the plot
            save_plot (bool): Whether to save the plot to a file
            output_dir (str): Directory to save the plot
        """
        if not hasattr(self, 'plot_data') or not self.plot_data:
            print(f"No plot data available for {self.comp}")
            return
            
        # Create output directory if it doesn't exist
        if save_plot:
            os.makedirs(output_dir, exist_ok=True)
            
        # Create a safe source string for filenames
        safe_source = ''.join(c if c.isalnum() else '_' for c in self.source.split(',')[0].strip())
        timestamp = datetime.now().strftime("%Y%m%d")
        base_filename = f'PDF_{self.comp}_{safe_source}' #_{timestamp}
        
        # Apply publication-style matplotlib settings to match tc_batch_cli.py
        plt.rcParams['font.family'] = 'Times New Roman'
        plt.rcParams['font.size'] = 14
        plt.rcParams['axes.labelsize'] = 14
        plt.rcParams['axes.labelweight'] = 'bold'
        plt.rcParams['axes.linewidth'] = 1.5
        plt.rcParams['xtick.labelsize'] = 14
        plt.rcParams['ytick.labelsize'] = 14
        plt.rcParams['xtick.direction'] = 'out'
        plt.rcParams['ytick.direction'] = 'out'
        plt.rcParams['xtick.major.width'] = 1.75
        plt.rcParams['ytick.major.width'] = 1.75
        # Legend and small text/annotation sizing (~0.85x of tick labels)
        _tick_size = 14
        _small_text = int(round(0.95 * 12))#_tick_size))
        plt.rcParams['legend.frameon'] = False
        plt.rcParams['legend.fontsize'] = _small_text
        # Ensure mathtext (e.g., $r_{i}$) uses Times New Roman
        plt.rcParams['mathtext.fontset'] = 'custom'
        plt.rcParams['mathtext.rm'] = 'Times New Roman'
        plt.rcParams['mathtext.it'] = 'Times New Roman:italic'
        plt.rcParams['mathtext.bf'] = 'Times New Roman:bold'

        plt.figure(figsize=(4.75, 4.25))
        ion_pair_colors = {}  # Dictionary to store colors for each ion pair

        # Count the number of cation-anion pairs
        ca_pair_count = sum(1 for pdf_data in self.ion_pairs.values() if pdf_data.type == "ca")

        for ion_pair, pdf_data in self.ion_pairs.items():
            # Always use the interpolated splines for plotting
            if ion_pair in self.weighted_splines:
                # Evaluate the spline at the interpolated x-range
                weighted_y = self.weighted_splines[ion_pair](self.x_new)
                spline, = plt.plot(self.x_new, weighted_y, label=f"{ion_pair}")
                spline_color = spline.get_color()
                ion_pair_colors[ion_pair] = spline_color  # Store the color for this ion pair

            # Access peak and minima
            peak = pdf_data.peak
            minima = pdf_data.minima

            # if peak[0] is not None:
            #     peakpoint, = plt.plot(peak[0], peak[1], 'ro', markersize=10)
            #     #plt.axvline(x=2*peak[0], color=peakpoint.get_color(),linestyle='dotted')
            #     plt.annotate(f"Peak: ({peak[0]:.2f}, {peak[1]:.2f})", 
            #                 (peak[0], peak[1]), 
            #                 xytext=(5, 5), 
            #                 textcoords='offset points')

            # if minima[0] is not None:
            #     plt.plot(minima[0], minima[1], 'bo', markersize=10)
            #     plt.annotate(f"Min: ({minima[0]:.2f}, {minima[1]:.2f})", 
            #                 (minima[0], minima[1]), 
            #                 xytext=(5, 5), 
            #                 textcoords='offset points')

        # In the plot_pdf method, update the S_i plotting section to:
        for data in self.plot_data:
            if 'ion_pair_ca_i' in data:
                ion_pair = data['ion_pair_ca_i']
                color = ion_pair_colors.get(ion_pair, 'black')
                
                # Get S_i data
                if 'S_i' in data and len(data['S_i']) == len(self.x_new):
                    S_i = data['S_i']
                elif 'S_i' in data and 'x_range' in data:
                    S_i = np.interp(self.x_new, data['x_range'], data['S_i'])
                else:
                    continue
                    
                # Find the steps in S_i
                step_indices = np.where(np.diff(S_i) != 0)[0] + 1
                step_indices = np.concatenate(([0], step_indices, [len(S_i)-1]))
                
                # Plot bars for each step
                for i in range(len(step_indices)-1):
                    start_idx = step_indices[i]
                    end_idx = step_indices[i+1] if (i < len(step_indices)-1) else len(self.x_new)-1
                    x_start = self.x_new[start_idx]
                    x_end = self.x_new[end_idx] if end_idx < len(self.x_new) else self.x_new[-1]
                    
                    # Calculate the midpoint between steps for bar width
                    if i < len(step_indices)-2:
                        next_x_start = self.x_new[step_indices[i+1]]
                        bar_width = (next_x_start - x_start)# / 2
                    else:
                        bar_width = (x_end - x_start)# / 2
                        
                    s_value = S_i[start_idx]
                    # New approach - scale between min_alpha and max_alpha
                    min_alpha = 0.1
                    max_alpha = 0.55
                    alpha = min_alpha + (max_alpha - min_alpha) * s_value
                    
                    # Plot the bar (centered on the step)
                    bar_x = (x_start + x_end) / 2
                    plt.bar(bar_x, s_value, width=bar_width, 
                        color=color, alpha=alpha*0.6, edgecolor='none', 
                        align='center', zorder=0)
                    
                    # # Plot the top line
                    # plt.hlines(s_value, x_start, x_end, colors=color, 
                    #         linewidth=1.5, alpha=0., zorder=0)

        for data in self.plot_data:
            # Check if 'ion_pair_ca_i' key exists
            if 'ion_pair_ca_i' in data:
                ion_pair = data['ion_pair_ca_i']
                color = ion_pair_colors.get(ion_pair, 'black')  # Default to black if not found
                
                # Ensure we're using the interpolated x-range for S_i
                if 'S_i' in data and len(data['S_i']) == len(self.x_new):
                    plt.plot(self.x_new, data['S_i'], label=f"S(r): {ion_pair}", color=color, linestyle='dotted')
                elif 'S_i' in data and 'x_range' in data:
                    # If S_i was calculated on a different x-range, interpolate it to self.x_new
                    interp_S_i = np.interp(self.x_new, data['x_range'], data['S_i'])
                    plt.plot(self.x_new, interp_S_i, label=f"S(r): {ion_pair}", color=color, linestyle='dotted')
                
                # Only plot the line if there is more than one cation-anion pair
                # if ca_pair_count > 1 and 'x_SCL_pair' in data:
                #     plt.axvline(x=data['x_SCL_pair'], color='gray', linestyle='--', 
                #                label=f"$\\lambda_{{{ion_pair}}}$ = {round(data['x_SCL_pair'], 2)}")

        # Use dark green for the average SCL line
        plt.axvline(x=data['avg_SCL'], color='g', linestyle='-.', label=f"$\\ell_{{\\mathrm{{sc}}}}$ = {round(data['avg_SCL'], 2)}")
        if self.gamma_bc > 0:
            plt.axvline(x=self.gamma_bc, color='k', linestyle='--', label=f"$\\ell_{{\\mathrm{{exp}}}}$ = {round(self.gamma_bc,2)}")

        current_date = datetime.now().strftime("%Y%m%d")
        save_title = f'PDF_{self.comp}_{current_date}.png'
        
        folder_name = "scl_plots"
        os.makedirs(os.path.join(os.getcwd(), folder_name), exist_ok=True)
        file_path = os.path.join(os.getcwd(), folder_name, save_title)

        plt.xlabel('r [Å]')
        plt.ylabel('g(r)')
        # plt.title(f'{format_composition_with_subscripts(self.comp)} ({self.temp}K)')
        # Calculate the x-range to start 0.5 Angstrom before the first nonzero data point
        x_min = min(self.x_new)
        x_max = max(self.x_new)

        # Find the first nonzero data point across all ion pairs
        first_nonzero_x = x_max  # Start with maximum as fallback
        for ion_pair, pdf_data in self.ion_pairs.items():
            if ion_pair in self.weighted_splines:
                weighted_y = self.weighted_splines[ion_pair](self.x_new)
                # Find first index where y > 0.01 (small threshold to avoid numerical noise)
                nonzero_indices = np.where(weighted_y > 0.01)[0]
                if len(nonzero_indices) > 0:
                    first_nonzero_x = min(first_nonzero_x, self.x_new[nonzero_indices[0]])

        # Set x-axis range to start 0.5 Angstrom before first nonzero data point
        x_start = max(x_min, first_nonzero_x - 0.75)
        x_pad = (x_max-x_min)*0.05
        plt.xlim(x_start, x_max+x_pad)

        ax = plt.gca()
        # Only create top axis for unary salts with a single endmember
        if len(self.composition) == 1:  # Check if it's a unary salt
            ca_pairs = [(ip, res) for ip, res in self.ion_pair_results.items() 
                    if res.get('type') == 'ca' and res.get('peak_x')]
            max_labels = 8  # limit labels per pair to avoid horizontal collisions
            # plt.figure(figsize=(4.75, 4.35))
            for idx, (ion_pair, result) in enumerate(ca_pairs):
                rep_peak = result.get('peak_x')
                if not rep_peak or rep_peak <= 0:
                    continue
                r_points_full = np.arange(0.0, (x_max + x_pad) + 0.5 * rep_peak, rep_peak)
                r_labels_full = ["" if i == 0 else rf"$r_{{{i}}}$" for i in range(len(r_points_full))]
                valid = (r_points_full >= x_start) & (r_points_full <= (x_max + x_pad))
                r_points = r_points_full[valid]
                r_labels = [lbl for lbl, v in zip(r_labels_full, valid) if v]
                # Thin labels to avoid horizontal overlap
                if len(r_points) > max_labels:
                    step = int(np.ceil(len(r_points) / max_labels))
                    r_points = r_points[::step]
                    r_labels = r_labels[::step]
                # Create a dedicated twin axis for this pair
                secax = ax.twiny()
                secax.set_xlim(ax.get_xlim())
                secax.set_xticks(r_points)
                secax.set_xticklabels(r_labels)
                # Color tick labels to match the ion pair curve
                color = ion_pair_colors.get(ion_pair, 'black')
                for lbl in secax.get_xticklabels():
                    lbl.set_color(color)
                # Styling: no axis label/spine, no tick lines, and offset pad to avoid overlap
                secax.set_xlabel("")
                if 'top' in secax.spines:
                    secax.spines['top'].set_visible(False)
                secax.tick_params(axis='x', which='major', length=5, width=1.25, colors=color, pad=6 + 12 * idx)

        # ax.legend(ncol=2, loc='upper right', bbox_to_anchor=(1.0, 1.0), facecolor='white', framealpha=1)
        legend = ax.legend(
            ncol=2,
            loc='upper right',
            bbox_to_anchor=(1.0, 1.0),
            facecolor='white',
            frameon=True,  # Remove the frame/border
            framealpha=0.75,
            edgecolor='none',  # Remove the border
            borderpad=0.5,  # Reduce padding around the legend
            borderaxespad=0.5,  # Reduce padding between border and axes
            handletextpad=0.5,  # Reduce space between legend line and text
            columnspacing=0.6,  # Reduce space between columns
            handlelength=1.5,   # Adjust the length of the legend lines
            labelspacing=0.3    # Reduce space between legend entries
        )
        
        # Get current y-axis limits
        ymin, ymax = plt.ylim()

        # Round up to nearest 0.1 for the maximum y-limit
        ymax_rounded = np.ceil(ymax * 10) / 10

        # Set the y-axis limits with the rounded max
        plt.ylim(ymin, ymax_rounded)

        # Format y-tick labels to show only one decimal place
        ax = plt.gca()
        ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%.1f'))

        #plt.grid(True)
        plt.tight_layout()
        
        # Save the plot if requested
        if save_plot:
            # Save as PNG
            plot_filename = os.path.join(output_dir, f'{base_filename}.png')
            plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
            print(f"Plot saved to {plot_filename}")
            
            # # Also save as PDF
            # pdf_filename = os.path.join(output_dir, f'{base_filename}.pdf')
            # plt.savefig(pdf_filename, bbox_inches='tight')
            # print(f"Plot saved to {pdf_filename}")
        
        if show_plot:
            plt.show()
        # plt.close()
        


class PDFAnalyzer:
    def __init__(self, save_plot_data=False):
        self.molten_salts = []
        self.save_plot_data = save_plot_data

    def add_molten_salt(self, pdf_file, comp, source, temp, gamma_bc):
        salt = MoltenSaltPDF(comp, pdf_file, source, temp, gamma_bc)
        self.molten_salts.append(salt)

    def analyze_all(self):
        for salt in self.molten_salts:
            salt.analyze_pdf()
            if self.save_plot_data:
                salt.save_plot_data()

    def plot_all(self):
        for salt in self.molten_salts:
            salt.plot_pdf(show_plot=False)

def main():
    # Set to True to save plot data for all salts
    save_all_plot_data = True  # Change to False to disable automatic saving
    
    analyzer = PDFAnalyzer(save_plot_data=save_all_plot_data)

    # Add molten salts to analyze
    analyzer.add_molten_salt('PDF_LiCl.csv',"1.0LiCl",'Walz, 2019', 878, 4.10511)   # Walz, 2019; 878K
    # analyzer.add_molten_salt('PDF_NaCl.csv',"1.0NaCl",'Walz, 2019', 1074, 4.48028)   # Walz, 2019; 1074.15K
    # analyzer.add_molten_salt('NaCl_Lu.csv', "1.0NaCl", 'Lu, 2021', 1200, 4.48028)
    # analyzer.add_molten_salt('PDF_KCl.csv',"1.0KCl",'Walz, 2019', 1043, 4.47675)   # Walz, 2019; 1043K

    # analyzer.add_molten_salt('PDF_LiF.csv',"1.0LiF",'Walz, 2019', 1121, 3.28553)   # Walz, 2019; 1121K
    # analyzer.add_molten_salt('PDF_NaF.csv',"1.0NaF",'Walz, 2019', 1266, 5.22361)   # Walz, 2019; 1266.15K
    analyzer.add_molten_salt('PDF_KF.csv',"1.0KF",'Walz, 2019', 1131, 4.63533)   # Walz, 2019; 1131.15K
    # analyzer.add_molten_salt('MgCl_Lu.csv', "1.0MgCl2", 'Lu, 2021', 1100, 4.76796)
    # analyzer.add_molten_salt('PDF_MgCl2.csv',"1.0MgCl2",'McGreevy, 1987', 998, 4.76796)   # McGreevy, 1987; 998K
    analyzer.add_molten_salt('PDF_MgCl2_Roy.csv',"1.0MgCl2",'Roy, 2021', 1073, 4.76796)   # Roy, 2021; 1073K
    # analyzer.add_molten_salt('PDF_CaCl2.csv',"1.0CaCl2",'McGreevy, 1987', 1093, 7.72598)   # McGreevy, 1987; 1093K
    # analyzer.add_molten_salt('CaCl_Bu.csv',"1.0CaCl2",'Bu, 2021', 1073, 7.72598)   # Bu, 2022; 1073K
    # analyzer.add_molten_salt('PDF_SrCl2.csv',"1.0SrCl2",'McGreevy, 1987', 1198, 0)   # McGreevy, 1987; 998K
    analyzer.add_molten_salt('PDF_NaCl-UCl3.csv',"0.64NaCl-0.36UCl3",'Andersson, 2022', 1250, 2.5393)   # PDF Andersson, 2022
    # analyzer.add_molten_salt('PDF_NaCl-UCl3.csv',"0.64NaCl-0.36UCl3",'ANL-Andersson, 2022', 1250, 2.5393)   # PDF Andersson, 2022
    # analyzer.add_molten_salt('PDF_FLiBe_Grizzi.csv',"0.5LiF-0.5BeF2",'Sun, 2024', 900, 0)   # PDF Sun, 900 K, 2024
    # analyzer.add_molten_salt('PDF_FLiBe_Langford.csv',"0.66LiF-0.34BeF2",'Langford, 2022', 973, 1.90187)   # PDF Langford, 2022 (Cylindrical)
    analyzer.add_molten_salt('PDF_FLiBe_Fayfar.csv',"0.66LiF-0.34BeF2",'Fayfar, 2023', 973, 1.90187)   # PDF Fayfar
    # analyzer.add_molten_salt('PDF_FLiNa.csv',"0.6LiF-0.4NaF",'Grizzi, 2024', 973, 2.63857)   # PDF Grizzi, 900 K, 2024
    analyzer.add_molten_salt('PDF_FLiNaK.csv',"0.465LiF-0.115NaF-0.42KF",'Frandsen, 2020', 873, 2.26059)   # PDF Frandsen, 940 K, 2020
    analyzer.add_molten_salt('PDF_FMgNaK.csv',"0.345NaF-0.59KF-0.065MgF2",'Solano, 2021', 1073, 3.92263)   # PDF Frandsen, 940 K, 2020
    # analyzer.add_molten_salt('PDF_38MgCl2-21NaCl-41KCl.csv',"0.38MgCl2-0.21NaCl-0.41KCl",'Jiang, 2024', 750, 0)
    # analyzer.add_molten_salt('PDF_45MgCl2-33NaCl-22KCl.csv',"0.45MgCl2-0.33NaCl-0.22KCl",'Jiang, 2024', 750, 0)
    analyzer.add_molten_salt('0.4903NaCl-0.5097CaCl2_Wei.csv', "0.4903NaCl-0.5097CaCl2", 'Wei, 2022', 1023, 3.76913)
    analyzer.add_molten_salt('0.535NaCl-0.15CaCl2-0.315MgCl2_Wei.csv', "0.535NaCl-0.15CaCl2-0.315MgCl2", 'Wei, 2022', 1023, 3.52027)

    # # # # # # No thermal conductivity measurements for these
    # analyzer.add_molten_salt('0.637LiCl-0.363KCl_Jiang.csv',"0.637LiCl-0.363KCl",'Jiang, 2024', 750, 0)# Need to fix CSV
    # analyzer.add_molten_salt('PDF_LiF-NaF-UF4.csv',"0.5454LiF-0.3636NaF-0.091UF4",'Grizzi, 2024', 1473, 0)   # PDF Grizzi, 2024
    # analyzer.add_molten_salt('PDF_NaCl-KCl-ZnCl2_1073.csv',"0.22NaCl-0.393KCl-0.387ZnCl2",'Xi, 2024', 1073, 0)      # Do have measurements for this one
    # analyzer.add_molten_salt('PDF_NaCl-KCl-ZnCl2_373.csv', "0.22NaCl-0.393KCl-0.387ZnCl2", 'Xi, 2024', 373, 0)
    # # analyzer.add_molten_salt('PDF_0.4NaF-0.186KF-0.414AlF3.csv', "0.4NaF-0.186KF-0.414AlF3", 'Zhang, 2024', 1123, 0)
    # # analyzer.add_molten_salt('PDF_CuCl-CuCl2.csv', "0.5CuCl-0.5CuCl2", 'Raskovalov, 2018', 835, 0)  #Multiple Cu oxidation states make this difficult to calculate
    # analyzer.add_molten_salt('PDF_LiCl-CaCl2.csv', "0.7LiCl-0.3CaCl2", 'Liang, 2024', 1073, 0)
    # analyzer.add_molten_salt('PDF_MgCl2_Roy.csv', "1.0MgCl2", 'Roy, 2021', 1073, 4.76796)
    # analyzer.add_molten_salt('PDF_ZnCl2_Roy.csv', "1.0ZnCl2", 'Roy, 2021', 1073, 0)
    # analyzer.add_molten_salt('0.718KCl-0.282CaCl2_Wei.csv', "0.718KCl-0.282CaCl2", 'Wei, 2022', 1300, 0)
    # analyzer.add_molten_salt('0.417NaCl-0.525CaCl2-0.058KCl_Wei.csv', "0.417NaCl-0.525CaCl2-0.058KCl", 'Wei, 2022', 1023, 0)
    # analyzer.add_molten_salt('0.5NaCl-0.5KCl_Manga.csv', "0.5NaCl-0.5KCl", 'Manga, 2013', 1100, 4.32778)
    # analyzer.add_molten_salt('0.5LiCl-0.5KCl_Jiang.csv', "0.5LiCl-0.5KCl", 'Jiang, 2016', 727, 0)
    # analyzer.add_molten_salt('0.637LiCl-0.363KCl_Jiang.csv', "0.637LiCl-0.363KCl", 'Jiang, 2016', 673, 0)
    # analyzer.add_molten_salt('LiCl-KCl-CeCl_Fuller.csv', "0.571LiCl-0.397KCl-0.032CeCl", 'Fuller, 2022', 773, 0)
    # analyzer.add_molten_salt('LiCl-KCl-EuCl_Fuller.csv', "0.571LiCl-0.397KCl-0.032EuCl", 'Fuller, 2022', 773, 0)
    # analyzer.add_molten_salt('LiCl-KCl-SmCl_Fuller.csv', "0.571LiCl-0.397KCl-0.032SmCl", 'Fuller, 2022', 773, 0)
    # analyzer.add_molten_salt('0.5UCl-0.5KCl_Andersson.csv',"0.5KCl-0.5UCl3",'Andersson, 2024', 1250, 0)   # PDF Andersson, 2024
    # analyzer.add_molten_salt('0.15UCl-0.85KCl_Andersson.csv',"0.85KCl-0.15UCl3",'Andersson, 2024', 1250, 0)   # PDF Andersson, 2022
    # analyzer.add_molten_salt('0.25UCl-0.75KCl_Andersson.csv',"0.75KCl-0.25UCl3",'Andersson, 2024', 1250, 0)   # PDF Andersson, 2022
    # analyzer.add_molten_salt('0.35UCl-0.65KCl_Andersson.csv',"0.65KCl-0.35UCl3",'Andersson, 2024', 1250, 0)   # PDF Andersson, 2022
    # analyzer.add_molten_salt('PDF_ThF_Dai.csv',"1.0ThF4",'Dai, 2015', 1633, 0)   # PDF Dai, 2024
    # analyzer.add_molten_salt('UF4_1357K_Ocadiz-Flores_2021.csv',"1.0UF4",'Ocadiz-Flores, 2021', 1357, 0)   # PDF Dai, 2024

    # Perform analysis
    analyzer.analyze_all()

    # Plot results
    analyzer.plot_all()

if __name__ == "__main__":
    main()