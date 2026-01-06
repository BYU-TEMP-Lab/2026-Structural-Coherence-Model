import os
import pandas as pd
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import re
from mendeleev import element
from scipy.signal import find_peaks
from datetime import datetime

def standardize_ion_pair(ion_pair):
    element1, element2 = ion_pair.split('-')
    el1 = element(element1)
    el2 = element(element2)
    el1_states = el1.oxistates
    el2_states = el2.oxistates

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
        self.peak,self.minima = self.find_first_peak_and_minimum(self.x, self.y)

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
        
    def find_first_peak_and_minimum(self, x, y):
        if self.type == "ca":
            # Get the weight for this ion pair
            weight = self.weights.get(self.ion_pair, 1.0)  # Default to 1.0 if not found
            
            # Weight the y values
            y_weighted = y * weight

            # Find the maximum value in y_weighted
            max_index = np.argmax(y_weighted)
            max_value = y_weighted[max_index]

            # Check if the maximum value is a peak
            peaks, _ = find_peaks(y_weighted)
            if max_index in peaks:
                first_peak_index = max_index
                first_peak_x = x[first_peak_index]
                first_peak_y = max_value
            else:
                first_peak_x, first_peak_y = None, None

            # Find the first minimum after the first peak
            if first_peak_x is not None:
                after_peak = y_weighted[first_peak_index:]
                minima, _ = find_peaks(-after_peak)
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
        self.composition, self.ion_counts = self.parse_composition(comp)
        self.weights = self.calculate_weights()
        self.ion_pairs = self.load_pdf_data(pdf_file)
        self.interpolated_splines, self.weighted_splines, self.x_new = self.interpolate_data()
        self.plot_data = []  # Store plot data here

    def load_pdf_data(self, pdf_file):
        folder_name = "RDF_Plots\RDF_CSV"
        current_dir = os.getcwd()
        rdf_plots_folder = os.path.join(current_dir, folder_name)
        file_path = os.path.join(rdf_plots_folder, pdf_file)
        
        df = pd.read_csv(file_path, header=None)
        ion_pairs = {}

        for i in range(0, df.shape[1], 2):
            ion_pair = standardize_ion_pair(df.iloc[0, i])  # Standardize the ion pair name
            x_values = df.iloc[2:, i].dropna().astype(float).values
            y_values = df.iloc[2:, i+1].dropna().astype(float).values

            # Extend x_values and y_values
            x_extended = np.linspace(0, x_values[0], num=int(x_values[0]) + 1)
            y_extended = np.zeros_like(x_extended)
            
            x_combined = np.concatenate((x_extended, x_values))
            y_combined = np.concatenate((y_extended, y_values))

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
        return composition, ion_counts

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

    def interpolate_data(self, num_pnts_interval=100):
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

        for ion_pair_ca_i, pdf_data_ca_i in self.ion_pairs.items():
            if pdf_data_ca_i.type == "ca":
                print(f"Cation-anion pair: {ion_pair_ca_i}")

                # Access peak and minima
                peak_i = pdf_data_ca_i.peak
                minima_i = pdf_data_ca_i.minima

                x_range_i = pdf_data_ca_i.x
                

                cation_i, anion_i = ion_pair_ca_i.split('-')
                y_combined = np.zeros(len(x_range_i))

                conc_ca_i = self.weights[ion_pair_ca_i]
                y_weighted_ca_i = self.weighted_splines[ion_pair_ca_i](x_range_i)
                c_ca_i = conc_ca_i * peak_i[1] * peak_i[0]**-1
                sum_c_ca_j = 0

                ion_pair_cc_i = '-'.join(sorted([cation_i, cation_i]))
                no_cc_spline = False
                try:
                    y_weighted_cc_i = self.weighted_splines[ion_pair_cc_i](x_range_i)
                except KeyError:
                    y_weighted_cc_i = np.linspace(self.weights[ion_pair_cc_i], self.weights[ion_pair_cc_i], len(x_range_i))
                    no_cc_spline = True
                    print("No cation-cation spline data available. Approximating based on concentration.")

                y_ideal = np.zeros(len(x_range_i))

                # Calculate peak properties, combined spline, etc.
                for ion_pair_j, pdf_data_j in self.ion_pairs.items():
                    ion1, ion2 = ion_pair_j.split('-')
                    y_weighted_j = self.weighted_splines[ion_pair_j](x_range_i)

                    if pdf_data_j.type == "ca":
                        # Access peak and minima
                        peak_j = pdf_data_j.peak
                        minima_j = pdf_data_j.minima
                        conc_ca_j = self.weights[ion_pair_j]
                        sum_c_ca_j += conc_ca_j * peak_j[1] * peak_j[0]**-1

                    if cation_i in [ion1, ion2]:
                        y_combined += y_weighted_j
                    
                # Create a spline for y_combined
                y_combined_spline = interp1d(x_range_i, y_combined, kind='linear', bounds_error=False, fill_value=0)

                transfer_points = select_elements_at_intervals(x_range_i, peak_i[0])
                transfer_points_indices = [np.where(x_range_i == point)[0][0] for point in transfer_points]

                print(f"Transfer points: {transfer_points}")

                b_ph = np.zeros(len(transfer_points))
                b_dif = np.zeros(len(transfer_points))
                b_sro = np.zeros(len(transfer_points))

                # Diffusion factor
                dif_f = ((peak_i[1]-minima_i[1])/peak_i[1])

                # Phonon transfer factor
                ph_f = (peak_i[1]/y_combined[transfer_points_indices[1]])

                for i in range(1,len(transfer_points)-1):
                    r_i = transfer_points_indices[i]
                    r_i_prev = transfer_points_indices[i-1]
                    if i % 2 != 0 and i == 1:
                        b_ph[i] += 1 - ph_f
                        b_dif[i] += 1 - dif_f
                        b_dif[i+1] -= b_dif[i]
                        #b_sro[i] = 1 - (y_weighted_cc_i[r_i]/y_weighted_ca_i[r_i_prev])
                    elif i % 2 == 0:
                        b_ph[i] += 1 - ph_f
                        b_dif[i] += 1 - dif_f
                        b_ph[i] += (1 - (y_weighted_cc_i[r_i]/y_weighted_ca_i[r_i_prev]) *ph_f *dif_f )
                        b_dif[i+1] -= b_dif[i]
                    elif i % 2 != 0 and i != 1:
                        b_ph[i] += 1 - ph_f
                        b_dif[i] += 1 - dif_f
                        b_ph[i] += (1 - (y_weighted_ca_i[r_i]/y_weighted_cc_i[r_i_prev]) *ph_f *dif_f )
                        b_dif[i+1] -= b_dif[i]

                    # Ensure b_ph, b_dif, and b_sro are within [0, 1]
                    b_ph[i] = min(max(b_ph[i], 0), 1)
                    b_dif[i] = min(max(b_dif[i], 0), 1)
                    b_sro[i] = min(max(b_sro[i], 0), 1)

                beta_i = []
                beta_i_integral = np.zeros(len(x_range_i))

                b_array = [b_ph,b_dif,b_sro]#, b_I] #b_SRO]

                P_i = np.zeros(len(x_range_i))
                
                for k in range(len(transfer_points_indices)):
                    b_i_sum = 0
                    for b_i in b_array:
                        if b_i[k] == 1:
                            b_i_sum += (float('inf'))
                        else:
                            b_i_sum += b_i[k] / (1 - b_i[k])
                    beta_i.append(b_i_sum)

                print(f'beta_i: {beta_i}')

                beta_x = []
                beta_y = [0]  # Initialize beta_y with a default value

                for i in range(len(x_range_i)):
                    beta_x.append(x_range_i[i])
                    for k in range(1, len(transfer_points_indices) - 1):
                        if transfer_points_indices[k] > i >= transfer_points_indices[k-1]:
                            beta_y.append(beta_i[k])
                            break  # Exit the loop once the correct beta_i is found

                    # Ensure beta_y is populated for each x_range_i
                    if len(beta_y) < len(beta_x):
                        if beta_y:  # Check if beta_y is not empty
                            beta_y.append(beta_y[-1])  # Repeat the last value to match the length
                        else:
                            beta_y.append(0)  # Append a default value if beta_y is empty

                    beta_i_integral[i] = np.trapz(beta_y[:i + 1], beta_x[:i + 1])
                    P_i[i] = np.exp(-beta_i_integral[i])
                
                # Debugging output
                print(f"Length of x_range_i: {len(x_range_i)}")
                print(f"Length of beta_x: {len(beta_x)}")
                print(f"Length of beta_y: {len(beta_y)}")
                print(f"Length of P_i: {len(P_i)}")

                P_exp = np.trapz(P_i * x_range_i, P_i) / np.trapz(x_range_i, P_i)
                index = np.argmin(np.abs(P_i - P_exp))
                x_SCL_pair = x_range_i[index]
                print(f'x_SCL_pair: {x_SCL_pair}')

                


                # Store plot data
                self.plot_data.append({
                    # 'x_range': x_range_i,
                    # 'y_ideal': y_ideal,
                    # 'y_combined': y_combined,
                    # 'beta_y': beta_y,
                    # 'beta_i': beta_i,
                    # 'beta_i_integral': beta_i_integral,
                    # 'P_i': P_i,
                    'x_range': x_range_i,
                    'x_SCL_pair': x_SCL_pair,
                    'P_i': P_i,
                    'ion_pair_ca_i': ion_pair_ca_i,
                    'MFP_bc': self.gamma_bc,
                })

                # print(f"lambda_P(r) - {ion_pair_ca_i} = {round(x_SCL_pair, 5)}")
        print(f"lambda_BC = {round(self.gamma_bc, 5)}")

    def plot_pdf(self):
        plt.figure(figsize=(12, 8))
        for ion_pair, pdf_data in self.ion_pairs.items():
            #plt.plot(pdf_data.x, pdf_data.y, label=f"{ion_pair} (raw)", alpha=0.5)
            
            if ion_pair in self.weighted_splines:
                weighted_y = self.weighted_splines[ion_pair](self.x_new)
                spline, = plt.plot(self.x_new, weighted_y, label=f"{ion_pair}")
                spline_color = spline.get_color()
                #plt.axhline(y=self.weights[ion_pair], color=spline_color, linestyle='--')

            # Access peak and minima
            peak = pdf_data.peak
            minima = pdf_data.minima

            if peak[0] is not None:
                peakpoint, = plt.plot(peak[0], peak[1], 'ro', markersize=10)
                plt.axvline(x=2*peak[0], color=peakpoint.get_color(),linestyle='dotted')
                plt.annotate(f"Peak: ({peak[0]:.2f}, {peak[1]:.2f})", 
                            (peak[0], peak[1]), 
                            xytext=(5, 5), 
                            textcoords='offset points')

            if minima[0] is not None:
                plt.plot(minima[0], minima[1], 'bo', markersize=10)
                plt.annotate(f"Min: ({minima[0]:.2f}, {minima[1]:.2f})", 
                            (minima[0], minima[1]), 
                            xytext=(5, 5), 
                            textcoords='offset points')
        #for data in self.plot_data:
            #plt.plot(data['x_range'], data['P_i'], label=f"P(r): {data['ion_pair_ca_i']}")
            # plt.plot(data['x_range'], data['y_ideal'], linestyle=':', label=fr"$g_{{Ideal}}(r)$: {data['ion_pair_ca_i']}")
            # plt.plot(data['x_range'], data['y_combined'], linestyle=':', label=fr"$g_{{Tot}}(r)$: {data['ion_pair_ca_i']}")
            #plt.plot(data['x_range'], data['beta_y'], linestyle='-.', label=fr"$\beta(r)$: {data['ion_pair_ca_i']}")
            # plt.plot(data['x_range'], data['beta_i_integral'], linestyle='-.', label=fr"$\beta_int(r)$: {data['ion_pair_ca_i']}")
            #   plt.plot(data['x_range'], data['P_i'], linestyle='-.', label=f"P(r): {data['ion_pair_ca_i']}")
            #plt.axvline(x=data['x_SCL_pair'], color='g', linestyle='--', label=f"$\lambda_{{P(r)}}$: {data['ion_pair_ca_i']} = {round(data['x_SCL_pair'], 2)}")

        if self.gamma_bc > 0:
            plt.axvline(x=self.gamma_bc, color='k', linestyle='--', label=f"$\lambda_{{BC}}$: = {round(self.gamma_bc,2)}")

        current_date = datetime.now().strftime("%Y%m%d")
        save_title = f'PDF_{self.comp}_{current_date}.png'
        
        folder_name = "V1_Plots_Blank"
        os.makedirs(os.path.join(os.getcwd(), folder_name), exist_ok=True)
        file_path = os.path.join(os.getcwd(), folder_name, save_title)

        plt.xlabel('r [Å]')
        plt.ylabel('g(r)')
        plt.title(f'Partial PDF for {self.comp} - {self.source}')
        plt.legend()
        #plt.grid(True)
        plt.savefig(file_path)
        plt.close()
        


class PDFAnalyzer:
    def __init__(self):
        self.molten_salts = []

    def add_molten_salt(self, pdf_file, comp, source, temp, gamma_bc):
        salt = MoltenSaltPDF(comp, pdf_file, source, temp, gamma_bc)
        self.molten_salts.append(salt)

    def analyze_all(self):
        for salt in self.molten_salts:
            salt.analyze_pdf()

    def plot_all(self):
        for salt in self.molten_salts:
            salt.plot_pdf()

def main():
    analyzer = PDFAnalyzer()

    # Add molten salts to analyze
    analyzer.add_molten_salt('PDF_LiCl.csv',"1.0LiCl",'Walz, 2019', 878, 3.893446)   # Walz, 2019; 878K
    analyzer.add_molten_salt('PDF_NaCl.csv',"1.0NaCl",'Walz, 2019', 1074, 4.352426)   # Walz, 2019; 1074.15K
    analyzer.add_molten_salt('PDF_KCl.csv',"1.0KCl",'Walz, 2019', 1043, 4.453593)   # Walz, 2019; 1043K
    analyzer.add_molten_salt('PDF_LiF.csv',"1.0LiF",'Walz, 2019', 1121, 3.305228)   # Walz, 2019; 1121K
    analyzer.add_molten_salt('PDF_NaF.csv',"1.0NaF",'Walz, 2019', 1266, 4.005403)   # Walz, 2019; 1266.15K
    analyzer.add_molten_salt('PDF_KF.csv',"1.0KF",'Walz, 2019', 1131, 4.225170)   # Walz, 2019; 1131.15K
    analyzer.add_molten_salt('PDF_MgCl2.csv',"1.0MgCl2",'McGreevy, 1987', 998, 4.049956)   # McGreevy, 1987; 998K
    analyzer.add_molten_salt('PDF_CaCl2.csv',"1.0CaCl2",'McGreevy, 1987', 1093, 3.372692)   # McGreevy, 1987; 1093K
    analyzer.add_molten_salt('PDF_SrCl2.csv',"1.0SrCl2",'McGreevy, 1987', 1198, 3.388477)   # McGreevy, 1987; 998K
    analyzer.add_molten_salt('PDF_NaCl-UCl3.csv',"0.64NaCl-0.36UCl3",'Andersson, 2022', 1250, 2.455254)   # PDF Andersson, 2022
    analyzer.add_molten_salt('PDF_NaCl-UCl3.csv',"0.64NaCl-0.36UCl3",'ANL-Andersson, 2022', 1250, 2.717874)   # PDF Andersson, 2022
    analyzer.add_molten_salt('PDF_FLiBe_Grizzi.csv',"0.5LiF-0.5BeF2",'Sun, 2024', 900, 1.922959)   # PDF Sun, 900 K, 2024
    analyzer.add_molten_salt('PDF_FLiBe_Langford.csv',"0.66LiF-0.34BeF2",'Langford, 2022', 973, 1.922959)   # PDF Langford, 2022 (Cylindrical)
    analyzer.add_molten_salt('PDF_FLiBe_Fayfar.csv',"0.66LiF-0.34BeF2",'Fayfar, 2023', 973, 1.922959)   # PDF Fayfar
    analyzer.add_molten_salt('PDF_FLiNa.csv',"0.6LiF-0.4NaF",'Grizzi, 2024', 973, 2.378058)   # PDF Grizzi, 900 K, 2024
    analyzer.add_molten_salt('PDF_FLiNaK.csv',"0.465LiF-0.115NaF-0.42KF",'Frandsen, 2020', 873, 2.421498)   # PDF Frandsen, 940 K, 2020
    analyzer.add_molten_salt('PDF_FMgNaK.csv',"0.345NaF-0.59KF-0.065MgF2",'Solano, 2021', 1073, 0)   # PDF Frandsen, 940 K, 2020
    analyzer.add_molten_salt('PDF_LiF-NaF-UF4.csv',"0.5454LiF-0.3636NaF-0.091UF4",'Grizzi, 2024', 1473, 0)   # PDF Grizzi, 2024
    analyzer.add_molten_salt('PDF_NaCl-KCl-ZnCl2_1073.csv',"0.22NaCl-0.393KCl-0.387ZnCl2",'Xi, 2024', 1073, 0)

    # Perform analysis
    analyzer.analyze_all()

    # Plot results
    analyzer.plot_all()

if __name__ == "__main__":
    main()