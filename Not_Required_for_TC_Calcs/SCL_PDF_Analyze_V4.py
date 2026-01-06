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
from scipy.ndimage import gaussian_filter1d
import csv

def standardize_ion_pair(ion_pair):
    """
    Standardize the ion pair name based on oxidation states:
    - Returns cation-anion as is, swaps if input is anion-cation.
    - For cation-cation or anion-anion, sorts alphabetically.
    """
    element1, element2 = ion_pair.split('-')
    el1 = element(element1)
    el2 = element(element2)
    el1_states = el1.oxistates
    el2_states = el2.oxistates

    if el1_states[0] > 0 and el2_states[0] < 0:
        return f"{element1}-{element2}"  # cation-anion
    elif el1_states[0] < 0 and el2_states[0] > 0:
        return f"{element2}-{element1}"  # anion-cation
    else:
        # For cation-cation or anion-anion, sort alphabetically
        return '-'.join(sorted([element1, element2]))

def interval_cut(splines, filtered_crossing_points, x_range, num_pnts_interval):
    """
    Interpolate splines over intervals defined by filtered crossing points.
    Each ion pair's spline is evaluated over each interval and stored.
    """
    for ion_pair in splines:
        splines[ion_pair]['intervals'] = []

    for index, (start, end) in enumerate(zip(filtered_crossing_points[:-1], filtered_crossing_points[1:])):
        x_interval = np.linspace(x_range[start], x_range[end], num_pnts_interval)
        for ion_pair, spline in splines.items():
            y_interpolated = spline['spline'](x_interval)
            splines[ion_pair]['intervals'].append({
                'interval_number': index,
                'x': x_interval,
                'y': y_interpolated
            })
    return splines

def extract_and_average(x, y, x_0, width):
    """
    Extract y-values within a range centered at x_0 and return their average.
    Returns None if no values are found.
    """
    x_min = x_0 - width / 2
    x_max = x_0 + width / 2
    mask = (x >= x_min) & (x <= x_max)
    y_extracted = y[mask]
    if len(y_extracted) == 0:
        return None
    return np.mean(y_extracted)

def select_elements_at_intervals(array, interval):
    """
    Select elements from a sorted array at approximately regular intervals.
    Returns a list of selected elements.
    """
    selected_elements = []
    current_index = 0
    while current_index < len(array):
        selected_elements.append(array[current_index])
        next_index = np.argmin(np.abs(array - (array[current_index] + interval)))
        if next_index <= current_index:
            break
        current_index = next_index
    return selected_elements

def P_factor(NI, KF, PH, DF):
    """
    Calculate phonon (P_ph) and diffusion (P_dif) factors based on input parameters.
    Returns tuple (P_ph, P_dif).
    """
    P_ph = NI * KF * PH
    P_dif = NI * KF * (1 - PH) * DF
    return P_ph, P_dif

class IonPairPDF:
    """
    Represents a pair of ions and their pair distribution function (PDF) data.
    Provides methods for classifying the ion pair type and extracting peak/minimum features.
    """
    def __init__(self, ion_pair, x_values, y_values, weights):
        """
        Initialize an IonPairPDF object.
        Args:
            ion_pair (str): Ion pair string (e.g., 'Na-Cl').
            x_values (array-like): x-values (distances).
            y_values (array-like): y-values (PDF values).
            weights (dict): Dictionary of weights for all ion pairs.
        """
        self.ion_pair = standardize_ion_pair(ion_pair)
        self.x = np.array(x_values)
        self.y = np.array(y_values)
        self.weights = weights
        self.spline = self._create_spline()
        self.type = self._type_ion()
        self.peak, self.minima = self.find_first_peak_and_minimum(self.x, self.y)

    def _create_spline(self):
        """
        Create a linear interpolation spline for the PDF data.
        Returns:
            interp1d: Linear interpolation function.
        """
        return interp1d(self.x, self.y, kind='linear', bounds_error=False, fill_value=0)

    def _type_ion(self):
        """
        Classify the ion pair type based on oxidation states:
        Returns:
            str: One of 'ca', 'ac', 'cc_sim', 'cc_diff', 'aa_sim', 'aa_diff'.
        """
        element1, element2 = self.ion_pair.split('-')
        el1 = element(element1)
        el2 = element(element2)
        el1_states = el1.oxistates
        el2_states = el2.oxistates
        if el1_states[0] > 0 and el2_states[0] < 0:
            return "ca"  # cation-anion
        elif el1_states[0] < 0 and el2_states[0] > 0:
            return "ac"  # anion-cation
        elif (el1_states[0] > 0 and el2_states[0] > 0) and el1 == el2:
            return "cc_sim"  # similar cation-cation
        elif (el1_states[0] > 0 and el2_states[0] > 0) and el1 != el2:
            return "cc_diff"  # different cation-cation
        elif (el1_states[0] < 0 and el2_states[0] < 0) and el1 == el2:
            return "aa_sim"  # similar anion-anion
        elif (el1_states[0] < 0 and el2_states[0] < 0) and el1 != el2:
            return "aa_diff"  # different anion-anion

    def find_first_peak_and_minimum(self, x, y, any_type=False):
        """
        Find the first significant peak and the first minimum after the peak in the PDF.
        For non-cation-anion pairs, returns (None, None).
        Args:
            x (np.ndarray): x-values (distances).
            y (np.ndarray): y-values (PDF values).
            any_type (bool): If True, process all types; otherwise only 'ca'.
        Returns:
            tuple: ((peak_x, peak_y), (min_x, min_y)) or ((None, None), (None, None))
        """
        if self.type == "ca" or any_type:
            weight = self.weights.get(self.ion_pair, 1.0)
            y_weighted = y * weight
            peaks, _ = find_peaks(y_weighted, height=0)
            sorted_peaks = sorted(peaks, key=lambda i: y_weighted[i], reverse=True)
            first_peak_index = None
            max_value = max(y_weighted) if len(y_weighted) > 0 else 0
            for peak in sorted_peaks:
                if y_weighted[peak] > 0.5 * max_value and peak < len(x) // 2:
                    first_peak_index = peak
                    break
            if first_peak_index is not None:
                first_peak_x = x[first_peak_index]
                first_peak_y = y_weighted[first_peak_index]
            else:
                first_peak_x, first_peak_y = None, None
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
            first_peak_x, first_peak_y = None, None
            first_min_x, first_min_y = None, None
        return (first_peak_x, first_peak_y), (first_min_x, first_min_y)


class MoltenSaltPDF:
    """
    Represents a molten salt system and provides methods to load PDF data, parse composition,
    calculate weights, and interpolate PDF data for all ion pairs in the system.
    """
    def __init__(self, comp, pdf_file, source, temp, gamma_bc):
        """
        Initialize a MoltenSaltPDF object.
        Args:
            comp (str): Composition string (e.g., '0.5NaCl-0.5KCl').
            pdf_file (str): Filename of the PDF data CSV.
            source (str): Data source identifier.
            temp (float): Temperature.
            gamma_bc (float): Additional parameter (purpose defined by user context).
        """
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
        """
        Load PDF data for all ion pairs from a CSV file.
        Args:
            pdf_file (str): Filename of the PDF data CSV.
        Returns:
            dict: Mapping of ion pair names to IonPairPDF objects.
        """
        folder_name = "RDF_Plots\RDF_CSV"
        current_dir = os.getcwd()
        rdf_plots_folder = os.path.join(current_dir, folder_name)
        file_path = os.path.join(rdf_plots_folder, pdf_file)
        df = pd.read_csv(file_path, header=None)
        ion_pairs = {}
        for i in range(0, df.shape[1], 2):
            ion_pair = standardize_ion_pair(df.iloc[0, i])
            x_values = df.iloc[2:, i].dropna().astype(float).values
            y_values = df.iloc[2:, i+1].dropna().astype(float).values
            x_extended = np.linspace(0, x_values[0], num=int(x_values[0]) + 1)
            y_extended = np.zeros_like(x_extended)
            x_combined = np.concatenate((x_extended, x_values))
            y_combined = np.concatenate((y_extended, y_values))
            ion_pairs[ion_pair] = IonPairPDF(ion_pair, x_combined, y_combined, self.weights)
        return ion_pairs

    def parse_composition(self, composition_str):
        """
        Parse a composition string into a dictionary of salt fractions and ion counts.
        Args:
            composition_str (str): Composition string (e.g., '0.5NaCl-0.5KCl').
        Returns:
            tuple: (composition dict, ion_counts dict)
        """
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
        """
        Calculate normalized weights for all ion pairs based on composition and ion counts.
        Returns:
            dict: Mapping of ion pair names to weights.
        """
        element_concentration = {}
        for salt, fraction in self.composition.items():
            for element, count in self.ion_counts[salt].items():
                if element in element_concentration:
                    element_concentration[element] += fraction * count
                else:
                    element_concentration[element] = fraction * count
        total_concentration = sum(element_concentration.values())
        rel_conc = {element: concentration / total_concentration for element, concentration in element_concentration.items()}
        weights = {}
        for ion1, conc1 in rel_conc.items():
            for ion2, conc2 in rel_conc.items():
                ion_pair = standardize_ion_pair(f"{ion1}-{ion2}")
                weights[ion_pair] = conc1 * conc2
        sumweights = sum(weights.values())
        weights = {key: value/sumweights for key, value in weights.items()}
        return weights

    def interpolate_data(self, num_pnts_interval=100):
        """
        Interpolate all ion pair PDFs onto a common x-range and create weighted splines.
        Args:
            num_pnts_interval (int): Number of points in the interpolation interval.
        Returns:
            tuple: (interpolated_splines, weighted_splines, x_new)
        """
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
        """
        Analyze the PDF data for all cation-anion pairs in the system.
        Calculates transport factors, disruption integrals, and writes summary results to CSV.
        Also stores intermediate results for plotting.
        """
        print(f"\n###  {self.comp}   #############")
        csv_data = [self.comp]
        weighted_SCLs = 0
        weights_pair = 0

        for ion_pair_ca, pdf_ca in self.ion_pairs.items():
            if pdf_ca.type != "ca":
                continue
            print(f"Cation-anion pair: {ion_pair_ca}")
            peak_ca, min_ca = pdf_ca.peak, pdf_ca.minima
            cation, anion = ion_pair_ca.split('-')

            cc_peak = self._get_cc_peak(cation)
            y_combined, y_weighted_ca_total, sum_c_ca, y_weighted_ca, y_weighted_cc = self._combine_weighted_pdfs(cation, pdf_ca.x)
            transfer_points, transfer_points_indices = self._calculate_transfer_points(pdf_ca.x, peak_ca[0])
            # Get the concentration of the current cation-anion pair and total cation-anion concentration
            conc_ca_i = self.weights.get(ion_pair_ca, 0)
            sum_conc_ca = sum(wt for pair, wt in self.weights.items() 
                            if any(p in pair for p in [cation, anion]))
            
            KF, NI, DF, cc_transfer, PH = self._calculate_transport_factors(
                peak_ca, min_ca, y_combined, transfer_points_indices, 
                cc_peak, conc_ca_i, sum_conc_ca
            )
            RTE = conc_ca_i / sum_conc_ca if sum_conc_ca > 0 else 0
            print(f"DF: {DF}")
            S_y = self._calculate_disruption_integral(
                y_weighted_ca, y_weighted_cc, y_combined, transfer_points, 
                transfer_points_indices, NI, KF, DF, PH, cc_transfer
            )
            # Use only the x-values corresponding to the transfer points for integration
            x_points = pdf_ca.x[transfer_points_indices[:len(S_y)]]
            x_SCL_pair = np.trapz(S_y, x_points)
            print(f'x_SCL_pair: {x_SCL_pair}')
            print(f'Pair_weight: {RTE}')
            weights_pair += RTE
            weighted_SCLs += x_SCL_pair * RTE
            csv_data.extend([ion_pair_ca, x_SCL_pair, peak_ca[0], peak_ca[1], min_ca[0], min_ca[1]])
            self.plot_data.append({
                'x_range': pdf_ca.x,
                'x_SCL_pair': x_SCL_pair,
                'S_i': S_y,
                'ion_pair_ca_i': ion_pair_ca,
                'MFP_bc': self.gamma_bc,
            })
        avg_SCL = weighted_SCLs / weights_pair if weights_pair else 0
        self.plot_data.append({'avg_SCL': round(avg_SCL, 5)})
        print(f"Average SCL: {round(avg_SCL, 5)}")
        print(f"lambda_BC = {round(self.gamma_bc, 5)}")
        self._write_csv_results(csv_data)

    def _get_cc_peak(self, cation):
        """Return the cation-cation peak for a given cation, or None if not found."""
        ion_pair_cc = f"{cation}-{cation}"
        if ion_pair_cc in self.ion_pairs:
            pdf_cc = self.ion_pairs[ion_pair_cc]
            cc_peak, _ = pdf_cc.find_first_peak_and_minimum(pdf_cc.x, pdf_cc.y, any_type=True)
            print(f"Cation-cation pair: {ion_pair_cc}")
            if cc_peak[0] is not None:
                print(f"First cc peak: x = {cc_peak[0]:.2f}, y = {cc_peak[1]:.2f}")
            else:
                print("No peak found for cation-cation pair")
            return cc_peak
        else:
            print(f"No data found for cation-cation pair: {ion_pair_cc}")
            return None

    def _combine_weighted_pdfs(self, cation, x_range):
        """Combine weighted PDFs for all relevant ion pairs."""
        y_combined = np.zeros(len(x_range))
        y_weighted_ca_total = np.zeros(len(x_range))
        sum_c_ca = 0
        y_weighted_ca = np.zeros(len(x_range))
        y_weighted_cc = np.zeros(len(x_range))
        for ion_pair, pdf in self.ion_pairs.items():
            ion1, ion2 = ion_pair.split('-')
            y_weighted = self.weighted_splines[ion_pair](x_range)
            if pdf.type == "ca":
                sum_c_ca += self.weights[ion_pair]
                y_weighted_ca_total += y_weighted
                if cation in [ion1, ion2]:
                    y_weighted_ca = y_weighted
            if cation in [ion1, ion2]:
                y_combined += y_weighted
                if ion1 == cation and ion2 == cation:
                    y_weighted_cc = y_weighted
        return y_combined, y_weighted_ca_total, sum_c_ca, y_weighted_ca, y_weighted_cc

    def _calculate_transfer_points(self, x_range, peak_x):
        """Calculate transfer points and their indices."""
        transfer_points = select_elements_at_intervals(x_range, peak_x)
        transfer_points_indices = [np.where(x_range == point)[0][0] for point in transfer_points]
        return transfer_points, transfer_points_indices

    def _calculate_transport_factors(self, peak, minima, y_combined, transfer_points_indices, cc_peak, conc_ca_i, sum_conc_ca):
        """Calculate transport factors KF, NI, DF, and cc_transfer flag.
        
        Optimized version using NumPy vectorization and pre-computing values.
        
        Args:
            peak: (x, y) tuple of the main peak
            minima: (x, y) tuple of the first minimum after the peak
            y_combined: Combined y-values from all relevant ion pairs
            transfer_points_indices: Indices of transfer points
            cc_peak: (x, y) tuple of the cation-cation peak (if any)
            conc_ca_i: Concentration of the current cation-anion pair
            sum_conc_ca: Sum of all cation-anion pair concentrations
            
        Returns:
            tuple: (KF, NI, DF, cc_transfer, PH)
        """
        # Pre-compute peak values to avoid multiple lookups
        peak_y = peak[1] if peak and peak[1] is not None else 0
        min_y = minima[1] if minima and minima[1] is not None else 0
        
        # Bond strength/diffusion factor (KF) - vectorized
        KF = abs(peak_y - min_y) / peak_y if peak_y > 0 else 0
        
        # Non-ideal recipient factor (NI) - vectorized
        y_combined_at_transfer = y_combined[transfer_points_indices[1]] if len(transfer_points_indices) > 1 else 0
        NI = peak_y / y_combined_at_transfer if y_combined_at_transfer > 0 else 0
        
        # Diffusion factor (DF) - only relevant if cation-cation transfer is possible
        DF = 0
        cc_transfer = False
        if (cc_peak is not None and cc_peak[0] and peak[0] and 
            cc_peak[0] > peak[0] and cc_peak[0] < transfer_points_indices[1]):
            DF = ((peak[0] / cc_peak[0]) ** 3 - (1 / 2) ** 3) / (1 - (1 / 2) ** 3)
            cc_transfer = True
            
        # Phonon factor (PH) - concentration-based as in V3
        PH = conc_ca_i / sum_conc_ca if sum_conc_ca > 0 else 0
        
        return KF, NI, DF, cc_transfer, PH

    def _calculate_disruption_integral(self, y_weighted_ca, y_weighted_cc, y_combined, transfer_points, 
                                    transfer_points_indices, NI, KF, DF, PH, cc_transfer):
        """Calculate the disruption integral S_y for transfer points.
        
        Args:
            y_weighted_ca: Weighted y-values for cation-anion pairs
            y_weighted_cc: Weighted y-values for cation-cation pairs
            y_combined: Combined y-values from all relevant ion pairs
            transfer_points: Array of transfer points (x-coordinates)
            transfer_points_indices: Indices of transfer points in the x-range
            NI: Non-ideal recipient factor
            KF: Bond strength/diffusion factor
            DF: Diffusion factor
            PH: Phonon factor
            cc_transfer: Boolean indicating if cation-cation transfer is active
            
        Returns:
            numpy.ndarray: Array of S_y values for each transfer point
        """
        # Initialize arrays for different disruption factors
        b_KF = np.ones(len(transfer_points))
        b_NI = np.ones(len(transfer_points))
        b_PC = np.ones(len(transfer_points))
        b_SRO = np.ones(len(transfer_points))
        
        # Calculate disruption factors for each transfer point
        for i in range(len(transfer_points)):
            if i == 0:  # Initial energy carrier density - cation
                pass
            else:
                r_i = transfer_points_indices[i]
                r_i_prev = transfer_points_indices[i-1]
                
                # Apply different disruption factors based on transfer point index
                if i % 2 == 0 and i == 2:  # First anion-to-cation transfer
                    b_KF[i] = 1 - KF
                    b_NI[i] = 1 - NI
                    b_PC[i] = 1 - PH
                elif i % 2 != 0 and i > 1:  # n cation-to-anion transfer (n > 1)
                    b_KF[i] = 1 - KF
                    b_NI[i] = 1 - (y_weighted_ca[r_i] / y_combined[r_i]) if y_combined[r_i] != 0 else 0
                    b_PC[i] = 1 - PH
                    b_SRO[i] = 1 - (y_weighted_ca[r_i] / y_weighted_cc[r_i_prev]) if y_weighted_cc[r_i_prev] != 0 else 0
                elif i % 2 == 0 and i > 2:  # n; anion-to-cation transfer (n > 1)
                    b_KF[i] = 1 - KF
                    b_NI[i] = 1 - (y_weighted_cc[r_i] / y_combined[r_i]) if y_combined[r_i] != 0 else 0
                    b_PC[i] = 1 - PH
                    b_SRO[i] = 1 - (y_weighted_cc[r_i] / y_weighted_ca[r_i_prev]) if y_weighted_ca[r_i_prev] != 0 else 0
        
        # Calculate beta_i and S_y
        S_y = np.zeros(len(transfer_points))
        beta_i_integral = 0
        
        for i in range(len(transfer_points) - 1):
            # Calculate combined disruption factor (product of all factors)
            b_disruption = b_KF[i] + b_NI[i] + b_PC[i] + b_SRO[i]
            
            # Ensure b_disruption is within [0, 1]
            b_disruption = max(0, min(1, b_disruption))
            
            # Calculate beta_i and update the integral
            beta_i = b_disruption / (1 - b_disruption) if b_disruption < 1 else 0
            beta_i_integral += beta_i * (transfer_points[i+1] - transfer_points[i])
            
            # Calculate S_y
            S_y[i] = np.exp(-beta_i_integral)
        
        return S_y

    def _write_csv_results(self, csv_data):
        """Write the results to a CSV file."""
        csv_filename = 'SCL_PDF_results.csv'
        file_exists = os.path.isfile(csv_filename)
        with open(csv_filename, mode='a', newline='') as file:
            writer = csv.writer(file)
            if not file_exists:
                header = [
                    'Composition', 'Pair1', 'x-SCL1', 'Peak-x1', 'Peak-y1', 'Min-x1', 'Min-y1',
                    'Pair2', 'x-SCL2', 'Peak-x2', 'Peak-y2', 'Min-x2', 'Min-y2',
                    'Pair3', 'x-SCL3', 'Peak-x3', 'Peak-y3', 'Min-x3', 'Min-y3',
                    'Pair4', 'x-SCL4', 'Peak-x4', 'Peak-y4', 'Min-x4', 'Min-y4',
                    'Pair5', 'x-SCL5', 'Peak-x5', 'Peak-y5', 'Min-x5', 'Min-y5'
                ]
                writer.writerow(header)
            writer.writerow(csv_data)

    def plot_pdf(self):
        plt.figure(figsize=(8, 6))
        ion_pair_colors = {}  # Dictionary to store colors for each ion pair

        # Count the number of cation-anion pairs
        ca_pair_count = sum(1 for pdf_data in self.ion_pairs.values() if pdf_data.type == "ca")


        for ion_pair, pdf_data in self.ion_pairs.items():
            #plt.plot(pdf_data.x, pdf_data.y, label=f"{ion_pair} (raw)", alpha=0.5)
            
            if ion_pair in self.weighted_splines:
                weighted_y = self.weighted_splines[ion_pair](self.x_new)
                spline, = plt.plot(self.x_new, weighted_y, label=f"{ion_pair}")
                spline_color = spline.get_color()
                ion_pair_colors[ion_pair] = spline_color  # Store the color for this ion pair
                #plt.axhline(y=self.weights[ion_pair], color=spline_color, linestyle='--')

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

        for data in self.plot_data:
            # Check if 'ion_pair_ca_i' key exists
            if 'ion_pair_ca_i' in data:
                ion_pair = data['ion_pair_ca_i']
                color = ion_pair_colors.get(ion_pair, 'black')  # Default to black if not found
                # Use transfer points for x-values when plotting S_i
                x_points = data['x_range'][:len(data['S_i'])]
                plt.plot(x_points, data['S_i'], label=f"S(r): {ion_pair}", color=color, linestyle='dotted',)
                # Only plot the line if there is more than one cation-anion pair
                # if ca_pair_count > 1:
                #     plt.axvline(x=data['x_SCL_pair'], color='gray', linestyle='--', label=f"$\\lambda_{{{ion_pair}}}$ = {round(data['x_SCL_pair'], 2)}")
            else:
                print("Warning: 'ion_pair_ca_i' key not found in data entry:", data)

        # Use dark green for the average SCL line
        plt.axvline(x=data['avg_SCL'], color='g', linestyle='-.', label=f"$\lambda_{{SCL}}$: = {round(data['avg_SCL'], 2)}")

        if self.gamma_bc > 0:
            plt.axvline(x=self.gamma_bc, color='k', linestyle='--', label=f"$\lambda_{{BC}}$: = {round(self.gamma_bc,2)}")

        current_date = datetime.now().strftime("%Y%m%d")
        save_title = f'PDF_{self.comp}_{current_date}.png'
        
        folder_name = "V4"
        os.makedirs(os.path.join(os.getcwd(), folder_name), exist_ok=True)
        file_path = os.path.join(os.getcwd(), folder_name, save_title)

        plt.xlabel('r [Å]')
        plt.ylabel('g(r)')
        plt.title(f'Partial PDF for {self.comp} - {self.source}')
        #plt.xlim(0,5)
        #plt.ylim(0,1.3)
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
    analyzer.add_molten_salt('PDF_FMgNaK.csv',"0.345NaF-0.59KF-0.065MgF2",'Solano, 2021', 1073, 3.95351295352575)   # PDF Frandsen, 940 K, 2020
    analyzer.add_molten_salt('PDF_LiF-NaF-UF4.csv',"0.5454LiF-0.3636NaF-0.091UF4",'Grizzi, 2024', 1473, 0)   # PDF Grizzi, 2024
    analyzer.add_molten_salt('PDF_NaCl-KCl-ZnCl2_1073.csv',"0.22NaCl-0.393KCl-0.387ZnCl2",'Xi, 2024', 1073, 0)      # Do have measurements for this one
    analyzer.add_molten_salt('PDF_38MgCl2-21NaCl-41KCl.csv',"0.38MgCl2-0.21NaCl-0.41KCl",'Jiang, 2024', 750, 3.93912648093083)
    analyzer.add_molten_salt('PDF_45MgCl2-33NaCl-22KCl.csv',"0.45MgCl2-0.33NaCl-0.22KCl",'Jiang, 2024', 750, 3.80599721669684)


    # Perform analysis
    analyzer.analyze_all()

    # Plot results
    analyzer.plot_all()

if __name__ == "__main__":
    main()