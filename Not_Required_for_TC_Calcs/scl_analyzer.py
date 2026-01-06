import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mendeleev import element
from scipy.interpolate import interp1d
from scipy.signal import find_peaks
from scipy.integrate import simps
from scipy.optimize import minimize_scalar
from matplotlib.gridspec import GridSpec

# Type aliases
ArrayLike = Union[np.ndarray, List[float]]
IonPair = str
PDFDataDict = Dict[IonPair, 'PDFData']

def get_script_dir() -> Path:
    """Get the directory where the script is located."""
    return Path(__file__).parent.absolute()

# Define paths
SCRIPT_DIR = get_script_dir()
DATA_DIR = SCRIPT_DIR / "RDF_Plots" / "RDF_CSV"
RESULTS_DIR = SCRIPT_DIR / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

@dataclass
class PDFData:
    """Class to store and analyze PDF data for an ion pair.
    
    Attributes:
        ion_pair: String identifier for the ion pair (e.g., 'Li-F')
        distances: Array of distance values (r) in Angstroms
        values: Array of PDF values (g(r))
        weight: Weight factor for this ion pair (default: 1.0)
    """
    ion_pair: str
    distances: np.ndarray
    values: np.ndarray
    weight: float = 1.0
    
    # Cached properties
    _pair_type: Optional[str] = field(init=False, default=None)
    _peak_info: Optional[dict] = field(init=False, default=None)
    _spline: Optional[interp1d] = field(init=False, default=None)
    
    def __post_init__(self):
        # Ensure numpy arrays and apply weight
        self.distances = np.asarray(self.distances, dtype=float)
        self.values = np.asarray(self.values, dtype=float) * self.weight
    
    @property
    def spline(self) -> interp1d:
        """Get or create a spline interpolation of the PDF data."""
        if self._spline is None:
            self._spline = interp1d(
                self.distances, 
                self.values, 
                kind='cubic',
                bounds_error=False,
                fill_value=0.0
            )
        return self._spline
    
    @property
    def pair_type(self) -> str:
        """Determine if the ion pair is cation-anion, anion-cation, or other."""
        if self._pair_type is None:
            elements = self.ion_pair.split('-')
            if len(elements) != 2:
                self._pair_type = 'other'
                return self._pair_type
                
            try:
                # Extract element symbols (handle cases like 'Na1+' -> 'Na')
                elem1 = ''.join(filter(str.isalpha, elements[0]))
                elem2 = ''.join(filter(str.isalpha, elements[1]))
                
                # Get element objects
                el1 = element(elem1.title())
                el2 = element(elem2.title())
                
                # Get electronegativity (prefer Pauling, fallback to Allen)
                en1 = getattr(el1, 'en_pauling', getattr(el1, 'en_allen', 0.0))
                en2 = getattr(el2, 'en_pauling', getattr(el2, 'en_allen', 0.0))
                
                # Determine ion types based on electronegativity difference
                en_diff = abs(en1 - en2)
                
                # Classify pair type
                if en_diff > 0.5:  # Threshold for ionic/polar bonds
                    self._pair_type = 'ac' if en1 > en2 else 'ca'  # anion-cation or cation-anion
                else:
                    self._pair_type = 'other'  # Not a cation-anion pair
                        
            except Exception as e:
                # Fallback to simple naming if mendeleev fails
                if elements[0].startswith(('Li', 'Na', 'K', 'Rb', 'Cs')) and elements[1].startswith(('F', 'Cl', 'Br', 'I')):
                    self._pair_type = 'ca'  # cation-anion
                elif elements[1].startswith(('Li', 'Na', 'K', 'Rb', 'Cs')) and elements[0].startswith(('F', 'Cl', 'Br', 'I')):
                    self._pair_type = 'ac'  # anion-cation
                else:
                    self._pair_type = 'other'  # Not a cation-anion pair
        
        return self._pair_type
    
    @property
    def peak_info(self) -> dict:
        """Find and cache peak information."""
        if self._peak_info is None:
            # Find peaks using scipy's find_peaks
            peaks, properties = find_peaks(
                self.values,
                height=0.1,  # Minimum peak height
                prominence=0.1,  # Minimum peak prominence
                distance=10  # Minimum distance between peaks
            )
            
            if len(peaks) > 0:
                # Get the first significant peak
                peak_idx = peaks[0]
                peak_dist = self.distances[peak_idx]
                peak_val = self.values[peak_idx]
                
                # Find the first minimum after the peak
                post_peak = self.values[peak_idx:]
                min_idx = np.argmin(post_peak) + peak_idx
                min_dist = self.distances[min_idx]
                min_val = self.values[min_idx]
                
                self._peak_info = {
                    'peak_distance': peak_dist,
                    'peak_value': peak_val,
                    'min_distance': min_dist,
                    'min_value': min_val
                }
            else:
                # Fallback values if no peaks found
                self._peak_info = {
                    'peak_distance': None,
                    'peak_value': None,
                    'min_distance': None,
                    'min_value': None
                }
                
        return self._peak_info
    
    def find_peaks_minima(self, prominence: float = 0.1):
        """Find the first peak and minimum in the PDF data.
        
        For cation-anion pairs, finds the first peak and first minimum after the peak.
        For other pair types, only finds the first peak.
        """
        from scipy.signal import find_peaks
        
        # Find peaks in the weighted PDF
        peaks, properties = find_peaks(
            self.values,
            distance=10,
            prominence=prominence
        )
        
        if len(peaks) > 0:
            self.peak_distance = self.distances[peaks[0]]
            self.peak_value = self.values[peaks[0]]
            
            # For cation-anion pairs, find the first minimum after the first peak
            if self.pair_type in ('ca', 'ac') and len(peaks) > 1:
                min_idx = np.argmin(self.values[peaks[0]:peaks[1]]) + peaks[0]
                self.min_distance = self.distances[min_idx]
                self.min_value = self.values[min_idx]

def load_pdf_data(csv_path: str, composition: Optional[Dict[str, float]] = None) -> Dict[str, PDFData]:
    """Load PDF data from a CSV file.
    
    Args:
        csv_path: Path to the CSV file containing PDF data
        composition: Optional dictionary of ion pair weights (e.g., {'Li-F': 0.5, 'Na-F': 0.5})
        
    Returns:
        Dictionary mapping ion pair names to PDFData objects
    """
    # Read the CSV file with two header rows
    df = pd.read_csv(csv_path, header=[0, 1])
    
    pdf_data = {}
    
    # Process each ion pair (every two columns)
    for i in range(0, len(df.columns), 2):
        if i + 1 >= len(df.columns):
            break
            
        # Get ion pair name from the first header row
        ion_pair = df.columns[i][0].strip()
        if not ion_pair:
            continue
        
        # Get weight from composition if provided, otherwise use 1.0
        weight = composition.get(ion_pair, 1.0) if composition else 1.0
        
        # Get the data columns
        dist_col = df.columns[i]
        pdf_col = df.columns[i + 1]
        
        # Get the numeric data, skipping any non-numeric rows
        pair_df = df[[dist_col, pdf_col]].copy()
        
        # Convert to numeric, forcing errors to NaN, then drop any rows with NaN
        pair_df = pair_df.apply(pd.to_numeric, errors='coerce').dropna()
        
        if len(pair_df) == 0:
            print(f"Warning: No valid numeric data found for {ion_pair}")
            continue
        
        # Sort by distance to ensure monotonic increasing x-values
        pair_df = pair_df.sort_values(by=dist_col)
        
        # Create PDFData object with weight
        pdf_data[ion_pair] = PDFData(
            ion_pair=ion_pair,
            distances=pair_df[dist_col].values,
            values=pair_df[pdf_col].values,
            weight=weight
        )
    
    return pdf_data

@dataclass
class SaltComposition:
    """Class to handle salt composition parsing and calculations."""
    formula: str
    components: List[Tuple[float, str]] = field(default_factory=list)
    
    def __post_init__(self):
        """Parse the composition formula."""
        if not self.formula:
            raise ValueError("Composition formula cannot be empty")
        
        # Split into components (e.g., "0.64NaCl-0.36UCl3" -> ["0.64NaCl", "0.36UCl3"])
        parts = self.formula.split('-')
        for part in parts:
            # Split into fraction and formula
            i = 0
            while i < len(part) and (part[i].isdigit() or part[i] in '.-'):
                i += 1
            if i == 0:
                raise ValueError(f"Invalid component format: {part}")
            
            fraction = float(part[:i])
            component_formula = part[i:]
            self.components.append((fraction, component_formula))
    
    def get_ions(self) -> Set[str]:
        """Get all unique ions in the composition."""
        ions = set()
        for _, formula in self.components:
            # Simple ion extraction - can be enhanced for more complex formulas
            if formula.endswith('Cl') or formula.endswith('F') or formula.endswith('Br') or formula.endswith('I'):
                cation = formula[:-2]  # Remove the anion
                ions.add(cation)
                ions.add(formula[-2:])  # Add the anion
            else:
                # Fallback for other formats
                parts = re.findall('[A-Z][^A-Z]*', formula)
                if len(parts) >= 2:
                    ions.update(parts[:2])  # Assume first two parts are cation and anion
        return ions

@dataclass
class TransferFactors:
    """Class to store transfer factors for an ion pair."""
    f_kf: float = 0.0  # Bond strength/diffusion factor
    f_ni: float = 0.0   # Non-ideal factor
    f_ph: np.ndarray = field(default_factory=lambda: np.array([]))  # Phonon transfer factor
    f_dt: np.ndarray = field(default_factory=lambda: np.array([]))  # Diffuson transfer factor
    transfer_points: np.ndarray = field(default_factory=lambda: np.array([]))  # r_tr_i, 2*r_tr_i, etc.

class SCLAnalyzer:
    """Main class for SCL analysis."""
    
    def __init__(self, pdf_data: Dict[str, PDFData], composition: SaltComposition, 
                 temperature: float, experimental_scl: Optional[float] = None,
                 author: str = "Unknown"):
        self.pdf_data = pdf_data
        self.composition = composition
        self.temperature = temperature
        self.experimental_scl = experimental_scl
        self.author = author
        self.transfer_factors: Dict[str, TransferFactors] = {}
        self.results = {}
        
    def analyze(self):
        """Run the complete SCL analysis pipeline."""
        # 1. Find peaks and minima for all ion pairs
        for pdf in self.pdf_data.values():
            pdf.find_peaks_minima()
            
        # 2. Calculate transfer factors
        self._calculate_transfer_factors()
        
        # 3. Calculate SCL for each ion pair
        self._calculate_scl()
        
        # 4. Calculate effective SCL
        self._calculate_effective_scl()
            
        return self.results
        
    def _calculate_transfer_factors(self):
        """Calculate all transfer factors (f_KF, f_NI, f_PH, f_DT).
        
        Follows V3.1's approach for calculating transfer factors while maintaining
        performance through vectorized operations and caching.
        """
        # Pre-calculate interpolation functions for all PDFs
        interp_funcs = {}
        for ion_pair, pdf in self.pdf_data.items():
            if len(pdf.distances) >= 2:
                interp_funcs[ion_pair] = interp1d(
                    pdf.distances,
                    pdf.values,
                    bounds_error=False,
                    fill_value=0.0
                )
        
        # First pass: calculate transfer points and initialize factors
        for ion_pair, pdf in self.pdf_data.items():
            if '-' not in ion_pair:
                continue
                
            tf = TransferFactors()
            
            # Calculate f_KF (Kubo factor) - peak to minimum ratio
            if pdf.peak_value is not None and pdf.min_value is not None:
                tf.f_kf = abs(pdf.peak_value - pdf.min_value) / pdf.peak_value
            
            # Set transfer points as multiples of peak distance
            if pdf.peak_distance is not None:
                # Use up to 5 coordination shells
                tf.transfer_points = np.linspace(
                    pdf.peak_distance * 0.8,  # Slightly before first peak
                    pdf.peak_distance * 5.0,  # Up to 5th coordination shell
                    num=50,  # Increased resolution for better integration
                    endpoint=True
                )
                
                # Initialize factor arrays
                tf.f_ni = np.zeros_like(tf.transfer_points)
                tf.f_ph = np.ones_like(tf.transfer_points)  # Will be updated later
                tf.f_dt = np.ones_like(tf.transfer_points)  # Will be updated later
                
                self.transfer_factors[ion_pair] = tf
        
        # Second pass: calculate f_NI (non-ideal factor) for each ion pair
        for ion_pair, tf in self.transfer_factors.items():
            if ion_pair not in self.pdf_data:
                continue
                
            pdf = self.pdf_data[ion_pair]
            cation, anion = ion_pair.split('-')
            
            # Get interpolation function for this PDF
            if ion_pair not in interp_funcs:
                continue
                
            interp_func = interp_funcs[ion_pair]
            
            # Calculate total PDF for this anion at each transfer point
            total_anion_pdf = np.zeros_like(tf.transfer_points)
            for other_pair, other_interp in interp_funcs.items():
                if other_pair.endswith(anion):  # All pairs sharing the same anion
                    other_pdf = self.pdf_data[other_pair]
                    total_anion_pdf += other_interp(tf.transfer_points) * other_pdf.weight
            
            # Calculate f_NI as ratio of this PDF to total anion PDF
            pair_pdf = interp_func(tf.transfer_points) * pdf.weight
            with np.errstate(divide='ignore', invalid='ignore'):
                tf.f_ni = np.where(
                    total_anion_pdf > 0,
                    pair_pdf / total_anion_pdf,
                    0.0
                )
            
            # Calculate phonon and diffuson factors
            self._calculate_phonon_diffuson_factors(ion_pair, tf)
            
            # Ensure all factors are within [0, 1]
            tf.f_ni = np.clip(tf.f_ni, 0.0, 1.0)
            tf.f_ph = np.clip(tf.f_ph, 0.0, 1.0)
            tf.f_dt = np.clip(tf.f_dt, 0.0, 1.0)

    def _calculate_non_ideal_factor(self, ion_pair: str, r: float) -> float:
        """Calculate the non-ideal factor f_NI at a specific distance r.
        
        Args:
            ion_pair: The ion pair identifier (e.g., 'Li-F')
            r: Distance at which to calculate the factor
            
        Returns:
            Non-ideal factor value at distance r
            
        Note: This method is kept for backward compatibility but is no longer used in V3.1.
        The non-ideal factor is now calculated directly in _calculate_transfer_factors.
        """
        return 1.0  # Default value for backward compatibility
        print(f"Warning: _calculate_non_ideal_factor called for {ion_pair} at r={r}")
        return 0.0

    def _calculate_phonon_diffuson_factors(self, ion_pair: str, tf: TransferFactors):
        """Calculate phonon and diffuson transfer factors using vectorized operations.
        
        This method implements the V3.1 approach for calculating:
        - f_PH: Phonon transfer factor based on coordination environment
        - f_DT: Diffuson transfer factor (set to 1.0 in V3.1)
        
        Args:
            ion_pair: The ion pair identifier (e.g., 'Li-F')
            tf: TransferFactors object to store results
        """
        # In V3.1, f_DT is set to 1.0 and not used in calculations
        tf.f_dt = np.ones_like(tf.transfer_points)
        
        # Get the current PDF data and check if it has valid peak info
        current_pdf = self.pdf_data[ion_pair]
        if not hasattr(current_pdf, 'peak_distance') or current_pdf.peak_distance is None:
            tf.f_ph = np.ones_like(tf.transfer_points)
            return
            
        # Pre-allocate arrays for vectorized operations
        weighted_sum = np.zeros_like(tf.transfer_points)
        weight_sum = 0.0
        
        # Calculate weights and contributions from all PDFs
        for pair, pdf in self.pdf_data.items():
            if pair == ion_pair:
                continue
                
            # Skip PDFs without valid peak information
            if not hasattr(pdf, 'peak_distance') or pdf.peak_distance is None:
                continue
                
            # Get weight from composition (default to 1.0 if not found)
            weight = 1.0
            for comp_weight, comp_formula in self.composition.components:
                if pair in comp_formula:
                    weight = comp_weight
                    break
            
            # Calculate coordination number ratio (r_ij / r_ii)
            r_ij = pdf.peak_distance
            r_ii = current_pdf.peak_distance
            coordination_ratio = r_ij / r_ii if r_ii > 0 else 1.0
            
            # Calculate effective weight considering coordination
            effective_weight = weight * coordination_ratio
            
            # Interpolate PDF values at transfer points
            interp_vals = np.interp(
                tf.transfer_points,
                pdf.distances,
                pdf.values,
                left=0.0,
                right=0.0
            )
            
            # Add to weighted sum
            weighted_sum += interp_vals * effective_weight
            weight_sum += effective_weight
        
        # Calculate f_PH: 1 - (weighted average PDF) / max(weighted average PDF)
        if weight_sum > 0 and np.max(weighted_sum) > 0:
            tf.f_ph = 1.0 - (weighted_sum / weight_sum) / np.max(weighted_sum / weight_sum)
        elif weight_sum == 0:
            tf.f_ph = np.ones_like(tf.transfer_points)
            
            # Apply any distance-dependent corrections if needed
            # (This is where V3.1 would apply any additional corrections)
            # For now, we keep it simple and just use 1.0
            pass

    def _calculate_scl(self):
        """Calculate SCL for each ion pair using V3.1's survival function approach.
        
        Implements the formula: SCL = ∫₀^∞ S(r) dr, where S(r) = exp(-∫₀^r β(r') dr')
        and β(r) = b(r) / (1 - b(r)), with b(r) = f_KF * f_NI(r) * f_PH(r).
        
        This implementation uses vectorized operations for better performance.
        """
        for ion_pair, tf in self.transfer_factors.items():
            if ion_pair not in self.pdf_data:
                continue
                
            # Get the PDF data for this ion pair
            pdf = self.pdf_data[ion_pair]
            
            # Calculate the integrand β(r) = b(r) / (1 - b(r))
            # Use np.clip to avoid division by zero and numerical instability
            b_r = np.clip(tf.f_kf * tf.f_ni * tf.f_ph, 1e-10, 1.0 - 1e-10)
            beta_r = b_r / (1.0 - b_r)
            
            # Calculate the cumulative integral of β(r) from 0 to r using trapezoidal rule
            # This is more accurate than the previous implementation
            cum_integral = np.zeros_like(tf.transfer_points)
            if len(tf.transfer_points) > 1:
                # Use cumulative trapezoidal integration
                cum_integral[1:] = np.cumulative_trapezoid(
                    beta_r, 
                    tf.transfer_points,
                    initial=0.0
                )
            
            # Calculate the survival function S(r) = exp(-∫₀^r β(r') dr')
            # Use np.exp(-x) for numerical stability
            survival = np.exp(-cum_integral)
            
            # Calculate SCL as the integral of the survival function
            # Use Simpson's rule for better accuracy with fewer points
            if len(tf.transfer_points) > 2:
                scl = simps(survival, tf.transfer_points)
            elif len(tf.transfer_points) > 1:
                scl = np.trapz(survival, tf.transfer_points)
            else:
                scl = 0.0
            
            # Store results with additional metadata
            self.results[ion_pair] = {
                'scl': scl,
                'survival': survival,
                'transfer_points': tf.transfer_points,
                'f_kf': tf.f_kf,
                'f_ni': tf.f_ni,
                'f_ph': tf.f_ph,
                'f_dt': tf.f_dt,
                'peak_distance': pdf.peak_distance,
                'min_distance': pdf.min_distance,
                'pair_type': pdf.pair_type
            }

    def _calculate_effective_scl(self):
        """Calculate the effective total SCL using composition and coordination weights.
        
        The effective SCL is calculated as a weighted average of individual ion pair SCLs,
        where weights consider both composition and coordination effects.
        
        The formula used is:
        SCL_eff = (∑_i w_i * c_i * SCL_i) / (∑_i w_i * c_i)
        where:
        - w_i is the composition weight
        - c_i is a coordination-based correction factor
        - SCL_i is the SCL for the i-th ion pair
        
        The method prioritizes cation-anion pairs but falls back to other pairs if needed.
        """
        if not self.results:
            self.results['effective_scl'] = 0.0
            return
            
        weighted_scl = 0.0
        total_weight = 0.0
        ion_pairs_used = []
        
        # First, collect all valid ion pairs with their weights and SCL values
        valid_pairs = []
        for ion_pair, result in self.results.items():
            if not isinstance(result, dict) or ion_pair == 'effective' or 'scl' not in result:
                continue
                
            if ion_pair not in self.pdf_data:
                continue
                
            pdf = self.pdf_data[ion_pair]
            if not hasattr(pdf, 'weight') or not hasattr(pdf, 'pair_type'):
                continue
                
            valid_pairs.append({
                'ion_pair': ion_pair,
                'scl': result['scl'],
                'weight': pdf.weight,
                'pair_type': pdf.pair_type,
                'peak_distance': result.get('peak_distance', 0)
            })
        
        if not valid_pairs:
            self.results['effective_scl'] = 0.0
            return
        
        # Separate cation-anion pairs from others
        ca_pairs = [p for p in valid_pairs if p['pair_type'] in ('ca', 'ac')]
        other_pairs = [p for p in valid_pairs if p['pair_type'] not in ('ca', 'ac')]
        
        # Use cation-anion pairs if available, otherwise fall back to other pairs
        pairs_to_use = ca_pairs if ca_pairs else other_pairs
        
        # Calculate total weight for normalization
        total_weight_sum = sum(p['weight'] for p in pairs_to_use)
        
        # Calculate weighted SCL
        for pair in pairs_to_use:
            if total_weight_sum > 0:
                weight = pair['weight'] / total_weight_sum
            else:
                weight = 1.0 / len(pairs_to_use)  # Equal weights if no weight info
            
            # Apply coordination-based correction factor
            # This could be enhanced with actual coordination number data
            coordination_factor = 1.0  # Default to no correction
            
            # Calculate weighted contribution
            weighted_contribution = pair['scl'] * weight * coordination_factor
            
            weighted_scl += weighted_contribution
            total_weight += weight * coordination_factor
            
            ion_pairs_used.append({
                'ion_pair': pair['ion_pair'],
                'scl': pair['scl'],
                'weight': weight,
                'pair_type': pair['pair_type'],
                'contribution': weighted_contribution
            })
        
        # Calculate effective SCL if we have valid pairs
        if total_weight > 0 and ion_pairs_used:
            effective_scl = weighted_scl / total_weight
            
            # Calculate the weighted variance for uncertainty estimation
            weighted_variance = 0.0
            for pair in ion_pairs_used:
                deviation = pair['scl'] - effective_scl
                weighted_variance += pair['weight'] * (deviation ** 2)
            
            # Normalize by total weight to get weighted variance
            if total_weight > 0:
                weighted_variance /= total_weight
            
            # Store the results
            self.results['effective'] = {
                'scl': effective_scl,
                'std_dev': np.sqrt(weighted_variance) if weighted_variance > 0 else 0.0,
                'ion_pairs_used': ion_pairs_used,
                'num_pairs': len(ion_pairs_used),
                'pair_types': list(set(p['pair_type'] for p in ion_pairs_used)),
                'timestamp': pd.Timestamp.now().isoformat()
            }
            
            # Also store as a top-level value for backward compatibility
            self.results['effective_scl'] = effective_scl
        else:
            # If no valid pairs were found
            self.results['effective'] = {
                'scl': 0.0,
                'std_dev': 0.0,
                'ion_pairs_used': [],
                'num_pairs': 0,
                'warning': 'No valid ion pairs found for effective SCL calculation',
                'timestamp': pd.Timestamp.now().isoformat()
            }
            self.results['effective_scl'] = 0.0
            variance = sum(
                weight * (scl - effective_scl) ** 2 
                for _, weight, scl in ion_pairs_used
            ) / total_weight if len(ion_pairs_used) > 1 else 0.0
            std_dev = np.sqrt(variance)
            
            # Store detailed results for debugging and analysis
            self.results['effective'] = {
                'scl': effective_scl,
                'std_dev': std_dev,
                'weight': total_weight,
                'ion_pairs': [
                    {
                        'ion_pair': pair,
                        'weight': float(weight),
                        'scl': float(scl),
                        'contribution': float(weight * scl / total_weight if total_weight > 0 else 0)
                    }
                    for pair, weight, scl in ion_pairs_used
                ],
                'method': 'weighted_average_composition',
                'temperature': self.temperature,
                'experimental_scl': self.experimental_scl,
                'num_ion_pairs': len(ion_pairs_used),
                'total_coordination': total_coordination,
                'total_ca_weight': total_ca_weight
            }
            
            # Print summary for debugging
            print(f"\nEffective SCL calculation (V3.1 method):")
            print(f"  Temperature: {self.temperature} K")
            if self.experimental_scl is not None:
                print(f"  Experimental SCL: {self.experimental_scl:.3f} Å")
            print(f"  Calculated SCL: {effective_scl:.3f} ± {std_dev:.3f} Å")
            print(f"  Number of ion pairs: {len(ion_pairs_used)}")
            print("  Ion pair contributions:")
            for pair, weight, scl in sorted(ion_pairs_used, key=lambda x: x[1], reverse=True):
                contribution = weight * scl / total_weight if total_weight > 0 else 0
                print(f"    {pair}: {scl:.3f} Å (weight: {weight:.3f}, contribution: {contribution:.3f} Å)")
            
            return effective_scl
        
        # If we couldn't calculate an effective SCL, return None
        # and remove any existing 'effective' entry to avoid confusion
        self.results.pop('effective', None)
        return None
        
    def plot_results(self, save_path: str = None, show_plot: bool = True):
        """Plot the PDFs with survival probabilities overlaid and mark the SCL values.
        
        Args:
            save_path: If provided, save the plot to this path
            show_plot: If True, display the plot
        """
        if not self.results or 'effective' not in self.results:
            print("No results to plot. Run analyze() first.")
            return
        
        # Calculate total weight of all cation-anion pairs for normalization
        total_ca_weight = sum(
            pdf.weight 
            for pair, pdf in self.pdf_data.items()
            if pdf.pair_type in ('ca', 'ac')
        )
            
        # Create figure with a single plot
        plt.figure(figsize=(12, 8))
        ax1 = plt.gca()
        
        # Plot all PDFs, weighted by composition
        for ion_pair, pdf in self.pdf_data.items():
            # Skip if no data or not a cation-anion pair
            if len(pdf.distances) == 0 or len(pdf.values) == 0 or pdf.pair_type not in ('ca', 'ac'):
                continue
                
            # Calculate weight for this ion pair
            weight = pdf.weight / total_ca_weight if total_ca_weight > 0 else 1.0
            
            # Plot the weighted PDF
            line, = ax1.plot(pdf.distances, pdf.values * weight,
                           linestyle='-',
                           alpha=0.7,
                           label=f"{ion_pair} PDF (x{weight:.2f})")
            
            # Mark peaks and minima with vertical lines
            if pdf.peak_distance is not None:
                ax1.axvline(pdf.peak_distance, color=line.get_color(), 
                          linestyle='--', alpha=0.5)
                if pdf.min_distance is not None:
                    ax1.axvline(pdf.min_distance, color=line.get_color(), 
                              linestyle=':', alpha=0.5)
        
        # Create second y-axis for survival probability
        ax2 = ax1.twinx()
        
        # Plot S(r) for each analyzed ion pair (only cation-anion pairs)
        for ion_pair, result in self.results.items():
            # Skip non-cation-anion pairs and the 'effective' entry
            if ion_pair == 'effective' or ion_pair not in self.pdf_data:
                continue
                
            if self.pdf_data[ion_pair].pair_type not in ('ca', 'ac'):
                continue
                
            # Check if we have valid data to plot
            if 'r_points' not in result or 'S_r' not in result:
                print(f"Warning: Missing S(r) data for {ion_pair} - skipping plot")
                continue
                
            # Get the weight for this ion pair
            weight = self.pdf_data[ion_pair].weight / total_ca_weight if total_ca_weight > 0 else 1.0
                
            # Get the S(r) data
            r_points = np.array(result['r_points'])
            S_r = np.array(result['S_r'])
            
            # Get color from PDF plot (matching the PDF line color)
            color = None
            for i, pair in enumerate(self.pdf_data.keys()):
                if pair == ion_pair:
                    color = ax1.lines[i].get_color()
                    break
            
            if color is None:
                color = f'C{len(ax2.lines)}'  # Fallback color
            
            # Plot the survival function S(r)
            ax2.plot(r_points, S_r, '-', color=color, linewidth=2,
                    alpha=0.8, label=f'{ion_pair} S(r)')
            
            # Mark the SCL point with a vertical line and label
            if 'scl' in result:
                scl = result['scl']
                # Add a vertical line at the SCL value
                ax2.axvline(scl, color=color, linestyle='--', alpha=0.6, linewidth=1.5)
                # Add a text label for the SCL value
                ax2.text(scl, 0.05, f'SCL={scl:.2f} Å', 
                        color=color, ha='left', va='bottom',
                        rotation=90, transform=ax2.get_xaxis_transform(),
                        bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1))
        
        # Add experimental SCL if available
        if self.experimental_scl is not None:
            exp_scl = self.experimental_scl
            ax2.axvline(exp_scl, color='black', linestyle='-', linewidth=2, alpha=0.8,
                       label=f'Exp. SCL: {exp_scl:.2f} Å')
        
        # Add 1/e line
        ax2.axhline(np.exp(-1), color='gray', linestyle=':', alpha=0.7, label='1/e')
        
        # Set plot titles and labels
        ax1.set_xlabel('Distance (Å)', fontsize=12)
        ax1.set_ylabel('Weighted PDF (arb. units)', fontsize=12)
        ax2.set_ylabel('Survival Probability S(r)', fontsize=12)
        
        # Set title with temperature and effective SCL
        title_parts = [f"SCL Analysis - {self.composition.formula}"]
        if 'effective' in self.results:
            scl = self.results['effective']['scl']
            std = self.results['effective'].get('std_dev', 0)
            title_parts.append(f"Effective SCL: {scl:.2f} ± {std:.2f} Å")
        if self.temperature is not None:
            title_parts.append(f"T = {self.temperature} K")
        
        plt.title(" | ".join(title_parts), fontsize=14, pad=20)
        
        # Add legends
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1, labels1, loc='upper right')
        ax2.legend(lines2, labels2, loc='upper left')
        
        # Adjust layout to prevent label cutoff
        plt.tight_layout()
        
        # Save the plot if requested
        if save_path:
            os.makedirs(os.path.dirname(os.path.abspath(save_path)), exist_ok=True)
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Plot saved to {save_path}")
        
        # Show the plot if requested
        if show_plot:
            plt.show()
        else:
            plt.close()
        
        return ax1, ax2

    def save_results(self, filepath: Union[str, Path]):
        """Save the analysis results to a CSV file with enhanced metadata and formatting.
        
        Args:
            filepath: Path to save the results CSV (str or Path)
            
        Returns:
            Path: The absolute path to the saved file
        """
        if not self.results:
            print("No results to save. Run analyze() first.")
            return None
            
        # Ensure filepath is a Path object
        output_path = Path(filepath).resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Prepare data for CSV export
        data = []
        
        # Add detailed results for each ion pair
        for ion_pair, result in self.results.items():
            if ion_pair == 'effective' or not isinstance(result, dict):
                continue
                
            row = {
                'ion_pair': ion_pair,
                'scl': result.get('scl', np.nan),
                'peak_distance': result.get('peak_distance', np.nan),
                'min_distance': result.get('min_distance', np.nan),
                'f_kf': result.get('f_kf', np.nan),
                'pair_type': result.get('pair_type', 'unknown'),
                'weight': self.pdf_data[ion_pair].weight if ion_pair in self.pdf_data else 1.0
            }
            
            # Add transfer factors if available
            for factor in ['f_ni', 'f_ph', 'f_dt']:
                if factor in result:
                    if isinstance(result[factor], (list, np.ndarray)) and len(result[factor]) > 0:
                        row[f'{factor}_mean'] = np.mean(result[factor])
                        row[f'{factor}_std'] = np.std(result[factor]) if len(result[factor]) > 1 else 0.0
                    else:
                        row[factor] = result[factor]
            
            data.append(row)
        
        # Create DataFrame and sort by SCL (descending)
        df = pd.DataFrame(data)
        if not df.empty and 'scl' in df.columns:
            df = df.sort_values('scl', ascending=False)
        
        # Add metadata as a comment at the top of the file
        metadata = [
            f"# SCL Analysis Results - {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}",
            f"# Composition: {getattr(self.composition, 'formula', 'N/A')}",
            f"# Temperature: {getattr(self, 'temperature', 'N/A')} K",
            f"# Author: {getattr(self, 'author', 'N/A')}",
            f"# Effective SCL: {self.results.get('effective_scl', 'N/A')} Å"
        ]
        
        # Write metadata and data to file
        with open(output_path, 'w', newline='') as f:
            # Write metadata as comments
            f.write('\n'.join(metadata) + '\n')
            # Write data
            df.to_csv(f, index=False, float_format='%.6f')
        
        print(f"Results saved to {output_path}")
        return output_path

    def run_analysis(self, plot_path: Optional[Union[str, Path]] = None, 
                    results_path: Optional[Union[str, Path]] = None, 
                    show_plot: bool = True,
                    progress_callback: Optional[callable] = None) -> dict:
        """Run the complete analysis pipeline with progress tracking and error handling.
        
        Args:
            plot_path: If provided, save the plot to this path (str or Path)
            results_path: If provided, save results to this CSV path (str or Path)
            show_plot: If True, display the plot
            progress_callback: Optional callback function for progress updates.
                             Signature: callback(progress: float, message: str)
            
        Returns:
            dict: Analysis results with the following structure:
                {
                    'ion_pairs': {
                        'pair1': {scl: float, peak_distance: float, ...},
                        'pair2': {...},
                        ...
                    },
                    'effective_scl': float,
                    'metadata': {
                        'temperature': float,
                        'composition': str,
                        'timestamp': str,
                        'version': str
                    },
                    'warnings': List[str]  # Any warnings generated during analysis
                }
                
        Raises:
            ValueError: If required data is missing or invalid
            RuntimeError: If analysis fails
        """
        import time
        start_time = time.time()
        warnings = []
        
        def update_progress(progress: float, message: str):
            if progress_callback and callable(progress_callback):
                try:
                    progress_callback(progress, message)
                except Exception as e:
                    warnings.append(f"Progress callback failed: {str(e)}")
        
        try:
            update_progress(0.1, "Starting analysis...")
            
            # Validate inputs
            if not hasattr(self, 'pdf_data') or not self.pdf_data:
                raise ValueError("No PDF data provided. Load data before running analysis.")
                
            if not hasattr(self, 'composition') or not self.composition:
                raise ValueError("No composition provided. Set composition before running analysis.")
            
            # Run the analysis steps with progress updates
            update_progress(0.2, "Finding peaks and minima...")
            for pdf in self.pdf_data.values():
                pdf.find_peaks_minima()
            
            update_progress(0.4, "Calculating transfer factors...")
            self._calculate_transfer_factors()
            
            update_progress(0.6, "Calculating SCL values...")
            self._calculate_scl()
            
            update_progress(0.8, "Calculating effective SCL...")
            self._calculate_effective_scl()
            
            # Save results if path is provided
            saved_path = None
            if results_path:
                update_progress(0.9, f"Saving results to {results_path}...")
                saved_path = self.save_results(results_path)
            
            # Generate and optionally save/show the plot
            if plot_path or show_plot:
                update_progress(0.95, "Generating plots...")
                self.plot_results(save_path=plot_path, show_plot=show_plot)
            
            # Add execution time to results
            exec_time = time.time() - start_time
            if '_metadata' not in self.results:
                self.results['_metadata'] = {}
            self.results['_metadata']['execution_time_seconds'] = exec_time
            
            update_progress(1.0, f"Analysis completed in {exec_time:.2f} seconds")
            
            # Add any warnings to results
            if warnings:
                if '_warnings' not in self.results:
                    self.results['_warnings'] = []
                self.results['_warnings'].extend(warnings)
            
            return self.results
            
        except Exception as e:
            error_msg = f"Analysis failed: {str(e)}"
            if progress_callback:
                progress_callback(-1, error_msg)
            raise RuntimeError(error_msg) from e

if __name__ == "__main__":
    # Example usage
    salt_name = "LiCl"  # Change this to match your CSV filename
    pdf_file = os.path.join(DATA_DIR, "PDF_LiCl.csv")
    composition = '1.0LiCl'  # Update this based on your salt
    temperature = 878
    experimental_value = 3.893446
    author = "Walz, 2019"
    
    # Generate output filenames with timestamp
    from datetime import datetime
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    plot_path = os.path.join(RESULTS_DIR, f"analysis_plot_{salt_name}_{timestamp}.png")
    results_path = os.path.join(RESULTS_DIR, f"analysis_results_{salt_name}_{timestamp}.csv")
    
    try:
        # Check if file exists
        if not os.path.exists(pdf_file):
            available_files = [f for f in os.listdir(DATA_DIR) if f.startswith('PDF_') and f.endswith('.csv')]
            raise FileNotFoundError(
                f"PDF file not found: {pdf_file}\n"
                f"Available PDF files in {DATA_DIR}:\n"
                + "\n".join(available_files) if available else "No PDF files found"
            )
        
        # Load PDF data
        print(f"Loading data from: {pdf_file}")
        pdf_data = load_pdf_data(pdf_file)
        
        # Create analyzer and run analysis
        analyzer = SCLAnalyzer(
            pdf_data=pdf_data,
            composition=SaltComposition(composition),
            temperature=temperature,
            experimental_scl=experimental_value,
            author=author
        )
        
        results = analyzer.run_analysis(
            plot_path=plot_path,
            results_path=results_path,
            show_plot=True
        )
        
        print(f"Analysis complete! Results saved to:")
        print(f"- Plot: {os.path.abspath(plot_path)}")
        print(f"- Results: {os.path.abspath(results_path)}")
        
    except Exception as e:
        print(f"Error during analysis: {str(e)}")
        raise