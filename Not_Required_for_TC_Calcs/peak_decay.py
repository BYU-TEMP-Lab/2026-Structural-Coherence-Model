import numpy as np
import pandas as pd
import glob, os
import re
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error

def get_y_at_x(x_values, y_values, target_x, tolerance=0.01):
    """Get the y-value at the target_x position using linear interpolation."""
    # Find the closest points for interpolation
    idx = np.searchsorted(x_values, target_x)
    if idx == 0:
        return y_values[0]
    if idx == len(x_values):
        return y_values[-1]
    
    x0, x1 = x_values[idx-1], x_values[idx]
    y0, y1 = y_values[idx-1], y_values[idx]
    
    # Linear interpolation
    if x1 - x0 < tolerance:
        return (y0 + y1) / 2
    
    return y0 + (target_x - x0) * (y1 - y0) / (x1 - x0)

# Manual peak positions for specific salts (in Angstroms)
MANUAL_PEAKS = {
    'CaCl2': [2.78, 5.66320254, 8.928837447, 11.7564512],    # Example values, replace with your actual peaks
    'KF': [2.481944469, 5.743995554, 9.226702399, 12.62528589],       # Example values, replace with your actual peaks
    'MgCl2': [2.442330743, 5.626370077, 8.587589082, 10.71217683],   # Example values, replace with your actual peaks
    'KCl': [2.985223581, 6.656822612, 10.84073779, 14.93072834],       # Example values, replace with your actual peaks
    'SrCl2': [2.9, 6.037339687, 10.204821, 13.87426139],     # Example values, replace with your actual peaks
    'RbF': [2.598862512, 6.030316518, 9.781366731, 13.33056671]       # Example values, replace with your actual peaks
}

# Predicted peak positions for other salt types (in Angstroms)
PREDICTED_PEAKS = {
    'NaCl': [2.8, 4.0, 5.5, 6.5],
    'LiF': [1.8, 4.5, 6.8, 9.5],
    'NaF': [2.3, 5.0, 7.9, 11.0],
    'LiCl': [2.6, 5.4, 8.8, 12.1],
    'CsF': [3.0, 6.5, 10.5, 14.2],
    'ZnCl2': [2.3, 5.5, 7.8, 11.0]
}


# --- Power law function ---
def power_law(n, a, b):
    return a * n**b

# --- Evaluate fit quality ---
def evaluate_fit(y_true, y_pred, n_params=2):
    """
    Evaluate model fit with multiple metrics.
    
    Args:
        y_true: Array of true values
        y_pred: Array of predicted values
        n_params: Number of parameters in the model (for adjusted R²)
        
    Returns:
        dict: Dictionary containing various fit metrics
    """
    n = len(y_true)
    if n <= n_params + 1:
        adj_r2 = float('nan')  # Not enough data points
    else:
        r2 = r2_score(y_true, y_pred)
        adj_r2 = 1 - (1 - r2) * (n - 1) / (n - n_params - 1)
    
    metrics = {
        'r2': r2_score(y_true, y_pred),
        'adj_r2': adj_r2,
        'rmse': np.sqrt(mean_squared_error(y_true, y_pred)),
        'mae': mean_absolute_error(y_true, y_pred),
        'max_abs_error': np.max(np.abs(y_true - y_pred)),
        'mean_rel_error': np.mean(np.abs((y_true - y_pred) / y_true)) * 100,
        'median_rel_error': np.median(np.abs((y_true - y_pred) / y_true)) * 100
    }
    
    return metrics

# --- Process one CSV file ---
def get_expected_peaks(salt_name, n_peaks=4):
    """Get expected peak positions based on salt type."""
    # First try to find in MANUAL_PEAKS
    for key in MANUAL_PEAKS:
        if key.lower() in salt_name.lower():
            print(f"{salt_name}: Using MANUAL peak positions")
            return np.array(MANUAL_PEAKS[key][:n_peaks])
    
    # If not in MANUAL_PEAKS, try PREDICTED_PEAKS
    for key in PREDICTED_PEAKS:
        if key.lower() in salt_name.lower():
            print(f"{salt_name}: Using predicted peak positions")
            return np.array(PREDICTED_PEAKS[key][:n_peaks])
    
    print(f"{salt_name}: No manual or predicted peaks found, using automatic detection")
    return None

def find_nearest_peaks(r, g, expected_positions, window=0.5):
    """Find actual peaks near expected positions."""
    peak_positions = []
    peak_heights = []
    
    for pos in expected_positions:
        # Find indices within window of expected position
        mask = (r >= pos - window) & (r <= pos + window)
        if not np.any(mask):
            continue
            
        # Find maximum in this window
        local_max_idx = np.argmax(g[mask])
        global_idx = np.where(mask)[0][local_max_idx]
        
        peak_positions.append(r[global_idx])
        peak_heights.append(g[global_idx])
    
    return np.array(peak_positions), np.array(peak_heights)

def process_pdf(csv_path, n_peaks=4, plot=False, plot_dir="plots"):
    with open(csv_path, "r") as f:
        salt_name = f.readline().strip()   # 1st line is salt name

    df = pd.read_csv(csv_path, skiprows=1, sep=None, engine="python")
    df.columns = df.columns.str.strip()
    
    # Find the first two numeric columns for X and Y data
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    x_col = numeric_cols[0]
    y_col = numeric_cols[1]
    
    r = df[x_col].values
    g = df[y_col].values
    
    # Check if we have manual peaks for this salt
    manual_peaks = None
    # Clean up the salt name by removing any non-alphanumeric characters
    clean_salt_name = re.sub(r'[^a-zA-Z0-9]', '', salt_name.lower())
    
    # First, try to match the full salt name (e.g., 'SrCl2' in 'SrCl2' or 'Sr-Cl')
    for key in MANUAL_PEAKS:
        # Check if the key is in the salt name (e.g., 'srcl2' in 'srcl')
        key_clean = key.lower().replace('2', '').replace('3', '')
        if key_clean in clean_salt_name or clean_salt_name in key_clean:
            manual_peaks = MANUAL_PEAKS[key][:n_peaks]
            print(f"{salt_name}: Found exact manual peaks for {key}")
            break
    
    # If no exact match, try matching by element parts (e.g., 'sr' and 'cl' in 'srcl')
    if manual_peaks is None:
        for key in MANUAL_PEAKS:
            # Split key into elements (e.g., 'SrCl2' -> ['sr', 'cl'])
            key_parts = re.findall('[A-Z][^A-Z]*', key)
            key_parts = [p.lower() for p in key_parts]
            
            # Check if all element parts are in the clean salt name
            if all(part.lower() in clean_salt_name for part in key_parts):
                manual_peaks = MANUAL_PEAKS[key][:n_peaks]
                print(f"{salt_name}: Found manual peaks by element match for {key}")
                break
    
    if manual_peaks is not None:
        print(f"{salt_name}: Using MANUAL peak positions")
        peak_positions = np.array(manual_peaks)
        peak_heights = np.array([get_y_at_x(r, g, x) for x in peak_positions])
        print(f"{salt_name}: Using manual peaks at {peak_positions}")
    else:
        print(f"{salt_name}: Using automatic peak detection")
        
        # Minimum peak distance (0.1 Å in data points) - very small to catch close peaks
        min_peak_distance = max(1, int(0.1 / (r[1] - r[0])))
        
        # Find peaks with extremely lenient parameters
        peaks, properties = find_peaks(g, 
                                     distance=min_peak_distance,
                                     height=np.percentile(g, 30),  # Very low height threshold (30th percentile)
                                     prominence=0.001,  # Extremely small prominence
                                     width=0,  # No minimum width requirement
                                     rel_height=0.5,
                                     wlen=None)  # No window length restriction
        
        if len(peaks) > 0:
            # Get peak positions and heights
            peak_positions = r[peaks]
            peak_heights = g[peaks]
            
            # Sort peaks by height (descending) and take top n_peaks
            sorted_indices = np.argsort(peak_heights)[::-1]
            peak_positions = peak_positions[sorted_indices][:n_peaks]
            peak_heights = peak_heights[sorted_indices][:n_peaks]
            
            # Sort peaks by position (ascending)
            sorted_pos_indices = np.argsort(peak_positions)
            peak_positions = peak_positions[sorted_pos_indices]
            peak_heights = peak_heights[sorted_pos_indices]
            
            print(f"{salt_name}: Found peaks at {peak_positions}")
            
            if len(peak_positions) < n_peaks:
                print(f"{salt_name}: Warning - Only found {len(peak_positions)} peaks out of requested {n_peaks}")
                return {
                    'file': os.path.basename(csv_path),
                    'salt': salt_name,
                    'status': f'Skipped (only {len(peak_positions)} peaks found)'
                }

    # Get peak heights from the PDF
    peak_heights = np.array([g[np.argmin(abs(r-p))] for p in peak_positions])
    
    # Define power law with baseline 1: y = a * x^(-b) + 1
    def power_law_with_baseline(x, a, b):
        return a * (x ** (-b)) + 1
    
    # Fit peak heights as a function of position
    try:
        popt, _ = curve_fit(power_law_with_baseline, 
                           peak_positions,  # x = peak positions
                           peak_heights,    # y = peak heights
                           p0=[max(peak_heights)-1, 1],  # Initial guess
                           bounds=([0, 0], [np.inf, np.inf]))
        a, b = popt
        
        # Generate smooth x values for plotting
        x_smooth = np.linspace(min(peak_positions), max(peak_positions), 100)
        
        # Predicted values from power law
        fit_vals = power_law_with_baseline(peak_positions, a, b)
        
        # Evaluate fit with enhanced metrics
        metrics = evaluate_fit(peak_heights, fit_vals, n_params=2)
        
    except Exception as e:
        print(f"{salt_name}: Error fitting power law - {str(e)}")
        return {
            'file': os.path.basename(csv_path),
            'salt': salt_name,
            'status': f'Error fitting power law: {str(e)}'
        }

    # --- Optional plotting ---
    if plot:
        os.makedirs(plot_dir, exist_ok=True)

        plt.figure(figsize=(6,4))
        plt.plot(r, g, label="PDF g(r)")
        # plt.scatter(peak_positions, [g[np.argmin(abs(r-p))] for p in peak_positions],
        #             color='red', zorder=5, label="Detected Peaks")

        # Generate smooth x values for the power law curve
        x_smooth = np.linspace(min(peak_positions), max(peak_positions), 100)
        
        # Generate the y-values using the height power law
        y_smooth = power_law_with_baseline(x_smooth, a, b)
        
        # Plot the actual peak points
        plt.plot(peak_positions, peak_heights, 'go', 
                markersize=6, alpha=0.7, 
                label='Peak Heights')
        
        # Plot the smooth power law curve
        plt.plot(x_smooth, y_smooth, 'g--', 
                label=f'Power Law Fit:\ny = {a:.2f}·x^(-{b:.2f}) + 1', 
                alpha=0.7, 
                linewidth=2)
        
        # # Plot vertical lines at the fitted peak positions
        # for n, fv in zip(n_vals, fit_vals):
        #     plt.axvline(x=fv, color='green', linestyle='--', alpha=0.5, linewidth=1)
            
        plt.xlabel('r (Å)')
        plt.ylabel('g(r)')
        plt.title(f'{os.path.basename(csv_path)}\nPower Law Fit: y = {a:.2f}·x^(-{b:.2f}) + 1')
        
        # Add text box with fit parameters and metrics
        textstr = '\n'.join((
            f'Fit Parameters:',
            f'a = {a:.4f}',
            f'b = {b:.4f}',
            '\nFit Quality:',
            f'R² = {metrics["r2"]:.4f}',
            f'Adj R² = {metrics["adj_r2"]:.4f}',
            f'RMSE = {metrics["rmse"]:.4f}',
            f'MAE = {metrics["mae"]:.4f}',
            f'Mean Rel Err = {metrics["mean_rel_error"]:.2f}%',
            f'Max Abs Err = {metrics["max_abs_error"]:.4f}'
        ))
        props = dict(boxstyle='round', facecolor='white', alpha=0.7, pad=0.5)
        plt.gca().text(0.02, 0.98, textstr, transform=plt.gca().transAxes,
                      fontsize=8, verticalalignment='top', bbox=props,
                      fontfamily='monospace')
        
        plt.legend()
        plt.grid(True, alpha=0.2)
        
        # Save the plot with higher resolution
        plot_path = os.path.join(plot_dir, f"{os.path.splitext(os.path.basename(csv_path))[0]}_fit.png")
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Plot saved to {plot_path}")

    # Create result dictionary with all metrics
    result = {
        "file": os.path.basename(csv_path),
        "salt": salt_name,
        "a": a,
        "b": b,
        "R2": metrics['r2'],
        "Adj_R2": metrics['adj_r2'],
        "RMSE": metrics['rmse'],
        "MAE": metrics['mae'],
        "Max_Abs_Error": metrics['max_abs_error'],
        "Mean_Rel_Error": metrics['mean_rel_error'],
        "Median_Rel_Error": metrics['median_rel_error'],
        "status": "OK"
    }
    
    # Add peak positions to the result
    for i, pos in enumerate(peak_positions, 1):
        result[f"peak_{i}"] = pos
        
    return result

# --- Run for all CSVs in a folder ---
def create_residuals_plot(results, n_peaks=4, plot_dir="plots"):
    """Create a residuals plot for all salts."""
    plt.figure(figsize=(12, 8))
    
    for result in results:
        if 'status' in result and result['status'] != 'OK':
            continue
            
        # Get peak positions and predicted values
        salt_name = result['salt'].strip(',')
        a, b = result['a'], result['b']
        n_vals = np.arange(1, n_peaks + 1)
        predicted = power_law(n_vals, a, b)
        
        # Get actual peak positions (from manual or detected)
        actual = []
        for i in range(n_peaks):
            actual.append(float(result.get(f'peak_{i+1}', np.nan)))
        actual = np.array(actual)
        
        # Calculate residuals
        residuals = actual - predicted
        
        # Plot residuals
        plt.plot(n_vals, residuals, 'o-', label=f"{salt_name} (a={a:.2f}, b={b:.2f})")
    
    plt.axhline(0, color='black', linestyle='--', alpha=0.5)
    plt.xlabel('Peak Index (n)')
    plt.ylabel('Residual (Actual - Predicted) (Å)')
    plt.title('Residuals of Power Law Fits for All Salts')
    plt.xticks(range(1, n_peaks + 1))
    plt.grid(True, alpha=0.2)
    
    # Adjust legend
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    plt.tight_layout()
    
    # Save the plot
    residuals_path = os.path.join(plot_dir, "residuals_summary.png")
    plt.savefig(residuals_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Residuals plot saved to {residuals_path}")

def process_folder(folder_path, n_peaks=4, output_csv="fit_summary.csv", plot=False):
    results = []
    files = glob.glob(os.path.join(folder_path, "*.csv"))
    
    for csv_file in files:
        result = process_pdf(csv_file, n_peaks=n_peaks, plot=plot)
        results.append(result)

    # Save summary table
    summary = pd.DataFrame(results)
    
    # Create residuals plot if we have results
    if plot and len(results) > 0:
        create_residuals_plot(results, n_peaks=n_peaks)
    
    summary.to_csv(output_csv, index=False)
    print(f"\nSummary saved to {output_csv}")
    return summary


# --- Example usage ---
results_df = process_folder("PeakDecayCSVs", n_peaks=4, plot=True)
print(results_df)
