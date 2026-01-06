import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D # For custom legend
from matplotlib.patches import ConnectionPatch
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def create_parity_plot(csv_filepath):
    """
    Generates a parity plot comparing SCM and KTM models against
    experimental data from the provided CSV file.
    """
    
    # --- 1. Load and Prepare Data ---
    try:
        # Read the CSV with headers in the first row (0-indexed)
        # Don't set any column as index to keep all data in columns
        df = pd.read_csv(csv_filepath, header=0)
    except FileNotFoundError:
        print(f"Error: The file '{csv_filepath}' was not found.")
        return
    except Exception as e:
        print(f"Error reading CSV: {e}")
        print("Please ensure the file is a valid CSV and the header is on the second row.")
        return

    # Clean data: drop rows with missing essential values
    # Note: Column names must match exactly with the CSV header
    required_cols = ['k_exp', 'k_KTM', 'k_SCM', 'salt type', 'salt abbreviated composition']
    
    # Check if all required columns exist in the dataframe
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f"Error: The following required columns are missing: {', '.join(missing_cols)}")
        print(f"Available columns: {', '.join(df.columns)}")
        return
        
    # Drop rows with missing values in any of the required columns
    df.dropna(subset=required_cols, inplace=True)

    # Ensure numeric columns are numeric
    for col in ['k_exp', 'k_KTM', 'k_SCM']:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    
    # Drop rows that became NaN during coercion
    df.dropna(subset=['k_exp', 'k_KTM', 'k_SCM'], inplace=True)

    # --- 2. Set up Plot ---
    # Apply publication-style matplotlib settings to match project
    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.rm'] = 'Times New Roman'
    plt.rcParams['font.size'] = 12
    plt.rcParams['axes.labelsize'] = 14
    plt.rcParams['axes.labelweight'] = 'bold'
    plt.rcParams['axes.linewidth'] = 1.5
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12
    plt.rcParams['xtick.direction'] = 'out'
    plt.rcParams['ytick.direction'] = 'out'
    plt.rcParams['xtick.major.width'] = 1.5
    plt.rcParams['ytick.major.width'] = 1.5
    plt.rcParams['legend.frameon'] = False
    fig, ax = plt.subplots(figsize=(8, 7))
    
    # Define mappings for plot styles
    # Models: KTM=blue squares; SCM=red circles
    model_color = {'KTM': 'C0', 'SCM': '#d62728'}
    model_marker = {'KTM': 's', 'SCM': 'o'}

    # --- 3. Plot Error Regions and Parity Line ---
    
    # Fixed axis limits per request
    limit_min, limit_max = 0.0, 1.5
    # Create x-values for the lines over fixed range
    x_line = np.array([limit_min, limit_max])

    # Plot 30% deviation region (back) in warm yellow
    ax.fill_between(x_line, x_line * 0.70, x_line * 1.30,
                    color='#FFF59D', alpha=0.7)
    
    # Plot 15% deviation region (front) in light green
    ax.fill_between(x_line, x_line * 0.85, x_line * 1.15,
                    color='#C5E1A5', alpha=0.8)
    
    # Plot parity line (no legend label)
    ax.plot(x_line, x_line, 'k-')
    
    # Percent labels at correct band positions within limits
    try:
        # Choose x so that y = (1+e)*x is within [0, limit_max]
        x15p = limit_max / 1.15
        y15p = 1.15 * x15p
        x28p = limit_max / 1.30
        y28p = 1.30 * x28p
        ax.text(x15p, y15p, '+15%', fontsize=10, ha='right', va='bottom')
        ax.text(x28p, y28p, '+28%', fontsize=10, ha='right', va='bottom')
        # Negative bands can simply use x at max
        ax.text(limit_max, 0.85 * limit_max, '-15%', fontsize=10, ha='right', va='top')
        ax.text(limit_max, 0.70 * limit_max, '-28%', fontsize=10, ha='right', va='top')
    except Exception:
        pass

    # --- 4. Plot all points in main plot (no labels) ---
    for _, row in df.iterrows():
        salt_type = str(row['salt type']).strip()
        k_exp = float(row['k_exp'])
        k_scm = float(row['k_SCM'])
        k_ktm = float(row['k_KTM'])
        
        # Plot KTM point
        if salt_type == 'Actinides':
            ax.plot(k_exp, k_ktm, marker='+', color=model_color['KTM'], markersize=8, linestyle='none')
        else:
            face_ktm = model_color['KTM'] if salt_type == 'Dissociating' else 'none'
            ax.plot(k_exp, k_ktm, marker=model_marker['KTM'],
                   markerfacecolor=face_ktm, markeredgecolor=model_color['KTM'],
                   markersize=8, linestyle='none')
        
        # Plot SCM point
        if salt_type == 'Actinides':
            ax.plot(k_exp, k_scm, marker='x', color=model_color['SCM'], markersize=8, linestyle='none')
        else:
            face_scm = model_color['SCM'] if salt_type == 'Dissociating' else 'none'
            ax.plot(k_exp, k_scm, marker=model_marker['SCM'],
                   markerfacecolor=face_scm, markeredgecolor=model_color['SCM'],
                   markersize=8, linestyle='none')

    # --- 5. Finalize Plot (Labels, Legend, Layout) ---
    
    # Set labels and title
    ax.set_xlabel('k$_{\mathrm{Exp.}}$ [W·m$^{-1}$·K$^{-1}$]')
    ax.set_ylabel('k$_{\mathrm{Predicted}}$ [W·m$^{-1}$·K$^{-1}$]')
    # ax.set_title('Parity Plot: Model vs. Experimental Thermal Conductivity', fontsize=14)
    
    # Set plot limits and aspect
    ax.set_xlim(limit_min, limit_max)
    ax.set_ylim(limit_min, limit_max)
    ax.set_aspect('equal', 'box')
    # No gridlines
    
    # --- 6. Table-like key (KTM/SCM vs salt type) ---
    # Create a custom legend-like table
    legend_elements = [
        # KTM column
        Line2D([0], [0], marker=model_marker['KTM'], color='w', 
               markerfacecolor=model_color['KTM'], markersize=8, 
               label='KTM', markeredgecolor=model_color['KTM']),
        Line2D([0], [0], marker=model_marker['SCM'], color='w',
               markerfacecolor=model_color['SCM'], markersize=8,
               label='SCM', markeredgecolor=model_color['SCM']),
        # Dissociating row
        Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', 
               markersize=8, label='Dissociating'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='gray',
               markersize=8, label=''),
        # Complexing row
        Line2D([0], [0], marker='o', color='gray', markerfacecolor='none',
               markersize=8, label='Complexing'),
        Line2D([0], [0], marker='o', color='gray', markerfacecolor='none',
               markersize=8, label=''),
        # Actinides row
        Line2D([0], [0], marker='+', color=model_color['KTM'], 
               markersize=8, label='Actinides', linewidth=1.5),
        Line2D([0], [0], marker='x', color=model_color['SCM'],
               markersize=8, label='', linewidth=1.5)
    ]
    
    # Create a single legend with two columns
    leg = ax.legend(handles=legend_elements, 
                   loc='upper left', 
                   ncol=2,
                   handlelength=1.5,
                   handletextpad=0.5,
                   columnspacing=0.5,
                   borderpad=0.5,
                   fontsize=9,
                   frameon=False)
    
    # Adjust the legend to look like a table
    for i, text in enumerate(leg.get_texts()):
        if i % 2 == 0:  # Left column (labels)
            text.set_ha('left')
            text.set_position((0, 0))
        else:  # Right column (empty)
            text.set_visible(False)
    
    # --- 7. Create subplots for main plot and zoomed plot ---
    # Create a figure with two subplots side by side
    fig = plt.figure(figsize=(12, 6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1], wspace=0.3)
    
    # Main plot (left)
    ax = fig.add_subplot(gs[0])
    
    # Plot all points in main plot (no labels)
    for _, row in df.iterrows():
        salt_type = str(row['salt type']).strip()
        k_exp = float(row['k_exp'])
        k_scm = float(row['k_SCM'])
        k_ktm = float(row['k_KTM'])
        
        if salt_type == 'Actinides':
            ax.plot(k_exp, k_ktm, marker='+', color=model_color['KTM'], markersize=8, linestyle='none')
            ax.plot(k_exp, k_scm, marker='x', color=model_color['SCM'], markersize=8, linestyle='none')
        else:
            face_ktm = model_color['KTM'] if salt_type == 'Dissociating' else 'none'
            face_scm = model_color['SCM'] if salt_type == 'Dissociating' else 'none'
            ax.plot(k_exp, k_ktm, marker=model_marker['KTM'], 
                   markerfacecolor=face_ktm, markeredgecolor=model_color['KTM'],
                   markersize=8, linestyle='none')
            ax.plot(k_exp, k_scm, marker=model_marker['SCM'],
                   markerfacecolor=face_scm, markeredgecolor=model_color['SCM'],
                   markersize=8, linestyle='none')
    
    # Add deviation bands and parity line to main plot
    ax.fill_between(x_line, x_line * 0.70, x_line * 1.30, color='#FFF59D', alpha=0.7)
    ax.fill_between(x_line, x_line * 0.85, x_line * 1.15, color='#C5E1A5', alpha=0.8)
    ax.plot(x_line, x_line, 'k-')
    
    # Set main plot labels and limits
    ax.set_xlabel('k$_{\mathrm{Exp.}}$ [W·m$^{-1}$·K$^{-1}$]')
    ax.set_ylabel('k$_{\mathrm{Predicted}}$ [W·m$^{-1}$·K$^{-1}$]')
    ax.set_xlim(limit_min, limit_max)
    ax.set_ylim(limit_min, limit_max)
    ax.set_aspect('equal', 'box')
    
    # Add zoom region rectangle
    zoom_xlim = (0.3, 0.7)
    zoom_ylim = (0.2, 0.8)
    x0, x1 = zoom_xlim
    y0, y1 = zoom_ylim
    ax.plot([x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0], 'k--')
    
    # Zoomed plot (right)
    ax_zoom = fig.add_subplot(gs[1])
    
    # Plot the same content in zoom plot but with labels
    ax_zoom.fill_between(x_line, x_line * 0.70, x_line * 1.30, color='#FFF59D', alpha=0.7)
    ax_zoom.fill_between(x_line, x_line * 0.85, x_line * 1.15, color='#C5E1A5', alpha=0.8)
    ax_zoom.plot(x_line, x_line, 'k-')
    
    # Plot points in zoom region with labels
    for _, row in df.iterrows():
        salt_type = str(row['salt type']).strip()
        k_exp = float(row['k_exp'])
        k_scm = float(row['k_SCM'])
        k_ktm = float(row['k_KTM'])
        
        # Only plot points that are in the zoom region
        if (zoom_xlim[0] <= k_exp <= zoom_xlim[1] and 
            zoom_ylim[0] <= min(k_scm, k_ktm) and max(k_scm, k_ktm) <= zoom_ylim[1]):
            
            if salt_type == 'Actinides':
                ax_zoom.plot(k_exp, k_ktm, marker='+', color=model_color['KTM'], markersize=8, linestyle='none')
                ax_zoom.plot(k_exp, k_scm, marker='x', color=model_color['SCM'], markersize=8, linestyle='none')
            else:
                face_ktm = model_color['KTM'] if salt_type == 'Dissociating' else 'none'
                face_scm = model_color['SCM'] if salt_type == 'Dissociating' else 'none'
                ax_zoom.plot(k_exp, k_ktm, marker=model_marker['KTM'], 
                           markerfacecolor=face_ktm, markeredgecolor=model_color['KTM'],
                           markersize=8, linestyle='none')
                ax_zoom.plot(k_exp, k_scm, marker=model_marker['SCM'],
                           markerfacecolor=face_scm, markeredgecolor=model_color['SCM'],
                           markersize=8, linestyle='none')
            
            # Add label
            ax_zoom.annotate(
                str(row['salt abbreviated composition']),
                xy=(k_exp, (k_scm + k_ktm)/2),
                xytext=(5, 0),
                textcoords='offset points',
                fontsize=9,
                color='0.3',
                va='center'
            )
    
    # Style zoom plot
    ax_zoom.set_xlim(*zoom_xlim)
    ax_zoom.set_ylim(*zoom_ylim)
    ax_zoom.set_aspect('equal', 'box')
    ax_zoom.set_xlabel('k$_{\mathrm{Exp.}}$ [W·m$^{-1}$·K$^{-1}$]')
    ax_zoom.set_ylabel('k$_{\mathrm{Predicted}}$ [W·m$^{-1}$·K$^{-1}$]')
    ax_zoom.set_title('Zoomed Region', fontsize=12, pad=10)
    
    # Add zoom indicator to main plot
    ax.text(0.05, 0.95, 'Zoomed region →', 
            transform=ax.transAxes, ha='left', va='top', 
            bbox=dict(facecolor='white', alpha=0.8, edgecolor='none', pad=4))
    
    plt.tight_layout()
    plt.show()

# --- Main execution ---
if __name__ == "__main__":
    # The user must place their CSV file in the same directory as this script,
    # or provide the correct relative/absolute path.
    csv_file_path = 'Not_Required_for_TC_Calcs/Ppts for Plots/SaltPlots.csv'
    create_parity_plot(csv_file_path)