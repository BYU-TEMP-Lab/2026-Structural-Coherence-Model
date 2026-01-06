import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

def create_parity_plot(df, output_dir='Comparison_Plot_Figs'):
    """
    Create a parity plot comparing experimental and predicted thermal conductivity values
    """

    # Apply publication-style matplotlib settings
    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.rm'] = 'Times New Roman'
    plt.rcParams['font.size'] = 12
    plt.rcParams['axes.labelsize'] = 14
    plt.rcParams['axes.labelweight'] = 'bold'
    plt.rcParams['axes.linewidth'] = 1.5
    plt.rcParams['xtick.labelsize'] = 10
    plt.rcParams['ytick.labelsize'] = 12
    plt.rcParams['xtick.direction'] = 'out'
    plt.rcParams['ytick.direction'] = 'out'
    plt.rcParams['xtick.major.width'] = 1.5
    plt.rcParams['ytick.major.width'] = 1.5
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['xtick.direction'] = 'out'
plt.rcParams['ytick.direction'] = 'out'
plt.rcParams['xtick.major.width'] = 1.5
plt.rcParams['ytick.major.width'] = 1.5
# Legend and small text/annotation sizing (~0.85x of tick labels)
_tick_size = 12
_small_text = int(round(0.85 * _tick_size))
plt.rcParams['legend.frameon'] = False
plt.rcParams['legend.fontsize'] = _small_text

fig, ax = plt.subplots(figsize=(6,6))

# 1:1 line and ±15% region
x = np.linspace(0, 1.5, 200)
ax.plot(x, x, 'k', lw=1)
ax.fill_between(x, 0.85*x, 1.15*x, color='green', alpha=0.15)

# === Marker styles ===
marker_style = {
    ('Simple', 'KT'): dict(marker='s', facecolors='blue', edgecolors='blue'),
    ('Complexing', 'KT'): dict(marker='s', facecolors='none', edgecolors='blue'),
    ('Actinide', 'KT'): dict(marker='+', color='blue'),
    ('Simple', 'SCM'): dict(marker='o', facecolors='red', edgecolors='red'),
    ('Complexing', 'SCM'): dict(marker='o', facecolors='none', edgecolors='red'),
    ('Actinide', 'SCM'): dict(marker='x', color='red'),
}

# === Plot data ===
for (cat, model), style in marker_style.items():
    subset = df[df['Type'] == cat]
    ax.scatter(subset['k_Exp'], subset[model], s=60, **style)

# === Annotate key salts ===
annotations = {
    '0.34BeF2-0.66LiF': 'FLiBe',
    '0.42KF-0.465LiF-0.115NaF': 'FLiNaK',
    '0.59KF-0.065MgF2-0.345NaF': 'FMgNaK',
    '0.63NaCl-0.37UCl3': 'NaCl-UCl3'
}
for salt, label in annotations.items():
    row = df[df['Salt Composition'] == salt]
    if not row.empty:
        xval, yval = row['k_Exp'].values[0], row['SCM'].values[0]
        ax.text(xval + 0.02, yval, label, fontsize=8, color='gray')
        ax.plot([xval, xval + 0.02], [yval, yval], color='gray', lw=0.8)
# === Custom Table Legend ===
# Create a custom table-like legend instead of standard legend
legend_x = 0.02  # Left position
legend_y = 0.98  # Top position
row_height = 0.04  # Vertical spacing between rows
col1_width = 0.35   # Wider first column for labels
col2_width = 0.15   # Second column for KT
col3_width = 0.15   # Third column for SCM

# Row 1: Headers
ax.text(legend_x, legend_y, '', fontsize=9, fontweight='bold')
ax.text(legend_x + col1_width, legend_y, 'KT', fontsize=9, fontweight='bold', ha='center')
ax.text(legend_x + col1_width + col2_width, legend_y, 'SCM', fontsize=9, fontweight='bold', ha='center')

# Row 2: Simple (Dissociating)
ax.text(legend_x, legend_y - row_height, 'Dissociating', fontsize=8, va='center')
# Use the marker style from the dictionary
simple_kt_style = marker_style[('Simple', 'KT')]
ax.scatter(legend_x + col1_width, legend_y - row_height, s=40, **simple_kt_style)
simple_scm_style = marker_style[('Simple', 'SCM')]
ax.scatter(legend_x + col1_width + col2_width, legend_y - row_height, s=40, **simple_scm_style)

# Row 3: Complexing
ax.text(legend_x, legend_y - 2*row_height, 'Complexing', fontsize=8, va='center')
complexing_kt_style = marker_style[('Complexing', 'KT')]
ax.scatter(legend_x + col1_width, legend_y - 2*row_height, s=40, **complexing_kt_style)
complexing_scm_style = marker_style[('Complexing', 'SCM')]
ax.scatter(legend_x + col1_width + col2_width, legend_y - 2*row_height, s=40, **complexing_scm_style)

# Row 4: Actinide
ax.text(legend_x, legend_y - 3*row_height, 'Actinides', fontsize=8, va='center')
actinide_kt_style = marker_style[('Actinide', 'KT')]
ax.scatter(legend_x + col1_width, legend_y - 3*row_height, s=50, **actinide_kt_style)
actinide_scm_style = marker_style[('Actinide', 'SCM')]
ax.scatter(legend_x + col1_width + col2_width, legend_y - 3*row_height, s=50, **actinide_scm_style)

# Row 5: MAE (%)
ax.text(legend_x, legend_y - 4*row_height, 'MAE (%)', fontsize=8, va='center')
ax.text(legend_x + col1_width, legend_y - 4*row_height, f'{mae_kt:.1f}', fontsize=8, ha='center', va='center')
ax.text(legend_x + col1_width + col2_width, legend_y - 4*row_height, f'{mae_scm:.1f}', fontsize=8, ha='center', va='center')

# === Labels and styling ===
ax.set_xlabel(r'$k_{\mathrm{Exp.}}\ [\mathrm{W\cdot m^{-1}\cdot K^{-1}}]$', fontsize=12, fontweight='bold')
ax.set_ylabel(r'$k_{\mathrm{Predicted}}\ [\mathrm{W\cdot m^{-1}\cdot K^{-1}}]$', fontsize=12, fontweight='bold')
ax.set_xlim(0, 1.5)
ax.set_ylim(0, 1.5)
ax.set_aspect('equal', adjustable='box')

# ±15% labels
ax.text(1.4, 1.32, '+15%', color='darkgreen', fontsize=9)
ax.text(1.4, 1.05, '-15%', color='darkgreen', fontsize=9)
ax.text(1.0, 0.05, r'$T = T_{\mathrm{Exp.\ Min}}$', color='gray', fontsize=8)

plt.tight_layout()
plt.show()
return fig, ax
