# At the top of your imports, add:
from matplotlib import rcParams
import matplotlib.pyplot as plt
import matplotlib.patheffects as patheffects
import numpy as np

# Set up Times New Roman font
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'Times New Roman'
plt.rcParams['mathtext.it'] = 'Times New Roman:italic'
plt.rcParams['mathtext.bf'] = 'Times New Roman:bold'

# Parameters
r = np.linspace(0, 20, 1000)  # Range of r values

def b_func(r,
           r_peak=1,
           b_inf=0.5,
           A=0.9,
           sigma=0.35,
           B=-0.2,
           k=0.9,
           omega=4,
           slope=9):
    """
    Function that stays near 0 before r_peak,
    rises steeply (adjustable slope),
    then decays with oscillations to b_inf.
    """
    # Smooth step (logistic) instead of Heaviside → avoids sharp jump
    smooth_step = 1 / (1 + np.exp(-slope * (r - r_peak)))

    peak_term = A * np.exp(-((r - r_peak)/sigma)**2)  # Gaussian peak
    decay_term = B * np.exp(-k*(r - r_peak)) * np.cos(omega*(r - r_peak))  # damped oscillations

    return smooth_step * (b_inf + peak_term + decay_term)


# Define the disruption factors b_i as a function of r
def get_b_values(r):
    return b_func(r)


# Calculate β(r) = b(r)/(1 - b(r))  (vectorized)
def calculate_beta(b_values):
    safe_b = np.where(b_values >= 1, 0.99999, b_values)
    return safe_b / (1 - safe_b)


# Calculate S(r) = exp(-∫β(r)dr)
def calculate_survival(r, b_values):
    beta = calculate_beta(b_values)
    integral = np.cumsum(beta) * (r[1] - r[0])  # trapezoidal could be better
    return np.exp(-integral), beta


# Calculate values
b_values = get_b_values(r)
S_r, beta_values = calculate_survival(r, b_values)

# Define colors
S_COLOR = '#0c4b8e'  # Darker blue
B_COLOR = '#cc0000'  # red

# Create the plot
fig, ax1 = plt.subplots(figsize=(6, 4))
fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.15)  # Add padding
ax1.tick_params(axis='both', which='major', labelsize=14)

# Plot S(r) on primary y-axis
ax1.set_xlabel('r [a.u.]', fontsize=14, fontweight='bold')
ax1.set_ylabel('S(r) [a.u.]', color=S_COLOR, fontsize=14, fontweight='bold')

# Plot S(r) and b(r) lines
ax1.plot(r, S_r, color=S_COLOR, linewidth=2.5, linestyle='-',
         path_effects=[patheffects.withStroke(linewidth=3.5, foreground='white', alpha=0.7)],
         label='S(r)')
ax1.plot(r, b_values, color=S_COLOR, linewidth=2, linestyle='--',  # Changed to S_COLOR
         path_effects=[patheffects.withStroke(linewidth=3, foreground='white', alpha=0.7)],
         label='$\sum$b$_i$(r)')
ax1.tick_params(axis='y', labelcolor=S_COLOR)

# Set y-limits to show the 10x relationship
s_max = 1.2
ax1.set_ylim(0, s_max)
ax1.set_xlim(0, 5)
# ax1.grid(True, linestyle=':', alpha=0.3, zorder=0)

# Create second y-axis for β(r)
ax2 = ax1.twinx()
ax2.tick_params(axis='both', which='major', labelsize=14)
ax2.set_ylabel('β(r) [10 a.u.]', color=B_COLOR, fontsize=14, fontweight='bold')
ax2.plot(r, beta_values, color=B_COLOR, linewidth=2, linestyle='-.',
         path_effects=[patheffects.withStroke(linewidth=3, foreground='white', alpha=0.7)],
         label='β(r)')
ax2.tick_params(axis='y', labelcolor=B_COLOR)
ax2.set_ylim(0, 12)  # 10x of s_max

# Add legend with better visibility
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
# Update the legend to remove the box and use the same color for S(r) and b(r)
legend = ax1.legend(lines1 + lines2, labels1 + labels2, 
                   loc='upper right', 
                   fontsize=14,
                   frameon=False)  # This removes the box
legend.get_frame().set_alpha(0.9)

plt.tight_layout()

# Calculate the area under S(r) to find l_scl
l_scl = np.trapz(S_r, r)  # Numerical integration using trapezoidal rule

# Find the index where r first exceeds l_scl
l_scl_idx = np.argmax(r >= l_scl)

if l_scl_idx > 0:
    # Add vertical line
    ax1.axvline(x=l_scl, color='black', linestyle=':', alpha=0.9, 
                linewidth=1.5, zorder=0)
    
    # Position the label inside the plot at 80% of the y-range
    y_pos = 0.8 * (ax1.get_ylim()[1] - ax1.get_ylim()[0]) + ax1.get_ylim()[0]
    
    # Add text inside the plot
    text = ax1.text(l_scl + 0.22, y_pos,  # Position inside plot
                   f'$l_{{\\mathit{{scl}}}} = {l_scl:.2f}$ a.u.',
                   color='black', 
                   fontsize=14,  # Slightly smaller font
                   fontfamily='Times New Roman',
                   bbox=dict(facecolor='white', alpha=0.8, edgecolor='none', pad=2),
                   zorder=10)
    
    # Add arrow from text to line
    ax1.annotate('', 
                xy=(l_scl, y_pos+0.001),  # Point to the line
                xytext=(l_scl + 0.2, y_pos+0.001),  # Start from text
                arrowprops=dict(facecolor='black', 
                              arrowstyle='-|>', 
                              linewidth=1,
                              shrinkA=0, 
                              shrinkB=0),
                zorder=10)

plt.show()
