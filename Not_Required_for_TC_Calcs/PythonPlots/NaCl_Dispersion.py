import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib.patheffects as patheffects

# Set the font to Times New Roman for all text
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman']
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'Times New Roman'
plt.rcParams['mathtext.it'] = 'Times New Roman:italic'
plt.rcParams['mathtext.bf'] = 'Times New Roman:bold'

# Set default font sizes
plt.rcParams['font.size'] = 12
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['legend.fontsize'] = 12

N_A = 6.02214076e23  # mol^-1

def plot_dispersion_with_lambda_axis(data_sources,
                                     v_s, v_inf,
                                     omega_D=None, omega_F=None,
                                     v_Na=None,
                                     title='Phonon Dispersion'):

    fig, ax = plt.subplots(figsize=(10, 6))

    # Lattice constant and zone boundary
    a = 2.667190256/10  # nm
    Q_zone = np.pi / a  # Zone boundary
    
    # Create extended Q range for plotting
    Q_extended = np.linspace(0, 2*Q_zone, 100)
    
    # Dispersion relations (only up to zone boundary)
    mask = Q_extended <= Q_zone
    Q_plot = Q_extended[mask]
    
    # Calculate frequencies up to zone boundary
    omega_ad = v_s * Q_plot * 1e9 / 1e12
    omega_inf = v_inf * Q_plot * 1e9 / 1e12

    # Plot each data source
    for i, source in enumerate(data_sources):
        marker = source.get('marker', 'o')
        color = source.get('color', f'C{i}')
        
        # Check what type of data is provided and convert to frequency (THz)
        if 'v_vals' in source:  # If velocity data is provided in m/s
            omega_vals = (source['v_vals'] * source['Q']) * 1e9 / 1e12  # Convert m/s * nm^-1 -> THz
        elif 'omega' in source:  # If frequency data is provided directly in THz
            omega_vals = source['omega']
        else:
            raise ValueError("Each data source must have either 'v_vals' (m/s) or 'omega' (THz)")
            
        # Set marker style and face color based on data source
        if 'Liquid Na' in source['label']:
            # Empty circles for Liquid Na data
            ax.scatter(source['Q'], omega_vals, 
                   marker='o', facecolor='none', edgecolor=color, s=40,
                   label=source['label'])
        else:
            # Filled markers for other data
            marker_style = 'o' if 'Exp' in source['label'] else 's'  # 'o' for Exp, 's' for MD
            ax.scatter(source['Q'], omega_vals, 
                   marker=marker_style, color=color, s=40,
                   label=source['label'])
    
    # Plot theoretical curves with improved visibility and grayscale compatibility
    # First plot lines with zorder=1 to ensure they stay behind data points
    ax.plot(Q_plot, omega_ad, '-.', color='#2ca02c', linewidth=2.5, zorder=0,
           label=f'v$_{{s}}$ = {v_s:.0f} m/s',
           path_effects=[patheffects.withStroke(linewidth=3.5, foreground='white', alpha=0.7)])
    
    ax.plot(Q_plot, omega_inf, '--', color='#1f77b4', linewidth=2.5, zorder=0,
           label=f'v$_{{∞}}$ = {v_inf:.0f} m/s',
           path_effects=[patheffects.withStroke(linewidth=3.5, foreground='white', alpha=0.7)])
           
    # Plot liquid Na speed of sound line if provided
    if v_Na is not None:
        omega_Na = v_Na * Q_plot * 1e9 / 1e12
        ax.plot(Q_plot, omega_Na, '-', color='#ff7f0e', linewidth=2.5, zorder=0,
               label=f'v$_{{Na}}$ = {v_Na:.0f} m/s',
               path_effects=[patheffects.withStroke(linewidth=3.5, foreground='white', alpha=0.7)])
    
    # Add a subtle grid for better readability
    #ax.grid(True, which='both', linestyle=':', linewidth=0.5, alpha=0.3, zorder=0)
    
    # Add vertical line at zone boundary
    ax.axvline(x=Q_zone, color='black', linestyle=':', alpha=0.5,
              label=r'$\lambda_{\text{min}}$ = 2$a_{\text{lattice}}$ = 0.533 nm')

    # # Horizontal lines if provided
    # if omega_D is not None:
    #     ax.axhline(omega_D, color='purple', linestyle=':', linewidth=2,
    #                label=r'$\omega_D$ (Debye)')
    # if omega_F is not None:
    #     ax.axhline(omega_F, color='brown', linestyle=':', linewidth=2,
    #                label=r'$\omega_F$ (Frenkel)')

    ax.set_xlabel(r'$\mathbf{Q}$ [nm$^{-1}$]', fontweight='bold', fontsize=12)
    ax.set_ylabel(r'$\mathbf{\omega}$ [THz]', fontweight='bold', fontsize=12)
    #ax.set_title(title, fontweight='bold', fontsize=12)
    # ax.grid(True)
    # Create legend without border
    ax.legend(frameon=False)
    
    # Find the overall data range from all sources
    all_q = np.concatenate([source['Q'] for source in data_sources])
    x_min, x_max = np.nanmin(all_q), np.nanmax(all_q)
    ax.set_xlim(0, 35)  # Set fixed x-limit
    
    # Collect all y-values to set appropriate y-limits
    all_omega = []
    for source in data_sources:
        if 'v_vals' in source:
            all_omega.extend((source['v_vals'] * source['Q']) * 1e9 / 1e12)
        elif 'omega' in source:
            all_omega.extend(source['omega'])
    
    # Set y-limits with some padding
    y_max = max(all_omega) * 1.1  # 10% padding
    ax.set_ylim(0, 40)# y_max)
    
    # Create a simple secondary x-axis with wavelength in Angstroms
    # Calculate the wavelength values for the current Q range
    q_values = np.linspace(0.1, x_max * 1.1, 1000)  # Avoid Q=0
    lam_values = 2 * np.pi / q_values * 10  # Convert to Angstroms
    
    # Create a second axes that shares the same y-axis
    secax = ax.twiny()
    
    # Set the limits to match the primary x-axis
    secax.set_xlim(ax.get_xlim())
    secax.set_ylim(ax.get_ylim())  # Match primary y-limits
    
    # Set the ticks at specific wavelength values (in nm)
    angstrom_ticks = [0.2, 0.3, 0.5, 1, 2, 5]
    q_ticks = 2 * np.pi / (np.array(angstrom_ticks))  # Convert to Q space
    
    # Filter ticks to be within the plot range
    valid_ticks = (q_ticks >= 0.1) & (q_ticks <= x_max * 1.1)
    if np.any(valid_ticks):
        secax.set_xticks(q_ticks[valid_ticks])
        secax.set_xticklabels([f'{x:.1f}' for x, valid in zip(angstrom_ticks, valid_ticks) if valid], fontsize=12)
    
    secax.set_xlabel(r'$\lambda$ [nm]', fontweight='bold', fontsize=12)
    secax.xaxis.set_label_position('top')
    
    # Adjust layout with padding
    plt.tight_layout(pad=2.0)
    plt.show()



def debye_and_frenkel(M_kg_per_mol, rho_kg_m3, v_L, v_T, eta_pa_s, theta_D=None):

    # Number density (formula units / m^3)
    n = (rho_kg_m3 / M_kg_per_mol) * N_A

    # Debye k-vector
    k_D = (6 * np.pi**2 * n)**(1/3)

    # Debye-average sound velocity (standard Debye average)
    v_m = ((1/3) * (1.0/v_L**3 + 2.0/v_T**3))**(-1/3)
    
    # Physical constants
    hbar = 1.0545718e-34  # Reduced Planck constant (J·s)
    kB = 1.380649e-23     # Boltzmann constant (J/K)
    
    if theta_D is not None:
        # Calculate Debye frequency from Debye temperature
        omega_D = (kB * theta_D) / hbar  # in rad/s
        f_D_THz = omega_D / (2 * np.pi) / 1e12  # Convert to THz
    else:
        # Calculate Debye frequency from sound velocity (fallback)
        f_D_THz = (v_m * k_D) / (2 * np.pi) / 1e12

    # Infinite-frequency shear modulus estimate
    G_inf = rho_kg_m3 * v_T**2

    # Frenkel frequency in THz (f_F = G_inf / (2*pi*eta))
    f_F_THz = (G_inf / eta_pa_s) / (2 * np.pi) / 1e12

    return {
        "n_m3": n,
        "k_D_m-1": k_D,
        "v_m_m_s": v_m,
        "f_D_THz": f_D_THz,
        "G_inf_Pa": G_inf,
        "f_F_THz": f_F_THz
    }


# === 1. DEFINE DATA SOURCES ===

# Sound velocities (m/s)
v_s = 1660 # Demmel, 2021; Adiabatic speed in m/s
v_inf = 4400 # Demmel, 2021; 2702.4785990040846  # High-frequency speed in m/s

# Data source 1: Demmel, 2004
Q_demmel2004 = np.array([2.2052845528455274, 3.478658536585365, 4.680894308943089, 5.840447154471544, 
                     7.04979674796748, 8.294715447154474, 9.496951219512194, 11.005081300813007,
                     12.20020325203252])
# Velocity in m/s
v_demmel2004 = np.array([2280.68, 2786.89, 3073.17, 2877.82, 2812.70, 2874.13, 2681.23, 2485.87, 1991.95])

# Data source 2: Example second source (replace with actual data)
Q_demmel2021 = np.array([2.2197140707298715, 3.4988713318284397, 4.683972911963883, 5.8690744920993225, 7.054176072234762, 8.314522197140704, 9.499623777276145, 11.004514672686234, 12.208427389014297, 13.487584650112865, 14.691497366440935, 15.838976674191123, 17.042889390519193, 18.265613243039883, 19.5071482317532, 22.197140707298722, 24.548532731376977])
v_demmel2021 = np.array([2593.4579439252343, 3275.7009345794395, 3733.6448598130846, 3560.747663551402, 3233.6448598130846, 3242.990654205607, 3121.4953271028053, 3000, 2700.9345794392525, 2401.869158878505, 2088.785046728972, 1887.8504672897193, 1696.2616822429904, 1700.934579439252, 1574.7663551401856, 1210.2803738317757, 1014.0186915887853])  # Example velocities in m/s

# Data source 3: 
Q_bryk2003_sound = np.array([1.958360442, 2.482108003, 3.529603123, 4.964216005, 8.58490566, 12.41054001, 17.19258295, 21.6330514, 25.77748861, 30.17241379, 34.38516591])
omega_bryk2003_sound = np.array([6.832627119, 8.73940678, 13.29449153, 18.96186441, 29.60805085, 35.38135593, 35.16949153, 37.4470339, 32.78601695, 33.26271186, 35.80508475])  # THz

# Data source 4:
Q_Na_Demmel2004 = np.array([1.5010162601626016, 2.2977642276422765, 3.9126016260162597, 4.702235772357723, 6.359756097560974, 7.113821138211381, 8.842479674796746, 9.496951219512194, 10.692073170731707, 11.296747967479675, 12.050813008130081])
v_Na_Demmel2004 = np.array([2544.8464163822528, 2915.904436860068, 2982.2525597269623, 3030.1706484641636, 2811.467576791809, 2854.4709897610924, 2527.6450511945395, 2450.238907849829, 2184.8464163822528, 2023.890784982935, 1877.6791808873722])

# # Data source 4: 
# Q_bryk2003_heat = np.array([1.981132075, 2.504879636, 3.50683149, 4.964216005, 8.607677293, 12.43331165, 17.19258295, 21.65582303, 25.80026025, 30.21795706, 34.40793754])
# omega_bryk2003_heat = np.array([0, 0.529661017, 0.953389831, 9.480932203, 12.5529661, 15.20127119, 19.01483051, 22.93432203, 25.15889831, 25.15889831, 26.90677966])  # THz

# Package data sources
data_sources = [
    {
        'Q': Q_demmel2004,
        'v_vals': v_demmel2004,  # Velocity in m/s
        'label': 'Molten NaCl - Demmel, 2004 (Exp)',
        'marker': 'o',
        'color': 'blue'
    },
    {
        'Q': Q_demmel2021,
        'v_vals': v_demmel2021,  # Velocity in m/s
        'label': 'Molten NaCl - Demmel, 2021 (Exp)',
        'marker': 's',
        'color': 'green'
    },
    {
        'Q': Q_bryk2003_sound,
        'omega': omega_bryk2003_sound,  # Frequency in THz
        'label': 'Molten NaCl - Bryk, 2003 (MD)',
        'marker': 'o',
        'color': 'red'
    },
    {
        'Q': Q_Na_Demmel2004,
        'v_vals': v_Na_Demmel2004,  # Frequency in THz
        'label': 'Liquid Na - Demmel, 2004 (Exp)',
        'marker': 's',
        'color': 'black'
    },
    # {
    #     'Q': Q_bryk2003_heat,
    #     'omega': omega_bryk2003_heat,  # Frequency in THz
    #     'label': 'Heat - Bryk et al. (2003)',
    #     'marker': 's',
    #     'color': 'green'
    # }
]

# Calculate Debye frequency (example calculation, adjust as needed)
M_NaCl = 0.05844   # kg/mol
rho_NaCl = 1.56E+03  # kg/m3
vL_NaCl = 1743.7  # m/s (longitudinal)
vT_NaCl = vL_NaCl/2   # m/s (transverse)
eta_NaCl = 0.001032139   # Pa·s (viscosity)
# theta_D = 178   # K

na = debye_and_frenkel(M_NaCl, rho_NaCl, vL_NaCl, vT_NaCl, eta_NaCl)
#na['f_D_THz'] = 6.41768

print("NaCl estimates:", na)

# Create the plot
plot_dispersion_with_lambda_axis(
    data_sources=data_sources,
    v_s=v_s,
    v_inf=v_inf,
    omega_D=na['f_D_THz'],
    omega_F=na['f_F_THz'],
    v_Na=2720,  # Speed of sound in liquid Na (m/s)
    title='Phonon Dispersion in NaCl (1170K)'
)
