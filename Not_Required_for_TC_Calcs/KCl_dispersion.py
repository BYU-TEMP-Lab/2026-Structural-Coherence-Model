import numpy as np
import matplotlib.pyplot as plt

N_A = 6.02214076e23  # mol^-1

def Q_to_lambda(Q):
    """Convert Q [nm^-1] to wavelength [nm].
    
    Parameters
    ----------
    Q : float or array-like
        Wavevector magnitude in nm^-1
        
    Returns
    -------
    float or array-like
        Wavelength in nm, or np.nan for Q ≤ 0
    """
    with np.errstate(divide='ignore', invalid='ignore'):
        lam = 2*np.pi / Q
        if isinstance(lam, np.ndarray):
            lam[Q <= 0] = np.nan
        elif Q <= 0:
            lam = np.nan
        return lam

def lambda_to_Q(lam):
    """Convert wavelength [nm] to Q [nm^-1].
    
    Parameters
    ----------
    lam : float or array-like
        Wavelength in nm
        
    Returns
    -------
    float or array-like
        Wavevector magnitude in nm^-1, or np.nan for λ ≤ 0
    """
    with np.errstate(divide='ignore', invalid='ignore'):
        Q = 2*np.pi / lam
        if isinstance(lam, np.ndarray):
            Q[lam <= 0] = np.nan
        elif lam <= 0:
            Q = np.nan
        return Q

def plot_dispersion_with_lambda_axis(Q_vals, y_vals,
                                     v_ad, v_inf,
                                     omega_D=None, omega_F=None,
                                     data_type="velocity",
                                     label="data",
                                     title='Phonon Dispersion'):
    """
    Plot ω(Q) with a secondary x-axis showing wavelength λ = 2π/Q.

    Parameters
    ----------
    Q_vals : array
        Momentum transfer values [nm^-1].
    y_vals : array
        Either velocity v(Q) [m/s] if data_type="velocity",
        or energy E(Q) [meV] if data_type="energy".
    v_ad : float
        Adiabatic sound velocity [m/s].
    v_inf : float
        High-frequency sound velocity [m/s].
    omega_D : float, optional
        Debye frequency [THz].
    omega_F : float, optional
        Frenkel frequency [THz].
    data_type : str
        "velocity" or "energy".
    label : str
        Label for the dataset.
    title : str
        Plot title.
    """

    # Convert input to ω(Q) in THz
    if data_type == "velocity":
        omega_vals = (y_vals * Q_vals) / 1000.0  # m/s * nm^-1 -> THz
    elif data_type == "energy":
        omega_vals = y_vals / 0.2418  # meV -> THz
    else:
        raise ValueError("data_type must be 'velocity' or 'energy'")

    # Lattice constant and zone boundary
    a = 2.659012/10  # nm
    Q_zone = np.pi / a  # Zone boundary
    
    # Create extended Q range for plotting
    Q_extended = np.linspace(0, 2*Q_zone, 100)
    
    # Folded dispersion relations
    omega_ad = np.where(Q_extended <= Q_zone,
                       v_ad * Q_extended / 1000.0,  # Original below zone boundary
                       v_ad * (2*Q_zone - Q_extended) / 1000.0)  # Folded above zone boundary
    
    omega_inf = np.where(Q_extended <= Q_zone,
                        v_inf * Q_extended / 1000.0,  # Original below zone boundary
                        v_inf * (2*Q_zone - Q_extended) / 1000.0)  # Folded above zone boundary

    fig, ax = plt.subplots(figsize=(8,5))
    ax.plot(Q_vals, omega_vals, 'o-', label=label)
    
    # Plot folded dispersion relations
    ax.plot(Q_extended, omega_ad, '--', color='C1', 
           label=f'Adiabatic speed: {v_ad} m/s')
    ax.plot(Q_extended, omega_inf, '--', color='C2',
           label=f'High-frequency speed: {v_inf} m/s')
    
    # Add vertical line at zone boundary
    ax.axvline(x=Q_zone, color='black', linestyle=':', alpha=0.5,
              label=rf'Zone boundary $\pi/a$')

    # Horizontal lines if provided
    if omega_D is not None:
        ax.axhline(omega_D, color='purple', linestyle=':', linewidth=2,
                   label=r'$\omega_D$ (Debye)')
    if omega_F is not None:
        ax.axhline(omega_F, color='brown', linestyle=':', linewidth=2,
                   label=r'$\omega_F$ (Frenkel)')

    ax.set_xlabel(r'$Q$ [nm$^{-1}$]')
    ax.set_ylabel(r'$\omega$ [THz]')
    ax.set_title(title)
    ax.grid(True)
    ax.legend()

    # Add vertical lines for λ = 2a
    a1 = 2.659012/10  # nm
    a2 = 3.93805/10   # nm
    
    # Convert λ = 2a to Q values
    Q1 = 2 * np.pi / (2 * a1)
    Q2 = 2 * np.pi / (2 * a2)
    
    # Add vertical lines with labels
    ax.axvline(x=Q1, color='gray', linestyle='--', alpha=0.7, 
              label=rf'$\lambda = 2a_1$ ({2*a1:.2E} nm)')
    ax.axvline(x=Q2, color='darkgray', linestyle='--', alpha=0.7,
              label=rf'$\lambda = 2a_2$ ({2*a2:.2E} nm)')
    
    # Update legend to include the new lines
    ax.legend()
    
    # Set reasonable axis limits
    x_min, x_max = np.nanmin(Q_vals), np.nanmax(Q_vals)
    ax.set_xlim(0, x_max * 1.1)  # Add 10% padding
    
    # Create a simple secondary axis with a fixed scale
    # Calculate the wavelength range that matches our Q range
    q_ticks = np.linspace(0, x_max, 5)
    q_ticks = q_ticks[q_ticks > 0]  # Remove zero to avoid division by zero
    lam_ticks = 2 * np.pi / q_ticks
    
    # Create the secondary axis with fixed ticks
    secax = ax.secondary_xaxis('top')
    secax.set_xticks(q_ticks)
    secax.set_xticklabels([f'{lam:.1f}' for lam in lam_ticks])
    secax.set_xlabel(r'Wavelength $\lambda$ [nm]')
    
    # Adjust layout with padding
    plt.tight_layout(pad=2.0)
    plt.show()



def debye_and_frenkel(M_kg_per_mol, rho_kg_m3, v_L, v_T, eta_pa_s):
    # number density (formula units / m^3)
    n = (rho_kg_m3 / M_kg_per_mol) * N_A

    # Debye k-vector
    k_D = (6 * np.pi**2 * n)**(1/3)

    # Debye-average sound velocity (standard Debye average)
    v_m = ((1/3) * (1.0/v_L**3 + 2.0/v_T**3))**(-1/3)

    # Debye frequency in THz (use f = v_m * k_D / (2*pi))
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


# === 1. INPUT DATA ===

###### NaCl, 1170K - Demmel, 2004 ######
# Q values in nm^-1 (from image)
Q_vals = np.array([2.2052845528455274, 3.478658536585365,4.680894308943089, 5.840447154471544, 7.04979674796748, 8.294715447154474,9.496951219512194,11.005081300813007,
12.20020325203252])

# Corresponding phase velocity values in m/s (estimated from image)
v_vals = np.array([2280.6825938566553,  2786.8941979522183, 3073.174061433447, 2877.815699658703, 2812.6962457337886,2874.1296928327647,2681.22866894198,2485.8703071672353,1991.9453924914676])

v_ad = 1749.6357634420633  # Adiabatic speed in m/s
v_inf = 2702.4785990040846  # High-frequency speed in m/s

# Molar masses:
M_NaCl = 0.05844   # kg/mol

# ASSUMPTIONS (example — change to real melt values)
rho_NaCl = 1503.99  # kg/m3

vL_NaCl = 1662.38   # m/s (longitudinal)
vT_NaCl = 831.19   # m/s (transverse)
eta_NaCl = 0.001032139   # Pa·s (viscosity)

na = debye_and_frenkel(M_NaCl, rho_NaCl, vL_NaCl, vT_NaCl, eta_NaCl)
na['f_D_THz'] = 6.41768
# na['f_F_THz']

print("NaCl estimates:", na)

plot_dispersion_with_lambda_axis(Q_vals, v_vals, v_ad, v_inf, omega_D=na['f_D_THz'], omega_F=na['f_F_THz'], data_type="velocity", label="NaCl")




###### KCl, 1070K - Demmel, 2008 ######
# energy transfer data (meV)
E_vals = np.array([3.708008505, 6.180014174, 8.15308292, 10.03543586, 12.03118356, 12.99503898, 13.33522325, 14.09496811, 14.38979447])

# Q values in nm^-1 (from image)
Q_vals = np.array([2.563995838, 3.995837669, 5.502601457, 6.934443288, 8.499479709, 9.898022893, 11.39646202, 12.89490114, 14.4516129])

# === 3. SOUND VELOCITY REFERENCE LINES ===
v_ad = 1600  # Adiabatic speed in m/s
v_inf = 2350  # High-frequency speed in m/s

# Molar masses:
M_KCl  = 0.07455   # kg/mol

# ASSUMPTIONS (example — change to real melt values)
rho_KCl  = 1511.98  # kg/m3

vL_KCl  = 1575.37   # m/s (longitudinal)
vT_KCl  = 787.68   # m/s (transverse)
eta_KCl = 0.001107867   # Pa·s (viscosity)

kc = debye_and_frenkel(M_KCl, rho_KCl, vL_KCl, vT_KCl, eta_KCl)
kc['f_D_THz'] = 4.959115
# kc['f_F_THz']

print("KCl  estimates:", kc)

plot_dispersion_with_lambda_axis(Q_vals, E_vals, v_ad, v_inf, omega_D=kc['f_D_THz'], omega_F=kc['f_F_THz'], data_type="energy", label="KCl")
