import numpy as np
import matplotlib.pyplot as plt

# ==========================
# Parameters and Initial Setup
# ==========================

# Physical constants
e = 1.602e-19         # elementary charge [C]
epsilon_0 = 8.854e-12 # vacuum permittivity [F/m]
k_B = 1.381e-23       # Boltzmann constant [J/K]
k_e = 1 / (4 * np.pi * epsilon_0)  # [N·m²/C²]
amu_to_kg = 1.66e-27  # [kg]

# Debye length (screening length) in meters
lambda_D = 4.58196E-11  # e.g., 2.5 Å Debye length

# Simulation parameters
Temp = 1000    # K, Molten salt temperature
N = 30         # Number of particles per chain
box_length = 100e-10  # Simulation box length in meters (100Å)
dt = 1e-15     # Time step (1 fs)
T = 10e-12    # Total simulation time (10 ps)
steps = int(T / dt)
output_freq = 100  # Output data every N steps

# Thermal bath parameters
thermostat_temp = Temp  # Temperature in K
thermostat_gamma = 0.1  # Coupling strength (1/ps)


# Sample setup: user-defined for each chain
# Masses (alternate m1, m2 for Chain A; uniform for Chain B)
mass_c = 39  # Mass of K+ (amu)
mass_a = 35  # Mass of Cl- (amu)
# Ensure Chain A alternates between cation and anion masses
mass_chain_A = np.array([(mass_c if i % 2 == 0 else mass_a) for i in range(N)])
# Chain B has uniform mass (all cations)
mass_chain_B = np.ones(N) * mass_c
radius_c = 1.38E-10
radius_a = 1.81E-10
r_sum_aa = radius_a + radius_a
r_sum_ac = radius_a + radius_c
r_sum_cc = radius_c + radius_c

# Charges in Coulombs
# Ensure charges alternate for Chain A
charges_A = np.array([(+1 if i % 2 == 0 else -1) for i in range(N)]) * e  # Chain A charges
charges_B = np.ones(N) * e  # Chain B (all cations)

# Rest separations (in Angstroms or normalized)
sep_chain_A = np.ones(N) * 2.96*1e-10
sep_chain_B = np.ones(N) * 4.38*1e-10

# Initial positions with periodic boundary conditions
def init_positions(sep, n):
    pos = np.cumsum(sep)
    # Center the chain in the box
    pos -= pos[n//2]
    pos += box_length/2
    return pos

x_A = init_positions(sep_chain_A, N)
x_B = init_positions(sep_chain_B, N)

v_A = np.zeros(N)
v_B = np.zeros(N)

# Impulse: initial velocity to first mass
v_thermal = np.sqrt(3*k_B*Temp/(np.mean([mass_c,mass_a])*amu_to_kg))
v_A[0] = v_thermal*0.01
v_B[0] = v_thermal*0.01

# ==========================
# Data Tracking
# ==========================

# Initialize arrays for tracking
ke_profile_A = []
ke_profile_B = []
time_points = []

# Reduced history to save memory
x_hist_A = np.zeros((steps//output_freq + 1, N))
x_hist_B = np.zeros((steps//output_freq + 1, N))
hist_idx = 0
prev_x_A = x_A.copy()
prev_x_B = x_B.copy()
total_distance_A = np.zeros(N)
total_distance_B = np.zeros(N)

ke_max_A = np.zeros(N)
ke_max_B = np.zeros(N)
pe_final_A = np.zeros(N)
pe_final_B = np.zeros(N)

# ==========================
# Simulation Loop
# ==========================

for t in range(steps):
    # Record positions at reduced frequency
    if t % output_freq == 0:
        x_hist_A[hist_idx] = x_A
        x_hist_B[hist_idx] = x_B
        hist_idx += 1

    # Yukawa-based acceleration for Chain A with distance cutoff
    a_A = np.zeros(N)
    min_distance_A = r_sum_ac  # Minimum distance cutoff
    for i in range(N):
        force = 0.0
        if i > 0:
            r = max(abs(x_A[i] - x_A[i-1]), min_distance_A  )  # Apply distance cutoff
            direction = np.sign(x_A[i] - x_A[i-1])  # Direction of force
            q1, q2 = charges_A[i], charges_A[i-1]
            # Calculate force magnitude with safe exponent
            exponent = min(r / lambda_D, 100)  # Cap the exponent to prevent overflow
            F = k_e * q1 * q2 / (r**2) * (1 + r / lambda_D) * np.exp(-exponent)
            force += direction * F
            
        if i < N - 1:
            r = max(abs(x_A[i+1] - x_A[i]), min_distance_A)  # Apply distance cutoff
            direction = np.sign(x_A[i+1] - x_A[i])  # Direction of force
            q1, q2 = charges_A[i], charges_A[i+1]
            # Calculate force magnitude with safe exponent
            exponent = min(r / lambda_D, 100)  # Cap the exponent to prevent overflow
            F = k_e * q1 * q2 / (r**2) * (1 + r / lambda_D) * np.exp(-exponent)
            force += direction * F
            
        a_A[i] = force / (mass_chain_A[i] * amu_to_kg)
    
    # Update velocities and positions using velocity Verlet with PBC
    v_half = v_A + 0.5 * a_A * dt
    x_A = (x_A + v_half * dt) % box_length  # Apply periodic boundary conditions
    
    # Apply thermostat (Nose-Hoover or Langevin)
    if t % 100 == 0:  # Apply thermostat every 100 steps
        ke = 0.5 * mass_chain_A * v_half**2
        current_temp = 2 * np.mean(ke) / k_B
        scale = np.sqrt(1 + thermostat_gamma * dt * (thermostat_temp/current_temp - 1))
        v_A = v_half * scale + 0.5 * a_A * dt
    else:
        v_A = v_half + 0.5 * a_A * dt
    
    # Update total distance traveled for Chain A
    total_distance_A += np.abs(x_A - prev_x_A)
    prev_x_A = x_A.copy()

    # Yukawa-based acceleration for Chain B with distance cutoff
    a_B = np.zeros(N)
    min_distance_B = r_sum_cc  # Minimum distance cutoff
    for i in range(N):
        force = 0.0
        if i > 0:
            r = max(abs(x_B[i] - x_B[i-1]), min_distance_B)  # Apply distance cutoff
            direction = np.sign(x_B[i] - x_B[i-1])  # Direction of force
            q1, q2 = charges_B[i], charges_B[i-1]
            # Calculate force magnitude with safe exponent
            exponent = min(r / lambda_D, 100)  # Cap the exponent to prevent overflow
            F = k_e * q1 * q2 / (r**2) * (1 + r / lambda_D) * np.exp(-exponent)
            force += direction * F
            
        if i < N - 1:
            r = max(abs(x_B[i+1] - x_B[i]), min_distance_B)  # Apply distance cutoff
            direction = np.sign(x_B[i+1] - x_B[i])  # Direction of force
            q1, q2 = charges_B[i], charges_B[i+1]
            # Calculate force magnitude with safe exponent
            exponent = min(r / lambda_D, 100)  # Cap the exponent to prevent overflow
            F = k_e * q1 * q2 / (r**2) * (1 + r / lambda_D) * np.exp(-exponent)
            force += direction * F
            
        a_B[i] = force / (mass_chain_B[i] * amu_to_kg)
    
    # Update velocities and positions using velocity Verlet with PBC
    v_half = v_B + 0.5 * a_B * dt
    x_B = (x_B + v_half * dt) % box_length  # Apply periodic boundary conditions
    
    # Apply thermostat
    if t % 100 == 0:  # Apply thermostat every 100 steps
        ke = 0.5 * mass_chain_B * v_half**2
        current_temp = 2 * np.mean(ke) / k_B
        scale = np.sqrt(1 + thermostat_gamma * dt * (thermostat_temp/current_temp - 1))
        v_B = v_half * scale + 0.5 * a_B * dt
    else:
        v_B = v_half + 0.5 * a_B * dt
    
    # Update total distance traveled for Chain B
    total_distance_B += np.abs(x_B - prev_x_B)
    prev_x_B = x_B.copy()

    # Update kinetic energy and track time evolution
    current_ke_A = 0.5 * mass_chain_A * v_A**2
    current_ke_B = 0.5 * mass_chain_B * v_B**2
    ke_max_A = np.maximum(ke_max_A, current_ke_A)
    ke_max_B = np.maximum(ke_max_B, current_ke_B)
    
    # Track energy propagation
    if t % output_freq == 0:
        # Record energy profile along the chain
        ke_profile_A.append(current_ke_A.copy())
        ke_profile_B.append(current_ke_B.copy())
        time_points.append(t * dt)

# ==========================
# Energy Analysis
# ==========================

# Calculate energy propagation speed
def analyze_energy_propagation(ke_profiles, time_points):
    """Analyze how energy propagates through the chain over time."""
    ke_profiles = np.array(ke_profiles)
    propagation_front = []
    
    for i, profile in enumerate(ke_profiles):
        # Find the front of the energy wave
        threshold = 0.1 * np.max(profile)  # 10% of max energy as threshold
        front_pos = np.argmax(profile > threshold)
        propagation_front.append((time_points[i], front_pos))
    
    return np.array(propagation_front)

# Analyze energy propagation
propagation_A = analyze_energy_propagation(ke_profile_A, time_points)
propagation_B = analyze_energy_propagation(ke_profile_B, time_points)

# Calculate propagation speeds
if len(propagation_A) > 1:
    speed_A = np.polyfit(propagation_A[:,0], propagation_A[:,1], 1)[0]  # particles/ps
    print(f"Energy propagation speed - Chain A: {speed_A:.2f} particles/ps")

if len(propagation_B) > 1:
    speed_B = np.polyfit(propagation_B[:,0], propagation_B[:,1], 1)[0]  # particles/ps
    print(f"Energy propagation speed - Chain B: {speed_B:.2f} particles/ps")

for i in range(1, N):
    r = abs(x_A[i] - x_A[i-1])
    pe_final_A[i] = k_e * charges_A[i] * charges_A[i-1] / r * np.exp(-r / lambda_D)

    r = abs(x_B[i] - x_B[i-1])
    pe_final_B[i] = k_e * charges_B[i] * charges_B[i-1] / r * np.exp(-r / lambda_D)

# ==========================
# Result Summaries
# ==========================

# Convert kinetic and potential energies to eV
ke_max_A /= e
ke_max_B /= e
pe_final_A /= e
pe_final_B /= e

# Convert distances to Å for plotting
x_hist_A_angstrom = x_hist_A * 1e10
x_hist_B_angstrom = x_hist_B * 1e10
total_distance_A_angstrom = total_distance_A * 1e10
total_distance_B_angstrom = total_distance_B * 1e10

# Calculate metrics for first 10 particles
particle_indices = np.arange(10)
displacement_A = x_hist_A_angstrom[-1, :10] - x_hist_A_angstrom[0, :10]  # Final - Initial position
displacement_B = x_hist_B_angstrom[-1, :10] - x_hist_B_angstrom[0, :10]  # Final - Initial position

# Calculate energy transfer efficiency (fraction of initial KE transferred to each particle)
initial_ke_A = 0.5 * mass_chain_A[0] * v_thermal**2
initial_ke_B = 0.5 * mass_chain_B[0] * v_thermal**2
eff_energy_A = (0.5 * mass_chain_A * v_A**2) / initial_ke_A
eff_energy_B = (0.5 * mass_chain_B * v_B**2) / initial_ke_B

# Plotting function
def plot_comparison_bar(data_A, data_B, ylabel, title):
    indices = np.arange(10)
    width = 0.35

    fig, ax = plt.subplots(figsize=(10, 5))
    bar1 = ax.bar(indices - width/2, data_A, width, label='Chain A (Cation-Anion)', color='skyblue')
    bar2 = ax.bar(indices + width/2, data_B, width, label='Chain B (Cation-Cation)', color='salmon')

    ax.set_xlabel('Mass Index')
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xticks(indices)
    ax.set_xticklabels([str(i) for i in range(10)])
    ax.legend()
    ax.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    return fig

# Create a figure with all metrics
plt.figure(figsize=(15, 10))

# 1. Potential Energy
plt.subplot(2, 2, 1)
width = 0.35
plt.bar(particle_indices - width/2, pe_final_A[:10], width, label='Chain A (Cation-Anion)', color='skyblue')
plt.bar(particle_indices + width/2, pe_final_B[:10], width, label='Chain B (Cation-Cation)', color='salmon')
plt.xlabel('Particle Index')
plt.ylabel('Potential Energy (eV)')
plt.title('Final Potential Energy')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.5)

# 2. Total distance traveled
plt.subplot(2, 2, 2)
plt.bar(particle_indices - width/2, total_distance_A_angstrom[:10], width, label='Chain A (Cation-Anion)', color='skyblue')
plt.bar(particle_indices + width/2, total_distance_B_angstrom[:10], width, label='Chain B (Cation-Cation)', color='salmon')
plt.xlabel('Particle Index')
plt.ylabel('Total Distance Traveled (Å)')
plt.title('Total Distance Traveled')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.5)

# 3. Max Kinetic Energy
plt.subplot(2, 2, 3)
plt.bar(particle_indices - width/2, ke_max_A[:10], width, label='Chain A (Cation-Anion)', color='skyblue')
plt.bar(particle_indices + width/2, ke_max_B[:10], width, label='Chain B (Cation-Cation)', color='salmon')
plt.xlabel('Particle Index')
plt.ylabel('Max Kinetic Energy (eV)')
plt.title('Maximum Kinetic Energy')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.5)

# 4. Energy Transfer Efficiency
plt.subplot(2, 2, 4)
plt.bar(particle_indices - width/2, eff_energy_A[:10], width, label='Chain A (Cation-Anion)', color='skyblue')
plt.bar(particle_indices + width/2, eff_energy_B[:10], width, label='Chain B (Cation-Cation)', color='salmon')
plt.xlabel('Particle Index')
plt.ylabel('Energy Transfer Efficiency')
plt.title('Energy Transfer Efficiency (Fraction of Initial KE)')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.5)

plt.tight_layout()

plt.show()

