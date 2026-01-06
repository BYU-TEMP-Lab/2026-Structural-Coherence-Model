#!/usr/bin/env python3
"""
Molten Salt Energy Transfer Simulation Framework
===============================================

This module implements molecular dynamics simulations and kinetic energy transfer
models for molten salt systems based on theoretical frameworks from literature.

Key Features:
- Molecular dynamics simulation with Lennard-Jones potentials
- Kinetic energy transfer pathway analysis
- Energy distribution calculations
- Visualization of energy transfer mechanisms

Author: SciSpace Research Agent
Date: 2025-08-05
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.integrate import odeint
from scipy.spatial.distance import cdist
import pandas as pd
from typing import Tuple, List, Dict, Optional
import warnings
warnings.filterwarnings('ignore')

class MoltenSaltEnergyTransfer:
    """
    Comprehensive framework for modeling energy transfer in molten salts
    """
    
    def __init__(self, n_particles: int = 1000, temperature: float = 1000.0, 
                 density: float = 2.5, box_size: float = 10.0):
        """
        Initialize the molten salt system
        
        Parameters:
        -----------
        n_particles : int
            Number of particles in the simulation
        temperature : float
            System temperature in Kelvin
        density : float
            System density in g/cm³
        box_size : float
            Simulation box size in Angstroms
        """
        self.n_particles = n_particles
        self.temperature = temperature
        self.density = density
        self.box_size = box_size
        self.kb = 8.617e-5  # Boltzmann constant in eV/K
        
        # Initialize particle properties
        self.positions = None
        self.velocities = None
        self.forces = None
        self.masses = None
        self.charges = None
        
        # Energy tracking
        self.kinetic_energies = []
        self.potential_energies = []
        self.total_energies = []
        self.energy_transfer_rates = []
        
        # Simulation parameters
        self.dt = 1e-15  # Time step in seconds
        self.epsilon = 0.1  # LJ potential well depth (eV)
        self.sigma = 2.5    # LJ potential characteristic length (Angstrom)
        
    def initialize_system(self) -> None:
        """Initialize particle positions, velocities, and properties"""
        
        # Random initial positions
        self.positions = np.random.uniform(0, self.box_size, (self.n_particles, 3))
        
        # Maxwell-Boltzmann velocity distribution
        v_thermal = np.sqrt(self.kb * self.temperature / 1.0)  # Assuming unit mass
        self.velocities = np.random.normal(0, v_thermal, (self.n_particles, 3))
        
        # Initialize masses (simplified: alternating cation/anion masses)
        self.masses = np.ones(self.n_particles)
        self.masses[::2] = 23.0  # Na+ mass (amu)
        self.masses[1::2] = 35.5  # Cl- mass (amu)
        
        # Initialize charges
        self.charges = np.ones(self.n_particles)
        self.charges[::2] = 1.0   # Cation charge
        self.charges[1::2] = -1.0  # Anion charge
        
        # Initialize forces
        self.forces = np.zeros((self.n_particles, 3))
        
    def lennard_jones_potential(self, r: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate Lennard-Jones potential and forces
        
        Parameters:
        -----------
        r : np.ndarray
            Distance matrix between particles
            
        Returns:
        --------
        potential : np.ndarray
            Potential energy matrix
        force_magnitude : np.ndarray
            Force magnitude matrix
        """
        # Avoid division by zero
        r_safe = np.where(r < 0.1, 0.1, r)
        
        # LJ potential: 4ε[(σ/r)^12 - (σ/r)^6]
        sigma_r = self.sigma / r_safe
        potential = 4 * self.epsilon * (sigma_r**12 - sigma_r**6)
        
        # Force magnitude: -dU/dr
        force_magnitude = 24 * self.epsilon * (2 * sigma_r**12 - sigma_r**6) / r_safe
        
        return potential, force_magnitude
        
    def coulomb_potential(self, r: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate Coulomb potential and forces
        
        Parameters:
        -----------
        r : np.ndarray
            Distance matrix between particles
            
        Returns:
        --------
        potential : np.ndarray
            Coulomb potential energy matrix
        force_magnitude : np.ndarray
            Coulomb force magnitude matrix
        """
        ke = 14.4  # Coulomb constant in eV·Å/e²
        
        # Avoid division by zero
        r_safe = np.where(r < 0.1, 0.1, r)
        
        # Coulomb potential: ke * qi * qj / r
        charge_products = np.outer(self.charges, self.charges)
        potential = ke * charge_products / r_safe
        
        # Force magnitude: ke * qi * qj / r²
        force_magnitude = ke * charge_products / (r_safe**2)
        
        return potential, force_magnitude
        
    def calculate_forces(self) -> None:
        """Calculate total forces on all particles"""
        
        # Reset forces
        self.forces.fill(0.0)
        
        # Calculate distance matrix
        distances = cdist(self.positions, self.positions)
        
        # Calculate LJ and Coulomb potentials/forces
        lj_potential, lj_force = self.lennard_jones_potential(distances)
        coulomb_potential, coulomb_force = self.coulomb_potential(distances)
        
        # Total force magnitude
        total_force_magnitude = lj_force + coulomb_force
        
        # Calculate force vectors
        for i in range(self.n_particles):
            for j in range(i + 1, self.n_particles):
                if distances[i, j] > 0.1:  # Avoid self-interaction
                    # Unit vector from i to j
                    r_vec = self.positions[j] - self.positions[i]
                    r_unit = r_vec / distances[i, j]
                    
                    # Force vector
                    force_vec = total_force_magnitude[i, j] * r_unit
                    
                    # Newton's third law
                    self.forces[i] += force_vec
                    self.forces[j] -= force_vec
                    
    def velocity_verlet_step(self) -> None:
        """Perform one velocity-Verlet integration step"""
        
        # Update positions
        self.positions += self.velocities * self.dt + 0.5 * self.forces / self.masses[:, np.newaxis] * self.dt**2
        
        # Apply periodic boundary conditions
        self.positions = self.positions % self.box_size
        
        # Store old forces
        old_forces = self.forces.copy()
        
        # Calculate new forces
        self.calculate_forces()
        
        # Update velocities
        self.velocities += 0.5 * (old_forces + self.forces) / self.masses[:, np.newaxis] * self.dt
        
    def calculate_kinetic_energy(self) -> float:
        """Calculate total kinetic energy of the system"""
        ke = 0.5 * np.sum(self.masses[:, np.newaxis] * self.velocities**2)
        return ke
        
    def calculate_potential_energy(self) -> float:
        """Calculate total potential energy of the system"""
        distances = cdist(self.positions, self.positions)
        lj_potential, _ = self.lennard_jones_potential(distances)
        coulomb_potential, _ = self.coulomb_potential(distances)
        
        # Sum upper triangle to avoid double counting
        total_potential = 0.5 * np.sum(np.triu(lj_potential + coulomb_potential, k=1))
        return total_potential
        
    def calculate_energy_transfer_rate(self, window_size: int = 10) -> float:
        """
        Calculate energy transfer rate based on kinetic energy fluctuations
        
        Parameters:
        -----------
        window_size : int
            Window size for calculating transfer rate
            
        Returns:
        --------
        transfer_rate : float
            Energy transfer rate in eV/ps
        """
        if len(self.kinetic_energies) < window_size:
            return 0.0
            
        recent_ke = np.array(self.kinetic_energies[-window_size:])
        transfer_rate = np.std(recent_ke) / (self.dt * 1e12)  # Convert to eV/ps
        return transfer_rate
        
    def run_simulation(self, n_steps: int = 10000, save_interval: int = 10) -> None:
        """
        Run molecular dynamics simulation
        
        Parameters:
        -----------
        n_steps : int
            Number of simulation steps
        save_interval : int
            Interval for saving energy data
        """
        
        print(f"Running MD simulation with {self.n_particles} particles for {n_steps} steps...")
        
        for step in range(n_steps):
            # Perform MD step
            self.velocity_verlet_step()
            
            # Save energy data
            if step % save_interval == 0:
                ke = self.calculate_kinetic_energy()
                pe = self.calculate_potential_energy()
                te = ke + pe
                transfer_rate = self.calculate_energy_transfer_rate()
                
                self.kinetic_energies.append(ke)
                self.potential_energies.append(pe)
                self.total_energies.append(te)
                self.energy_transfer_rates.append(transfer_rate)
                
            # Progress indicator
            if step % 1000 == 0:
                print(f"Step {step}/{n_steps} completed")
                
        print("Simulation completed!")
        
    def analyze_energy_pathways(self) -> Dict[str, np.ndarray]:
        """
        Analyze energy transfer pathways in the system
        
        Returns:
        --------
        pathways : Dict[str, np.ndarray]
            Dictionary containing energy pathway analysis
        """
        
        # Calculate velocity autocorrelation function
        velocities_flat = self.velocities.flatten()
        autocorr = np.correlate(velocities_flat, velocities_flat, mode='full')
        autocorr = autocorr[autocorr.size // 2:]
        autocorr = autocorr / autocorr[0]  # Normalize
        
        # Calculate radial distribution function
        distances = cdist(self.positions, self.positions)
        hist, bin_edges = np.histogram(distances.flatten(), bins=50, range=(0, self.box_size/2))
        r_values = (bin_edges[:-1] + bin_edges[1:]) / 2
        
        # Energy flux analysis
        energy_flux = np.gradient(self.kinetic_energies)
        
        pathways = {
            'velocity_autocorr': autocorr[:100],  # First 100 points
            'radial_distribution': hist,
            'r_values': r_values,
            'energy_flux': np.array(energy_flux),
            'transfer_rates': np.array(self.energy_transfer_rates)
        }
        
        return pathways

def create_energy_transfer_visualizations(sim_results: MoltenSaltEnergyTransfer, 
                                        pathways: Dict[str, np.ndarray],
                                        save_path: str = '/home/sandbox/') -> None:
    """
    Create comprehensive visualizations of energy transfer results
    
    Parameters:
    -----------
    sim_results : MoltenSaltEnergyTransfer
        Simulation results object
    pathways : Dict[str, np.ndarray]
        Energy pathway analysis results
    save_path : str
        Path to save visualization files
    """
    
    # Set up the plotting style
    plt.style.use('seaborn-v0_8')
    sns.set_palette("husl")
    
    # Create figure with subplots
    fig = plt.figure(figsize=(20, 16))
    
    # 1. Energy evolution over time
    ax1 = plt.subplot(3, 3, 1)
    time_steps = np.arange(len(sim_results.kinetic_energies)) * 10  # Assuming save_interval=10
    plt.plot(time_steps, sim_results.kinetic_energies, label='Kinetic Energy', alpha=0.8)
    plt.plot(time_steps, sim_results.potential_energies, label='Potential Energy', alpha=0.8)
    plt.plot(time_steps, sim_results.total_energies, label='Total Energy', alpha=0.8)
    plt.xlabel('Time Step')
    plt.ylabel('Energy (eV)')
    plt.title('Energy Evolution During Simulation')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # 2. Energy transfer rates
    ax2 = plt.subplot(3, 3, 2)
    plt.plot(time_steps, sim_results.energy_transfer_rates, color='red', alpha=0.7)
    plt.xlabel('Time Step')
    plt.ylabel('Transfer Rate (eV/ps)')
    plt.title('Energy Transfer Rate Evolution')
    plt.grid(True, alpha=0.3)
    
    # 3. Velocity autocorrelation function
    ax3 = plt.subplot(3, 3, 3)
    plt.plot(pathways['velocity_autocorr'], color='green', linewidth=2)
    plt.xlabel('Time Lag')
    plt.ylabel('Autocorrelation')
    plt.title('Velocity Autocorrelation Function')
    plt.grid(True, alpha=0.3)
    
    # 4. Radial distribution function
    ax4 = plt.subplot(3, 3, 4)
    plt.plot(pathways['r_values'], pathways['radial_distribution'], color='blue', linewidth=2)
    plt.xlabel('Distance (Å)')
    plt.ylabel('g(r)')
    plt.title('Radial Distribution Function')
    plt.grid(True, alpha=0.3)
    
    # 5. Energy flux analysis
    ax5 = plt.subplot(3, 3, 5)
    plt.plot(time_steps, pathways['energy_flux'], color='orange', alpha=0.7)
    plt.xlabel('Time Step')
    plt.ylabel('Energy Flux (eV/step)')
    plt.title('Energy Flux Analysis')
    plt.grid(True, alpha=0.3)
    
    # 6. Kinetic energy distribution
    ax6 = plt.subplot(3, 3, 6)
    plt.hist(sim_results.kinetic_energies, bins=30, alpha=0.7, color='purple', density=True)
    plt.xlabel('Kinetic Energy (eV)')
    plt.ylabel('Probability Density')
    plt.title('Kinetic Energy Distribution')
    plt.grid(True, alpha=0.3)
    
    # 7. Phase space visualization (2D projection)
    ax7 = plt.subplot(3, 3, 7)
    ke_array = np.array(sim_results.kinetic_energies)
    pe_array = np.array(sim_results.potential_energies)
    plt.scatter(ke_array, pe_array, alpha=0.5, s=10, c=time_steps, cmap='viridis')
    plt.xlabel('Kinetic Energy (eV)')
    plt.ylabel('Potential Energy (eV)')
    plt.title('Phase Space Evolution')
    plt.colorbar(label='Time Step')
    plt.grid(True, alpha=0.3)
    
    # 8. Energy transfer rate histogram
    ax8 = plt.subplot(3, 3, 8)
    plt.hist(pathways['transfer_rates'], bins=25, alpha=0.7, color='red', density=True)
    plt.xlabel('Transfer Rate (eV/ps)')
    plt.ylabel('Probability Density')
    plt.title('Energy Transfer Rate Distribution')
    plt.grid(True, alpha=0.3)
    
    # 9. 3D particle positions (final configuration)
    ax9 = plt.subplot(3, 3, 9, projection='3d')
    cations = sim_results.positions[::2]  # Every other particle (cations)
    anions = sim_results.positions[1::2]  # Remaining particles (anions)
    
    ax9.scatter(cations[:, 0], cations[:, 1], cations[:, 2], 
                c='red', s=20, alpha=0.6, label='Cations')
    ax9.scatter(anions[:, 0], anions[:, 1], anions[:, 2], 
                c='blue', s=20, alpha=0.6, label='Anions')
    ax9.set_xlabel('X (Å)')
    ax9.set_ylabel('Y (Å)')
    ax9.set_zlabel('Z (Å)')
    ax9.set_title('Final Particle Configuration')
    ax9.legend()
    
    plt.tight_layout()
    plt.savefig(f'{save_path}molten_salt_energy_transfer_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create additional detailed plots
    create_detailed_energy_analysis(sim_results, pathways, save_path)

def create_detailed_energy_analysis(sim_results: MoltenSaltEnergyTransfer, 
                                   pathways: Dict[str, np.ndarray],
                                   save_path: str) -> None:
    """Create detailed energy analysis plots"""
    
    # Energy correlation analysis
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
    
    # 1. Energy autocorrelation
    ke_autocorr = np.correlate(sim_results.kinetic_energies, sim_results.kinetic_energies, mode='full')
    ke_autocorr = ke_autocorr[ke_autocorr.size // 2:]
    ke_autocorr = ke_autocorr / ke_autocorr[0]
    
    ax1.plot(ke_autocorr[:50], linewidth=2, color='blue')
    ax1.set_xlabel('Time Lag')
    ax1.set_ylabel('KE Autocorrelation')
    ax1.set_title('Kinetic Energy Autocorrelation')
    ax1.grid(True, alpha=0.3)
    
    # 2. Energy transfer efficiency
    efficiency = np.abs(pathways['energy_flux']) / np.array(sim_results.kinetic_energies[1:])
    ax2.plot(efficiency, color='green', alpha=0.7)
    ax2.set_xlabel('Time Step')
    ax2.set_ylabel('Transfer Efficiency')
    ax2.set_title('Energy Transfer Efficiency')
    ax2.grid(True, alpha=0.3)
    
    # 3. Temperature evolution
    temperature = np.array(sim_results.kinetic_energies) / (1.5 * sim_results.kb * sim_results.n_particles)
    ax3.plot(temperature, color='red', linewidth=2)
    ax3.axhline(y=sim_results.temperature, color='black', linestyle='--', label='Target T')
    ax3.set_xlabel('Time Step')
    ax3.set_ylabel('Temperature (K)')
    ax3.set_title('System Temperature Evolution')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # 4. Energy variance analysis
    window_size = 50
    energy_variance = []
    for i in range(window_size, len(sim_results.total_energies)):
        window_data = sim_results.total_energies[i-window_size:i]
        energy_variance.append(np.var(window_data))
    
    ax4.plot(energy_variance, color='purple', alpha=0.8)
    ax4.set_xlabel('Time Step')
    ax4.set_ylabel('Energy Variance')
    ax4.set_title('Total Energy Variance (Rolling Window)')
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{save_path}detailed_energy_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()

def main():
    """Main execution function for molten salt energy transfer simulation"""
    
    print("=== Molten Salt Energy Transfer Simulation ===")
    print("Initializing simulation parameters...")
    
    # Initialize simulation
    sim = MoltenSaltEnergyTransfer(n_particles=500, temperature=1000.0, 
                                  density=2.5, box_size=15.0)
    
    # Initialize system
    sim.initialize_system()
    print(f"System initialized with {sim.n_particles} particles")
    print(f"Temperature: {sim.temperature} K")
    print(f"Density: {sim.density} g/cm³")
    print(f"Box size: {sim.box_size} Å")
    
    # Run simulation
    sim.run_simulation(n_steps=5000, save_interval=10)
    
    # Analyze energy pathways
    print("Analyzing energy transfer pathways...")
    pathways = sim.analyze_energy_pathways()
    
    # Create visualizations
    print("Creating visualizations...")
    create_energy_transfer_visualizations(sim, pathways)
    
    # Save simulation data
    print("Saving simulation data...")
    simulation_data = {
        'kinetic_energies': sim.kinetic_energies,
        'potential_energies': sim.potential_energies,
        'total_energies': sim.total_energies,
        'energy_transfer_rates': sim.energy_transfer_rates,
        'final_positions': sim.positions.tolist(),
        'final_velocities': sim.velocities.tolist()
    }
    
    import json
    with open('/home/sandbox/simulation_data.json', 'w') as f:
        json.dump(simulation_data, f, indent=2)
    
    # Generate summary statistics
    print("\n=== Simulation Results Summary ===")
    print(f"Average kinetic energy: {np.mean(sim.kinetic_energies):.3f} eV")
    print(f"Average potential energy: {np.mean(sim.potential_energies):.3f} eV")
    print(f"Average total energy: {np.mean(sim.total_energies):.3f} eV")
    print(f"Average energy transfer rate: {np.mean(sim.energy_transfer_rates):.3f} eV/ps")
    print(f"Energy conservation (std/mean): {np.std(sim.total_energies)/np.mean(sim.total_energies):.6f}")
    
    print("\nSimulation completed successfully!")
    print("Results saved to /home/sandbox/")

if __name__ == "__main__":
    main()