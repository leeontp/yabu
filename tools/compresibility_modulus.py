import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

# --- Constants ---
k_B = 1.380649e-23  # Boltzmann constant in J/K

# --- User Inputs ---
print("\n=== Membrane Compressibility Analysis ===")
gro_file = input("Enter structure file name (.gro or .pdb): ").strip()
xtc_file = input("Enter trajectory file name (.xtc or .trr): ").strip()
lipid_resnames = input("Enter lipid resnames (space-separated, e.g., POPC POPE CHOL): ").strip().split()
n_residues = int(input("Enter number of residues in ONE leaflet: "))
t_start_ns = float(input("Enter starting time for analysis (ns): "))
T_kelvin = float(input("Enter simulation temperature (K): "))
output_name = input("Enter base name for output files (CSV and PNG): ").strip()

# Convert starting time to ps
t_start_ps = t_start_ns * 1000

# --- Load Trajectory ---
u = mda.Universe(gro_file, xtc_file)

# --- Select all atoms of the lipid residues ---
lipid_atoms = u.select_atoms(f"resname {' '.join(lipid_resnames)}")
z_center = lipid_atoms.positions[:, 2].mean()

# --- Initialize arrays ---
times = []      # Time in ns
areas = []      # Area per lipid (proxy) in Å²/residue

# --- Process Trajectory ---
for ts in u.trajectory:
    if ts.time < t_start_ps:
        continue
    
    upper_leaflet = lipid_atoms.positions[:, 2] > z_center
    lower_leaflet = lipid_atoms.positions[:, 2] < z_center
    
    if not upper_leaflet.any() or not lower_leaflet.any():
        continue
    
    upper_z = lipid_atoms.positions[upper_leaflet, 2].mean()
    lower_z = lipid_atoms.positions[lower_leaflet, 2].mean()
    
    width = upper_z - lower_z
    area = (width ** 2) / n_residues
    times.append(ts.time / 1000)  # ps -> ns
    areas.append(area)

# Convert to numpy arrays
areas = np.array(areas)
times = np.array(times)

# --- Calculate compressibility modulus K_A per frame ---
areas_m2 = areas * 1e-20
K_A_frames = k_B * T_kelvin * areas_m2 / np.var(areas_m2, ddof=1)

# --- Calculate mean and std ---
K_A_mean = np.mean(K_A_frames)
K_A_std = np.std(K_A_frames, ddof=1)

# --- Save to CSV ---
csv_file = f"{output_name}.csv"
with open(csv_file, "w") as f:
    f.write("Time (ns),K_A (J/m²)\n")
    for t, k in zip(times, K_A_frames):
        f.write(f"{t:.4f},{k:.4e}\n")
    f.write(f"Mean ± Std,{K_A_mean:.4e} ± {K_A_std:.4e}\n")

# --- Plot ---
png_file = f"{output_name}.png"
plt.figure(figsize=(8,5))
plt.plot(times, K_A_frames, color="teal", label="K_A per frame")
plt.axhline(K_A_mean, color="red", linestyle="--", label="Mean K_A")
plt.xlabel("Time (ns)")
plt.ylabel("Compressibility Modulus K_A (J/m²)")
plt.title("Membrane Compressibility Modulus")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(png_file, dpi=300)

# --- Print Results ---
print(f"\nData saved to '{csv_file}' and plot saved to '{png_file}'!")
print(f"Average compressibility modulus K_A: {K_A_mean:.4e} J/m²")
print(f"Standard deviation: {K_A_std:.4e} J/m²")