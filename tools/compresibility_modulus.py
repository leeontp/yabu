import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

# --- Constants ---
k_B = 1.380649e-23  # Boltzmann constant (J/K)

# --- User Inputs ---
print("\n=== Membrane Area per Lipid & Compressibility Analysis ===")
gro_file = input("Enter structure file name (.gro or .pdb): ").strip()
xtc_file = input("Enter trajectory file name (.xtc or .trr): ").strip()
lipid_resnames = input("Enter lipid resnames (space-separated, e.g., POPC POPE CHOL): ").strip().split()
n_residues = int(input("Enter number of residues in ONE leaflet: "))
t_start_ns = float(input("Enter starting time for analysis (ns): "))
T_kelvin = float(input("Enter simulation temperature (K): "))
output_name = input("Enter base name for output files (CSV and PNG): ").strip()

# Convert starting time to ps
t_start_ps = t_start_ns * 1000

# --- Load trajectory ---
u = mda.Universe(gro_file, xtc_file)

# --- Initialize arrays ---
times = []     # time in ns
areas = []     # area per lipid in Å²

# --- Process trajectory ---
for ts in u.trajectory:
    if ts.time < t_start_ps:
        continue
    
    # Get box lengths
    Lx, Ly, Lz, _, _, _ = ts.dimensions
    
    # Area per lipid (Å² / lipid)
    area = (Lx * Ly) / n_residues
    
    times.append(ts.time / 1000.0)  # convert ps → ns
    areas.append(area)

# Convert to numpy arrays
areas = np.array(areas)
times = np.array(times)

# --- Convert to SI units ---
areas_m2 = areas * 1e-20   # Å² → m²

# --- Calculate statistics ---
A_mean = np.mean(areas_m2)
A_var  = np.var(areas_m2, ddof=1)

# Compressibility modulus in N/m
K_A = (k_B * T_kelvin * A_mean) / A_var

# --- Save results to CSV ---
csv_file = f"{output_name}.csv"
with open(csv_file, "w") as f:
    f.write("Time (ns),Area per lipid (Å^2)\n")
    for t, a in zip(times, areas):
        f.write(f"{t:.4f},{a:.4f}\n")
    f.write(f"\nMean area per lipid (Å^2),{np.mean(areas):.4f}\n")
    f.write(f"Variance (Å^4),{np.var(areas, ddof=1):.4f}\n")
    f.write(f"Compressibility modulus K_A (N/m),{K_A:.4e}\n")

# --- Plot area per lipid ---
png_file = f"{output_name}.png"
plt.figure(figsize=(8,5))
plt.plot(times, areas, color="teal", label="Area per lipid")
plt.axhline(np.mean(areas), color="red", linestyle="--", label="Mean")
plt.xlabel("Time (ns)")
plt.ylabel("Area per lipid (Å²)")
plt.title("Membrane Area per Lipid")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(png_file, dpi=300)

# --- Print results ---
print(f"\nData saved to '{csv_file}' and plot saved to '{png_file}'!")
print(f"Mean area per lipid: {np.mean(areas):.2f} Å²")
print(f"Compressibility modulus K_A: {K_A:.4e} N/m")