import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

# --- User Inputs ---
print("\n=== Lipid per Area Analysis ===")
gro_file = input("Enter structure file name (.gro or .pdb): ").strip()
xtc_file = input("Enter trajectory file name (.xtc or .trr): ").strip()
lipid_resnames = input("Enter lipid resnames (space-separated, e.g., POPC POPE CHOL): ").strip().split()
n_residues = int(input("Enter number of residues in ONE leaflet: "))
t_start_ns = float(input("Enter starting time for analysis (ns): "))
output_name = input("Enter base name for output files (CSV and PNG): ").strip()

# Convert starting time to ps
t_start_ps = t_start_ns * 1000

# --- Load trajectory ---
u = mda.Universe(gro_file, xtc_file)

# --- Initialize arrays ---
times = []          # time in ns
lipid_per_area = [] # lipid density (lipids / Å²)

# --- Process trajectory ---
for ts in u.trajectory:
    if ts.time < t_start_ps:
        continue
    
    # Get box lengths
    Lx, Ly, Lz, _, _, _ = ts.dimensions
    
    # Lipids per area (lipids / Å²)
    value = n_residues / (Lx * Ly)
    
    times.append(ts.time / 1000.0)  # ps → ns
    lipid_per_area.append(value)

# Convert to numpy arrays
lipid_per_area = np.array(lipid_per_area)
times = np.array(times)

# --- Statistics ---
mean_density = np.mean(lipid_per_area)
std_density = np.std(lipid_per_area, ddof=1)

# --- Save to CSV ---
csv_file = f"{output_name}.csv"
with open(csv_file, "w") as f:
    f.write("Time (ns),Lipids per area (lipids/Å^2)\n")
    for t, val in zip(times, lipid_per_area):
        f.write(f"{t:.4f},{val:.6f}\n")
    f.write(f"\nMean density (lipids/Å^2),{mean_density:.6f}\n")
    f.write(f"Std deviation (lipids/Å^2),{std_density:.6f}\n")
    f.write(f"Mean density (lipids/nm^2),{mean_density*100:.4f}\n")
    f.write(f"Std deviation (lipids/nm^2),{std_density*100:.4f}\n")

# --- Plot ---
png_file = f"{output_name}.png"
plt.figure(figsize=(8,5))
plt.plot(times, lipid_per_area*100, color="teal", label="Lipids per nm²")
plt.axhline(mean_density*100, color="red", linestyle="--", label="Mean")
plt.xlabel("Time (ns)")
plt.ylabel("Lipids per area (lipids/nm²)")
plt.title("Lipid Density in Membrane")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(png_file, dpi=300)

# --- Print results ---
print(f"\nData saved to '{csv_file}' and plot saved to '{png_file}'!")
print(f"Mean lipid density: {mean_density*100:.2f} lipids/nm²")
print(f"Standard deviation: {std_density*100:.2f} lipids/nm²")

