import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

# --- User Inputs ---
print("\n=== Membrane Width Analysis ===")
gro_file = input("Enter structure file name (.gro or .pdb): ").strip()
xtc_file = input("Enter trajectory file name (.xtc or .trr): ").strip()
lipid_resnames = input("Enter lipid resnames (space-separated, e.g., POPC POPE CHOL): ").strip().split()
n_residues = int(input("Enter number of residues in ONE leaflet: "))
t_start_ns = float(input("Enter starting time for analysis (ns): "))
output_name = input("Enter base name for output files (CSV and PNG): ").strip()

# Convert starting time to ps
t_start_ps = t_start_ns * 1000

# --- Load Trajectory ---
try:
    u = mda.Universe(gro_file, xtc_file)
except Exception as e:
    print(f"\nError loading files: {e}")
    exit()

# --- Select all atoms of the lipid residues ---
try:
    lipid_atoms = u.select_atoms(f"resname {' '.join(lipid_resnames)}")
    if lipid_atoms.n_atoms == 0:
        print(f"\nError: No atoms found with resnames {lipid_resnames}!")
        exit()
    
    z_center = lipid_atoms.positions[:, 2].mean()
    print(f"\nAutomatically calculated membrane center (Z-value): {z_center:.2f} Å")
    print(f"Using user-defined residues per leaflet: {n_residues}")
    print(f"Analysis will start from {t_start_ns:.2f} ns (frame times >= {t_start_ps:.1f} ps)")

except Exception as e:
    print(f"\nError selecting lipid atoms: {e}")
    exit()

# --- Initialize arrays ---
times = []      # Time in ps
widths = []     # Width in Å

# --- Process Trajectory ---
print("\nProcessing frames...")
for ts in u.trajectory:
    if ts.time < t_start_ps:
        continue  # Skip frames before starting time
    
    upper_leaflet = lipid_atoms.positions[:, 2] > z_center
    lower_leaflet = lipid_atoms.positions[:, 2] < z_center
    
    if not upper_leaflet.any() or not lower_leaflet.any():
        print(f"Warning: Skipping frame {ts.frame} due to missing atoms in leaflets.")
        continue
    
    upper_z = lipid_atoms.positions[upper_leaflet, 2].mean()
    lower_z = lipid_atoms.positions[lower_leaflet, 2].mean()
    
    widths.append(upper_z - lower_z)
    times.append(ts.time)

# Convert to numpy arrays
times = np.array(times) / 1000  # Convert ps to ns
widths = np.array(widths)
areas = (widths ** 2) / n_residues  # Area proxy: width² / user-defined n_residues

# --- Save to CSV with statistics ---
csv_file = f"{output_name}.csv"
with open(csv_file, "w") as f:
    f.write("Time (ns),Width (Å),Area (Å²/residue)\n")
    for t, w, a in zip(times, widths, areas):
        f.write(f"{t:.4f},{w:.4f},{a:.4f}\n")
    f.write(f"Mean ± Std,{widths.mean():.4f} ± {widths.std():.4f},{areas.mean():.4f} ± {areas.std():.4f}\n")

print(f"\nData saved to '{csv_file}' with mean and std at the end!")

# --- Plot Results ---
png_file = f"{output_name}.png"
plt.figure(figsize=(10, 5))
plt.subplot(1, 2, 1)
plt.plot(times, widths, color="teal", label="Width")
plt.xlabel("Time (ns)")
plt.ylabel("Width (Å)")
plt.title("Membrane Width")
plt.grid(True, alpha=0.3)

plt.subplot(1, 2, 2)
plt.plot(times, areas, color="purple", label="Area")
plt.xlabel("Time (ns)")
plt.ylabel("Area (Å²/residue)")
plt.title("Membrane Area Proxy")
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(png_file, dpi=300)
print(f"Plot saved as '{png_file}'!")

# --- Statistics ---
print("\n=== Results ===")
print(f"Average width: {widths.mean():.2f} ± {widths.std():.2f} Å")
print(f"Average area: {areas.mean():.2f} ± {areas.std():.2f} Å²/residue")

