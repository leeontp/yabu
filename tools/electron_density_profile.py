import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# ------------------------------
# User input
# ------------------------------
gro_file = input("Enter the .gro file: ")
traj_file = input("Enter the trajectory file (.xtc/.trr): ")
residues_input = input("Enter the residue names separated by commas (e.g., POPC,DPPC): ")
output_name = input("Enter the base name for output files (without extension): ")
n_bins = int(input("Enter the number of bins along z (e.g., 100): "))

residues = [r.strip() for r in residues_input.split(",")]

# ------------------------------
# Load the simulation
# ------------------------------
u = mda.Universe(gro_file, traj_file)

# ------------------------------
# Define number of electrons per atom type
electron_dict = {
    'H': 1,
    'C': 6,
    'N': 7,
    'O': 8,
    'P': 15,
    'S': 16
}

# ------------------------------
# Prepare bins
z_min = 0.0
z_max = u.dimensions[2]  # z limit
bin_edges = np.linspace(z_min, z_max, n_bins+1)
bin_centers = 0.5*(bin_edges[:-1] + bin_edges[1:])

# ------------------------------
# Dictionary to store all profiles
all_profiles = {}

# ------------------------------
# Function to calculate electron density of an atom group
def calculate_density(atomgroup):
    density_sum = np.zeros(n_bins)
    n_frames = 0
    for ts in u.trajectory:
        n_frames += 1
        z_coords = []
        electron_counts = []
        for atom in atomgroup:
            elem = atom.name[0].upper()
            electron_counts.append(electron_dict.get(elem, 0))
            z_coords.append(atom.position[2])
        hist, _ = np.histogram(z_coords, bins=bin_edges, weights=electron_counts)
        density_sum += hist
    density_avg = density_sum / n_frames
    area = u.dimensions[0] * u.dimensions[1]
    bin_width = bin_edges[1] - bin_edges[0]
    density_avg = density_avg / (area * bin_width)
    return density_avg

# ------------------------------
# Calculate density for each residue
for resname in residues:
    print(f"Processing residue: {resname}")
    selection = u.select_atoms(f"resname {resname}")
    if len(selection) == 0:
        print(f"Warning: No atoms found for residue {resname}. Skipping.")
        continue
    density_avg = calculate_density(selection)
    all_profiles[resname] = density_avg

# ------------------------------
# Calculate density for the whole system
print("Processing density for the whole system")
all_density = calculate_density(u.atoms)
all_profiles['ALL'] = all_density

# ------------------------------
# Create DataFrame and save CSV
df_all = pd.DataFrame({'z (Å)': bin_centers})
for name, profile in all_profiles.items():
    df_all[name] = profile

csv_filename = f"{output_name}_all_residues.csv"
df_all.to_csv(csv_filename, index=False)
print(f"CSV with all profiles saved: {csv_filename}")

# ------------------------------
# Plot all profiles together with double y-axis
fig, ax1 = plt.subplots(figsize=(8,5))

# Left axis -> full system profile
ax1.plot(bin_centers, all_profiles['ALL'], color='black', linewidth=2, label='ALL')
ax1.set_xlabel('z (Å)')
ax1.set_ylabel('Electron Density (System) [e/Å³]', color='black')
ax1.tick_params(axis='y', labelcolor='black')

# Adjust limits of left axis
ax1.set_ylim(min(all_profiles['ALL'])*0.95, max(all_profiles['ALL'])*1.05)

# Right axis -> residue profiles
ax2 = ax1.twinx()
for name, profile in all_profiles.items():
    if name == 'ALL':
        continue
    ax2.plot(bin_centers, profile, label=name, linestyle='--')

ax2.set_ylabel('Electron Density (Residues) [e/Å³]', color='blue')
ax2.tick_params(axis='y', labelcolor='blue')

# Adjust limits of right axis
residue_profiles = [profile for name, profile in all_profiles.items() if name != 'ALL']
if residue_profiles:  # avoid error if no residues
    all_res_vals = np.concatenate(residue_profiles)
    ax2.set_ylim(all_res_vals.min()*0.95, all_res_vals.max()*1.05)

# Combine legends
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='best')

plt.title('Electron Density Profiles (System vs Residues)')
plt.tight_layout()
png_filename = f"{output_name}_all_residues.png"
plt.savefig(png_filename, dpi=300)
plt.close()

print(f"Plot with all profiles saved: {png_filename}")
print("Final calculation completed for all residues and the total system.")
