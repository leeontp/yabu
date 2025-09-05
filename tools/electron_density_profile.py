import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse

# ------------------------------
# Electron dictionary
electron_dict = {
    'H': 1,
    'C': 6,
    'N': 7,
    'O': 8,
    'P': 15,
    'S': 16
}

# ------------------------------
# Function to calculate electron density (optimized: no full memory storage)
def calculate_density(atomgroup, u, bin_edges):
    density_sum = np.zeros(len(bin_edges) - 1)
    n_frames = 0

    for ts in u.trajectory:
        z_coords = []
        electron_counts = []
        for atom in atomgroup:
            elem = atom.name[0].upper()
            electron_counts.append(electron_dict.get(elem, 0))
            z_coords.append(atom.position[2])

        hist, _ = np.histogram(z_coords, bins=bin_edges, weights=electron_counts)
        density_sum += hist
        n_frames += 1

    # Average over frames
    density_avg = density_sum / n_frames

    # Normalize per area and bin width
    area = u.dimensions[0] * u.dimensions[1]
    bin_width = bin_edges[1] - bin_edges[0]
    density_avg = density_avg / (area * bin_width)

    return density_avg

# ------------------------------
# Main
def main():
    parser = argparse.ArgumentParser(description="Calculate electron density profiles")
    parser.add_argument("-g", "--gro", help="Input .gro file")
    parser.add_argument("-t", "--traj", help="Input trajectory file (.xtc/.trr)")
    parser.add_argument("-r", "--residues", help="Residues separated by commas (e.g., POPC,DPPC)")
    parser.add_argument("-o", "--output", help="Base name for output files")
    parser.add_argument("-n", "--nbins", type=int, help="Number of bins along z")
    parser.add_argument("-help", action="store_true", help="Show usage and exit")

    args = parser.parse_args()

    if args.help:
        print("""
Usage:
    python script.py -g system.gro -t traj.xtc -r POPC,DPPC -o output -n 200

Flags:
    -g   Input .gro file
    -t   Input trajectory file (.xtc/.trr)
    -r   Residue names separated by commas (e.g., POPC,DPPC)
    -o   Base name for output files
    -n   Number of bins along z
""")
        return

    # Interactive mode if no args
    if not any([args.gro, args.traj, args.residues, args.output, args.nbins]):
        gro_file = input("Enter the .gro file: ")
        traj_file = input("Enter the trajectory file (.xtc/.trr): ")
        residues_input = input("Enter the residue names separated by commas (e.g., POPC,DPPC): ")
        output_name = input("Enter the base name for output files (without extension): ")
        n_bins = int(input("Enter the number of bins along z (e.g., 100): "))
    else:
        gro_file = args.gro
        traj_file = args.traj
        residues_input = args.residues
        output_name = args.output
        n_bins = args.nbins

    residues = [r.strip() for r in residues_input.split(",")]

    # Load simulation
    u = mda.Universe(gro_file, traj_file)

    # Prepare bins
    z_min = 0.0
    z_max = u.dimensions[2]  # z box length
    bin_edges = np.linspace(z_min, z_max, n_bins + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    # Store profiles
    all_profiles = {}

    # System density
    print("Processing system density...")
    density_all = calculate_density(u.atoms, u, bin_edges)
    all_profiles["ALL"] = density_all

    # Residues density
    for resname in residues:
        print(f"Processing residue {resname}...")
        selection = u.select_atoms(f"resname {resname}")
        if len(selection) == 0:
            print(f"Warning: No atoms found for {resname}, skipping.")
            continue
        density_res = calculate_density(selection, u, bin_edges)
        all_profiles[resname] = density_res

    # Save CSV
    df_all = pd.DataFrame({'z (Å)': bin_centers})
    for name, profile in all_profiles.items():
        df_all[name] = profile
    csv_filename = f"{output_name}_profiles.csv"
    df_all.to_csv(csv_filename, index=False)
    print(f"CSV saved: {csv_filename}")

    # Plot
    fig, ax1 = plt.subplots(figsize=(8, 5))

    ax1.plot(bin_centers, all_profiles["ALL"], color="black", linewidth=2, label="ALL")
    ax1.set_xlabel("z (Å)")
    ax1.set_ylabel("Electron Density (System) [e/Å³]", color="black")
    ax1.tick_params(axis="y", labelcolor="black")

    # Residue profiles in right axis
    ax2 = ax1.twinx()
    for name, profile in all_profiles.items():
        if name == "ALL":
            continue
        ax2.plot(bin_centers, profile, linestyle="--", label=name)

    ax2.set_ylabel("Electron Density (Residues) [e/Å³]", color="blue")
    ax2.tick_params(axis="y", labelcolor="blue")

    # Combine legends
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="best")

    plt.title("Electron Density Profiles (System vs Residues)")
    plt.tight_layout()
    png_filename = f"{output_name}_profiles.png"
    plt.savefig(png_filename, dpi=300)
    plt.close()
    print(f"Plot saved: {png_filename}")


if __name__ == "__main__":
    main()
