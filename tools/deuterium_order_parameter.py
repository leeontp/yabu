import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import argparse
import csv
import sys
import os

def calculate_scd(universe, residue_name, carbon_hydrogen_pairs):
    """Calculate deuterium order parameter S_CD for a given residue and C-H pairs"""
    scd_values = []

    for carbon, hydrogen in carbon_hydrogen_pairs:
        carbon_atoms = universe.select_atoms(f"resname {residue_name} and name {carbon}")
        hydrogen_atoms = universe.select_atoms(f"resname {residue_name} and name {hydrogen}")

        if len(carbon_atoms) == 0 or len(hydrogen_atoms) == 0:
            print(f"⚠️ Pair not found: {carbon}-{hydrogen} in {residue_name}")
            continue

        carbon = carbon_atoms[0]
        hydrogen = hydrogen_atoms[0]

        cos2 = []
        for ts in universe.trajectory:
            vec = hydrogen.position - carbon.position
            norm = np.linalg.norm(vec)
            if norm == 0:
                continue
            vec /= norm
            cos_theta = vec[2]  # z-component
            cos2.append(cos_theta ** 2)

        if cos2:
            avg_cos2 = np.mean(cos2)
            scd = -0.5 * (3 * avg_cos2 - 1)
            scd_values.append(scd)

    return scd_values


def main():
    parser = argparse.ArgumentParser(
        description="Calculate Deuterium Order Parameter (S_CD) for lipid membranes",
        add_help=False
    )

    parser.add_argument("-g", "--gro", type=str, help="Input .gro file")
    parser.add_argument("-t", "--traj", type=str, help="Input trajectory file (.xtc/.trr)")
    parser.add_argument("-i", "--input", type=str, help="Input file with residues and C-H pairs")
    parser.add_argument("-o", "--output", type=str, help="Output CSV file name")
    parser.add_argument("--help", action="store_true", help="Show this help message and exit")

    args = parser.parse_args()

    if args.help:
        print("""
Usage:
  python deuterium_order_parameter.py -g system.gro -t traj.xtc -i input_pairs.txt -o results.csv

Flags:
  -g   Input .gro file
  -t   Input trajectory file (.xtc/.trr)
  -i   Input file with residues and C-H pairs
  -o   Output CSV file name

If no flags are provided, the script will ask interactively.

📌 Input file template (-i):
Each line: RESNAME C1-H1 C2-H2 C3-H3 ...
Example:
POPC C2-H11 C3-H12 C4-H13
DPPC C2-H21 C3-H22 C4-H23
""")
        sys.exit()

    # Interactive mode if no arguments are given
    if not any([args.gro, args.traj, args.input, args.output]):
        gro_file = input("Enter the .gro file: ")
        traj_file = input("Enter the trajectory file (.xtc/.trr): ")
        output_csv = input("Enter the output CSV file name: ")
        input_file = input("Enter the input file with residues and C-H pairs: ")
    else:
        gro_file = args.gro
        traj_file = args.traj
        output_csv = args.output
        input_file = args.input

    # Load trajectory
    u = mda.Universe(gro_file, traj_file)

    # Read residue and pairs
    residues_data = {}
    with open(input_file, "r") as f:
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue
            resname = parts[0]
            pairs = [p.split("-") for p in parts[1:]]
            residues_data[resname] = pairs

    # Calculate and save results
    results = []
    for resname, pairs in residues_data.items():
        scd_values = calculate_scd(u, resname, pairs)
        avg_scd = np.mean(scd_values) if scd_values else None
        if avg_scd is not None:
            for idx, scd in enumerate(scd_values, start=1):
                results.append([resname, idx, scd])

    # Save to CSV
    with open(output_csv+".csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Residue", "Carbon", "S_CD"])
        writer.writerows(results)

    print(f"✅ Results saved to {output_csv}")

    # Plot results
    plt.figure(figsize=(8, 6))
    markers = ["o", "s", "D", "^", "v", "<", ">", "p", "h", "x", "*"]
    for i, resname in enumerate(residues_data.keys()):
        x = [r[1] for r in results if r[0] == resname]
        y = [r[2] for r in results if r[0] == resname]
        plt.plot(x, y, marker=markers[i % len(markers)], label=resname, linestyle="--")

    plt.xlabel("Carbon number", fontsize=12)
    plt.ylabel("S_CD", fontsize=12)
    plt.title("Deuterium Order Parameters", fontsize=14)
    plt.xticks(range(1, max(r[1] for r in results) + 1))  # integers en el eje X
    plt.legend()
    plt.tight_layout()

    # Save figure instead of show
    fig_name = os.path.splitext(output_csv)[0] + ".png"
    plt.savefig(fig_name, dpi=300)
    print(f"📊 Plot saved to {fig_name}")


if __name__ == "__main__":
    main()
