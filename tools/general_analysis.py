#!/usr/bin/env python3
"""
system_summary_interactive.py

Interactive script that generates a molecular system summary
and saves it to a CSV file.
"""

import csv
from collections import Counter
import MDAnalysis as mda

def main():
    print("=== Molecular System Summary ===")

    # Ask for input files
    top = input("Enter topology file (.tpr, .gro, .pdb, etc.): ").strip()
    traj = input("Enter trajectory file (optional, press Enter if not applicable): ").strip()
    
    # Output file name
    output_name = input("Enter output file name (without extension) [system_summary]: ").strip()
    if output_name == "":
        output_name = "system_summary"
    out_csv = f"{output_name}.csv"

    # Load universe
    if traj:
        u = mda.Universe(top, traj)
    else:
        u = mda.Universe(top)

    # Count residues
    res_counts = Counter([res.resname for res in u.residues])

    # Prepare table
    rows = []
    for res in sorted(res_counts):
        res_obj = u.select_atoms(f"resname {res}")
        n_atoms = len(res_obj)
        n_residues = res_counts[res]
        rows.append({"Resname": res, "Residues": n_residues, "Atoms": n_atoms})

    # Calculate global characteristics
    total_atoms = len(u.atoms)
    total_residues = len(u.residues)
    total_segments = len(u.segments)
    box = u.dimensions if u.dimensions is not None else ["NA"]*6
    coords_min = u.atoms.positions.min(axis=0)
    coords_max = u.atoms.positions.max(axis=0)

    # Write CSV
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["=== Residue Table ==="])
        writer.writerow(["Resname", "Residues", "Atoms"])
        for row in rows:
            writer.writerow([row["Resname"], row["Residues"], row["Atoms"]])

        writer.writerow([])
        writer.writerow(["=== Global System Properties ==="])
        writer.writerow(["Total atoms", total_atoms])
        writer.writerow(["Total residues", total_residues])
        writer.writerow(["Total segments", total_segments])
        writer.writerow(["Box dimensions [Å] (lx, ly, lz, alpha, beta, gamma)"] + list(box))
        writer.writerow(["Bounding box min [Å]"] + list(coords_min))
        writer.writerow(["Bounding box max [Å]"] + list(coords_max))

    print(f"\n✅ Summary saved to {out_csv}")

if __name__ == "__main__":
    main()
    