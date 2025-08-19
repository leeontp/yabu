#!/usr/bin/env python3
"""
tilt_angle_analysis_interactive.py

Interactive script to calculate tilt angles of a residue with respect to the Z-axis,
plot a histogram with Gaussian PDF, and save results to CSV and PNG files.
"""

import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import csv

def run_tilt_angle():
    print("=== Tilt Angle Analysis ===")

    # Input files
    tpr_file = input("Enter the topology file (.tpr): ").strip()
    xtc_file = input("Enter the trajectory file (.xtc): ").strip()

    u = mda.Universe(tpr_file, xtc_file)

    # Select residue
    residue_name = input("Enter residue name: ").strip()
    residues = u.select_atoms(f"resname {residue_name}").residues
    if len(residues) == 0:
        raise ValueError(f"No residues found with name '{residue_name}'")

    residue = residues[0]

    # Show atoms in first residue
    print(f"\nFound {len(residues)} residues with name '{residue_name}'.")
    print(f"Atoms in the first residue (resid {residue.resid}):\n")
    for atom in residue.atoms:
        print(f"  {atom.index:5d}  {atom.name:6s}  {atom.type:6s}")

    # Reference atoms
    atom1_name = input("\nEnter first reference atom name: ").strip()
    atom2_name = input("Enter second reference atom name: ").strip()

    atom1 = residue.atoms.select_atoms(f"name {atom1_name}")[0]
    atom2 = residue.atoms.select_atoms(f"name {atom2_name}")[0]

    vectors = []
    for ts in u.trajectory:
        pos1 = atom1.position
        pos2 = atom2.position
        vectors.append(pos2 - pos1)
    vectors = np.array(vectors)

    # --- Compute tilt angles with respect to Z-axis ---
    z_axis = np.array([0, 0, 1])
    angles = []
    for v in vectors:
        cos_theta = np.dot(v, z_axis) / (np.linalg.norm(v) * np.linalg.norm(z_axis))
        theta = np.degrees(np.arccos(cos_theta))
        if theta > 90:
            theta = 180 - theta
        angles.append(theta)
    angles = np.array(angles)

    # --- Ask for output file names ---
    output_base = input("\nEnter base name for output files (without extension) [tilt_angle]: ").strip()
    if output_base == "":
        output_base = "tilt_angle"
    out_csv = f"{output_base}.csv"
    out_png = f"{output_base}.png"

    # --- Histogram + Gaussian PDF ---
    plt.hist(angles, bins=50, density=True, alpha=0.6, color="skyblue", edgecolor="black", label="Histogram")

    mu, sigma = np.mean(angles), np.std(angles)
    x_vals = np.linspace(0, 90, 500)
    pdf_vals = norm.pdf(x_vals, mu, sigma)

    max_idx = np.argmax(pdf_vals)
    most_probable_angle = x_vals[max_idx]
    max_probability = pdf_vals[max_idx]

    plt.plot(x_vals, pdf_vals, 'r-', lw=2, label=f"Normal PDF (μ={mu:.2f}, σ={sigma:.2f})")
    plt.axvline(most_probable_angle, color="k", linestyle="--", label=f"Most probable: {most_probable_angle:.2f}°")
    plt.xlabel("Tilt angle (degrees)")
    plt.ylabel("Probability")
    plt.title(f"Tilt Angle Distribution (0°–90°) for {residue_name}")
    plt.legend()
    plt.savefig(out_png, dpi=300)
    plt.close()

    # --- Save angles and statistics to CSV ---
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Angle (deg)"])
        for a in angles:
            writer.writerow([a])
        writer.writerow([])
        writer.writerow(["Mean (deg)", "Std dev (deg)"])
        writer.writerow([mu, sigma])
        writer.writerow([])
        writer.writerow(["Most probable angle (deg)", "Max probability"])
        writer.writerow([most_probable_angle, max_probability])

    print("\n✅ Analysis finished.")
    print(f"Figure saved as '{out_png}'")
    print(f"Angles and statistics saved in '{out_csv}'")

if __name__ == "__main__":
    run_tilt_angle()

