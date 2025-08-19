#!/usr/bin/env python3
"""
membrane_thickness_auto.py

Calculate membrane thickness by measuring the distance between phosphorus atoms
of upper and lower leaflets. Automatically separates top and bottom leaflet.
"""

import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import csv

def run_membrane_thickness():
    print("=== Membrane Thickness Calculation ===")

    # Input files
    gro_file = input("Enter the structure file (.gro): ").strip()
    traj_file = input("Enter the trajectory file (.xtc/.trr): ").strip()
    
    u = mda.Universe(gro_file, traj_file)

    # Residue selection (same for top and bottom)
    resname = input("Enter residue name for the lipid (e.g., POPC): ").strip()
    atom_name = input("Enter atom name for phosphorus (usually 'P'): ").strip()

    # Time selection
    start_time_ns = float(input("Enter start time in ns: ").strip())

    # Output files
    output_base = input("Enter base name for output files (without extension) [membrane_thickness]: ").strip()
    if output_base == "":
        output_base = "membrane_thickness"
    out_csv = f"{output_base}.csv"
    out_png = f"{output_base}.png"

    # Select phosphorus atoms
    phosphorus_atoms = u.select_atoms(f"resname {resname} and name {atom_name}")
    if len(phosphorus_atoms) == 0:
        raise ValueError("No phosphorus atoms found for the specified residue.")

    thickness_list = []
    time_list = []

    for ts in u.trajectory:
        time_ns = ts.time / 1000.0  # ps -> ns
        if time_ns < start_time_ns:
            continue

        # Separar hoja superior e inferior por coordenada Z
        z_coords = phosphorus_atoms.positions[:, 2]
        z_median = np.median(z_coords)
        upper_atoms = z_coords[z_coords >= z_median]
        lower_atoms = z_coords[z_coords < z_median]

        # Calcular posición promedio de cada hoja
        z_upper = np.mean(upper_atoms)
        z_lower = np.mean(lower_atoms)
        thickness = abs(z_upper - z_lower)

        thickness_list.append(thickness)
        time_list.append(time_ns)

    thickness_array = np.array(thickness_list)
    time_array = np.array(time_list)

    # Statistics
    mean_thickness = np.mean(thickness_array)
    std_thickness = np.std(thickness_array)

    # --- Save CSV ---
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Time (ns)", "Membrane thickness (nm)"])
        for t, d in zip(time_array, thickness_array):
            writer.writerow([t, d])
        writer.writerow([])
        writer.writerow(["Mean thickness (nm)", "Std dev (nm)"])
        writer.writerow([mean_thickness, std_thickness])

    # --- Plot thickness over time ---
    plt.figure(figsize=(8,5))
    plt.plot(time_array, thickness_array, '-o', markersize=3, color='skyblue', label='Membrane thickness')
    plt.axhline(mean_thickness, color='r', linestyle='--', label=f"Mean: {mean_thickness:.3f} nm")
    plt.fill_between(time_array, mean_thickness-std_thickness, mean_thickness+std_thickness, color='r', alpha=0.2, label=f"Std dev: {std_thickness:.3f} nm")
    plt.xlabel("Time (ns)")
    plt.ylabel("Membrane thickness (nm)")
    plt.title("Membrane Thickness Over Time")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

    print("\n✅ Membrane thickness calculation finished.")
    print(f"CSV data saved as '{out_csv}'")
    print(f"Thickness plot saved as '{out_png}'")
    print(f"Average thickness: {mean_thickness:.3f} nm ± {std_thickness:.3f} nm")

if __name__ == "__main__":
    run_membrane_thickness()