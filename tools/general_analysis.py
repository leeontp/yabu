#!/usr/bin/env python3
"""
system_summary_interactive.py

Script interactivo que genera un resumen del sistema molecular
y lo guarda en un archivo CSV.
"""

import csv
from collections import Counter
import MDAnalysis as mda

def main():
    print("=== Resumen de sistema molecular ===")

    # Pedir archivos al usuario
    top = input("Ingrese el archivo de topología (.tpr, .gro, .pdb, etc.): ").strip()
    traj = input("Ingrese el archivo de trayectoria (opcional, presione Enter si no aplica): ").strip()
    out = input("Ingrese el nombre del archivo de salida (.csv) [system_summary.csv]: ").strip()
    if out == "":
        out = "system_summary.csv"

    # Cargar universo
    if traj:
        u = mda.Universe(top, traj)
    else:
        u = mda.Universe(top)

    # Contar residuos
    res_counts = Counter([res.resname for res in u.residues])

    # Preparar tabla
    rows = []
    for res in sorted(res_counts):
        res_obj = u.select_atoms(f"resname {res}")
        n_atoms = len(res_obj)
        n_residues = res_counts[res]
        rows.append({"Resname": res, "Residues": n_residues, "Atoms": n_atoms})

    # Calcular características globales
    total_atoms = len(u.atoms)
    total_residues = len(u.residues)
    total_segments = len(u.segments)
    box = u.dimensions if u.dimensions is not None else ["NA"]*6
    coords_min = u.atoms.positions.min(axis=0)
    coords_max = u.atoms.positions.max(axis=0)

    # Escribir CSV
    with open(out, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["=== Tabla de residuos ==="])
        writer.writerow(["Resname", "Residues", "Atoms"])
        for row in rows:
            writer.writerow([row["Resname"], row["Residues"], row["Atoms"]])

        writer.writerow([])
        writer.writerow(["=== Características generales ==="])
        writer.writerow(["Total atoms", total_atoms])
        writer.writerow(["Total residues", total_residues])
        writer.writerow(["Total segments", total_segments])
        writer.writerow(["Box dimensions [Å] (lx, ly, lz, alpha, beta, gamma)"] + list(box))
        writer.writerow(["Bounding box min [Å]"] + list(coords_min))
        writer.writerow(["Bounding box max [Å]"] + list(coords_max))

    print(f"\n✅ Resumen guardado en {out}")

if __name__ == "__main__":
    main()

