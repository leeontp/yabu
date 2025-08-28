import MDAnalysis as mda
import numpy as np
import pandas as pd

def calculate_distances(gro_file, traj_file=None, ref_atoms=None, dist_atoms=None, output_csv="distances.csv"):
    """
    Calculate distances between two groups of atoms and save the results in a tidy .csv file
    
    Parameters:
    - gro_file: topology file (.gro)
    - traj_file: trajectory file (optional, e.g. .xtc, .trr)
    - ref_atoms: atom selection for the reference group (str, MDAnalysis selection)
    - dist_atoms: atom selection for the distance group (str, MDAnalysis selection)
    - output_csv: name of the output .csv file
    """
    
    # Load the system
    if traj_file:
        u = mda.Universe(gro_file, traj_file)
    else:
        u = mda.Universe(gro_file)
    
    # Select atoms
    ref_group = u.select_atoms(ref_atoms)
    dist_group = u.select_atoms(dist_atoms)

    if len(ref_group) == 0 or len(dist_group) == 0:
        raise ValueError("One of the atom selections is empty. Please check your inputs.")

    # Prepare list for results
    records = []

    # Calculate distances for each frame
    for ts in u.trajectory:
        for ref in ref_group:
            for dist in dist_group:
                d_ang = np.linalg.norm(ref.position - dist.position)  # in Å
                d_nm = d_ang / 10.0  # convert to nm
                atom_pair = f"{ref.id}-{ref.name} ↔ {dist.id}-{dist.name}"
                records.append([ts.frame, atom_pair, d_ang, d_nm])

    # Create dataframe
    distances_df = pd.DataFrame(records, columns=["Frame", "Atom Pair", "Distance (Å)", "Distance (nm)"])

    # Save to .csv
    distances_df.to_csv(output_csv, index=False)
    print(f"\nDistance file saved as {output_csv}")


if __name__ == "__main__":
    print("=== Script to calculate distances between atoms with MDAnalysis ===\n")
    
    gro_file = input("Enter the .gro file name: ").strip()
    traj_file = input("Enter the trajectory file (optional, leave empty if none): ").strip()
    traj_file = traj_file if traj_file else None
    ref_atoms = input("Enter the atom selection for the reference group (e.g., resid 1 and name CA): ").strip()
    dist_atoms = input("Enter the atom selection for the distance group (e.g., resname SOL and name OW): ").strip()
    output_csv = input("Enter the output .csv file name: ").strip()
    if not output_csv.endswith(".csv"):
        output_csv += ".csv"

    calculate_distances(
        gro_file=gro_file,
        traj_file=traj_file,
        ref_atoms=ref_atoms,
        dist_atoms=dist_atoms,
        output_csv=output_csv
    )