import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import sys
import time
import os

fox_confused = r"""
   |\_/|  
   (o.O)?  YABU is confused...
   > ^ < 
"""

fox_frames = [
r"""
   |\_/|  
   (o.o)  YABU
   > ^ < 
""",
r"""
   |\_/|  
   ( -.-) zZ
   > ^ < 
""",
r"""
   |\_/|  
   (o.o) /
   > ^ < 
"""
]

print(r"""
==============================
   🦊  Tilt Angle Tool (YABU)
==============================
""")

try:
    tpr = input("Enter topology file (.tpr/.gro/.pdb): ")
    xtc = input("Enter trajectory file (.xtc/.trr): ")

    if not os.path.exists(tpr) or not os.path.exists(xtc):
        raise FileNotFoundError("One or both input files do not exist.")

    u = mda.Universe(tpr, xtc)

    res_input = input("Enter the residue (e.g., 'CYS', '45'): ")

    if res_input.isdigit():
        selection = f"resid {res_input}"
    else:
        selection = f"resname {res_input}"

    res_atoms = u.select_atoms(selection)

    if len(res_atoms) == 0:
        raise ValueError("No atoms found for the given residue.")

    print("Atoms in this selection:", [a.name for a in res_atoms])

    atom1_name = input("Enter first atom name (e.g., 'CA'): ")
    atom2_name = input("Enter second atom name (e.g., 'CB'): ")

    atom1 = res_atoms.select_atoms(f"name {atom1_name}")
    atom2 = res_atoms.select_atoms(f"name {atom2_name}")

    if len(atom1) == 0 or len(atom2) == 0:
        raise ValueError("One or both atom names not found in the residue.")

    vectors = []
    frame_count = len(u.trajectory)

    print("\nCalculating vectors:")
    for i, ts in enumerate(u.trajectory):
        vec = atom2.positions[0] - atom1.positions[0]
        vectors.append(vec)
        fox = fox_frames[i % len(fox_frames)]
        sys.stdout.write(f"\r{fox}\nFrame {i+1}/{frame_count}")
        sys.stdout.flush()
        time.sleep(0.05)

    vectors = np.array(vectors)

    print("\n\n🦊 Done! YABU finished calculating vectors.\n")

    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111, projection='3d')

    origin = np.zeros((len(vectors), 3))
    ax.quiver(origin[:,0], origin[:,1], origin[:,2],
              vectors[:,0], vectors[:,1], vectors[:,2],
              length=1.0, normalize=True, color="blue", alpha=0.6)

    ax.set_title(f"Vectors from {atom1_name} to {atom2_name} in {res_input}")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    plt.show()

except Exception as e:
    print(fox_confused)
    print(f"ERROR: {e}\n")


