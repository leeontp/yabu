import MDAnalysis as mda
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import itertools

# ------------------------------
# User Input
# ------------------------------
gro_file = input("Enter the .gro file: ")
traj_file = input("Enter the trajectory file (.xtc/.trr): ")

u = mda.Universe(gro_file, traj_file)

# Normal to the membrane (z-axis in most MD setups)
membrane_normal = np.array([0, 0, 1])

# Store all results
all_results = []

# ------------------------------
# Interactive Selection
# ------------------------------
while True:
    resid_input = input("Enter the residue ID you want to analyze: ")
    resid = int(resid_input)
    res = u.residues[resid]

    print(f"\nAtoms in residue {res.resname}-{resid}:")
    for atom in res.atoms:
        print(f"{atom.name} ({atom.element})")

    pairs = []
    while True:
        C_name = input("Enter the Carbon atom name (or press Enter to finish this residue): ")
        if C_name.strip() == "":
            break
        H_name = input(f"Enter the Hydrogen atom name bound to {C_name}: ")
        pairs.append((C_name, H_name))

    if not pairs:
        print("No pairs selected, skipping this residue.")
    else:
        # ------------------------------
        # Calculate S_CD
        # ------------------------------
        residue_results = {C: [] for C, H in pairs}

        for ts in u.trajectory:
            for C_name, H_name in pairs:
                try:
                    C_atom = res.atoms.select_atoms(f"name {C_name}")[0]
                    H_atom = res.atoms.select_atoms(f"name {H_name}")[0]

                    CH_vector = H_atom.position - C_atom.position
                    CH_vector /= np.linalg.norm(CH_vector)

                    cos_theta = np.dot(CH_vector, membrane_normal)
                    cos2 = cos_theta**2

                    S_cd = 0.5 * (3 * cos2 - 1)
                    residue_results[C_name].append(S_cd)
                except IndexError:
                    print(f"Warning: atom {C_name} or {H_name} not found in residue {res.resname}-{resid}")

        # Average results
        for C_name in residue_results:
            avg_Scd = np.mean(residue_results[C_name])
            all_results.append({
                "Residue": f"{res.resname}-{resid}",
                "Carbon": C_name,
                "S_CD": avg_Scd
            })

    more = input("Do you want to analyze another residue? (y/n): ")
    if more.lower() != "y":
        break

# ------------------------------
# Save results
# ------------------------------
df = pd.DataFrame(all_results)
df.to_csv("deuterium_order_parameters.csv", index=False)
print("\nResults saved to deuterium_order_parameters.csv")

# ------------------------------
# Plot results
# ------------------------------
markers = itertools.cycle(["o", "s", "^", "D", "v", "x", "*", "p", "h"])
plt.figure(figsize=(8,6))

for resid, group in df.groupby("Residue"):
    x = np.arange(1, len(group)+1)
    y = group["S_CD"].values
    plt.plot(x, y, marker=next(markers), label=resid, linestyle='-', linewidth=1.2)

plt.xlabel("Carbon index")
plt.ylabel("S_CD")
plt.title("Deuterium Order Parameters")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()