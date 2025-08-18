import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import csv

def run_tilt_angle():
    print("Working...")

    # Archivos de entrada
    tpr_file = input("Enter the .tpr file: ")
    xtc_file = input("Enter the .xtc file: ")

    u = mda.Universe(tpr_file, xtc_file)

    # Selección de residuo
    residue_name = input("Enter residue name: ")

    residues = u.select_atoms(f"resname {residue_name}").residues
    if len(residues) == 0:
        raise ValueError(f"No residues found with resname {residue_name}")

    residue = residues[0]

    # Mostrar información de los átomos disponibles en ese residuo
    print(f"\nFound {len(residues)} residues with name '{residue_name}'.")
    print(f"Showing atoms of the first residue (resid {residue.resid}):\n")
    for atom in residue.atoms:
        print(f"  {atom.index:5d}  {atom.name:6s}  {atom.type:6s}")

    # Pedir átomos de referencia
    atom1_name = input("\nEnter first reference atom name: ")
    atom2_name = input("Enter second reference atom name: ")

    atom1 = residue.atoms.select_atoms(f"name {atom1_name}")[0]
    atom2 = residue.atoms.select_atoms(f"name {atom2_name}")[0]

    vectors = []
    for ts in u.trajectory:
        pos1 = atom1.position
        pos2 = atom2.position
        vector = pos2 - pos1
        vectors.append(vector)

    vectors = np.array(vectors)

    # --- Calcular ángulos de inclinación respecto al eje Z ---
    angles = []
    z_axis = np.array([0, 0, 1])

    for v in vectors:
        cos_theta = np.dot(v, z_axis) / (np.linalg.norm(v) * np.linalg.norm(z_axis))
        theta = np.degrees(np.arccos(cos_theta))

        # Plegar ángulos para que estén en [0°, 90°]
        if theta > 90:
            theta = 180 - theta

        angles.append(theta)

    angles = np.array(angles)

    # --- Histograma + PDF Gaussiana ---
    plt.hist(angles, bins=50, density=True, alpha=0.6, color="skyblue", edgecolor="black", label="Histogram")

    # Calcular PDF normal
    mu, sigma = np.mean(angles), np.std(angles)
    x_vals = np.linspace(0, 90, 500)
    pdf_vals = norm.pdf(x_vals, mu, sigma)

    # Encontrar el ángulo más probable (máximo de la PDF)
    max_idx = np.argmax(pdf_vals)
    most_probable_angle = x_vals[max_idx]
    max_probability = pdf_vals[max_idx]

    # Dibujar PDF
    plt.plot(x_vals, pdf_vals, 'r-', lw=2, label=f"Normal PDF (μ={mu:.2f}, σ={sigma:.2f})")
    plt.axvline(most_probable_angle, color="k", linestyle="--", label=f"Most probable: {most_probable_angle:.2f}°")
    plt.xlabel("Tilt angle (degrees)")
    plt.ylabel("Probability")
    plt.title(f"Tilt Angle Distribution (0°–90°) for {residue_name}")
    plt.legend()
    plt.savefig("angles_hist.png", dpi=300)
    plt.close()

    # --- Guardar resultados en CSV ---
    with open("angles_data.csv", "w", newline="") as f:
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
    print("Figures saved as 'vectors.png' and 'angles_hist.png'")
    print("Angles and statistics saved in 'angles_data.csv'")

# --- Ejecutar automáticamente si se llama desde terminal ---
if __name__ == "__main__":
    run_tilt_angle()

