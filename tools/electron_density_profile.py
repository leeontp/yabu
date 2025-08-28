import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# ------------------------------
# Input del usuario
# ------------------------------
gro_file = input("Ingrese el archivo .gro: ")
traj_file = input("Ingrese el archivo de trayectoria (.xtc/.trr): ")
residues_input = input("Ingrese los nombres de los residuos separados por coma (ej: POPC,DPPC): ")
output_name = input("Ingrese el nombre base para archivos de salida (sin extensión): ")
n_bins = int(input("Ingrese el número de bins a lo largo de z (ej: 100): "))

residues = [r.strip() for r in residues_input.split(",")]

# ------------------------------
# Cargar la simulación
# ------------------------------
u = mda.Universe(gro_file, traj_file)

# ------------------------------
# Definir números de electrones por tipo de átomo
electron_dict = {
    'H': 1,
    'C': 6,
    'N': 7,
    'O': 8,
    'P': 15,
    'S': 16
}

# ------------------------------
# Preparar bins
z_min = 0.0
z_max = u.dimensions[2]  # límite en z
bin_edges = np.linspace(z_min, z_max, n_bins+1)
bin_centers = 0.5*(bin_edges[:-1] + bin_edges[1:])

# ------------------------------
# Diccionario para almacenar todos los perfiles
all_profiles = {}

# ------------------------------
# Función para calcular densidad electrónica de un grupo de átomos
def calculate_density(atomgroup):
    density_sum = np.zeros(n_bins)
    n_frames = 0
    for ts in u.trajectory:
        n_frames += 1
        z_coords = []
        electron_counts = []
        for atom in atomgroup:
            elem = atom.name[0].upper()
            electron_counts.append(electron_dict.get(elem, 0))
            z_coords.append(atom.position[2])
        hist, _ = np.histogram(z_coords, bins=bin_edges, weights=electron_counts)
        density_sum += hist
    density_avg = density_sum / n_frames
    area = u.dimensions[0] * u.dimensions[1]
    bin_width = bin_edges[1] - bin_edges[0]
    density_avg = density_avg / (area * bin_width)
    return density_avg

# ------------------------------
# Calcular densidad para cada residuo
for resname in residues:
    print(f"Procesando residuo: {resname}")
    selection = u.select_atoms(f"resname {resname}")
    if len(selection) == 0:
        print(f"Advertencia: No se encontraron átomos para el residuo {resname}. Se omite.")
        continue
    density_avg = calculate_density(selection)
    all_profiles[resname] = density_avg

# ------------------------------
# Calcular densidad para todo el sistema
print("Procesando densidad de todo el sistema")
all_density = calculate_density(u.atoms)
all_profiles['ALL'] = all_density

# ------------------------------
# Crear DataFrame y guardar CSV
df_all = pd.DataFrame({'z (Å)': bin_centers})
for name, profile in all_profiles.items():
    df_all[name] = profile

csv_filename = f"{output_name}_all_residues.csv"
df_all.to_csv(csv_filename, index=False)
print(f"CSV con todos los perfiles guardado: {csv_filename}")

# ------------------------------
# Graficar todos los perfiles juntos con doble eje
fig, ax1 = plt.subplots(figsize=(8,5))

# Eje izquierdo -> perfil del sistema completo
ax1.plot(bin_centers, all_profiles['ALL'], color='black', linewidth=2, label='ALL')
ax1.set_xlabel('z (Å)')
ax1.set_ylabel('Electron Density (System) [e/Å³]', color='black')
ax1.tick_params(axis='y', labelcolor='black')

# Ajustar límites del eje izquierdo
ax1.set_ylim(min(all_profiles['ALL'])*0.95, max(all_profiles['ALL'])*1.05)

# Eje derecho -> perfiles de los residuos
ax2 = ax1.twinx()
for name, profile in all_profiles.items():
    if name == 'ALL':
        continue
    ax2.plot(bin_centers, profile, label=name, linestyle='--')

ax2.set_ylabel('Electron Density (Residues) [e/Å³]', color='blue')
ax2.tick_params(axis='y', labelcolor='blue')

# Ajustar límites del eje derecho
residue_profiles = [profile for name, profile in all_profiles.items() if name != 'ALL']
if residue_profiles:  # evitar error si no hay residuos
    all_res_vals = np.concatenate(residue_profiles)
    ax2.set_ylim(all_res_vals.min()*0.95, all_res_vals.max()*1.05)

# Combinar leyendas
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='best')

plt.title('Electron Density Profiles (System vs Residues)')
plt.tight_layout()
png_filename = f"{output_name}_all_residues.png"
plt.savefig(png_filename, dpi=300)
plt.close()

print(f"Gráfico con todos los perfiles guardado: {png_filename}")
print("Cálculo finalizado para todos los residuos y el sistema total.")
