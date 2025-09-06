#!/usr/bin/env python3
"""
electron_density_weighted.py

Perfiles de densidad electrónica a lo largo de Z, centrados en la membrana (COM de P* si existe),
memoria-optimizado: acumula histogramas por frame y pondera cada átomo por su número de electrones.

Uso ejemplo:
python electron_density_weighted.py -g system.gro -x traj.xtc -r POPC,CHOL --bins 300 --sigma 2 --plot

Opciones:
--sigma : sigma en *bins* (como gaussian_filter1d)
--sigma_nm : alternativa más intuitiva, sigma en nm (convertido internamente)
--center-by : "P" (por defecto) o "all"
"""
import argparse
import MDAnalysis as mda
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt
import sys

# Números de electrones por elemento (solo los más relevantes)
ELEC = {
    "H": 1, "C": 6, "N": 7, "O": 8, "P": 15, "S": 16,
    "CL": 17, "NA": 11, "K": 19
}

def infer_element(atom):
    """Intentar obtener elemento; si no, inferir del nombre del átomo"""
    el = getattr(atom, "element", None)
    if el and isinstance(el, str) and el.strip():
        return el.strip().upper()
    # inferir desde atom.name
    name = atom.name.strip().upper()
    letters = "".join([c for c in name if c.isalpha()])
    if not letters:
        return None
    if len(letters) >= 2 and letters[:2] in ELEC:
        return letters[:2]
    return letters[0]

def build_electron_array(universe):
    elems = []
    unknown = set()
    for atom in universe.atoms:
        el = infer_element(atom)
        if el is None:
            elems.append(0)
        else:
            val = ELEC.get(el.upper(), None)
            if val is None:
                val = ELEC.get(el[0].upper(), 0)
                unknown.add(el)
            elems.append(val)
    if unknown:
        print(f"⚠️ Warning: elementos no conocidos (asumidos según primera letra o 0): {sorted(list(unknown))}", file=sys.stderr)
    return np.array(elems, dtype=float)

def sigma_nm_to_bins(sigma_nm, z_min, z_max, bins):
    bin_width = (z_max - z_min) / float(bins)
    if bin_width <= 0:
        return 1.0
    return max(1.0, sigma_nm / bin_width)

def main():
    p = argparse.ArgumentParser(description="Electron density profile (weighted by electrons), memory-optimized")
    p.add_argument("-g", "--gro", required=True, help=".gro topology")
    p.add_argument("-x", "--xtc", required=True, help=".xtc trajectory")
    p.add_argument("-r", "--resnames", required=False, default="", help="Residue names comma separated (ej: POPC,CHOL)")
    p.add_argument("-o", "--output", default="electron_density.csv", help="CSV output")
    p.add_argument("--plot", action="store_true", help="Guardar figura electron_density.png")
    p.add_argument("--bins", type=int, default=300, help="Número de bins en Z (resolución)")
    p.add_argument("--sigma", type=float, default=None, help="Sigma en bins para suavizado")
    p.add_argument("--sigma_nm", type=float, default=None, help="Sigma en nm (sobreescribe --sigma si se da)")
    p.add_argument("--center-by", choices=["P", "all"], default="P",
                   help='Cómo calcular el centro de la membrana por frame: "P" (COM de átomos P* si existen) o "all"')
    args = p.parse_args()

    resnames = [r.strip() for r in args.resnames.split(",") if r.strip()]
    u = mda.Universe(args.gro, args.xtc)
    n_atoms = len(u.atoms)
    print(f"Cargado sistema: {n_atoms} átomos, {len(u.trajectory)} frames")

    electron_per_atom = build_electron_array(u)

    res_indices = {}
    for r in resnames:
        sel = u.select_atoms(f"resname {r}")
        res_indices[r] = sel.indices if len(sel) > 0 else np.array([], dtype=int)

    # ---- Selección para centrar ----
    if args.center_by == "P":
        selP = u.select_atoms("name P*")
        if len(selP) == 0:
            print("⚠️ No se encontraron átomos con name P*; se usará COM de todos los átomos.", file=sys.stderr)
            center_selection = None
        else:
            center_selection = selP.indices
    else:
        center_selection = None

    # ---- Primera pasada: rango Z ----
    z_min_global = np.inf
    z_max_global = -np.inf
    frames_count = 0
    for ts in u.trajectory:
        frames_count += 1
        if center_selection is None:
            z_center = u.atoms.center_of_mass()[2]
        else:
            z_center = u.atoms[center_selection].center_of_mass()[2]
        z_coords = u.atoms.positions[:, 2] - z_center
        z_min_global = min(z_min_global, z_coords.min())
        z_max_global = max(z_max_global, z_coords.max())
    if frames_count == 0:
        raise RuntimeError("La trayectoria no tiene frames o no se pudo leer.")
    print(f"Rango Z centrado: {z_min_global:.3f} .. {z_max_global:.3f} Å")

    bins = args.bins
    edges = np.linspace(z_min_global, z_max_global, bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])

    if args.sigma_nm is not None:
        sigma_bins = sigma_nm_to_bins(args.sigma_nm, z_min_global, z_max_global, bins)
        print(f"Sigma pedido: {args.sigma_nm} nm -> {sigma_bins:.2f} bins")
    elif args.sigma is not None:
        sigma_bins = args.sigma
    else:
        sigma_bins = 2.0
    print(f"Usando sigma = {sigma_bins} bins")

    histograms = {r: np.zeros(bins, dtype=float) for r in resnames}
    histograms["SYSTEM"] = np.zeros(bins, dtype=float)
    counts = 0

    for ts in u.trajectory:
        if center_selection is None:
            z_center = u.atoms.center_of_mass()[2]
        else:
            z_center = u.atoms[center_selection].center_of_mass()[2]
        z_coords = u.atoms.positions[:, 2] - z_center

        hist_sys, _ = np.histogram(z_coords, bins=edges, weights=electron_per_atom)
        histograms["SYSTEM"] += hist_sys

        for r in resnames:
            idx = res_indices.get(r, None)
            if idx is None or len(idx) == 0:
                continue
            z_r = z_coords[idx]
            w_r = electron_per_atom[idx]
            hist_r, _ = np.histogram(z_r, bins=edges, weights=w_r)
            histograms[r] += hist_r

        counts += 1

    bin_width = centers[1] - centers[0] if bins > 1 else (z_max_global - z_min_global)
    density_profiles = {}
    for k, h in histograms.items():
        density = (h / counts) / bin_width
        density_sm = gaussian_filter1d(density, sigma=sigma_bins)
        density_profiles[k] = density_sm

    df = pd.DataFrame({"z (Å)": centers})
    for k, v in density_profiles.items():
        df[k] = v

    extrema = []
    for k, v in density_profiles.items():
        i_max = int(np.argmax(v))
        i_min = int(np.argmin(v))
        extrema.append([k, centers[i_max], float(v[i_max]), centers[i_min], float(v[i_min])])
    extrema_df = pd.DataFrame(extrema, columns=["Profile", "z_max (Å)", "max", "z_min (Å)", "min"])

    with open(args.output, "w") as f:
        df.to_csv(f, index=False)
        f.write("\n# Maxima and minima (Profile, z_max (Å), max, z_min (Å), min)\n")
        extrema_df.to_csv(f, index=False)
    print(f"Guardado: {args.output}")

    if args.plot:
        fig, ax1 = plt.subplots(figsize=(7, 4))
        for r in resnames:
            if r in density_profiles:
                ax1.plot(centers, density_profiles[r], label=r)
        ax1.set_xlabel("Z position (Å)")
        ax1.set_ylabel("Residue electron density (e⁻ / Å)")
        ax1.legend(loc="upper left")

        ax2 = ax1.twinx()
        ax2.plot(centers, density_profiles["SYSTEM"], linestyle="--", label="SYSTEM")
        ax2.set_ylabel("System electron density (electrons / Å)")
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper right", fontsize="small")
        fig.tight_layout()
        plt.savefig("electron_density.png", dpi=300)
        print("Guardado plot: electron_density.png")

if __name__ == "__main__":
    main()

