# main.py

import os
import sys
import shutil
import numpy as np
from ase.io import write
from ase.constraints import FixAtoms

from config import ELEMENT, REPEAT, LATTICE_CONSTANT, CAVITY_BUFFER, DISPLACEMENT_THRESHOLD
from crystal_utils import create_crystal, carve_elliptical_cavity
from molecule_utils import load_and_center_molecule, compute_cavity_radii
from relaxation import relax_with_lj
from output_utils import print_cell_dimensions, print_crystal_info

def main():
    # Parse args
    if len(sys.argv) > 1:
        molecule_filename = sys.argv[1]
    else:
        molecule_filename = 'water.xyz'  # Default molecule

    # Optional manual cavity center
    if len(sys.argv) > 4:
        try:
            cavity_center = np.array([float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])])
            print(f"Using manual cavity center: {cavity_center}")
        except ValueError:
            print("Invalid cavity center coordinates provided; using crystal center instead.")
            cavity_center = None
    else:
        cavity_center = None

    # Setup output dir
    BASE_OUTPUT_DIR = "output"
    molecule_name = os.path.splitext(os.path.basename(molecule_filename))[0]
    OUTPUT_DIR = os.path.join(BASE_OUTPUT_DIR, molecule_name)

    if os.path.exists(OUTPUT_DIR):
        shutil.rmtree(OUTPUT_DIR)
    os.makedirs(OUTPUT_DIR)

    # Create crystal
    crystal = create_crystal(ELEMENT, REPEAT, LATTICE_CONSTANT)

    # Print and save cell dimensions and cavity suggestions
    print_cell_dimensions(crystal, output_dir=OUTPUT_DIR)
    print_crystal_info(crystal, label="Original", output_dir=OUTPUT_DIR)

    n_atoms_original = len(crystal)
    print(f"Number of atoms in original crystal: {n_atoms_original}")
    with open(os.path.join(OUTPUT_DIR, "info.txt"), "w") as f:
        f.write(f"Number of atoms in original crystal: {n_atoms_original}\n")

    # Determine cavity center
    if cavity_center is None:
        cavity_center = crystal.get_cell().sum(axis=0) / 2

    # Load and center molecule
    molecule = load_and_center_molecule(molecule_filename, cavity_center)

    # Compute cavity radii (molecule size + buffer)
    radius_x, radius_y, radius_z = compute_cavity_radii(molecule, CAVITY_BUFFER)

    # Carve cavity in crystal
    crystal_cavity = carve_elliptical_cavity(crystal, cavity_center, radius_x, radius_y, radius_z)
    print_crystal_info(crystal_cavity, label="After Cavity", output_dir=OUTPUT_DIR)

    n_atoms_after_cavity = len(crystal_cavity)
    print(f"Number of atoms after cavity carving: {n_atoms_after_cavity}")
    with open(os.path.join(OUTPUT_DIR, "info.txt"), "a") as f:
        f.write(f"Number of atoms after cavity carving: {n_atoms_after_cavity}\n")

    # Combine molecule + cavity crystal
    combined = crystal_cavity + molecule

    # Freeze molecule atoms (last atoms)
    num_mol_atoms = len(molecule)
    molecule_indices = list(range(len(combined) - num_mol_atoms, len(combined)))
    combined.set_constraint(FixAtoms(indices=molecule_indices))

    n_atoms_after_molecule = len(combined)
    print(f"Number of atoms after adding molecule: {n_atoms_after_molecule}")
    with open(os.path.join(OUTPUT_DIR, "info.txt"), "a") as f:
        f.write(f"Number of atoms after adding molecule: {n_atoms_after_molecule}\n")

    original_positions = combined.get_positions().copy()

    # Save structure before relaxation
    write(os.path.join(OUTPUT_DIR, "before_lj.xyz"), combined)

    # Relax structure with custom LJ potential
    relaxed, displacements, num_displaced = relax_with_lj(
        combined, original_positions, DISPLACEMENT_THRESHOLD, OUTPUT_DIR
    )
    print(f"[LJ] Displaced atoms > {DISPLACEMENT_THRESHOLD} Å: {num_displaced}")
    with open(os.path.join(OUTPUT_DIR, "info.txt"), "a") as f:
        f.write(f"[LJ] Displaced atoms > {DISPLACEMENT_THRESHOLD} Å: {num_displaced}\n")

    # Tag displaced atoms for ASE visualization
    tags = np.zeros(len(relaxed), dtype=int)
    tags[displacements > DISPLACEMENT_THRESHOLD] = 1
    relaxed.set_tags(tags)

    # Save relaxed structure
    write(os.path.join(OUTPUT_DIR, "after_lj.xyz"), relaxed)

if __name__ == "__main__":
    main()

