import os
import shutil
import sys
import numpy as np
from ase.build import bulk
from ase.io import read, write
from ase.constraints import FixAtoms
from ase.calculators.lj import LennardJones
from ase.optimize import BFGS

# --- CONFIGURATION ---
ELEMENT = 'C'  # Carbon crystal
REPEAT = (3, 3, 3)
LATTICE_CONSTANT = 3.57  # Angstrom for Carbon diamond cubic approx.
CAVITY_BUFFER = 2.0  # extra buffer around molecule (Å)
DISPLACEMENT_THRESHOLD = 0.1  # Å

# --- CRYSTAL UTILS ---

def create_crystal(element, repeat, lattice_constant):
    crystal = bulk(name=element, crystalstructure='fcc', a=lattice_constant)
    crystal = crystal.repeat(repeat)
    return crystal

def carve_elliptical_cavity(crystal, center, radius_x, radius_y, radius_z):
    new_atoms = crystal.copy()
    positions = new_atoms.get_positions()
    mask = []
    for pos in positions:
        val = ((pos[0] - center[0]) / radius_x) ** 2 + \
              ((pos[1] - center[1]) / radius_y) ** 2 + \
              ((pos[2] - center[2]) / radius_z) ** 2
        mask.append(val > 1.0)
    return new_atoms[mask]

def print_cell_dimensions(crystal, output_dir=None):
    cell = crystal.get_cell()
    max_x = np.linalg.norm(cell[0])
    max_y = np.linalg.norm(cell[1])
    max_z = np.linalg.norm(cell[2])

    msg = (
        f"Crystal cell dimensions (Å):\n"
        f"  x: {max_x:.3f}\n"
        f"  y: {max_y:.3f}\n"
        f"  z: {max_z:.3f}\n"
        f"Suggested cavity centers:\n"
        f"  center:        {[max_x/2, max_y/2, max_z/2]}\n"
        f"  bottom-right:  {[max_x*0.9, max_y*0.1, max_z/2]}\n"
        f"  top-left:      {[max_x*0.1, max_y*0.9, max_z/2]}\n"
    )
    print(msg)
    if output_dir:
        with open(f"{output_dir}/cell_dimensions.txt", "w") as f:
            f.write(msg)

# --- MOLECULE UTILS ---

def load_and_center_molecule(filename, center):
    mol = read(filename)
    mol_positions = mol.get_positions()
    mol_center = mol_positions.mean(axis=0)
    shift = center - mol_center
    mol.translate(shift)
    return mol

def compute_cavity_radii(molecule, buffer):
    positions = molecule.get_positions()
    min_pos = positions.min(axis=0)
    max_pos = positions.max(axis=0)
    size = max_pos - min_pos
    radius_x = size[0]/2 + buffer
    radius_y = size[1]/2 + buffer
    radius_z = size[2]/2 + buffer
    return radius_x, radius_y, radius_z

# --- LJ RELAXATION ---

def relax_with_lj(atoms, original_positions, displacement_threshold, output_dir):
    lj = LennardJones()
    atoms.set_calculator(lj)

    dyn = BFGS(atoms, logfile=os.path.join(output_dir, "relax.log"))
    dyn.run(fmax=0.01, steps=100)

    new_positions = atoms.get_positions()
    displacements = np.linalg.norm(new_positions - original_positions, axis=1)
    num_displaced = np.sum(displacements > displacement_threshold)

    with open(os.path.join(output_dir, "lj_displacements.txt"), "w") as f:
        for i, disp in enumerate(displacements):
            f.write(f"{i} {disp:.4f}\n")

    return atoms, displacements, num_displaced

# --- OUTPUT UTILS ---

def print_crystal_info(crystal, label="", output_dir=None):
    cell = crystal.get_cell()
    size = crystal.get_positions().ptp(axis=0)
    info = (f"{label} crystal info:\n"
            f"  Number of atoms: {len(crystal)}\n"
            f"  Cell vectors:\n{cell}\n"
            f"  Atomic positions range (Å):\n  x: {size[0]:.3f}, y: {size[1]:.3f}, z: {size[2]:.3f}\n")
    print(info)
    if output_dir:
        with open(os.path.join(output_dir, f"{label.lower().replace(' ', '_')}_info.txt"), "w") as f:
            f.write(info)

# --- MAIN SCRIPT ---

def main():
    if len(sys.argv) > 1:
        molecule_filename = sys.argv[1]
    else:
        molecule_filename = 'water.xyz'  # default molecule

    # Optional manual cavity center coordinates: python main.py water.xyz x y z
    if len(sys.argv) > 4:
        try:
            cavity_center = np.array([float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])])
            print(f"Using manual cavity center: {cavity_center}")
        except ValueError:
            print("Invalid cavity center coordinates provided; using crystal center instead.")
            cavity_center = None
    else:
        cavity_center = None

    BASE_OUTPUT_DIR = "output"
    molecule_name = os.path.splitext(os.path.basename(molecule_filename))[0]
    OUTPUT_DIR = os.path.join(BASE_OUTPUT_DIR, molecule_name)

    if os.path.exists(OUTPUT_DIR):
        shutil.rmtree(OUTPUT_DIR)
    os.makedirs(OUTPUT_DIR)

    # Create crystal
    crystal = create_crystal(ELEMENT, REPEAT, LATTICE_CONSTANT)

    # Print and save cell dimensions + cavity suggestions
    print_cell_dimensions(crystal, output_dir=OUTPUT_DIR)
    print_crystal_info(crystal, label="Original", output_dir=OUTPUT_DIR)

    n_atoms_original = len(crystal)
    print(f"Number of atoms in original crystal: {n_atoms_original}")
    with open(os.path.join(OUTPUT_DIR, "info.txt"), "w") as f:
        f.write(f"Number of atoms in original crystal: {n_atoms_original}\n")

    # Choose cavity center
    if cavity_center is None:
        cavity_center = crystal.get_cell().sum(axis=0) / 2

    # Load molecule and center it at cavity center
    molecule = load_and_center_molecule(molecule_filename, cavity_center)

    # Compute cavity radii based on molecule size + buffer
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

    # Save before relaxation
    write(os.path.join(OUTPUT_DIR, "before_lj.xyz"), combined)

    # Relax with Lennard-Jones potential
    relaxed, displacements, num_displaced = relax_with_lj(
        combined, original_positions, DISPLACEMENT_THRESHOLD, OUTPUT_DIR
    )
    print(f"[LJ] Displaced atoms > {DISPLACEMENT_THRESHOLD} Å: {num_displaced}")
    with open(os.path.join(OUTPUT_DIR, "info.txt"), "a") as f:
        f.write(f"[LJ] Displaced atoms > {DISPLACEMENT_THRESHOLD} Å: {num_displaced}\n")

    # Tag displaced atoms for ASE GUI visualization
    tags = np.zeros(len(relaxed), dtype=int)
    tags[displacements > DISPLACEMENT_THRESHOLD] = 1
    relaxed.set_tags(tags)

    # Save after relaxation
    write(os.path.join(OUTPUT_DIR, "after_lj.xyz"), relaxed)

if __name__ == "__main__":
    main()

