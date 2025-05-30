# molecule_utils.py

import numpy as np
from ase.io import read

def load_and_center_molecule(filename, center):
    mol = read(filename)
    mol_center = mol.get_positions().mean(axis=0)
    mol.translate(center - mol_center)
    return mol

def compute_cavity_radii(molecule, buffer):
    positions = molecule.get_positions()
    min_pos = positions.min(axis=0)
    max_pos = positions.max(axis=0)
    size = max_pos - min_pos
    return size[0]/2 + buffer, size[1]/2 + buffer, size[2]/2 + buffer

