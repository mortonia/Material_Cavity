# crystal_utils.py

import numpy as np
from ase.build import bulk

def create_crystal(element, repeat, lattice_constant):
    crystal = bulk(name=element, crystalstructure='fcc', a=lattice_constant)
    return crystal.repeat(repeat)

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

