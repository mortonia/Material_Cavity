# lj_calculator.py

from ase.calculators.calculator import Calculator, all_changes
import numpy as np

class CustomLennardJones(Calculator):
    implemented_properties = ['energy', 'forces']

    def __init__(self, lj_params, rc=6.0):
        super().__init__()
        self.lj_params = lj_params
        self.rc = rc

    def calculate(self, atoms=None, properties=['energy', 'forces'], system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        positions = atoms.get_positions()
        symbols = atoms.get_chemical_symbols()
        n_atoms = len(atoms)

        energy = 0.0
        forces = np.zeros((n_atoms, 3))

        for i in range(n_atoms):
            for j in range(i+1, n_atoms):
                r_vec = positions[j] - positions[i]
                r = np.linalg.norm(r_vec)

                if r > self.rc:
                    continue

                elem_i = symbols[i]
                elem_j = symbols[j]
                sigma_ij = 0.5 * (self.lj_params[elem_i]['sigma'] + self.lj_params[elem_j]['sigma'])
                epsilon_ij = np.sqrt(self.lj_params[elem_i]['epsilon'] * self.lj_params[elem_j]['epsilon'])

                sr6 = (sigma_ij / r) ** 6
                sr12 = sr6 ** 2
                energy += 4 * epsilon_ij * (sr12 - sr6)

                force_mag = 24 * epsilon_ij / r * (2 * sr12 - sr6)
                force_vec = force_mag * r_vec / r
                forces[i] -= force_vec
                forces[j] += force_vec

        self.results['energy'] = energy
        self.results['forces'] = forces

