# config.py

ELEMENT = 'C'
REPEAT = (5, 5, 5)
LATTICE_CONSTANT = 3.57  # Å
CAVITY_BUFFER = 2.0
DISPLACEMENT_THRESHOLD = 0.1

# Lennard-Jones parameters per element (eV, Å)
LJ_PARAMS = {
    'C': {'epsilon': 0.00284, 'sigma': 3.40},
    'H': {'epsilon': 0.00236, 'sigma': 2.96},
    'O': {'epsilon': 0.00674, 'sigma': 3.00},
    'N': {'epsilon': 0.00739, 'sigma': 3.20},
    'Ar': {'epsilon': 0.0103,  'sigma': 3.40},
}

