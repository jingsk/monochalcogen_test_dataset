import numpy as np
from ase.io import write
from ase import Atoms


def add_inversion(atoms, inv_c):
    inv_c = np.array(inv_c)
    scaled_pos = 2* inv_c - atoms.get_scaled_positions()
    inverted_atoms = Atoms(
    atoms.get_chemical_formula(),
    scaled_positions=scaled_pos,
    cell = atoms.cell,
    pbc =True
    )
    return atoms + inverted_atoms


a = 4.1534
b = a
c = 10.9821
h = 0.25

scaled_pos = [
    [0.25, 0, 0.25 - h / 2], # M
    [0.75, 0.5, 0.25 + h / 2], # M
    [0.25, 0, 0.25 + h / 2], # X
    [0.75, 0.5, 0.25 - h / 2], # X
]


atoms = Atoms(
    'Ge2Se2',
    scaled_positions=scaled_pos,
    cell = [a,b,c],
    pbc =True
    )

full_layers = add_inversion(atoms, [0.5, 0.5, 0.5])
write('pristine_GeSe.vasp', full_layers, format='vasp', sort=True)
