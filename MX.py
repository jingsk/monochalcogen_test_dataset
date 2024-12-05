import numpy as np
from ase.io import write
from ase import Atoms
from numpy.typing import NDArray

class MX:
    def __init__(
            self, 
            M: str, 
            X: str,
            cell: NDArray= [4.1534,4.1534,10.9821],
            h = 0.25,
            dh = 0): #For GeSe dh = 0.02279
        self.M = M
        self.X = X
        self.cell = cell
        #self.ab_ratio = ab_ratio
        self._h = h
        self._dh = dh
        self._atoms = self.update_atoms()
    
    @property
    def h(self):
        return self._h
    
    @property
    def dh(self):
        return self._dh
    
    # @property
    # def atoms(self):
    #     return self._atoms
    
    def update_atoms(self):
        atoms = Atoms(
            f'{self.M}2{self.X}2',
            scaled_positions=self.update_scaled_positions(),
            cell = self.cell,
            pbc =True
        )
        return atoms
        
    def update_scaled_positions(self):
        scaled_pos = [
            [0.25, 0, 0.25 - self.h / 2], # M
            [0.75, 0.5, 0.25 + self.h / 2], # M
            [0.25, 0, 0.25 + self.h / 2 - self.dh], # X
            [0.75, 0.5, 0.25 - self.h / 2 + self.dh], # X
        ]
        return np.array(scaled_pos)
    
    def add_inversion(self, inv_c = [0.5, 0.5, 0.5]):
        inv_c = np.array(inv_c)
        scaled_pos = 2* inv_c - self._atoms.get_scaled_positions()
        inverted_atoms = Atoms(
            self._atoms.get_chemical_formula(),
            scaled_positions=scaled_pos,
            cell = self.cell,
            pbc =True
        )
        return self._atoms + inverted_atoms
    
    def toatoms(self):
        return self.add_inversion()