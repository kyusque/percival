from rdkit.Chem import rdchem
from typing import List
from percival.analysis.domain.value_objects import ElectronicStructure, Energy


class Conformer:
    mol_with_structure: rdchem.RWMol
    electronic_structure: ElectronicStructure

    def __init__(self, mol: rdchem.RWMol, electronic_struture: ElectronicStructure):
        self.mol_with_structure = mol
        self.electronic_structure = electronic_struture


class Conformers:
    conformer_list: List[Conformer]

    def __init__(self, conformer_list: List[Conformer]):
        self.conformer_list = conformer_list

    def to_list(self):
        return self.conformer_list

    @property
    def minimum_energy_conformer(self) -> Conformer:
        return sorted(self.conformer_list, key=lambda x: x.electronic_structure.energy.hartree)[0]

    @property
    def minimum_energy(self) -> Energy:
        return self.minimum_energy_conformer.electronic_structure.energy
