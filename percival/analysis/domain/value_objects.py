from typing import List


class Smiles:
    value: str

    def __init__(self, smiles: str):
        """

        :type smiles: str
        """
        self.value: str = smiles


class Energy:
    hartree_energy: float

    # https://www.hpc.co.jp/chem/software/gaussian/help/keywords/constants/

    def __init__(self, hartree_energy: float):
        self.hartree_energy = hartree_energy

    @property
    def hartree(self):
        return self.hartree_energy

    @property
    def eV(self):
        return self.hartree_energy * 27.2114

    @property
    def kcal_mol(self):
        return self.hartree_energy * 627.5095

    @property
    def kJ_mol(self):
        return self.kcal_mol * 4.184


class MolecularOrbital:
    energy: Energy
    occupied: bool

    def __init__(self, energy: Energy, occupied: bool):
        self.energy = energy
        self.occupied = occupied


class MolecularOrbitals:
    molecular_orbitals: List[MolecularOrbital]

    def __init__(self, molecular_orbitals: List[MolecularOrbital]):
        self.molecular_orbitals = molecular_orbitals

    def homo(self, minus: int = 0):
        occupied_orbitals = [mo for mo in self.molecular_orbitals if mo.occupied]
        occupied_orbitals.sort(key=lambda x: x.energy.hartree, reverse=True)
        return occupied_orbitals[minus]

    def lumo(self, plus: int = 0):
        unoccupied_orbitals = [mo for mo in self.molecular_orbitals if not mo.occupied]
        unoccupied_orbitals.sort(key=lambda x: x.energy.hartree)
        return unoccupied_orbitals[plus]


class ElectronicStructure:
    program: str
    basis: str
    method: str
    energy: Energy
    molecular_orbitals: MolecularOrbitals

    def __init__(self, program: str, basis: str, method: str, energy: Energy,
                 molecular_orbitals: MolecularOrbitals = None):
        self.program = program
        self.basis = basis
        self.method = method
        self.energy = energy
        self.molecular_orbitals = molecular_orbitals
