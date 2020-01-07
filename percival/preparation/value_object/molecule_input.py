from abc import abstractmethod
from rdkit import Chem
from rdkit.Chem import rdchem


class MoleculeInput(object):

    @abstractmethod
    def generate_mol(self) -> rdchem.Mol:
        raise NotImplementedError()


class Smiles(MoleculeInput):
    value: str

    def __init__(self, smiles: str):
        self.value: str = smiles

    def generate_mol(self) -> rdchem.Mol:
        mol: rdchem.RWMol = Chem.MolFromSmiles(self.value)
        mol = Chem.AddHs(mol)
        return mol

