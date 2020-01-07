from rdkit.Chem import rdchem, AllChem
from rdkit import Chem

from .conformer import Conformer, Conformers
from percival.preparation.value_object.input import Smiles


class Molecule:
    mol: rdchem.RWMol
    title: str

    def __init__(self, mol: rdchem.RWMol, title: str):
        self.mol = mol
        self.title = title

    def generate_conformers(self, num_confs=100, prune_rms_thresh=2) -> Conformers:
        conf_ids = AllChem.EmbedMultipleConfs(self.mol, numConfs=num_confs, pruneRmsThresh=prune_rms_thresh)
        prop = AllChem.MMFFGetMoleculeProperties(self.mol)

        conformers = []
        for conf_id in conf_ids:
            temp_mol: rdchem.Mol = Chem.Mol(self.mol)
            temp_mol.RemoveAllConformers()
            temp_mol.AddConformer(self.mol.GetConformer(conf_id))

            mmff = AllChem.MMFFGetMoleculeForceField(self.mol, prop, confId=conf_id)
            mmff.Minimize()
            energy = mmff.CalcEnergy()

            temp_mol.SetProp("_Name", "Conformer_" + str(conf_id) + " at " + str(
                energy) + " MMFF energy generated from " + self.title)
            temp_mol.SetProp("MMFF_Energy", str(energy))

            conformer = Conformer(temp_mol)
            conformers.append(conformer)

        return Conformers(conformers)


class MoleculeFactory:
    def create_molecule_from_smiles(smiles: Smiles, title = "no title") -> Molecule:
        mol: rdchem.RWMol = Chem.MolFromSmiles(smiles.value)
        mol = Chem.AddHs(mol)
        return Molecule(mol, title)
