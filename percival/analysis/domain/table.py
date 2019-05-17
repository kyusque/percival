from .entry import EntryRepository
import pandas as pd
from typing import List
from rdkit.Chem.rdchem import Mol
from rdkit import Chem


class Table:
    entries: EntryRepository

    def __init__(self, entries: EntryRepository):
        self.entries = entries

    def to_dataframe(self):
        name = pd.Series([entry.directory_path.split("/")[-1] for entry in self.entries.to_list()], dtype=str)
        lumo = pd.Series(
            [entry.conformers.minimum_energy_conformer.electronic_structure.molecular_orbitals.lumo().energy.eV for
             entry in
             self.entries.to_list()], dtype=float)
        lumo1 = pd.Series(
            [entry.conformers.minimum_energy_conformer.electronic_structure.molecular_orbitals.lumo(1).energy.eV for
             entry in
             self.entries.to_list()], dtype=float)
        homo = pd.Series(
            [entry.conformers.minimum_energy_conformer.electronic_structure.molecular_orbitals.homo().energy.eV for
             entry in
             self.entries.to_list()], dtype=float)
        smiles = pd.Series([Chem.MolToSmiles(entry.conformers.minimum_energy_conformer.mol_with_structure) for entry in
                            self.entries.to_list()], dtype=str)

        data = pd.DataFrame({'name': name, 'lumo+1(eV)': lumo1, 'lumo(eV)': lumo, "homo(eV)": homo, "smiles": smiles})
        return data

    def to_csv(self):
        return self.to_dataframe().to_csv()

    def get_entry_representative_mols(self) -> List[Mol]:
        return [entry.minimum_energy_conformer.mol_with_structure for entry in self.entries.to_list()]

    def get_entry_smiles(self) -> List[str]:
        return [Chem.MolToSmiles(mol) for mol in self.get_entry_representative_mols()]
