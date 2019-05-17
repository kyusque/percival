from rdkit.Chem import rdchem
from rdkit import Chem
import io
from typing import List


class Conformer:
    mol_with_structure: rdchem.RWMol

    def __init__(self, mol: rdchem.RWMol):
        self.mol_with_structure = mol

    def to_sdf(self):
        output = io.StringIO()
        sdf_writer = Chem.SDWriter(output)
        sdf_writer.write(self.mol_with_structure)
        sdf_writer.close()
        return output.getvalue().strip()


class Conformers:
    conformer_list: List[Conformer]

    def __init__(self, conformer_list: List[Conformer]):
        self.conformer_list = conformer_list

    def to_list(self):
        return self.conformer_list

    def to_sdf(self):
        output = io.StringIO()
        sdf_writer = Chem.SDWriter(output)
        for conformer in self.conformer_list:
            sdf_writer.write(conformer.mol_with_structure)
        sdf_writer.close()
        return output.getvalue().strip()
