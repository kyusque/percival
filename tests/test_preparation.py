# -*- coding: utf-8 -*-
from percival.preparation.domain.value_objects import Smiles
from percival.preparation.domain.molecule import MoleculeFactory, Molecule
from percival.preparation.service.input_generator import InputGenerator

import unittest

class PreparationTest(unittest.TestCase):
    """Basic test cases."""

    def test_mol_factory(self):
        smiles: Smiles = Smiles("CCC")
        molecule: Molecule = MoleculeFactory.create_molecule_from_smiles(smiles)
        res = molecule.mol.GetNumAtoms()
        self.assertEqual(11, res)

    def test_mol_conformer(self):
        smiles: Smiles = Smiles("C")
        molecule: Molecule = MoleculeFactory.create_molecule_from_smiles(smiles)
        confs = molecule.generate_conformers()
        lines = confs.to_sdf().split("\n")
        self.assertEqual("$$$$", lines[-1])

    def test_InputGenerator_gjf(self):
        smiles: Smiles = Smiles("C")
        molecule: Molecule = MoleculeFactory.create_molecule_from_smiles(smiles)
        confs = molecule.generate_conformers()
        lines = InputGenerator.generate_gjf(confs.to_list()[0], "test").split("\n")
        self.assertEqual("%chk=test", lines[0])


if __name__ == '__main__':
    unittest.main()
