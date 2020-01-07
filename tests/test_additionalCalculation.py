from unittest import TestCase
from percival.preparation.domain.value_objects import Smiles
from percival.preparation.domain.molecule import MoleculeFactory, Molecule
from percival.preparation.service.input_generator import InputGenerator


class TestAdditionalCalculation(TestCase):
    def test_generate_gaussian_input(self):
        smiles = Smiles("C")
        molecule: Molecule = MoleculeFactory.create_molecule_from_smiles(smiles)
        confs = molecule.generate_conformers()
        lines = InputGenerator.generate_gjf(confs.to_list()[0], "test").split("\n")
        self.assertEqual("%chk=test", lines[0])
