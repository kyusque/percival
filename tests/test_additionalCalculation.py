from unittest import TestCase
from percival.preparation.value_object.molecule_input import Smiles
from percival.preparation.domain.molecule import MoleculeFactory, Molecule
from percival.preparation.service.input_generator import InputGenerator


class TestAdditionalCalculation(TestCase):
    def test_generate_gaussian_input(self):
        smiles: Smiles = Smiles("C")
        molecule: Molecule = MoleculeFactory.create_molecule(smiles)
        confs = molecule.generate_conformers()

        generated_lines = InputGenerator.generate_gjf(confs.to_list()[0], "test").split("\n")
        self.assertEqual("%chk=test", generated_lines[0])
