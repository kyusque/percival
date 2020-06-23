import pybel
import re
import tempfile
import os
from percival.preparation.domain.conformer import Conformer, Conformers
from percival.preparation.domain.molecule import Molecule, MoleculeFactory
from percival.preparation.domain.value_objects import Smiles

class InputGenerator:
    @staticmethod
    def generate_gjf_directory(smiles: Smiles, directory_name:str) -> None:
        os.makedirs(directory_name, exist_ok=True)

        molecule: Molecule = MoleculeFactory.create_molecule_from_smiles(smiles, directory_name)
        confs = molecule.generate_conformers()
        for i, conformer in enumerate(confs.to_list()):
            with open('{}/conf{}.gjf'.format(directory_name, i), 'w') as outfile:
                outfile.write(InputGenerator.generate_gjf(conformer, 'conf{}.gjf'))
                outfile.close()


    @staticmethod
    def generate_gjf(conformer: Conformer, chk_filename: str) -> str:
        with tempfile.TemporaryDirectory() as tmpdir:
            temp_file = os.path.join(tmpdir, 'temp.sdf')
            with open(temp_file, 'w') as fp:
                fp.write(conformer.to_sdf())
            mols = [mol for mol in pybel.readfile("sdf", temp_file)]
            mol = mols[0]

        result = mol.write(format='gjf')
        result = re.sub("#Put Keywords Here, check Charge and Multiplicity.",
                        "%chk=" + chk_filename + "\n%nprocshared=4\n%mem=260MB\n#opt hf/3-21g",
                        result)
        result = result + "--Link1--\n%chk=" + chk_filename + \
                 "\n%nprocshared=4\n%mem=260MB\n#opt b3lyp/6-31g(d) Geom=AllCheck Guess=Read\n\n"
        return result
