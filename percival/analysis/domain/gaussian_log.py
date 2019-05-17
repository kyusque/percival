from .value_objects import MolecularOrbital, MolecularOrbitals, ElectronicStructure, Energy
from .conformer import Conformer

from typing import List
import pandas as pd
import pybel
import os
import tempfile
from rdkit.Chem import AllChem


class GaussianLog:
    path: str
    lines_before_entering_link1: List[str]
    list_link1: List[List[str]]

    def __init__(self, path: str, lines_before_entering_link1: List[str], list_link1: List[List[str]]):
        self.path = path
        self.lines_before_entering_link1 = lines_before_entering_link1
        self.list_link1 = list_link1

    def get_molecular_orbitals(self):
        lines = iter(self.list_link1[-1])
        mo_eigen_list = []
        for line in lines:
            if line.startswith(' The electronic state is '):
                break
        for line in lines:
            if line.find('eigenvalues --') == -1:
                break
            block = line.split(' eigenvalues --  ')
            for item in block[1].split():
                mo_info = block[0].split()
                mo_eigen_list.append([mo_info[0], mo_info[1], float(item)])

        data = pd.DataFrame(mo_eigen_list)
        data.columns = ['type', 'electron', 'eigen']
        mos = [MolecularOrbital(Energy(eigen), electron == "occ.") for electron, eigen in
               zip(data["electron"], data["eigen"])]
        return MolecularOrbitals(mos)

    def get_conformer(self) -> Conformer:
        with tempfile.TemporaryDirectory() as tempdir:
            last_job_log = os.path.join(tempdir, "last.log")
            with open(last_job_log, "w") as file:
                for line in self.lines_before_entering_link1:
                    file.write(line)
                for line in self.list_link1[-1]:
                    file.write(line)

            mols = [mol for mol in pybel.readfile('g09', last_job_log)]
            mol = mols[0]
            data = mol.data

            electronic_structure = ElectronicStructure(data["program"], data["basis"], data["method"],
                                                       Energy(mol.energy),
                                                       self.get_molecular_orbitals())

            temp_file = os.path.join(tempdir, "temp.sdf")
            with open(temp_file, "w") as fp:
                fp.write(mol.write("sdf"))
            sdf_reader = AllChem.SDMolSupplier(temp_file)
            mols = [mol for mol in sdf_reader]
            mol = mols[0]
        return Conformer(mol, electronic_structure)

    @staticmethod
    def parse_file(path: str):
        chunks = []

        with open(path, "r") as file:
            i = 0
            chunks.append([])
            for line in file:
                if "Entering Link 1" in line:
                    i += 1
                    chunks.append([])
                chunks[i].append(line)

        return GaussianLog(path, chunks.pop(0), chunks)

    @staticmethod
    def parse_files(paths: List[str]):
        result = [GaussianLog.parse_file(path) for path in paths]
        return GaussianLogs(result)


class GaussianLogs:
    gaussian_logs: List[GaussianLog]

    def __init__(self, gaussian_logs: List[GaussianLog]):
        self.gaussian_logs = gaussian_logs

    def to_list(self):
        return self.gaussian_logs
