from typing import List
import tempfile
import os
import pybel


class GaussianLog:
    __path: str
    __lines_before_entering_link1: List[str]
    __list_link1: List[List[str]]

    def __init__(self, path: str, lines_before_entering_link1: List[str], list_link1: List[List[str]]):
        self.__path = path
        self.lines_before_entering_link1 = lines_before_entering_link1
        self.list_link1 = list_link1

    @property
    def path(self) -> str:
        return self.__path

    def get_mol(self) -> pybel.Molecule:
        with tempfile.TemporaryDirectory() as tempdir:
            last_job_log = os.path.join(tempdir, "last.log")
            with open(last_job_log, "w") as file:
                for line in self.lines_before_entering_link1:
                    file.write(line)
                for line in self.list_link1[-1]:
                    file.write(line)

            mols = [mol for mol in pybel.readfile('g09', last_job_log)]
            mol = mols[0]
        return mol

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


class CalculationMethod:
    __level: str

    def __init__(self, level: str):
        self.__level = level

    @property
    def level(self):
        return self.__level

