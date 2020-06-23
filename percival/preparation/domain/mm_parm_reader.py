import abc
import os
import re

class MM_param(metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def print(self):
        pass

    @abc.abstractmethod
    def id(self):
        pass

class VDW(MM_param):
    atom: str
    radius: float
    well_depth: float

    def __init__(self, atom: str, radius: float, well_depth: float):
        self.atom = atom.strip().upper()
        self.radius = radius
        self.well_depth = well_depth

    def print(self):
        return f"VDW {self.atom} {self.radius} {self.well_depth}"

    def id(self):
        return self.atom

class Stretch(MM_param):
    atom1: str
    atom2: str

    force_constant: float
    length: float

    def __init__(self, atom1, atom2, force_constant, length):
        self.atom1 = atom1.strip().upper()
        self.atom2 = atom2.strip().upper()
        self.force_constant = force_constant
        self.length = length
        self.atom1, self.atom2 = sorted([self.atom1, self.atom2])

    def print(self):
        return f"HrmStr1 {self.atom1} {self.atom2} {self.force_constant} {self.length}"

    def id(self):
        return f"{self.atom1}-{self.atom2}"

class Bend(MM_param):
    atom1: str
    atom2: str
    atom3: str
    force_constant: float
    angle: float

    def __init__(self, atom1, atom2, atom3, force_constant, angle):
        self.atom1 = atom1.strip().upper()
        self.atom2 = atom2.strip().upper()
        self.atom3 = atom3.strip().upper()
        self.force_constant = force_constant
        self.angle = angle
        atoms = [self.atom1, self.atom2, self.atom3]
        atoms = atoms if atoms[0] <= atoms[2] else list(reversed(atoms))
        self.atom1, self.atom2, self.atom3 = atoms

    def print(self):
        return f"HrmBnd1 {self.atom1} {self.atom2} {self.atom3} {self.force_constant} {self.angle}"

    def id(self):
        return f"{self.atom1}-{self.atom2}-{self.atom3}"


class ParmReader():
    def __init__(self, path):
        self.file = open(path, "r")
        for line in self.file:
            if len(line) < 5:
                break
        self.regexes = {
            "vdw": re.compile(
                "^\s*([a-z].)\s+([0-9]+?\.[0-9]+?)\s+([0-9]+?\.[0-9]+?)\s+"
            ),
            "str": re.compile(
            "^\s*([a-z].)-([a-z].)\s+([0-9]+?\.[0-9]+?)\s+([0-9]+?\.[0-9]+?)\s+"
            ),
            "bend": re.compile(
            "^\s*([a-z].)-([a-z].)-([a-z].)\s+([0-9]+?\.[0-9]+?)\s+([0-9]+?\.[0-9]+?)\s+"
            )
        }
        self.dispatch = {
            "vdw": lambda l: VDW(l[0], l[1], l[2]),
            "str": lambda l: Stretch(l[0], l[1], l[2], l[3]),
            "bend": lambda l: Bend(l[0], l[1], l[2], l[3], l[4])
        }

    def __iter__(self):
        return self

    def __next__(self):
        for line in self.file:
            i = 0
            for k, v in self.regexes.items():
                i += 1
                res = v.search(line)
                if res:
                    res = res.groups()
                    res = self.dispatch[k](res)
                    return [res.id(), res.print()]
        raise StopIteration


if __name__ == '__main__':
    for i in ParmReader(os.path.dirname(__file__) + "/../service/src/gaff2.dat"):
        print(i)
