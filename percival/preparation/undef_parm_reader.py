from percival.preparation.domain.mm_parm_reader import ParmReader
import argparse
from importlib.resources import path
from . import src
import re


class UndefParmReader:
    def __init__(self, path):
        self.file = open(path, "r")
        self.atom_set = set()
        self.atom_set_iter = None
        self.memo = set()
        self.regexps = {
            'str': re.compile("Bondstretch undefined between atoms\s+"
                              "[0-9]+\s+[0-9]+\s+"
                              "(.+?)\s*-(.+?)\s+"),
            'bend': re.compile("Angle bend\s+undefined between atoms\s+"
                               "[0-9]+\s+[0-9]+\s+[0-9]+\s+"
                               "(.+?)\s*-(.+?)\s*-(.+?)\s+")}
        self.sort_dispatcher = {
            "str": lambda l: sorted(l),
            "bend": lambda l: l if l[0] <= l[2] else list(reversed(l))
        }

    def __iter__(self):
        return self

    def __next__(self):
        for line in self.file:
            for key, regexp in self.regexps.items():
                res = regexp.search(line)
                if res:
                    atoms = res.groups()
                    atoms = self.sort_dispatcher[key](atoms)
                    parm_id = "-".join(atoms)
                    # Suppress duplicate results
                    if parm_id in self.memo:
                        continue
                    # For VDW parameters
                    for atom in atoms:
                        self.atom_set.add(atom)
                    self.memo.add(res)
                    return parm_id
        # Return VDW keys
        if not self.atom_set_iter:
            self.atom_set_iter = iter(self.atom_set)
        for i in self.atom_set_iter:
            return i
        raise StopIteration


if __name__ == '__main__':
    with path(src, "gaff2.dat") as file:
        data = dict([i for i in ParmReader(file)])
    p = argparse.ArgumentParser()
    p.add_argument("path")
    p.add_argument("--prefix")
    a = p.parse_args()
    reader = UndefParmReader(a.path)
    for key in reader:
        try:
            print(data[key.replace(a.prefix, "")].print(a.prefix))
        except:
            print(key)
