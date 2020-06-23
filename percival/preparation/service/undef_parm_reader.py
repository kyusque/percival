from ..domain.mm_parm_reader import ParmReader, VDW, Bend, Stretch
import argparse
import os
import re

class UndefParmReader():
    def __init__(self, path):
        self.file = open(path, "r")
        self.atom_set = set()
        self.atom_set_iter = None
        self.results = set()
        self.regex_str = re.compile(
            "Bondstretch undefined between atoms\s+[0-9]+\s+[0-9]+\s+(.+?)\s*-(.+?)\s+"
        )
        self.regex_bend = re.compile(
            "Angle bend\s+undefined between atoms\s+[0-9]+\s+[0-9]+\s+[0-9]+\s+(.+?)\s*-(.+?)\s*-(.+?)\s+"
        )

    def __iter__(self):
        return self

    def __next__(self):
        for line in self.file:
            res = self.regex_str.search(line)
            if res is not None:
                res = sorted(res.groups())
                for atom in res:
                    self.atom_set.add(atom)
                res = "{}-{}".format(res[0], res[1])
                if res in self.results:
                    continue
                self.results.add(res)
                return res
            res = self.regex_bend.search(line)
            if res is not None:
                res = res.groups()
                for atom in res:
                    self.atom_set.add(atom)
                res = res if res[0] <= res[2] else list(reversed(res))
                res = "{}-{}-{}".format(res[0], res[1], res[2])
                if res in self.results:
                    continue
                self.results.add(res)
                return res
        if self.atom_set_iter is None:
            self.atom_set_iter = iter(self.atom_set)

        for i in self.atom_set_iter:
            return i

        raise StopIteration

if __name__ == '__main__':
    data = dict([i for i in ParmReader(os.path.dirname(__file__) + "/src/gaff.dat")])
    p = argparse.ArgumentParser()
    p.add_argument("path")
    p.add_argument("--prefix")
    a = p.parse_args()
    reader = UndefParmReader(a.path)
    for a in reader:
        try:
            print(data[a])
        except:
            pass

