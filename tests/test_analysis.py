# -*- coding: utf-8 -*-
from percival.analysis.domain.table import Table
from percival.analysis.domain.gaussian_log import GaussianLog
from percival.analysis.domain.entry import EntryRepository

import os

import unittest


class AnalysisTest(unittest.TestCase):
    """Basic test cases."""

    def test_log_parser(self):
        file_path = os.path.join(os.path.dirname(__file__), "resources/entries/entry1/entry1_confs_0.out")
        gaussian_log = GaussianLog.parse_file(file_path)
        conformer = gaussian_log.get_conformer()

        self.assertEqual("b3lyp", conformer.electronic_structure.method)
        self.assertEqual(-604627.161725122, conformer.electronic_structure.energy.hartree)
        self.assertEqual(-0.20125, conformer.electronic_structure.molecular_orbitals.homo().energy.hartree)
        self.assertEqual(-0.24548, conformer.electronic_structure.molecular_orbitals.homo(1).energy.hartree)
        self.assertEqual(0.00585, conformer.electronic_structure.molecular_orbitals.lumo().energy.hartree)
        self.assertEqual(0.0135, conformer.electronic_structure.molecular_orbitals.lumo(1).energy.hartree)

    def test_log_parser_for_multiple_path(self):
        file_paths = [os.path.join(os.path.dirname(__file__), "resources/entries/entry1/entry1_confs_{}.out".format(i))
                      for i in range(2)]
        gaussian_logs = GaussianLog.parse_files(file_paths)
        conformers = EntryRepository.read_gaussian_logs(gaussian_logs)

        self.assertEqual("b3lyp", conformers.to_list()[0].electronic_structure.method)
        self.assertEqual(-604627.161725122, conformers.to_list()[0].electronic_structure.energy.hartree)
        self.assertEqual("b3lyp", conformers.to_list()[1].electronic_structure.method)
        self.assertEqual(-604626.8425123101, conformers.to_list()[1].electronic_structure.energy.hartree)
        # minimum in two files
        self.assertEqual(-604627.161725122, conformers.minimum_energy.hartree)

    def test_create_entry(self):
        dir_path = os.path.join(os.path.dirname(__file__), "resources/entries/entry1")
        entry = EntryRepository.create_entry(dir_path)
        # minimum in a directory
        self.assertEqual(-604627.161725122, entry.minimum_energy.hartree)

    def test_create_entries(self):
        dir_paths = [os.path.join(os.path.dirname(__file__), "resources/entries/entry{}".format(i)) for i in
                     range(1, 4)]
        entries = EntryRepository.get_entries(dir_paths)
        # minimum in all entries
        self.assertEqual(-604627.161725122, entries.minimum_energy.hartree)

    def test_table(self):
        dir_paths = [os.path.join(os.path.dirname(__file__), "resources/entries/entry{}".format(i)) for i in
                     range(1, 4)]
        entries = EntryRepository.get_entries(dir_paths)
        table = Table(entries)
        data = table.to_dataframe()
        self.assertEqual('Nc1nc2c(ncn2[C@H]2C[C@H](O)[C@@H](CO)O2)c(=O)[nH]1', data["smiles"][0])


if __name__ == '__main__':
    unittest.main()
