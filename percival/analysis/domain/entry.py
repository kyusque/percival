from typing import List
import glob
import re

from ..domain.conformer import Conformer, Conformers
from percival.analysis.domain.gaussian_log import GaussianLog, GaussianLogs, Energy


class Entry:
    conformers: Conformers
    directory_path: str

    def __init__(self, conformers: Conformers, directory_path: str):
        self.conformers = conformers
        self.directory_path = directory_path

    @property
    def minimum_energy(self) -> Energy:
        return self.conformers.minimum_energy

    @property
    def minimum_energy_conformer(self) -> Conformer:
        return self.conformers.minimum_energy_conformer


class EntryRepository:
    entries: List[Entry]

    def __init__(self, entries: List[Entry]):
        self.entries = entries

    def to_list(self) -> List[Entry]:
        return self.entries

    @staticmethod
    def read_gaussian_logs(gaussian_logs: GaussianLogs) -> Conformers:
        conformers = [log.get_conformer() for log in gaussian_logs.to_list()]
        return Conformers(conformers)

    @staticmethod
    def get_gaussian_logs_from_directory(dir_path: str) -> GaussianLogs:
        file_candidates = [path for path in glob.glob("{}/*".format(dir_path)) if re.search(r"\.(log|out)$", path)]
        return GaussianLog.parse_files(file_candidates)

    @staticmethod
    def create_entry(directory_path: str) -> Entry:
        gaussian_logs = EntryRepository.get_gaussian_logs_from_directory(directory_path)
        conformers = EntryRepository.read_gaussian_logs(gaussian_logs)
        return Entry(conformers, directory_path)

    @property
    def minimum_energy(self) -> Energy:
        self.to_list().sort(key=lambda x: x.minimum_energy.hartree)
        return self.to_list()[0].minimum_energy

    @staticmethod
    def get_entries(dir_paths: List[str]):
        entries = [EntryRepository.create_entry(path) for path in dir_paths]
        return EntryRepository(entries)
