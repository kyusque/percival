import os

from percival.preparation.service.additional_calculation import AdditionalCalculation


class AdditionalCalculationDirectory:
    name: str
    additional_calculation: AdditionalCalculation

    def __init__(self, name: str, additional_calculation: AdditionalCalculation):
        self.name = name
        self.additional_calculation = additional_calculation

    @property
    def directory_path(self) -> str:
        parent_directory = os.path.dirname(self.additional_calculation.gaussian_log.path)
        return os.path.join(parent_directory, self.name)

    def create_directory(self) -> None:
        os.makedirs(self.directory_path, exist_ok=True)