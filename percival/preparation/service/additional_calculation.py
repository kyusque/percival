from percival.additional_calculation.domain.value_objects import GaussianLog, CalculationMethod
import os
import re


class AdditionalCalculation:
    gaussian_log: GaussianLog
    calculation_method: CalculationMethod

    def __init__(self, gaussian_log: GaussianLog, calculation_method: CalculationMethod):
        self.gaussian_log = gaussian_log
        self.calculation_method = calculation_method

    def generate_gaussian_input(self, filename: str) -> None:
        mol = self.gaussian_log.get_mol()
        result = mol.write(format='gjf')
        result = re.sub("#Put Keywords Here, check Charge and Multiplicity.",
                        "%chk=" + os.path.splitext(os.path.basename(filename))[
                            0] + ".chk\n%nprocshared=4\n%mem=260MB\n#" + self.calculation_method.level,
                        result)
        with open(filename, "w") as file:
            file.write(result)


