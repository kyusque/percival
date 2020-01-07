class CalculationMethod:
    __level: str

    def __init__(self, level: str):
        self.__level = level

    @property
    def level(self):
        return self.__level