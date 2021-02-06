import os
import yaml
from typing import Optional
import samplematic


class Settings:

    def __init__(self, filename: Optional[str] = None) -> None:
        self.filename = self.get_defaults_filename()
        self.data: dict = self.load_data(self.filename)
        self.get = self.data.get
        self.__getitem__ = self.data.__getitem__
        if filename is not None:
            self.data.update(self.load_data(filename))
            self.filename = filename

    @staticmethod
    def get_defaults_filename():
        samplematic_home = os.getenv("SAMPLEMATIC_HOME")
        path = samplematic_home if samplematic_home is not None else os.path.join(samplematic.__path__[0], "data")
        return os.path.join(path, "defaults", "config.yml")

    @staticmethod
    def load_data(filename: str) -> dict:
        with open(filename, 'r') as file:
            d = yaml.load(file, Loader=yaml.FullLoader)
        return d if d is not None else {}
