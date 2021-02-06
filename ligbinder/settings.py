import os
import yaml
from typing import Optional
import sys
import site
import logging


logger = logging.getLogger(__file__)

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
        ligbinder_home = os.getenv("LIGBINDER_HOME")
        home_candidates = [ligbinder_home]
        home_candidates += [os.path.join(sys.prefix, "ligbinder")]
        home_candidates += [os.path.join(path, "ligbinder") for path in site.PREFIXES]
        home_candidates = [home for home in home_candidates if home is not None and os.path.exists(home)]
        configs = [os.path.join(path, "default_config.yml") for path in home_candidates]
        config = [config for config in configs if os.path.exists(config)][0]
        logger.info(f"Using {config} for default settings")
        return config

    @staticmethod
    def load_data(filename: str) -> dict:
        with open(filename, 'r') as file:
            d = yaml.load(file, Loader=yaml.FullLoader)
        return d if d is not None else {}
