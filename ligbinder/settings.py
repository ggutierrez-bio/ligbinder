import os
import yaml
from typing import Optional
import sys
import site
import logging
from collections.abc import Mapping


logger = logging.getLogger(__name__)


class Settings:
    def __init__(self, filename: Optional[str] = None) -> None:
        super().__init__()
        self.filename = self.get_defaults_filename()
        self.data: dict = self.load_data(self.filename)
        self.get = self.data.get
        if filename is not None:
            self.update_settings_with_file(filename)

    def __getitem__(self, k):
        return self.data[k]

    def __len__(self):
        return len(self.data)

    def __contains__(self, k):
        return k in self.data

    def update_settings_with_file(self, filename: str):
        self.merge(self.load_data(filename))
        self.filename = filename

    @staticmethod
    def get_defaults_filename():
        ligbinder_home = os.getenv("LIGBINDER_HOME")
        paths = [sys.prefix, os.path.join(sys.prefix, "local")]
        paths += site.PREFIXES
        paths += ["/usr", "/usr/local"]
        home_candidates = [ligbinder_home]
        home_candidates += [os.path.join(path, "ligbinder") for path in paths] + ["./ligbinder/data"]
        logger.info(f"Searching for ligbinder home in: {home_candidates}")
        home_candidates = [
            home
            for home in home_candidates
            if home is not None and os.path.exists(home)
        ]
        configs = [os.path.join(path, "default_config.yml") for path in home_candidates]
        logger.info(f"found configs: {configs}")
        config = next(iter([config for config in configs if os.path.exists(config)]), "config.yml")
        logger.info(f"Using {config} for default settings")
        return config

    @staticmethod
    def load_data(filename: str) -> dict:
        with open(filename, "r") as file:
            d = yaml.load(file, Loader=yaml.FullLoader)
        return d if d is not None else {}

    def merge(self, data):
        def _merge(dict1, dict2):
            for key, value in dict2.items():
                if (
                    key in dict1
                    and isinstance(dict1[key], dict)
                    and isinstance(value, Mapping)
                ):
                    _merge(dict1[key], value)
                elif (
                    key in dict1
                    and isinstance(dict1[key], list)
                    and isinstance(value, list)
                ):
                    dict1[key] += value
                else:
                    dict1[key] = dict2[key]

        _merge(self.data, data)


SETTINGS = Settings()
