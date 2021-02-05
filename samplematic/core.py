import os
from typing import Optional
from samplematic.settings import Settings
from samplematic.tree import Tree


class BigChemExplorer():
    def __init__(self, path: str = ".", config_file: Optional[str] = None) -> None:
        self.path = path
        self.settings = Settings(self.get_config_file(config_file))

    def get_config_file(self, config_file: Optional[str] = None) -> Optional[str]:
        local_default_config_file = os.path.join(self.path, "data", "config.yml")
        if config_file is not None and os.path.exists(config_file):
            return config_file
        elif os.path.exists(local_default_config_file):
            return local_default_config_file
        return None


def get_tree(path: str = ".") -> Tree:
    return Tree(path)


def should_keep_running(path: str = ".") -> bool:
    # tree = get_tree(path)
    return False


def run_samplematic(crd: Optional[str] = None, top: Optional[str] = None, path: str = "."):

    while should_keep_running(path):
        break
