import os
from typing import Optional
from parmed.amber import AmberParm
from parmed.tools.actions import HMassRepartition
from samplematic.settings import Settings
from samplematic.tree import Node, Tree
from samplematic.md import AmberMDEngine


class LigBinder():
    def __init__(self, path: str = ".", config_file: Optional[str] = None) -> None:
        self.path = path
        self.settings = Settings(self.get_config_file(config_file))
        self.tree = Tree(self.path, **self.settings.get("tree"))

    def get_config_file(self, config_file: Optional[str] = None) -> Optional[str]:
        local_default_config_file = os.path.join(self.path, "config.yml")
        if config_file is not None and os.path.exists(config_file):
            return config_file
        elif os.path.exists(local_default_config_file):
            return local_default_config_file
        return None

    def run(self):
        self.setup_hmr()
        if len(self.tree.nodes) == 0:
            self.tree.create_root_node(**self.settings.get("data_files"))
        while(not self.tree.has_converged() and self.tree.can_grow()):
            node: Node = self.tree.create_node_from_candidate()
            engine = AmberMDEngine(node.path, **self.settings["md"])
            engine.run()
            node.calc_node_rmsd()
    
    def setup_hmr(self):
        if not self.settings["md"]["use_hmr"]:
            return
        top_file = self.settings["data_files"]["top_file"]
        parm = AmberParm(top_file)
        HMassRepartition(parm).execute()
        parm.write_parm(top_file)
