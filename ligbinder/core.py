from ligbinder.report import Reporter
import os
from typing import Optional
import logging
from parmed.amber import AmberParm
from parmed.tools.actions import HMassRepartition
from ligbinder.settings import SETTINGS
from ligbinder.tree import Node, Tree
from ligbinder.md import AmberMDEngine
from ligbinder import VERSION

logger = logging.getLogger(__name__)


class LigBinder:
    def __init__(self, path: str = ".", config_file: Optional[str] = None) -> None:
        self.path = path
        if config_file is not None:
            SETTINGS.update_settings_with_file(self.get_config_file(config_file))
        self.tree = Tree(self.path, **SETTINGS["tree"])

    def get_config_file(self, config_file: Optional[str] = None) -> Optional[str]:
        local_default_config_file = os.path.join(self.path, "config.yml")
        if config_file is not None and os.path.exists(config_file):
            return config_file
        elif os.path.exists(local_default_config_file):
            return local_default_config_file
        return None

    def run(self):
        self.log_initial_info()
        metric_name = "nrmsd" if self.tree.use_normalized_rmsd else "rmsd"
        self.setup_hmr()
        if len(self.tree.nodes) == 0:
            logger.info("No root node found. Instantiating...")
            self.tree.create_root_node(**SETTINGS["data_files"])
        while not self.tree.has_converged() and self.tree.can_grow():
            node: Node = self.tree.create_node_from_candidate()
            logger.info("New node chosen for expansion.")
            logger.info(f"\tdepth: {node.depth}/{self.tree.max_depth}")
            parent_node = self.tree.nodes[node.parent_id]
            parent_metric = self.tree.get_metric(parent_node)
            logger.info(f"\t{metric_name}: {parent_metric:.3f}")
            engine = AmberMDEngine(node.path, **SETTINGS["md"])
            engine.run()
            node.calc_node_rmsd()
            node_metric = self.tree.get_metric(node)
            if node_metric < parent_metric:
                logger.info(
                    f"Node {node.node_id} improved {metric_name} by {parent_metric - node_metric:.3f}!"
                    f" current {metric_name}: {node_metric:.3f}"
                )
        logger.info("Exploration finished.")
        Reporter(self.tree).compile_results()

    def setup_hmr(self):
        if not SETTINGS["md"]["use_hmr"]:
            return
        top_file = SETTINGS["data_files"]["top_file"]
        logger.info(f"Applying HMR on topology file {top_file}")
        parm = AmberParm(top_file)
        HMassRepartition(parm).execute()
        parm.write_parm(top_file)
        logger.info("HMR applied")

    def log_initial_info(self):
        logger.info(f"Running ligbinder v{VERSION}")
        logger.info(f"Current directory: {os.path.abspath(os.getcwd())}")
        logger.info("Full settings for this run:")
        logger.info(SETTINGS.data)
