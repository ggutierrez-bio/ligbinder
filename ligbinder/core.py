import os
from typing import List, Optional
import logging
import pytraj
import yaml
from parmed.amber import AmberParm
from parmed.tools.actions import HMassRepartition
from ligbinder.settings import SETTINGS
from ligbinder.tree import Node, Tree
from ligbinder.md import AmberMDEngine


logger = logging.getLogger(__file__)


class LigBinder:
    def __init__(self, path: str = ".", config_file: Optional[str] = None) -> None:
        self.path = path
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
        self.setup_hmr()
        if len(self.tree.nodes) == 0:
            logger.info("No root node found. Instantiating...")
            self.tree.create_root_node(**SETTINGS["data_files"])
        while not self.tree.has_converged() and self.tree.can_grow():
            node: Node = self.tree.create_node_from_candidate()
            logger.info(f"New node chosen for expansion. Current depth: {node.depth}")
            engine = AmberMDEngine(node.path, **SETTINGS["md"])
            engine.run()
            node.calc_node_rmsd()
            parent_rmsd = self.tree.nodes[node.parent_id].rmsd
            if node.rmsd < parent_rmsd:
                logger.info(f"Node {node.rmsd} improved rmsd by {parent_rmsd - node.rmsd}! current rmsd: {node.rmsd}")
        logger.info("Exploration finished.")
        self.compile_results()

    def compile_results(self):

        path = self.tree.path
        report_dir = os.path.join(path, SETTINGS["results"]["report_dir"]) 

        def _create_report_dir():
            os.makedirs(report_dir, exist_ok=True)

        def _concat_trajectory(indices: List[int]):
            # get filenames
            traj_files = [os.path.join(path, f"node_{index}", SETTINGS["md"]["trj_file"]) for index in indices]
            top_file = os.path.join(path, SETTINGS["data_files"]["top_file"])
            ref_file = os.path.join(path, SETTINGS["data_files"]["ref_file"])
            full_traj_file = os.path.join(report_dir, SETTINGS["results"]["trj_file"])
            
            # load, align write
            traj = pytraj.iterload(traj_files, top=top_file)
            ref = pytraj.load(ref_file, top=top_file)
            mask = SETTINGS["system"]["protein_mask"]
            pytraj.rmsd(traj, mask=mask, ref=ref)
            pytraj.write_traj(full_traj_file, traj)

        def _write_node_list_file(indices: List[int]):
            node_list_file = os.path.join(report_dir, SETTINGS["results"]["idx_file"])
            with open(node_list_file, 'w') as idx_file:
                idx_file.writelines(indices)

        def _write_rmsd_file(indices: List[int]):
            rmsd_file = os.path.join(report_dir, SETTINGS["results"]["rms_file"])
            rmsds = [self.tree.nodes[index].rmsd for index in indices]
            with open(rmsd_file, 'w') as rms_file:
                rms_file.writelines(rmsds)

        def _write_stats(indices: List[int]):
            stats_filename = os.path.join(report_dir, SETTINGS["results"]["stats_file"])
            report = {
                "converged": self.tree.has_converged(),
                "total_nodes": len(self.tree.nodes),
                "max_depth": max([node.depth for node in self.tree.nodes.values()]),
                "best_rmsd": min([node.rmsd for node in self.tree.nodes.values()]),
            }
            with open(stats_filename, 'w') as stats_file:
                yaml.dump(report, stats_file)

        node_ids = self.tree.get_solution_path()
        if self.tree.has_converged():
            logger.warning("SUCCESS: LIGAND BOUND!!!")

            _create_report_dir()
            _concat_trajectory(node_ids)
            _write_node_list_file(node_ids)
            _write_rmsd_file(node_ids)
        else:
            logger.warning("FAILURE: UNABLE TO BIND")

        logger.info("writing report")
        _write_stats(node_ids)

    def setup_hmr(self):
        if not SETTINGS["md"]["use_hmr"]:
            return
        top_file = SETTINGS["data_files"]["top_file"]
        logger.info(f"Applying HMR on topology file {top_file}")
        parm = AmberParm(top_file)
        HMassRepartition(parm).execute()
        parm.write_parm(top_file)
        logger.info(f"HMR applied")
