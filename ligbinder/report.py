import os
import logging
from typing import List
import pytraj
import yaml
from ligbinder.settings import SETTINGS
from ligbinder.tree import Tree
import math

logger = logging.getLogger(__name__)


class Reporter:
    def __init__(self, tree: Tree) -> None:
        self.tree = tree
        self.path = tree.path
        self.report_dir = os.path.join(self.path, SETTINGS["results"]["report_dir"])

    def _create_report_dir(self):
        os.makedirs(self.report_dir, exist_ok=True)

    def _concat_trajectory(self, indices: List[int]):
        # get filenames
        traj_files = [
            os.path.join(self.path, f"node_{index}", SETTINGS["md"]["trj_file"])
            for index in indices if index != 0
        ]
        top_file = os.path.join(self.path, SETTINGS["data_files"]["top_file"])
        ref_file = os.path.join(self.path, SETTINGS["data_files"]["ref_file"])
        ref_top_file = os.path.join(self.path, SETTINGS["data_files"]["ref_top_file"])
        full_traj_file = os.path.join(self.report_dir, SETTINGS["results"]["trj_file"])

        # load, align write
        load_mask = SETTINGS["system"]["load_mask"]
        traj = pytraj.iterload(traj_files, top=top_file, mask=load_mask)
        ref = pytraj.load(ref_file, top=ref_top_file, mask=load_mask)
        mask = SETTINGS["system"]["protein_mask"]
        pytraj.rmsd(traj, mask=mask, ref=ref)
        pytraj.write_traj(full_traj_file, traj)

    def _write_node_list_file(self, indices: List[int]):
        node_list_file = os.path.join(self.report_dir, SETTINGS["results"]["idx_file"])
        with open(node_list_file, "w") as idx_file:
            idx_file.write("\n".join(map(str, indices))+"\n")

    def _write_rmsd_file(self, indices: List[int]):
        rmsd_file = os.path.join(self.report_dir, SETTINGS["results"]["rms_file"])
        rmsds = [self.tree.nodes[index].rmsd for index in indices]
        with open(rmsd_file, "w") as rms_file:
            rms_file.write("\n".join(map(str, rmsds))+"\n")

    def _write_stats(self):
        stats_filename = os.path.join(
            self.report_dir, SETTINGS["results"]["stats_file"]
        )

        best_node = self.tree.get_best_node()
        report = {
            "converged": self.tree.has_converged(),
            "total_nodes": len(self.tree.nodes),
            "max_depth": max([node.depth for node in self.tree.nodes.values()]),
            "best_node": {
                "node_id": best_node.node_id,
                "rmsd": best_node.rmsd,
                "nrmsd": best_node.nrmsd,
                "pBP": -math.log10(self.tree.get_biasing_power(best_node)),
                "path": self.tree.get_path_to_node(best_node),
            }
        }
        with open(stats_filename, "w") as stats_file:
            yaml.dump(report, stats_file)

    def compile_results(self):
        node_ids = self.tree.get_solution_path()
        self._create_report_dir()
        if self.tree.has_converged():
            logger.warning("SUCCESS: LIGAND BOUND!!!")
            if SETTINGS["results"]["join_trajectories"]:
                self._concat_trajectory(node_ids)
            self._write_node_list_file(node_ids)
            self._write_rmsd_file(node_ids)
        else:
            logger.warning("FAILURE: UNABLE TO BIND")

        logger.info("writing report")
        self._write_stats()
