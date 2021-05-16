from glob import glob
import os
import yaml
from typing import Dict, List, Optional
import pytraj
import logging
import math

from ligbinder.settings import SETTINGS

logger = logging.getLogger(__name__)


class Node:
    def __init__(
        self, path: str, id: int, parent_id: int, rmsd: float = None, nrmsd: float = None, depth: int = None
    ) -> None:
        self.path = path
        self.node_id = id
        self.parent_id = parent_id
        self.rmsd = rmsd
        self.nrmsd = nrmsd
        self.depth = depth
        self.children: List[Node] = []

    @staticmethod
    def load_node(path: str) -> "Node":
        logger.info(f"loading node from {path}")
        filename = f"{path}/.info.yml"
        config: dict = Node.load_node_info(filename)
        id = Node.get_id_from_path(path)
        parent_id = config.get("parent_id", 0)
        rmsd = config.get("rmsd")
        return Node(path, id, parent_id, rmsd)

    @staticmethod
    def get_id_from_path(path: str):
        sep = os.path.sep
        node_rel_path = path.rstrip(sep).split(sep)[-1]
        return int(node_rel_path[len("node_"):])

    @staticmethod
    def load_node_info(filename) -> dict:
        with open(filename, "r") as cfg_file:
            return yaml.load(cfg_file, Loader=yaml.FullLoader)

    def write_node_info(self) -> None:
        filename = os.path.join(self.path, ".info.yml")
        with open(filename, "w") as file:
            yaml.dump({"parent_id": self.parent_id, "rmsd": self.rmsd, "nrmsd": self.nrmsd}, file)

    def get_relative_file(self, filename) -> str:
        return os.path.relpath(os.path.join(self.path, filename))

    def calc_node_rmsd(self) -> float:
        logger.info(f"calculating rmsd for node {self.node_id}")
        crd_file = self.get_relative_file(SETTINGS["md"]["rst_file"])
        top_file = self.get_relative_file(SETTINGS["md"]["top_file"])
        ref_file = self.get_relative_file(SETTINGS["md"]["ref_file"])
        if not all(os.path.exists(file) for file in [crd_file, top_file, ref_file]):
            logger.warning(
                f"Unable to calculate rmsd for node {self.node_id}. Some files are missing"
            )
        load_mask = SETTINGS["system"]["load_mask"]
        traj = pytraj.load(crd_file, top=top_file, mask=load_mask)
        # this one is causing some trouble in CI. need to inspect why, though.
        # traj = pytraj.autoimage(traj)
        ref_top_file = self.get_relative_file(os.path.join(os.path.pardir, SETTINGS["data_files"]["ref_top_file"]))
        ref = pytraj.load(ref_file, top=ref_top_file, mask=load_mask)
        align_mask = SETTINGS["system"]["protein_mask"]
        pytraj.rmsd(traj, mask=align_mask, ref=ref)
        ligand_mask = SETTINGS["system"]["ligand_mask"]
        rmsd = pytraj.rmsd_nofit(traj, mask=ligand_mask, ref=ref)
        self.rmsd = float(rmsd[0])
        self.nrmsd = self.rmsd / traj[ligand_mask].top.n_atoms
        self.write_node_info()
        return self.rmsd


class Tree:
    def __init__(
        self,
        path: str = ".",
        max_children: int = 5,
        max_depth: int = 15,
        tolerance: float = 0.5,
        max_nodes: int = 500,
        min_relative_improvement: float = 0.1,
        min_absolute_improvement: float = 0.1,
        max_steps_back: int = 5,
        use_normalized_rmsd: bool = True,
    ) -> None:
        self.path = path
        self.max_children = max_children
        self.max_depth = max_depth
        self.tolerance = tolerance
        self.use_normalized_rmsd = use_normalized_rmsd
        self.max_nodes = max_nodes
        self.min_relative_improvement = min_relative_improvement
        self.min_absolute_improvement = min_absolute_improvement
        self.max_steps_back = max_steps_back
        self.nodes: Dict[int, Node] = {}
        self.load_nodes()
        self.current_max_depth = max([node.depth for node in self.nodes.values()] + [0])

    def load_nodes(self) -> Dict[int, Node]:
        node_paths = glob(f"{self.path}/node_*/")
        nodes: List[Node] = [Node.load_node(node_path) for node_path in node_paths]
        self.nodes = {node.node_id: node for node in nodes}
        [
            self.nodes[node.parent_id].children.append(node)
            for node in nodes
            if node.parent_id is not None
        ]
        self.set_depths()
        return self.nodes

    def get_next_id(self) -> int:
        return max(list(self.nodes), default=0) + 1

    def get_metric(self, node: Node) -> float:
        return node.nrmsd if self.use_normalized_rmsd else node.rmsd

    def create_node(
        self,
        parent_id: int,
        ref_file: Optional[str] = None,
        top_file: Optional[str] = None,
    ) -> Node:
        id = self.get_next_id()
        node_path = os.path.join(self.path, f"node_{id}")
        os.mkdir(node_path)
        parent_path = os.path.join(os.path.pardir, f"node_{parent_id}")
        top_file = (
            top_file
            if top_file is not None
            else os.path.join(self.path, SETTINGS["data_files"]["top_file"])
        )
        rel_top_file = os.path.relpath(top_file, node_path)
        os.symlink(rel_top_file, os.path.join(node_path, SETTINGS["md"]["top_file"]))

        ref_file = (
            ref_file
            if ref_file is not None
            else os.path.join(self.path, SETTINGS["data_files"]["ref_file"])
        )
        rel_ref_file = os.path.relpath(ref_file, node_path)
        os.symlink(rel_ref_file, os.path.join(node_path, SETTINGS["md"]["ref_file"]))

        initial_frame = os.path.join(parent_path, SETTINGS["md"]["rst_file"])
        os.symlink(initial_frame, os.path.join(node_path, SETTINGS["md"]["crd_file"]))

        parent_node = self.nodes[parent_id]
        node = Node(node_path, id, parent_id, depth=parent_node.depth + 1)
        node.write_node_info()
        self.nodes[id] = node
        parent_node.children.append(node)
        if node.depth > self.current_max_depth:
            self.current_max_depth = node.depth
        return node

    def create_root_node(
        self,
        crd_file: Optional[str] = None,
        top_file: Optional[str] = None,
        ref_file: Optional[str] = None,
        ref_top_file: Optional[str] = None
    ) -> Node:
        node_id = 0
        node_path = os.path.relpath(os.path.join(self.path, f"node_{node_id}"))
        os.mkdir(node_path)

        crd_file = (
            crd_file
            if crd_file is not None
            else os.path.join(self.path, SETTINGS["data_files"]["crd_file"])
        )
        rel_crd_file = os.path.relpath(crd_file, node_path)
        os.symlink(rel_crd_file, os.path.join(node_path, SETTINGS["md"]["rst_file"]))

        top_file = (
            top_file
            if top_file is not None
            else os.path.join(self.path, SETTINGS["data_files"]["top_file"])
        )
        rel_top_file = os.path.relpath(top_file, node_path)
        os.symlink(rel_top_file, os.path.join(node_path, SETTINGS["md"]["top_file"]))

        ref_file = (
            ref_file
            if ref_file is not None
            else os.path.join(self.path, SETTINGS["data_files"]["ref_file"])
        )
        rel_ref_file = os.path.relpath(ref_file, node_path)
        os.symlink(rel_ref_file, os.path.join(node_path, SETTINGS["md"]["ref_file"]))

        node = Node(node_path, node_id, None)
        node.calc_node_rmsd()
        node.depth = 0
        self.nodes[0] = node
        return node

    def create_node_from_candidate(self) -> Node:
        parent = self.choose_candidate_node()
        node = self.create_node(
            parent.node_id,
            SETTINGS["data_files"]["ref_file"],
            SETTINGS["data_files"]["top_file"],
        )
        return node

    def choose_candidate_node(self) -> Node:
        candidate: Node = sorted(
            [node for node in self.nodes.values() if self.is_expandable(node)],
            key=self.get_metric,
        )[0]
        logger.info(f"candidate selected: {candidate.node_id}")
        return candidate

    def set_depths(self):
        if len(self.nodes) > 0:
            self.set_node_depth(self.nodes[0], 0)

    def set_node_depth(self, node: Node, current_depth: int) -> None:
        node.depth = current_depth
        for child in node.children:
            self.set_node_depth(child, current_depth + 1)

    def can_grow(self):
        can_grow = (
            len(self.nodes) < self.max_nodes
            and any(self.is_expandable(node) for node in list(self.nodes.values()))
        )
        if not can_grow:
            logger.warning("tree can't grow any bigger")
        return can_grow

    def is_expandable(self, node: Node) -> bool:
        parent_metric = (
            self.get_metric(self.nodes[node.parent_id]) if node.parent_id is not None else math.inf
        )
        node_depth_is_acceptable = self.max_depth > node.depth > (self.current_max_depth - self.max_steps_back)
        children_count_is_acceptable = len(node.children) < self.max_children
        min_acceptable_metric = max(
            parent_metric - self.min_absolute_improvement,
            parent_metric * (1 - self.min_relative_improvement)
        )
        metric_improvement_is_enough = self.get_metric(node) < min_acceptable_metric
        is_expandable = all([
            node_depth_is_acceptable, children_count_is_acceptable, metric_improvement_is_enough
        ])
        if not is_expandable:
            logger.debug(f"node {node.node_id} can't be expanded")
        return is_expandable

    def has_converged(self) -> bool:
        has_converged = any(
            self.get_metric(node) <= self.tolerance for node in list(self.nodes.values())
        )
        msg = "Tree converged " if has_converged else "Tree didn't converge "
        msg += f"after exploring {len(self.nodes)}/{self.max_nodes}"
        logger.info(msg)
        return has_converged

    def get_best_node(self) -> Node:
        return sorted([node for node in self.nodes.values()], key=self.get_metric)[0]

    def get_solution_path(self) -> List[int]:
        node_ids = []
        if self.has_converged():
            node = self.get_best_node()
            node_ids = self.get_path_to_node(node)
        return node_ids

    def get_path_to_node(self, node: Node) -> List[int]:
        node_ids = []
        node_ids.append(node.node_id)
        while node.parent_id is not None:
            node_ids.append(node.parent_id)
            node = self.nodes[node.parent_id]
        node_ids.reverse()
        return node_ids

    def get_biasing_power(self, node: Node) -> float:
        node_ids = self.get_path_to_node(node)
        biasing_power = 1
        for node_id in node_ids[:-1]:
            parent_node = self.nodes[node_id]
            sibling_nodes = [child for child in parent_node.children if child.node_id not in node_ids]
            weights = sum(map(self.get_node_weight, sibling_nodes)) + 1
            biasing_power *= 1 / weights
        return biasing_power

    def get_node_weight(self, node: Node) -> float:
        if len(node.children) == 0:
            return 1
        return sum(map(self.get_node_weight, node.children))
