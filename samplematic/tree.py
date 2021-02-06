from glob import glob
import os
import yaml
from typing import Dict
import pytraj
import logging

from samplematic.settings import Settings

logger = logging.getLogger(__file__)


class Node:

    def __init__(self, path: str, id: int, parent_id: int, rmsd: float = None, depth: int = None) -> None:
        self.path = path
        self.id = id
        self.parent_id = parent_id
        self.rmsd = rmsd
        self.depth = depth
        self.children = []

    @staticmethod
    def load_node(path: str) -> "Node":
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
        with open(filename, 'r') as cfg_file:
            return yaml.load(cfg_file, Loader=yaml.FullLoader)

    def write_node_info(self) -> None:
        filename = os.path.join(self.path, ".info.yml")
        with open(filename, 'w') as file:
            yaml.dump({
                "parent_id": self.parent_id,
                "rmsd": self.rmsd
            }, file)

    def get_relative_file(self, filename) -> str:
        return os.path.relpath(os.path.join(self.path, filename))

    def calc_node_rmsd(self) -> float:
        crd_file = self.get_relative_file("final.rst7")
        top_file = self.get_relative_file("top.prmtop")
        ref_file = self.get_relative_file("ref.crd")
        if not all(os.path.exists(file) for file in [crd_file, top_file, ref_file]):
            logger.warning(f"Unable to calculate rmsd for node {self.id}. Some files are missing")
        traj = pytraj.autoimage(pytraj.load(crd_file, top=top_file))
        ref = pytraj.load(ref_file, top=top_file)
        pytraj.rmsd(traj, mask="@CA", ref=ref)
        self.rmsd = pytraj.rmsd_nofit(traj, mask=":LIG", ref=ref)[0]
        self.write_node_info()
        return self.rmsd


class Tree():

    def __init__(self, path: str = ".", max_children: int = 5, max_depth: int = 15, tolerance: float = 0.5, settings = None) -> None:
        self.path = path
        self.max_children = max_children
        self.max_depth = max_depth
        self.tolerance = tolerance
        self.nodes: Dict[int, Node] = {}
        self.settings = settings if settings is not None else Settings()
        self.load_nodes()

    def load_nodes(self) -> Dict[int, Node]:
        node_paths = glob(f"{self.path}/node_*/")
        nodes = [Node(node_path) for node_path in node_paths]
        self.nodes = {node.id: node for node in nodes}
        [self.nodes[node.parent_id].children.append(node) for node in nodes if node.parent_id is not None]
        self.set_depths()
        return self.nodes

    def get_next_id(self) -> int:
        return max(list(self.nodes), default=0) + 1

    def create_node(self, id: int, parent_id: int, ref_file, top_file) -> Node:
        node_path = os.path.join(self.path, f"node_{id}")
        os.mkdir(node_path)
        parent_path = os.path.relpath(os.path.join(os.path.pardir, f"node_{parent_id}"), start=node_path)
        initial_frame = os.path.join(parent_path, "final.rst7")
        os.link(initial_frame, os.path.join(node_path, "initial.rst7"))
        os.link(os.path.relpath(ref_file, node_path), os.path.join(node_path, "ref.crd"))
        os.link(os.path.relpath(top_file, node_path), os.path.join(node_path, "top.prmtop"))
        parent_node = self.nodes[parent_id]
        node = Node(node_path, id, parent_id, depth=parent_node.depth + 1)
        node.write_node_info()
        return node

    def create_root_node(
        self,
        crd_file: str = "./data/pose.rst7",
        top_file: str = "./data/top.prmtop",
        ref_file: str = "./data/ref.crd"
    ) -> Node:
        node_path = os.path.relpath(os.path.join(self.path, f"node_{id}"))
        os.mkdir(node_path)
        rel_crd_file = os.path.relpath(crd_file, node_path)
        rel_top_file = os.path.relpath(top_file, node_path)
        rel_ref_file = os.path.relpath(ref_file, node_path)
        os.link(rel_crd_file, os.path.join(node_path, self.settings["md"]["rst_file"]))
        os.link(rel_top_file, os.path.join(node_path, self.settings["md"]["top_file"]))
        os.link(rel_ref_file, os.path.join(node_path, self.settings["md"]["ref_file"]))
        node = Node(node_path, id, None)
        node.calc_node_rmsd()
        node.depth = 0
        return node

    def create_node_from_candidate(self):
        new_id = self.get_next_id()
        parent = self.choose_candidate_node()
        node = self.create_node(
            new_id,
            parent.id,
            self.settings["data_files"]["ref_file"],
            self.settings["datafiles"]["top_file"]
        )
        
    def choose_candidate_node(self) -> Node:
        return [node for node in self.nodes.values() if self.is_expandable(node)].sort(key=lambda n: n.rmsd)[0]

    def set_depths(self):
        if len(self.nodes) > 0:
            self.set_node_depth(self.nodes[0], 0)

    def set_node_depth(self, node: Node, current_depth: int) -> None:
        node.depth = current_depth
        for child in node.children:
            self.set_node_depth(child, current_depth + 1)

    def can_grow(self):
        return any(self.is_expandable(node) for node in list(self.nodes.values()))

    def is_expandable(self, node: Node) -> bool:
        return not (
            node.depth >= self.max_depth
            or len(node.children) >= self.max_children
            or node.rmsd >= self.nodes[node.parent_id].rmsd
        )

    def has_converged(self) -> bool:
        return any(node.rmsd <= self.tolerance for node in list(self.nodes.values()))
