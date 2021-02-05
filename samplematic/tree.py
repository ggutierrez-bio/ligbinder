from glob import glob
import os
import yaml
from typing import Dict
import pytraj
import logging


logger = logging.getLogger(__file__)


class Node:

    def __init__(self, path: str, id: int, parent_id: int, rmsd: float = None) -> None:
        self.path = path
        self.id = id
        self.parent_id = parent_id
        self.rmsd = None
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
    def create_node(path: str, id: int, parent_id: int, ref_file: str, top_file: str) -> "Node":
        node_path = os.path.join(path, f"node_{id}")
        os.mkdir(os.path.join(path, f"node_{id}"))
        parent_path = os.path.join(os.path.pardir, f"node_{parent_id}")
        initial_frame = os.path.join(parent_path, "final.rst7")
        os.link(initial_frame, os.path.join(node_path, "initial.rst7"))
        os.link(os.path.relpath(ref_file, node_path), os.path.join(node_path, "ref.crd"))
        os.link(os.path.relpath(top_file, node_path), os.path.join(node_path, "top.prmtop"))
        node = Node(path, id, parent_id)
        node.write_node_info()

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

    def __init__(self, path: str = ".", max_children: int = 5, max_depth: int = 15, tolerance: float = 0.5) -> None:
        self.path = path
        self.max_children = max_children
        self.max_depth = max_depth
        self.tolerance = tolerance
        self.nodes: Dict[int, Node] = {}
        self.load_nodes()

    def load_nodes(self) -> Dict[int, Node]:
        node_paths = glob(f"{self.path}/node*/")
        nodes = [Node(node_path) for node_path in node_paths]
        [self.nodes[node.parent_id].children.append(node) for node in nodes]

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

    def create_root_node(self, crd_file, top_file, ref_file) -> Node:
        node_path = os.path.relpath(os.path.join(self.path, f"node_{id}"))
        os.mkdir(node_path)
        rel_crd_file = os.path.relpath(crd_file, node_path)
        rel_top_file = os.path.relpath(top_file, node_path)
        rel_ref_file = os.path.relpath(ref_file, node_path)
        os.link(rel_crd_file, os.path.join(node_path, "final.rst7"))
        os.link(rel_top_file, os.path.join(node_path, "top.prmtop"))
        os.link(rel_ref_file, os.path.join(node_path, "ref.crd"))
        node = Node(node_path, id, None)
        return node
