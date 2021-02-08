import os
from distutils.dir_util import copy_tree

from ligbinder.tree import Tree

import pytest
from .conftest import FPATH


def copy_input_files_to_tmpdir(tmpdir):
    project_folder = os.path.join(FPATH, "project")
    copy_tree(str(project_folder), str(tmpdir))


@pytest.fixture
def basic_tree(tmpdir):
    copy_input_files_to_tmpdir(tmpdir)
    tree = Tree(tmpdir)
    return tree


def test_tree_create_root_node(basic_tree: Tree):
    node = basic_tree.create_root_node()
    assert node.rmsd == pytest.approx(16.4291511)
    assert node.parent_id is None
    assert node.depth == 0
    assert os.path.exists(os.path.join(node.path, ".info.yml"))
    assert len(basic_tree.nodes) == 1


def test_tree_get_next_id(basic_tree):
    assert basic_tree.get_next_id() == 1
    basic_tree.create_root_node()
    assert basic_tree.get_next_id() == 1
    basic_tree.create_node(0)
    assert basic_tree.get_next_id() == 2


def test_tree_create_node(basic_tree: Tree):
    basic_tree.create_root_node()
    node = basic_tree.create_node(0)
    assert len(basic_tree.nodes) == 2
    assert node.parent_id == 0
    assert node.rmsd is None
    assert node.depth == 1
    assert len(basic_tree.nodes[0].children) == 1
    assert os.path.exists(os.path.join(node.path, "top.prmtop"))
    assert os.path.exists(os.path.join(node.path, "initial.rst7"))
    assert os.path.exists(os.path.join(node.path, "ref.crd"))


def test_tree_has_converged(basic_tree):
    basic_tree.create_root_node()
    assert basic_tree.has_converged() is False
    basic_tree.nodes[0].rmsd = 0.2
    assert basic_tree.has_converged() is True


def test_tree_load_nodes(basic_tree):
    basic_tree.create_root_node()
    basic_tree.create_node(0)
    basic_tree.create_node(0)
    tree = Tree(basic_tree.path)
    assert len(tree.nodes) == len(basic_tree.nodes)
    assert len(tree.nodes[0].children) == len(basic_tree.nodes[0].children)
    for id in tree.nodes:
        assert tree.nodes[id].depth == basic_tree.nodes[id].depth
        assert tree.nodes[id].parent_id == basic_tree.nodes[id].parent_id


def test_tree_expandability(basic_tree: Tree):
    basic_tree.max_children = 2
    basic_tree.max_depth = 1
    basic_tree.create_root_node()
    node = basic_tree.create_node(0)
    node.rmsd = basic_tree.nodes[0].rmsd + 1
    assert basic_tree.can_grow()
    node = basic_tree.create_node_from_candidate()
    node.rmsd = basic_tree.nodes[0].rmsd + 1
    assert basic_tree.can_grow() is False
