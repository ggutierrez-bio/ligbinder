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
    assert node.rmsd == pytest.approx(17.29049097)
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


@pytest.mark.parametrize(
    "use_normalized_rmsd",
    [
        pytest.param(False, id="rmsd"),
        pytest.param(True, id="normalized rmsd")
    ]
)
def test_tree_has_converged(use_normalized_rmsd, basic_tree: Tree):
    basic_tree.use_normalized_rmsd = use_normalized_rmsd
    basic_tree.create_root_node()
    assert basic_tree.has_converged() is False
    basic_tree.nodes[0].nrmsd = basic_tree.nodes[0].rmsd = 0.2
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


@pytest.mark.parametrize(
    "use_normalized_rmsd",
    [
        pytest.param(False, id="rmsd"),
        pytest.param(True, id="normalized rmsd")
    ]
)
def test_tree_expandability(use_normalized_rmsd, basic_tree: Tree):
    basic_tree.use_normalized_rmsd = use_normalized_rmsd
    basic_tree.max_children = 2
    basic_tree.max_depth = 1
    basic_tree.create_root_node()
    node = basic_tree.create_node(0)
    node.rmsd = basic_tree.nodes[0].rmsd + 1
    node.nrmsd = basic_tree.nodes[0].nrmsd + 1
    basic_tree.max_nodes = 1
    assert basic_tree.can_grow() is False
    basic_tree.max_nodes = 500
    assert basic_tree.can_grow()
    node = basic_tree.create_node_from_candidate()
    node.rmsd = basic_tree.nodes[0].rmsd + 1
    node.nrmsd = basic_tree.nodes[0].nrmsd + 1
    assert basic_tree.can_grow() is False


@pytest.mark.parametrize(
    "use_normalized_rmsd",
    [
        pytest.param(False, id="rmsd"),
        pytest.param(True, id="normalized rmsd")
    ]
)
def test_node_expandability(basic_tree: Tree, use_normalized_rmsd):
    basic_tree.use_normalized_rmsd = use_normalized_rmsd
    basic_tree.create_root_node()
    parent_node = basic_tree.nodes[0]
    parent_metric = basic_tree.get_metric(parent_node)
    # deactivate relative improvement check
    basic_tree.min_relative_improvement = 1

    node = basic_tree.create_node(0)
    node_metric = parent_metric - 0.5 * basic_tree.min_absolute_improvement
    node.nrmsd = node.rmsd = node_metric
    assert basic_tree.is_expandable(node) is False

    node_metric = parent_metric - 1.5 * basic_tree.min_absolute_improvement
    node.nrmsd = node.rmsd = node_metric
    assert basic_tree.is_expandable(node) is True

    # deactivate absolute improvement check and reactivate relative
    basic_tree.min_absolute_improvement = parent_metric
    basic_tree.min_relative_improvement = 0.1

    node_metric = parent_metric * (1 - 0.5 * basic_tree.min_relative_improvement)
    node.nrmsd = node.rmsd = node_metric
    assert basic_tree.is_expandable(node) is False

    node_metric = parent_metric * (1 - 1.5 * basic_tree.min_relative_improvement)
    node.nrmsd = node.rmsd = node_metric
    assert basic_tree.is_expandable(node) is True


def test_node_weight(basic_tree: Tree):
    root_node = basic_tree.create_root_node()
    assert basic_tree.get_node_weight(root_node) == 1

    node = basic_tree.create_node(0)
    assert basic_tree.get_node_weight(root_node) == 1
    assert basic_tree.get_node_weight(node) == 1

    basic_tree.create_node(0)
    last_node = basic_tree.create_node(0)
    assert basic_tree.get_node_weight(root_node) == 3
    assert basic_tree.get_node_weight(node) == 1

    basic_tree.create_node(last_node.node_id)
    basic_tree.create_node(last_node.node_id)
    assert basic_tree.get_node_weight(root_node) == 4
    assert basic_tree.get_node_weight(node) == 1
    assert basic_tree.get_node_weight(last_node) == 2
