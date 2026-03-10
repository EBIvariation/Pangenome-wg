"""Tests for canonical serialization."""

import json

from pangenome_id.canonicalize import serialize, serialize_path, serialize_topology
from pangenome_id.model import AbstractGraph, Edge, Node, Path, Step


def _make_graph():
    nodes = [
        Node(id="aaaa"),
        Node(id="bbbb"),
    ]
    edges = [
        Edge(node_a="aaaa", orient_a="+", node_b="bbbb", orient_b="+"),
    ]
    steps = [Step(node_id="aaaa", orient="+"), Step(node_id="bbbb", orient="+")]
    paths = [Path(name="p1", steps=steps)]
    return AbstractGraph(nodes=nodes, edges=edges, paths=paths)


def test_serialize_is_bytes():
    g = _make_graph()
    result = serialize(g)
    assert isinstance(result, bytes)
    doc = json.loads(result)
    assert "graph_topology" in doc
    assert "names" in doc
    assert "paths" in doc


def test_node_order_does_not_affect_output():
    """Reordering nodes in the AGM must not change canonical bytes."""
    g1 = _make_graph()
    g2 = _make_graph()
    g2.nodes = list(reversed(g2.nodes))
    assert serialize(g1) == serialize(g2)


def test_edge_order_does_not_affect_output():
    """Reordering edges must not change canonical bytes."""
    nodes = [Node(id="aaaa"), Node(id="bbbb"), Node(id="cccc")]
    e1 = Edge(node_a="aaaa", orient_a="+", node_b="bbbb", orient_b="+")
    e2 = Edge(node_a="aaaa", orient_a="+", node_b="cccc", orient_b="+")
    path = Path(name="p", steps=[])

    g1 = AbstractGraph(nodes=nodes, edges=[e1, e2], paths=[path])
    g2 = AbstractGraph(nodes=nodes, edges=[e2, e1], paths=[path])
    assert serialize(g1) == serialize(g2)


def test_path_step_order_is_preserved():
    """Reversing path steps must produce different canonical bytes."""
    nodes = [Node(id="aaaa"), Node(id="bbbb")]
    edge = Edge(node_a="aaaa", orient_a="+", node_b="bbbb", orient_b="+")

    steps_fwd = [Step(node_id="aaaa", orient="+"), Step(node_id="bbbb", orient="+")]
    steps_rev = [Step(node_id="bbbb", orient="+"), Step(node_id="aaaa", orient="+")]

    g1 = AbstractGraph(nodes=nodes, edges=[edge], paths=[Path(name="p", steps=steps_fwd)])
    g2 = AbstractGraph(nodes=nodes, edges=[edge], paths=[Path(name="p", steps=steps_rev)])
    assert serialize(g1) != serialize(g2)


def test_overlap_policy_embedded():
    """Different overlap policies produce different bytes."""
    g1 = _make_graph()
    g1.overlap_policy = "discard"
    g2 = _make_graph()
    g2.overlap_policy = "length_only"
    assert serialize(g1) != serialize(g2)


def test_serialize_topology_excludes_paths():
    """Adding paths must not change serialize_topology() output."""
    nodes = [Node(id="aaaa"), Node(id="bbbb")]
    edge = Edge(node_a="aaaa", orient_a="+", node_b="bbbb", orient_b="+")
    steps = [Step(node_id="aaaa", orient="+"), Step(node_id="bbbb", orient="+")]

    g1 = AbstractGraph(nodes=nodes, edges=[edge], paths=[])
    g2 = AbstractGraph(nodes=nodes, edges=[edge], paths=[Path(name="p1", steps=steps)])
    assert serialize_topology(g1) == serialize_topology(g2)


def test_serialize_path_excludes_name():
    """Two paths with different names but identical steps have equal serialize_path() output."""
    steps = [Step(node_id="aaaa", orient="+"), Step(node_id="bbbb", orient="+")]
    p1 = Path(name="alpha", steps=steps)
    p2 = Path(name="beta", steps=steps)
    assert serialize_path(p1) == serialize_path(p2)


def test_serialize_paths_sorted_by_digest():
    """serialize(g)['paths'] is sorted lexicographically."""
    nodes = [Node(id="aaaa"), Node(id="bbbb"), Node(id="cccc")]
    steps_ab = [Step(node_id="aaaa", orient="+"), Step(node_id="bbbb", orient="+")]
    steps_ac = [Step(node_id="aaaa", orient="+"), Step(node_id="cccc", orient="+")]
    paths = [Path(name="p2", steps=steps_ac), Path(name="p1", steps=steps_ab)]
    g = AbstractGraph(nodes=nodes, edges=[], paths=paths)

    doc = json.loads(serialize(g))
    assert doc["paths"] == sorted(doc["paths"])
