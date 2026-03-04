"""Tests for canonical serialization."""

from pangenome_id.canonicalize import serialize
from pangenome_id.model import AbstractGraph, Edge, Node, Path, Step


def _make_graph():
    nodes = [
        Node(id="aaaa", sequence="ACGT"),
        Node(id="bbbb", sequence="TTGC"),
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
    assert len(result) > 0


def test_node_order_does_not_affect_output():
    """Reordering nodes in the AGM must not change canonical bytes."""
    g1 = _make_graph()
    g2 = _make_graph()
    g2.nodes = list(reversed(g2.nodes))
    assert serialize(g1) == serialize(g2)


def test_edge_order_does_not_affect_output():
    """Reordering edges must not change canonical bytes."""
    nodes = [Node(id="aaaa", sequence="ACGT"), Node(id="bbbb", sequence="TTGC"), Node(id="cccc", sequence="GGGG")]
    e1 = Edge(node_a="aaaa", orient_a="+", node_b="bbbb", orient_b="+")
    e2 = Edge(node_a="aaaa", orient_a="+", node_b="cccc", orient_b="+")
    path = Path(name="p", steps=[])

    g1 = AbstractGraph(nodes=nodes, edges=[e1, e2], paths=[path])
    g2 = AbstractGraph(nodes=nodes, edges=[e2, e1], paths=[path])
    assert serialize(g1) == serialize(g2)


def test_path_step_order_is_preserved():
    """Reversing path steps must produce different canonical bytes."""
    nodes = [Node(id="aaaa", sequence="ACGT"), Node(id="bbbb", sequence="TTGC")]
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
