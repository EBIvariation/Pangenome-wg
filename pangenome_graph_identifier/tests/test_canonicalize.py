"""Tests for canonical serialization."""

import json

from pangenome_id.canonicalize import (
    reverse_complement,
    serialize,
    serialize_path,
    serialize_topology,
)
from pangenome_id.hasher import sha512t24u
from pangenome_id.model import AbstractGraph, Edge, Node, Path, Step
from pangenome_id.parsers.gfa1 import GFA1Parser


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
    assert "sequences" in doc


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


# --- reverse_complement tests ---

def test_reverse_complement_simple():
    assert reverse_complement("ACGT") == "ACGT"


def test_reverse_complement_asymmetric():
    assert reverse_complement("AACC") == "GGTT"


def test_reverse_complement_single_base():
    assert reverse_complement("A") == "T"
    assert reverse_complement("C") == "G"
    assert reverse_complement("G") == "C"
    assert reverse_complement("T") == "A"


def test_reverse_complement_non_acgt_passthrough():
    # Non-ACGT chars pass through untranslated but are reversed
    assert reverse_complement("NNA") == "TNN"


# --- sequence_digest_for_path tests ---

def test_sequence_digest_correct_value():
    """Parse a known GFA and verify the sequence digest matches sha512t24u of the concatenated sequence."""
    gfa = "H\tVN:Z:1.0\nS\ts1\tACGT\nS\ts2\tTTGC\nL\ts1\t+\ts2\t+\t0M\nP\tp1\ts1+,s2+\t*\n"
    g = GFA1Parser().parse_string(gfa)
    doc = json.loads(serialize(g))
    expected = sha512t24u("ACGTTTGC".encode("ascii"))
    assert doc["sequences"][0] == expected


def test_sequence_digest_with_reverse_complement():
    """'-' oriented steps use reverse complement of the node sequence."""
    gfa = "H\tVN:Z:1.0\nS\ts1\tACGT\nS\ts2\tTTGC\nL\ts1\t+\ts2\t-\t0M\nP\tp1\ts1+,s2-\t*\n"
    g = GFA1Parser().parse_string(gfa)
    doc = json.loads(serialize(g))
    expected = sha512t24u(("ACGT" + reverse_complement("TTGC")).encode("ascii"))
    assert doc["sequences"][0] == expected


def test_sequence_digest_null_for_star_nodes():
    """Paths containing '*' wildcard nodes get null in 'sequences'."""
    import warnings
    gfa = "H\tVN:Z:1.0\nS\ts1\t*\nP\tp1\ts1+\t*\n"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        g = GFA1Parser().parse_string(gfa)
    doc = json.loads(serialize(g))
    assert doc["sequences"][0] is None


def test_sequences_parallel_to_paths_and_names():
    """'sequences', 'paths', and 'names' arrays have the same length."""
    gfa = "H\tVN:Z:1.0\nS\ts1\tACGT\nS\ts2\tTTGC\nL\ts1\t+\ts2\t+\t0M\nP\tp1\ts1+,s2+\t*\nP\tp2\ts2+,s1+\t*\n"
    g = GFA1Parser().parse_string(gfa)
    doc = json.loads(serialize(g))
    assert len(doc["sequences"]) == len(doc["paths"]) == len(doc["names"]) == 2


def test_sequences_empty_path():
    """Empty path (zero steps) gets sha512t24u of empty string."""
    nodes = [Node(id="aaaa", sequence="ACGT")]
    g = AbstractGraph(nodes=nodes, edges=[], paths=[Path(name="empty", steps=[])])
    doc = json.loads(serialize(g))
    assert doc["sequences"][0] == sha512t24u(b"")


def test_sequences_no_paths():
    """Graph with no paths has empty 'sequences' array."""
    g = AbstractGraph(nodes=[Node(id="aaaa", sequence="ACGT")], edges=[], paths=[])
    doc = json.loads(serialize(g))
    assert doc["sequences"] == []
