"""Tests for GFA v1 parser."""

import hashlib

from pangenome_id.parsers.gfa1 import GFA1Parser

MINIMAL_GFA1 = """\
H\tVN:Z:1.0
S\ts1\tACGT
S\ts2\tTTGC
L\ts1\t+\ts2\t+\t0M
P\tp1\ts1+,s2+\t*
"""


def _node_id(seq):
    return hashlib.sha256(seq.encode("ascii")).hexdigest()[:16]


def test_node_count():
    g = GFA1Parser().parse_string(MINIMAL_GFA1)
    assert len(g.nodes) == 2


def test_node_ids_are_sequence_derived():
    g = GFA1Parser().parse_string(MINIMAL_GFA1)
    ids = {n.id for n in g.nodes}
    assert _node_id("ACGT") in ids
    assert _node_id("TTGC") in ids


def test_node_sequences():
    g = GFA1Parser().parse_string(MINIMAL_GFA1)
    seqs = {n.sequence for n in g.nodes}
    assert "ACGT" in seqs
    assert "TTGC" in seqs


def test_edge_count():
    g = GFA1Parser().parse_string(MINIMAL_GFA1)
    assert len(g.edges) == 1


def test_edge_connects_correct_nodes():
    g = GFA1Parser().parse_string(MINIMAL_GFA1)
    edge = g.edges[0]
    expected_a = _node_id("ACGT")
    expected_b = _node_id("TTGC")
    # Edge may be in canonical order
    assert {edge.node_a, edge.node_b} == {expected_a, expected_b}


def test_path_count():
    g = GFA1Parser().parse_string(MINIMAL_GFA1)
    assert len(g.paths) == 1


def test_path_name():
    g = GFA1Parser().parse_string(MINIMAL_GFA1)
    assert g.paths[0].name == "p1"


def test_path_step_count():
    g = GFA1Parser().parse_string(MINIMAL_GFA1)
    assert len(g.paths[0].steps) == 2


def test_path_step_order():
    g = GFA1Parser().parse_string(MINIMAL_GFA1)
    steps = g.paths[0].steps
    assert steps[0].node_id == _node_id("ACGT")
    assert steps[1].node_id == _node_id("TTGC")


def test_path_step_orientations():
    g = GFA1Parser().parse_string(MINIMAL_GFA1)
    steps = g.paths[0].steps
    assert steps[0].orient == "+"
    assert steps[1].orient == "+"


def test_overlap_policy_discard():
    g = GFA1Parser(overlap_policy="discard").parse_string(MINIMAL_GFA1)
    assert g.edges[0].overlap is None


def test_overlap_policy_length_only():
    g = GFA1Parser(overlap_policy="length_only").parse_string(MINIMAL_GFA1)
    # "0M" → total length 0
    assert g.edges[0].overlap == "0"


def test_comments_skipped():
    gfa = "# comment\n" + MINIMAL_GFA1
    g = GFA1Parser().parse_string(gfa)
    assert len(g.nodes) == 2


def test_is_reference_false():
    g = GFA1Parser().parse_string(MINIMAL_GFA1)
    assert g.paths[0].is_reference is False
