"""Tests for GFA v2 parser."""

from pangenome_id.hasher import sha512t24u
from pangenome_id.parsers.gfa2 import GFA2Parser

MINIMAL_GFA2 = """\
H\tVN:Z:2.0
S\tnode_a\t4\tACGT
S\tnode_b\t4\tTTGC
E\t*\tnode_a+\tnode_b+\t4\t4$\t0\t0$\t*
O\tpath1\tnode_a+ node_b+
"""


def _node_id(seq):
    return sha512t24u(seq.encode("ascii"))


def test_node_count():
    g = GFA2Parser().parse_string(MINIMAL_GFA2)
    assert len(g.nodes) == 2


def test_node_ids_are_sequence_derived():
    g = GFA2Parser().parse_string(MINIMAL_GFA2)
    ids = {n.id for n in g.nodes}
    assert _node_id("ACGT") in ids
    assert _node_id("TTGC") in ids


def test_edge_count():
    g = GFA2Parser().parse_string(MINIMAL_GFA2)
    assert len(g.edges) == 1


def test_edge_connects_correct_nodes():
    g = GFA2Parser().parse_string(MINIMAL_GFA2)
    edge = g.edges[0]
    expected_a = _node_id("ACGT")
    expected_b = _node_id("TTGC")
    assert {edge.node_a, edge.node_b} == {expected_a, expected_b}


def test_path_count():
    g = GFA2Parser().parse_string(MINIMAL_GFA2)
    assert len(g.paths) == 1


def test_path_name():
    g = GFA2Parser().parse_string(MINIMAL_GFA2)
    assert g.paths[0].name == "path1"


def test_path_step_count():
    g = GFA2Parser().parse_string(MINIMAL_GFA2)
    assert len(g.paths[0].steps) == 2


def test_path_step_order():
    g = GFA2Parser().parse_string(MINIMAL_GFA2)
    steps = g.paths[0].steps
    assert steps[0].node_id == _node_id("ACGT")
    assert steps[1].node_id == _node_id("TTGC")


def test_header_skipped():
    g = GFA2Parser().parse_string(MINIMAL_GFA2)
    # H line shouldn't add a node
    assert len(g.nodes) == 2


def test_fragment_lines_skipped():
    gfa = MINIMAL_GFA2 + "F\tseg\tread\t0\t10\t0\t10\t+\t*\n"
    g = GFA2Parser().parse_string(gfa)
    assert len(g.nodes) == 2


def test_unordered_group_skipped():
    gfa = MINIMAL_GFA2 + "U\tset1\tnode_a node_b\n"
    g = GFA2Parser().parse_string(gfa)
    assert len(g.paths) == 1


def test_overlap_policy_discard():
    g = GFA2Parser(overlap_policy="discard").parse_string(MINIMAL_GFA2)
    assert g.edges[0].overlap is None


def test_overlap_policy_length_only():
    # beg1=4, end1=4$ → overlap = 4-4 = 0
    g = GFA2Parser(overlap_policy="length_only").parse_string(MINIMAL_GFA2)
    assert g.edges[0].overlap == "0"


def test_parse_gzipped_file():
    """parse() transparently decompresses .gfa.gz files."""
    import gzip
    import os
    import tempfile

    with tempfile.NamedTemporaryFile(suffix=".gfa.gz", delete=False) as f:
        tmp = f.name
    try:
        with gzip.open(tmp, "wt") as gz:
            gz.write(MINIMAL_GFA2)
        g = GFA2Parser().parse(tmp)
        assert len(g.nodes) == 2
        assert len(g.paths) == 1
    finally:
        os.unlink(tmp)
