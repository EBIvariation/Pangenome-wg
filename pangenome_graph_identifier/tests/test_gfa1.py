"""Tests for GFA v1 parser."""

from pangenome_id.hasher import sha512t24u
from pangenome_id.parsers.gfa1 import GFA1Parser

MINIMAL_GFA1 = """\
H\tVN:Z:1.0
S\ts1\tACGT
S\ts2\tTTGC
L\ts1\t+\ts2\t+\t0M
P\tp1\ts1+,s2+\t*
"""


def _node_id(seq):
    return sha512t24u(seq.encode("ascii"))


def test_node_count():
    g = GFA1Parser().parse_string(MINIMAL_GFA1)
    assert len(g.nodes) == 2


def test_node_ids_are_sequence_derived():
    g = GFA1Parser().parse_string(MINIMAL_GFA1)
    ids = {n.id for n in g.nodes}
    assert _node_id("ACGT") in ids
    assert _node_id("TTGC") in ids



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


# --- GFA 1.2 Walk ('W') lines ---

WALK_GFA = """\
H\tVN:Z:1.2
S\ts1\tACGT
S\ts2\tTTGC
W\tNA12878\t1\tchr1\t0\t8\t>s1>s2
"""


def test_walk_parsed_as_path():
    g = GFA1Parser().parse_string(WALK_GFA)
    assert len(g.paths) == 1


def test_walk_path_name_uses_pansn():
    g = GFA1Parser().parse_string(WALK_GFA)
    assert g.paths[0].name == "NA12878#1#chr1"


def test_walk_steps_empty():
    # W paths are processed into pre-computed digests; steps are not expanded
    g = GFA1Parser().parse_string(WALK_GFA)
    assert g.paths[0].steps == []


def test_walk_topology_digest_set():
    g = GFA1Parser().parse_string(WALK_GFA)
    assert g.paths[0]._topology_digest is not None


def test_walk_sequence_digest_set():
    g = GFA1Parser().parse_string(WALK_GFA)
    assert g.paths[0]._sequence_digest is not None


def test_walk_reverse_orientation_sequence_digest():
    # Reverse orientation uses reverse complement; digest should differ from forward
    gfa_fwd = "H\tVN:Z:1.2\nS\ts1\tACGT\nS\ts2\tTTGC\nW\tsample\t0\tchr1\t0\t8\t>s1>s2\n"
    gfa_rev = "H\tVN:Z:1.2\nS\ts1\tACGT\nS\ts2\tTTGC\nW\tsample\t0\tchr1\t0\t8\t>s1<s2\n"
    g_fwd = GFA1Parser().parse_string(gfa_fwd)
    g_rev = GFA1Parser().parse_string(gfa_rev)
    assert g_fwd.paths[0]._sequence_digest != g_rev.paths[0]._sequence_digest


def test_walk_topology_digest_matches_rfc8785():
    """Streaming hash must match sha512t24u(rfc8785.dumps(serialize_path(p))) for a W path."""
    import rfc8785
    from pangenome_id.canonicalize import serialize_path, sha512t24u
    from pangenome_id.model import Step

    g = GFA1Parser().parse_string(WALK_GFA)
    walk_path = g.paths[0]

    # Reconstruct expected path with actual steps for reference hash
    id_s1 = _node_id("ACGT")
    id_s2 = _node_id("TTGC")
    from pangenome_id.model import Path
    reference_path = Path(name=walk_path.name, steps=[Step(id_s1, "+"), Step(id_s2, "+")], is_reference=False)
    expected_digest = sha512t24u(serialize_path(reference_path))

    assert walk_path._topology_digest == expected_digest


# --- GFA 1.2 Jump ('J') lines ---

JUMP_GFA = """\
H\tVN:Z:1.2
S\ts1\tACGT
S\ts2\tTTGC
J\ts1\t+\ts2\t+\t100
"""


def test_jump_parsed_as_edge():
    g = GFA1Parser().parse_string(JUMP_GFA)
    assert len(g.edges) == 1


def test_jump_connects_correct_nodes():
    g = GFA1Parser().parse_string(JUMP_GFA)
    edge = g.edges[0]
    assert {edge.node_a, edge.node_b} == {_node_id("ACGT"), _node_id("TTGC")}


def test_jump_overlap_discard():
    g = GFA1Parser(overlap_policy="discard").parse_string(JUMP_GFA)
    assert g.edges[0].overlap is None


def test_jump_overlap_length_only():
    g = GFA1Parser(overlap_policy="length_only").parse_string(JUMP_GFA)
    assert g.edges[0].overlap == "100"


def test_jump_overlap_full_cigar():
    g = GFA1Parser(overlap_policy="full_cigar").parse_string(JUMP_GFA)
    assert g.edges[0].overlap == "100"


# --- Node.sequence population ---

def test_node_sequence_populated():
    g = GFA1Parser().parse_string(MINIMAL_GFA1)
    seqs = {n.id: n.sequence for n in g.nodes}
    assert seqs[_node_id("ACGT")] == "ACGT"
    assert seqs[_node_id("TTGC")] == "TTGC"


def test_node_sequence_none_for_star():
    import warnings
    gfa = "H\tVN:Z:1.0\nS\ts1\t*\n"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        g = GFA1Parser().parse_string(gfa)
    assert g.nodes[0].sequence is None


def test_node_sequence_uppercase_no_u_to_t():
    """normalize_sequence no longer converts U to T."""
    gfa = "H\tVN:Z:1.0\nS\ts1\tacgu\n"
    g = GFA1Parser().parse_string(gfa)
    assert g.nodes[0].sequence == "ACGU"


def test_parse_gzipped_file():
    """parse() transparently decompresses .gfa.gz files."""
    import gzip
    import os
    import tempfile

    with tempfile.NamedTemporaryFile(suffix=".gfa.gz", delete=False) as f:
        tmp = f.name
    try:
        with gzip.open(tmp, "wt") as gz:
            gz.write(MINIMAL_GFA1)
        g = GFA1Parser().parse(tmp)
        assert len(g.nodes) == 2
        assert len(g.paths) == 1
    finally:
        os.unlink(tmp)
