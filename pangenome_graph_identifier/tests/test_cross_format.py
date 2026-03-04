"""Cross-format identity tests — the same graph in GFA v1 and GFA v2 must hash identically."""

from pangenome_id.hasher import identify_from_string

GFA1_GRAPH = """\
H\tVN:Z:1.0
S\tseq1\tACGT
S\tseq2\tTTGC
L\tseq1\t+\tseq2\t+\t0M
P\tpath1\tseq1+,seq2+\t*
"""

GFA2_GRAPH = """\
H\tVN:Z:2.0
S\tnode_a\t4\tACGT
S\tnode_b\t4\tTTGC
E\t*\tnode_a+\tnode_b+\t4\t4$\t0\t0$\t*
O\tpath1\tnode_a+ node_b+
"""


def test_cross_format_identity():
    id1 = identify_from_string(GFA1_GRAPH, format="gfa1")
    id2 = identify_from_string(GFA2_GRAPH, format="gfa2")
    assert id1 == id2, f"GFA1={id1!r}, GFA2={id2!r}"


def test_cross_format_identity_ga4gh():
    id1 = identify_from_string(GFA1_GRAPH, format="gfa1", style="ga4gh")
    id2 = identify_from_string(GFA2_GRAPH, format="gfa2", style="ga4gh")
    assert id1 == id2
    assert id1.startswith("ga4gh:pg.")


# ── Collision sanity checks ──────────────────────────────────────────────────

GFA1_EXTRA_NODE = """\
H\tVN:Z:1.0
S\tseq1\tACGT
S\tseq2\tTTGC
S\tseq3\tGGGG
L\tseq1\t+\tseq2\t+\t0M
P\tpath1\tseq1+,seq2+\t*
"""


def test_extra_node_different_hash():
    base_id = identify_from_string(GFA1_GRAPH, format="gfa1")
    extra_id = identify_from_string(GFA1_EXTRA_NODE, format="gfa1")
    assert base_id != extra_id


GFA1_REVERSED_STEP = """\
H\tVN:Z:1.0
S\tseq1\tACGT
S\tseq2\tTTGC
L\tseq1\t+\tseq2\t+\t0M
P\tpath1\tseq2+,seq1+\t*
"""


def test_reversed_step_different_hash():
    base_id = identify_from_string(GFA1_GRAPH, format="gfa1")
    rev_id = identify_from_string(GFA1_REVERSED_STEP, format="gfa1")
    assert base_id != rev_id


def test_overlap_policy_different_hash():
    id_discard = identify_from_string(GFA1_GRAPH, format="gfa1", overlap_policy="discard")
    id_length = identify_from_string(GFA1_GRAPH, format="gfa1", overlap_policy="length_only")
    assert id_discard != id_length
