"""Canonical serialization of AbstractGraph to deterministic bytes."""

import struct

from pangenome_id.model import AbstractGraph

FORMAT_VERSION = 0x02

OVERLAP_POLICY_BYTES = {
    "discard": 0x00,
    "length_only": 0x01,
    "full_cigar": 0x02,
}

# Delimiter bytes — never appear in node ids (URL-safe base64), orientations
# ('+'/'-'), overlap strings (CIGAR/digits), or path names.
_SEP = b"\x00"  # field separator within a record
_REC = b"\x01"  # record separator
_SEC = b"\x02"  # section separator (nodes / edges / paths)


def serialize(graph: AbstractGraph) -> bytes:
    """Produce a deterministic binary representation of an AbstractGraph.

    The layout is:
        [FORMAT_VERSION 1B] [OVERLAP_POLICY 1B]
        for each node (sorted by id):
            id \\x01
        \\x02
        for each edge (sorted by node_a, orient_a, node_b, orient_b):
            node_a \\x00 orient_a \\x00 node_b \\x00 orient_b \\x00 overlap \\x01
        \\x02
        for each path (sorted by name), step order preserved:
            name \\x00 is_reference \\x00 step_id \\x00 orient [ \\x00 step_id \\x00 orient ...] \\x01
        \\x02

    Delimiter bytes:
        \\x00  field separator within a record
        \\x01  record separator
        \\x02  section separator

    Sorting nodes and edges makes the output independent of insertion order.
    Path step order is intentionally preserved — it encodes biological coordinates.
    The overlap policy byte is embedded so graphs hashed under different policies
    never collide even when the graph topology is identical.
    """
    buf = bytearray()

    # Header
    buf += struct.pack("<B", FORMAT_VERSION)
    buf += struct.pack("<B", OVERLAP_POLICY_BYTES[graph.overlap_policy])

    # Nodes — sorted by node id
    for node in sorted(graph.nodes, key=lambda n: n.id):
        buf += node.id.encode("ascii") + _REC
    buf += _SEC

    # Edges — sorted by (node_a, orient_a, node_b, orient_b)
    for edge in sorted(graph.edges, key=lambda e: (e.node_a, e.orient_a, e.node_b, e.orient_b)):
        buf += edge.node_a.encode("ascii") + _SEP
        buf += edge.orient_a.encode("ascii") + _SEP
        buf += edge.node_b.encode("ascii") + _SEP
        buf += edge.orient_b.encode("ascii") + _SEP
        buf += (edge.overlap or "").encode("ascii") + _REC
    buf += _SEC

    # Paths — sorted by name; step order preserved
    for path in sorted(graph.paths, key=lambda p: p.name):
        buf += path.name.encode("utf-8") + _SEP
        buf += (b"1" if path.is_reference else b"0")
        for step in path.steps:
            buf += _SEP + step.node_id.encode("ascii") + _SEP + step.orient.encode("ascii")
        buf += _REC
    buf += _SEC

    return bytes(buf)
