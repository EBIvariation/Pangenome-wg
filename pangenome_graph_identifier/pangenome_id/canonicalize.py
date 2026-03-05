"""Canonical serialization of AbstractGraph to deterministic bytes."""

import struct

from pangenome_id.model import AbstractGraph

FORMAT_VERSION = 0x01

OVERLAP_POLICY_BYTES = {
    "discard": 0x00,
    "length_only": 0x01,
    "full_cigar": 0x02,
}

ORIENT_BYTES = {"+": 0x00, "-": 0x01}


def _write_len_prefixed(buf: bytearray, data: bytes, len_format: str) -> None:
    """Append a length-prefixed field to buf.

    Writes the byte-length of data using struct format len_format (little-endian),
    followed by the raw bytes of data. This avoids null-terminator ambiguity when
    fields may themselves contain zero bytes.
    """
    buf += struct.pack(len_format, len(data))
    buf += data


def serialize(graph: AbstractGraph) -> bytes:
    """Produce a deterministic binary representation of an AbstractGraph.

    The layout is:
        [FORMAT_VERSION 1B] [OVERLAP_POLICY 1B]
        [N_NODES 8B] for each node sorted by id:
            [ID_LEN 2B] [ID] [SEQ_LEN 4B] [SEQ]
        [N_EDGES 8B] for each edge sorted by (node_a, orient_a, node_b, orient_b):
            [NODE_A_LEN 2B] [NODE_A] [ORIENT_A 1B]
            [NODE_B_LEN 2B] [NODE_B] [ORIENT_B 1B]
            [OVERLAP_LEN 2B] [OVERLAP]
        [N_PATHS 8B] for each path sorted by name:
            [NAME_LEN 2B] [NAME] [IS_REFERENCE 1B] [N_STEPS 4B]
            for each step (order preserved):
                [NODE_ID_LEN 2B] [NODE_ID] [ORIENT 1B]

    Sorting nodes and edges makes the output independent of insertion order.
    Path step order is intentionally preserved — it encodes biological coordinates.
    The overlap policy byte is embedded so graphs hashed under different policies
    never collide even when the graph topology is identical.

    All integers are little-endian. Orientation: '+' → 0x00, '-' → 0x01.
    """
    buf = bytearray()

    # Header
    buf += struct.pack("<B", FORMAT_VERSION)
    buf += struct.pack("<B", OVERLAP_POLICY_BYTES[graph.overlap_policy])

    # Nodes — sorted by node id
    sorted_nodes = sorted(graph.nodes, key=lambda n: n.id)
    buf += struct.pack("<Q", len(sorted_nodes))
    for node in sorted_nodes:
        id_bytes = node.id.encode("utf-8")
        seq_bytes = node.sequence.encode("ascii")
        _write_len_prefixed(buf, id_bytes, "<H")
        _write_len_prefixed(buf, seq_bytes, "<I")

    # Edges — sorted by (node_a, orient_a, node_b, orient_b)
    sorted_edges = sorted(
        graph.edges, key=lambda e: (e.node_a, e.orient_a, e.node_b, e.orient_b)
    )
    buf += struct.pack("<Q", len(sorted_edges))
    for edge in sorted_edges:
        _write_len_prefixed(buf, edge.node_a.encode("utf-8"), "<H")
        buf += struct.pack("<B", ORIENT_BYTES[edge.orient_a])
        _write_len_prefixed(buf, edge.node_b.encode("utf-8"), "<H")
        buf += struct.pack("<B", ORIENT_BYTES[edge.orient_b])
        overlap_bytes = (edge.overlap or "").encode("utf-8")
        _write_len_prefixed(buf, overlap_bytes, "<H")

    # Paths — sorted by name; step order preserved
    sorted_paths = sorted(graph.paths, key=lambda p: p.name)
    buf += struct.pack("<Q", len(sorted_paths))
    for path in sorted_paths:
        _write_len_prefixed(buf, path.name.encode("utf-8"), "<H")
        buf += struct.pack("<B", 0x01 if path.is_reference else 0x00)
        buf += struct.pack("<I", len(path.steps))
        for step in path.steps:
            _write_len_prefixed(buf, step.node_id.encode("utf-8"), "<H")
            buf += struct.pack("<B", ORIENT_BYTES[step.orient])

    return bytes(buf)
