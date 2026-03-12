"""Canonical JSON serialization of AbstractGraph using RFC-8785."""

import base64
import hashlib
import json

import rfc8785

from pangenome_id.model import AbstractGraph, Path

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence (uppercase ACGT; other chars pass through)."""
    return seq.translate(_COMPLEMENT)[::-1]


def sequence_digest_for_path(path: Path, seq_by_id: dict[str, str | None]) -> str | None:
    """Compute refget-compatible sequence digest for a path.

    Returns None if any node has a missing sequence (wildcard '*' node).
    Returns sha512t24u of the empty string for an empty path.
    """
    parts = []
    for step in path.steps:
        seq = seq_by_id.get(step.node_id)
        if seq is None:
            return None
        parts.append(reverse_complement(seq) if step.orient == "-" else seq)
    return sha512t24u("".join(parts).encode("ascii"))


def sha512t24u(data: bytes) -> str:
    """GA4GH sha512t24u: first 24 bytes of SHA-512, URL-safe base64."""
    return base64.urlsafe_b64encode(hashlib.sha512(data).digest()[:24]).decode("ascii")


def serialize_topology(graph: AbstractGraph) -> bytes:
    """RFC-8785 canonical bytes for the graph topology (nodes + edges).

    Structure:
      {
        "edges": [{"from": ..., "from_orient": ..., "overlap": ..., "to": ..., "to_orient": ...}, ...],
        "nodes": [...],
        "overlap_policy": "discard" | "length_only" | "full_cigar"
      }

    Nodes are sorted lexicographically by id. Edges are sorted by
    (from, from_orient, to, to_orient). Path information is excluded.
    """
    doc = {
        "edges": [
            {
                "from": e.node_a,
                "from_orient": e.orient_a,
                "overlap": e.overlap,
                "to": e.node_b,
                "to_orient": e.orient_b,
            }
            for e in sorted(graph.edges, key=lambda e: (e.node_a, e.orient_a, e.node_b, e.orient_b))
        ],
        "nodes": sorted(n.id for n in graph.nodes),
        "overlap_policy": graph.overlap_policy,
    }
    return rfc8785.dumps(doc)


def serialize_path(path: Path) -> bytes:
    """RFC-8785 canonical bytes for a single path (name excluded).

    Structure:
      {
        "is_reference": false,
        "steps": [{"id": "<node_id>", "orient": "+"}, ...]
      }

    Steps are in path traversal order (preserved, not sorted).
    The path name is not included — it is linked externally in the graph document.
    """
    doc = {
        "is_reference": path.is_reference,
        "steps": [{"id": s.node_id, "orient": s.orient} for s in path.steps],
    }
    return rfc8785.dumps(doc)


def serialize(graph: AbstractGraph) -> bytes:
    """RFC-8785 canonical JSON document for the complete AbstractGraph.

    Structure:
      {
        "graph_topology": sha512t24u digest of topology (nodes + edges),
        "names":  path names, in the order determined by "paths",
        "paths":  sha512t24u digests, sorted lexicographically by digest value
      }

    Path names and digests are parallel arrays ordered by digest value,
    making the document order independent of path name ordering.
    The returned bytes can be parsed as JSON directly.
    """
    topology_digest = sha512t24u(serialize_topology(graph))
    seq_by_id = {n.id: n.sequence for n in graph.nodes}
    triples = sorted(
        [(sha512t24u(serialize_path(p)), p.name, sequence_digest_for_path(p, seq_by_id))
         for p in graph.paths]
    )  # sort by path topology digest (first element)
    doc = {
        "graph_topology": topology_digest,
        "names": [name for _, name, _ in triples],
        "paths": [digest for digest, _, _ in triples],
        "sequences": [seq for _, _, seq in triples],
    }
    return rfc8785.dumps(doc)
