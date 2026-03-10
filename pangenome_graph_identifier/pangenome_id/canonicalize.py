"""Canonical JSON serialization of AbstractGraph using RFC-8785."""

import base64
import hashlib
import json

import rfc8785

from pangenome_id.model import AbstractGraph, Path


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
    pairs = sorted(
        [(sha512t24u(serialize_path(p)), p.name) for p in graph.paths]
    )  # sort by digest (first element)
    doc = {
        "graph_topology": topology_digest,
        "names": [name for _, name in pairs],
        "paths": [digest for digest, _ in pairs],
    }
    return rfc8785.dumps(doc)
