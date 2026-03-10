"""Identify a pangenome graph by hashing its canonical JSON serialization."""

import base64
import hashlib
import json

from pangenome_id.canonicalize import serialize
from pangenome_id.model import AbstractGraph


def sha512t24u(data: bytes) -> str:
    """GA4GH sha512t24u digest: first 24 bytes of SHA-512, URL-safe base64."""
    digest = hashlib.sha512(data).digest()
    return base64.urlsafe_b64encode(digest[:24]).decode("ascii")


def identify_graph(graph: AbstractGraph) -> dict:
    """Return the canonical JSON document for a graph as a Python dict."""
    return json.loads(serialize(graph))


def identify_from_string(
    text: str,
    format: str,
    overlap_policy: str = "discard",
) -> dict:
    """Parse a GFA string and return its identifier dict. Useful in tests."""
    from pangenome_id.parsers.gfa1 import GFA1Parser
    from pangenome_id.parsers.gfa2 import GFA2Parser

    if format == "gfa1":
        parser = GFA1Parser(overlap_policy=overlap_policy)
        graph = parser.parse_string(text)
    elif format == "gfa2":
        parser = GFA2Parser(overlap_policy=overlap_policy)
        graph = parser.parse_string(text)
    else:
        raise ValueError(f"Unknown format: {format!r}")
    return identify_graph(graph)
