"""Hash canonical bytes into a graph identifier."""

import base64
import hashlib

from pangenome_id.canonicalize import serialize
from pangenome_id.model import AbstractGraph


def sha512t24u(data: bytes) -> str:
    """GA4GH sha512t24u digest: first 24 bytes of SHA-512, URL-safe base64."""
    digest = hashlib.sha512(data).digest()
    return base64.urlsafe_b64encode(digest[:24]).decode("ascii")


def compute_identifier(canonical_bytes: bytes) -> str:
    """Hash canonical bytes into a GA4GH-style graph identifier.

    Returns a string of the form "ga4gh:pg.<32 url-safe base64 chars>",
    derived from the sha512t24u digest of the canonical bytes.
    """
    return f"ga4gh:pg.{sha512t24u(canonical_bytes)}"


def identify_graph(graph: AbstractGraph) -> str:
    return compute_identifier(serialize(graph))


def identify_from_string(
    text: str,
    format: str,
    overlap_policy: str = "discard",
) -> str:
    """Parse a GFA string and return its identifier. Useful in tests."""
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
