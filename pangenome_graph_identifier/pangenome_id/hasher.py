"""Hash canonical bytes into a graph identifier."""

import base64
import hashlib

from pangenome_id.canonicalize import serialize
from pangenome_id.model import AbstractGraph


def compute_identifier(canonical_bytes: bytes, style: str = "hex") -> str:
    """Hash canonical bytes into a graph identifier string.

    Args:
        canonical_bytes: Output of canonicalize.serialize().
        style: "hex" returns the full 64-character SHA-512 hexdigest.
               "ga4gh" returns a compact namespaced identifier of the form
               "ga4gh:pg.<24 url-safe base64 chars>" derived from the first
               18 bytes of the digest (truncated after stripping padding).

    Raises:
        ValueError: If style is not "hex" or "ga4gh".
    """
    digest = hashlib.sha512(canonical_bytes).digest()
    if style == "hex":
        return digest.hex()
    elif style == "ga4gh":
        b64 = base64.urlsafe_b64encode(digest).decode("ascii").rstrip("=")
        return f"ga4gh:pg.{b64[:24]}"
    else:
        raise ValueError(f"Unknown style: {style!r}")


def identify_graph(graph: AbstractGraph, style: str = "hex") -> str:
    return compute_identifier(serialize(graph), style)


def identify_from_string(
    text: str,
    format: str,
    overlap_policy: str = "discard",
    style: str = "hex",
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
    return identify_graph(graph, style)
