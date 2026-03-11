"""Abstract base parser."""

import re
import warnings
from abc import ABC, abstractmethod

from pangenome_id.hasher import sha512t24u
from pangenome_id.model import AbstractGraph, Edge, Node


class BaseParser(ABC):
    def __init__(self, overlap_policy: str = "discard"):
        self.overlap_policy = overlap_policy

    @abstractmethod
    def parse(self, filepath: str) -> AbstractGraph:
        ...

    def normalize_sequence(self, seq: str) -> str:
        return seq.strip().upper()

    def flip_orient(self, orient: str) -> str:
        return "-" if orient == "+" else "+"

    def node_id_from_sequence(self, seq: str, fallback_name: str) -> str:
        """Return node_id derived from sequence.

        If seq is '*', use fallback_name as node_id.
        """
        return self.node_from_sequence(seq, fallback_name).id

    def node_from_sequence(self, seq: str, fallback_name: str) -> Node:
        """Return a Node with id and sequence derived from seq.

        If seq is '*', id is fallback_name and sequence is None.
        """
        if seq == "*":
            warnings.warn(
                f"Node '{fallback_name}' has no sequence; sequence-derived identity unavailable.",
                stacklevel=3,
            )
            return Node(id=fallback_name, sequence=None)
        normalized = self.normalize_sequence(seq)
        return Node(id=sha512t24u(normalized.encode("ascii")), sequence=normalized)

    def canonical_edge(
        self,
        node_a: str,
        orient_a: str,
        node_b: str,
        orient_b: str,
        overlap: str | None = None,
    ) -> Edge:
        """Return the edge in canonical (lexicographically smallest) form.

        An undirected overlap edge (A+, B+) and its reverse complement (B-, A-) are
        biologically the same adjacency. To ensure a unique representation regardless
        of which endpoint a parser encounters first, we always store the tuple that
        sorts lexicographically smaller between the forward form (node_a, orient_a,
        node_b, orient_b) and the reverse-complement form (node_b, ~orient_b, node_a,
        ~orient_a). Comparison is purely lexicographic on the 4-tuple.
        """
        fwd = (node_a, orient_a, node_b, orient_b)
        rev = (node_b, self.flip_orient(orient_b), node_a, self.flip_orient(orient_a))
        a, oa, b, ob = min(fwd, rev)
        return Edge(node_a=a, orient_a=oa, node_b=b, orient_b=ob, overlap=overlap)

    def _resolve_jump_distance(self, distance: str) -> str | None:
        """Apply overlap_policy to a Jump line distance field.

        Jump distances are plain integers (or '*'), not CIGAR strings.
        Both 'length_only' and 'full_cigar' return the value verbatim.
        """
        if self.overlap_policy == "discard":
            return None
        if distance == "*" or not distance:
            return None
        return distance

    def _resolve_overlap(self, cigar: str) -> str | None:
        """Apply overlap_policy to a CIGAR string. Returns processed overlap or None.

        - "discard"     → always returns None; overlaps are ignored entirely.
        - "full_cigar"  → returns the CIGAR string verbatim.
        - "length_only" → sums all consumed lengths in the CIGAR (e.g. "5M2I3M" → "10")
                          and returns the total as a decimal string.

        A CIGAR of "*" (absent) is treated as no overlap regardless of policy.
        """
        if self.overlap_policy == "discard":
            return None
        if cigar == "*" or not cigar:
            return None
        if self.overlap_policy == "full_cigar":
            return cigar
        if self.overlap_policy == "length_only":
            total = sum(
                int(n) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)
            )
            return str(total)
        return None
