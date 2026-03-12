"""Abstract Graph Model dataclasses."""

from dataclasses import dataclass, field


@dataclass
class Node:
    id: str  # stable identifier derived from sequence (or source name if sequence is "*")
    sequence: str | None = None  # normalized sequence; None for "*" wildcard nodes


@dataclass
class Edge:
    node_a: str      # node id
    orient_a: str    # "+" or "-"
    node_b: str      # node id
    orient_b: str    # "+" or "-"
    overlap: str | None = None  # populated based on overlap_policy


@dataclass
class Step:
    node_id: str
    orient: str  # "+" or "-"


@dataclass
class Path:
    name: str
    steps: list[Step]
    is_reference: bool = False
    _topology_digest: str | None = field(default=None, compare=False, repr=False)
    _sequence_digest: str | None = field(default=None, compare=False, repr=False)


@dataclass
class AbstractGraph:
    nodes: list[Node]
    edges: list[Edge]
    paths: list[Path]
    overlap_policy: str = "discard"  # "discard" | "length_only" | "full_cigar"
