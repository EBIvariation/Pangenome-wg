"""GFA v1 parser."""

from __future__ import annotations

from pangenome_id.model import AbstractGraph, Node, Path, Step
from pangenome_id.parsers.base import BaseParser


class GFA1Parser(BaseParser):
    def parse(self, filepath: str) -> AbstractGraph:
        with open(filepath, "r") as fh:
            return self._parse_lines(fh)

    def parse_string(self, text: str) -> AbstractGraph:
        return self._parse_lines(iter(text.splitlines()))

    def _parse_lines(self, lines) -> AbstractGraph:
        name_to_id: dict[str, str] = {}
        nodes: list[Node] = []
        raw_edges: list[tuple] = []  # (from_name, from_orient, to_name, to_orient, cigar)
        raw_paths: list[tuple] = []  # (path_name, segments_str)

        for line in lines:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            record_type = fields[0]

            if record_type == "S":
                name = fields[1]
                seq = fields[2] if len(fields) > 2 else "*"
                node_id, normalized = self.node_id_from_sequence(seq, name)
                name_to_id[name] = node_id
                nodes.append(Node(id=node_id, sequence=normalized))

            elif record_type == "L":
                raw_edges.append((fields[1], fields[2], fields[3], fields[4], fields[5] if len(fields) > 5 else "*"))

            elif record_type == "P":
                raw_paths.append((fields[1], fields[2]))

        edges_seen: set = set()
        edges = []
        for from_name, from_orient, to_name, to_orient, cigar in raw_edges:
            node_a = name_to_id[from_name]
            node_b = name_to_id[to_name]
            overlap = self._resolve_overlap(cigar)
            edge = self.canonical_edge(node_a, from_orient, node_b, to_orient, overlap)
            key = (edge.node_a, edge.orient_a, edge.node_b, edge.orient_b, edge.overlap)
            if key not in edges_seen:
                edges_seen.add(key)
                edges.append(edge)

        paths = []
        for path_name, segments_str in raw_paths:
            steps = []
            for seg in segments_str.split(","):
                seg = seg.strip()
                orient = seg[-1]
                seg_name = seg[:-1]
                node_id = name_to_id[seg_name]
                steps.append(Step(node_id=node_id, orient=orient))
            paths.append(Path(name=path_name, steps=steps, is_reference=False))

        return AbstractGraph(
            nodes=nodes,
            edges=edges,
            paths=paths,
            overlap_policy=self.overlap_policy,
        )
