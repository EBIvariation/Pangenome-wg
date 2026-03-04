"""GFA v2 parser."""

from __future__ import annotations

from pangenome_id.model import AbstractGraph, Node, Path, Step
from pangenome_id.parsers.base import BaseParser


class GFA2Parser(BaseParser):
    def parse(self, filepath: str) -> AbstractGraph:
        with open(filepath, "r") as fh:
            return self._parse_lines(fh)

    def parse_string(self, text: str) -> AbstractGraph:
        return self._parse_lines(iter(text.splitlines()))

    def _parse_lines(self, lines) -> AbstractGraph:
        name_to_id: dict[str, str] = {}
        nodes: list[Node] = []
        raw_edges: list[tuple] = []
        raw_paths: list[tuple] = []

        for line in lines:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            record_type = fields[0]

            if record_type == "S":
                # S <id> <len> <seq> [tags...]
                seg_id = fields[1]
                seq = fields[3] if len(fields) > 3 else "*"
                node_id, normalized = self.node_id_from_sequence(seq, seg_id)
                name_to_id[seg_id] = node_id
                nodes.append(Node(id=node_id, sequence=normalized))

            elif record_type == "E":
                # E <edge_id> <sid1><orient> <sid2><orient> <beg1> <end1> <beg2> <end2> <align>
                # Discard edge_id (fields[1])
                sid1_field = fields[2]
                sid2_field = fields[3]
                beg1 = fields[4]
                end1 = fields[5]
                align = fields[8] if len(fields) > 8 else "*"

                orient1 = sid1_field[-1]
                seg1 = sid1_field[:-1]
                orient2 = sid2_field[-1]
                seg2 = sid2_field[:-1]

                overlap = self._resolve_overlap_gfa2(beg1, end1, align)
                raw_edges.append((seg1, orient1, seg2, orient2, overlap))

            elif record_type == "O":
                # O <group_id> <ref1><orient> <ref2><orient> ...
                group_id = fields[1]
                refs_str = fields[2] if len(fields) > 2 else ""
                raw_paths.append((group_id, refs_str))

            # Skip F, U, H, and everything else

        edges_seen: set = set()
        edges = []
        for seg1, orient1, seg2, orient2, overlap in raw_edges:
            node_a = name_to_id[seg1]
            node_b = name_to_id[seg2]
            edge = self.canonical_edge(node_a, orient1, node_b, orient2, overlap)
            key = (edge.node_a, edge.orient_a, edge.node_b, edge.orient_b, edge.overlap)
            if key not in edges_seen:
                edges_seen.add(key)
                edges.append(edge)

        paths = []
        for group_id, refs_str in raw_paths:
            steps = []
            for ref in refs_str.split():
                ref = ref.strip()
                if not ref:
                    continue
                orient = ref[-1]
                seg_name = ref[:-1]
                node_id = name_to_id[seg_name]
                steps.append(Step(node_id=node_id, orient=orient))
            paths.append(Path(name=group_id, steps=steps, is_reference=False))

        return AbstractGraph(
            nodes=nodes,
            edges=edges,
            paths=paths,
            overlap_policy=self.overlap_policy,
        )

    def _resolve_overlap_gfa2(self, beg1: str, end1: str, align: str) -> str | None:
        """Compute overlap from GFA2 position fields, respecting overlap_policy."""
        if self.overlap_policy == "discard":
            return None
        if self.overlap_policy == "full_cigar":
            return align if align != "*" else None
        if self.overlap_policy == "length_only":
            # overlap_length = end1 - beg1; strip trailing '$' from position strings
            try:
                b = int(beg1.rstrip("$"))
                e = int(end1.rstrip("$"))
                return str(e - b)
            except ValueError:
                return None
        return None
