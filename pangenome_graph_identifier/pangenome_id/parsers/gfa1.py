"""GFA v1 / v1.2 parser."""

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
        raw_edges: list[tuple] = []  # (from_name, from_orient, to_name, to_orient, overlap_field, line_type)
        raw_paths: list[tuple] = []  # (path_name, raw_str, fmt)  fmt = "P" | "W"

        for line in lines:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            record_type = fields[0]

            if record_type == "S":
                name = fields[1]
                seq = fields[2] if len(fields) > 2 else "*"
                node = self.node_from_sequence(seq, name)
                name_to_id[name] = node.id
                nodes.append(node)

            elif record_type == "L":
                # L <from> <from_orient> <to> <to_orient> <cigar>
                raw_edges.append((fields[1], fields[2], fields[3], fields[4],
                                   fields[5] if len(fields) > 5 else "*", "L"))

            elif record_type == "J":
                # J <from> <from_orient> <to> <to_orient> <distance> [tags]
                raw_edges.append((fields[1], fields[2], fields[3], fields[4],
                                   fields[5] if len(fields) > 5 else "*", "J"))

            elif record_type == "P":
                # P <path_name> <seg_names> <cigars>
                raw_paths.append((fields[1], fields[2], "P"))

            elif record_type == "W":
                # W <sample_id> <hap_index> <seq_id> <seq_start> <seq_end> <walk>
                path_name = f"{fields[1]}#{fields[2]}#{fields[3]}"
                walk = fields[6] if len(fields) > 6 else ""
                raw_paths.append((path_name, walk, "W"))

        edges_seen: set = set()
        edges = []
        for from_name, from_orient, to_name, to_orient, overlap_field, line_type in raw_edges:
            node_a = name_to_id[from_name]
            node_b = name_to_id[to_name]
            if line_type == "J":
                overlap = self._resolve_jump_distance(overlap_field)
            else:
                overlap = self._resolve_overlap(overlap_field)
            edge = self.canonical_edge(node_a, from_orient, node_b, to_orient, overlap)
            key = (edge.node_a, edge.orient_a, edge.node_b, edge.orient_b, edge.overlap)
            if key not in edges_seen:
                edges_seen.add(key)
                edges.append(edge)

        paths = []
        for path_name, raw_str, fmt in raw_paths:
            steps = []
            if fmt == "P":
                for seg in raw_str.split(","):
                    seg = seg.strip()
                    orient = seg[-1]
                    seg_name = seg[:-1]
                    steps.append(Step(node_id=name_to_id[seg_name], orient=orient))
            else:  # "W"
                i = 0
                while i < len(raw_str):
                    if raw_str[i] in "><":
                        orient = "+" if raw_str[i] == ">" else "-"
                        j = i + 1
                        while j < len(raw_str) and raw_str[j] not in "><":
                            j += 1
                        steps.append(Step(node_id=name_to_id[raw_str[i + 1:j]], orient=orient))
                        i = j
                    else:
                        i += 1
            paths.append(Path(name=path_name, steps=steps, is_reference=False))

        return AbstractGraph(
            nodes=nodes,
            edges=edges,
            paths=paths,
            overlap_policy=self.overlap_policy,
        )
