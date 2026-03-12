"""GFA v1 / v1.2 parser."""

import base64
import hashlib

from pangenome_id.canonicalize import reverse_complement
from pangenome_id.model import AbstractGraph, Node, Path, Step
from pangenome_id.parsers.base import BaseParser, _open_file


class GFA1Parser(BaseParser):
    def parse(self, filepath: str) -> AbstractGraph:
        with _open_file(filepath) as fh:
            return self._parse_lines(fh)

    def parse_string(self, text: str) -> AbstractGraph:
        return self._parse_lines(iter(text.splitlines()))

    def _parse_lines(self, lines) -> AbstractGraph:
        name_to_id: dict[str, str] = {}
        seq_by_name: dict[str, str | None] = {}
        nodes: list[Node] = []
        raw_edges: list[tuple] = []  # (from_name, from_orient, to_name, to_orient, overlap_field, line_type)
        raw_paths: list[tuple] = []  # (path_name, seg_str) — P lines only
        paths: list[Path] = []

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
                seq_by_name[name] = node.sequence
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
                raw_paths.append((fields[1], fields[2]))

            elif record_type == "W":
                # W <sample_id> <hap_index> <seq_id> <seq_start> <seq_end> <walk>
                path_name = f"{fields[1]}#{fields[2]}#{fields[3]}"
                walk = fields[6] if len(fields) > 6 else ""
                t_h = hashlib.sha512()
                t_h.update(b'{"is_reference":false,"steps":[')
                s_h = hashlib.sha512()
                seq_ok = True
                first = True
                i = 0
                while i < len(walk):
                    if walk[i] in "><":
                        orient = "+" if walk[i] == ">" else "-"
                        j = i + 1
                        while j < len(walk) and walk[j] not in "><":
                            j += 1
                        node_name = walk[i + 1:j]
                        node_id = name_to_id[node_name]
                        if not first:
                            t_h.update(b',')
                        t_h.update(b'{"id":"')
                        t_h.update(node_id.encode())
                        t_h.update(b'","orient":"')
                        t_h.update(orient.encode())
                        t_h.update(b'"}')
                        first = False
                        if seq_ok:
                            seq = seq_by_name.get(node_name)
                            if seq is None:
                                seq_ok = False
                            else:
                                s_h.update((reverse_complement(seq) if orient == "-" else seq).encode("ascii"))
                        i = j
                    else:
                        i += 1
                t_h.update(b']}')
                topo_d = base64.urlsafe_b64encode(t_h.digest()[:24]).decode("ascii")
                seq_d = base64.urlsafe_b64encode(s_h.digest()[:24]).decode("ascii") if seq_ok else None
                paths.append(Path(name=path_name, steps=[], is_reference=False,
                                  _topology_digest=topo_d, _sequence_digest=seq_d))

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

        for path_name, raw_str in raw_paths:
            steps = []
            for seg in raw_str.split(","):
                seg = seg.strip()
                orient = seg[-1]
                seg_name = seg[:-1]
                steps.append(Step(node_id=name_to_id[seg_name], orient=orient))
            paths.append(Path(name=path_name, steps=steps, is_reference=False))

        return AbstractGraph(
            nodes=nodes,
            edges=edges,
            paths=paths,
            overlap_policy=self.overlap_policy,
        )
