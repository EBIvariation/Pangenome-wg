"""vg PackedGraph (.pg) parser — requires libbdsg (pip install libbdsg)."""
from pangenome_id.model import AbstractGraph, Node, Edge, Path, Step
from pangenome_id.parsers.base import BaseParser


class PackedGraphParser(BaseParser):
    def parse(self, filepath: str) -> AbstractGraph:
        try:
            from bdsg.bdsg import PackedGraph
        except ImportError:
            raise ImportError(
                "libbdsg is required to parse PackedGraph files. "
                "Install it with: pip install libbdsg"
            )
        g = PackedGraph()
        g.deserialize(filepath)
        return self._extract(g)

    def _extract(self, g) -> AbstractGraph:
        # 1. Nodes: map bdsg int ID → sha512t24u node_id
        int_to_id: dict[int, str] = {}
        nodes: list[Node] = []

        def _collect_node(h):
            seq = g.get_sequence(h)
            nid = self.node_id_from_sequence(seq, fallback_name=str(g.get_id(h)))
            int_to_id[g.get_id(h)] = nid
            nodes.append(Node(id=nid))
            return True

        g.for_each_handle(_collect_node)

        # 2. Edges: follow both directions per node; canonical_edge + edges_seen deduplicate
        edges_seen: set = set()
        edges: list[Edge] = []

        def _follow(h, go_left):
            orient_a = "-" if g.get_is_reverse(h) else "+"
            nid_a = int_to_id[g.get_id(h)]

            def _add_edge(nbr):
                orient_b = "-" if g.get_is_reverse(nbr) else "+"
                nid_b = int_to_id[g.get_id(nbr)]
                # When going left, nbr is the FROM node; swap to get correct edge direction
                if go_left:
                    edge = self.canonical_edge(nid_b, orient_b, nid_a, orient_a)
                else:
                    edge = self.canonical_edge(nid_a, orient_a, nid_b, orient_b)
                key = (edge.node_a, edge.orient_a, edge.node_b, edge.orient_b, edge.overlap)
                if key not in edges_seen:
                    edges_seen.add(key)
                    edges.append(edge)
                return True

            g.follow_edges(h, go_left, _add_edge)
            return True

        g.for_each_handle(lambda h: _follow(h, False) and _follow(h, True))

        # 3. Paths
        paths: list[Path] = []

        def _collect_path(ph):
            steps: list[Step] = []

            def _collect_step(sh):
                nh = g.get_handle_of_step(sh)
                orient = "-" if g.get_is_reverse(nh) else "+"
                steps.append(Step(node_id=int_to_id[g.get_id(nh)], orient=orient))
                return True

            g.for_each_step_in_path(ph, _collect_step)
            paths.append(Path(name=g.get_path_name(ph), steps=steps))
            return True

        g.for_each_path_handle(_collect_path)

        return AbstractGraph(nodes=nodes, edges=edges, paths=paths,
                             overlap_policy=self.overlap_policy)
