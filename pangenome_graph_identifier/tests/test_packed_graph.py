"""Tests for vg PackedGraph parser (no real .pg file or libbdsg needed)."""

from unittest.mock import MagicMock, patch
import pytest

from pangenome_id.hasher import sha512t24u
from pangenome_id.parsers.packed_graph import PackedGraphParser


def _node_id(seq: str) -> str:
    return sha512t24u(seq.upper().encode("ascii"))


def _make_mock_graph(node_specs, edge_specs, path_specs):
    """
    Build a mock bdsg-style graph.

    node_specs: list of (int_id, sequence)
    edge_specs: list of (int_id_a, orient_a, int_id_b, orient_b)
                where orient is True=reverse, False=forward
    path_specs: list of (path_name, [(int_id, is_reverse), ...])
    """
    g = MagicMock()

    # Build handle mocks keyed by int id
    handles: dict[int, MagicMock] = {}
    for int_id, seq in node_specs:
        h = MagicMock()
        h.configure_mock(name=f"handle_{int_id}")
        handles[int_id] = h

    def _get_id(h):
        for iid, hh in handles.items():
            if hh is h:
                return iid
        raise KeyError(h)

    def _get_sequence(h):
        iid = _get_id(h)
        return dict(node_specs)[iid]

    def _get_is_reverse(h):
        return False  # for_each_handle yields forward-strand handles only

    g.get_id.side_effect = _get_id
    g.get_sequence.side_effect = _get_sequence
    g.get_is_reverse.side_effect = _get_is_reverse

    def _for_each_handle(fn):
        for h in handles.values():
            fn(h)

    g.for_each_handle.side_effect = _for_each_handle

    # Build edge adjacency: int_id -> (go_left=False list, go_left=True list)
    adjacency: dict[int, tuple[list, list]] = {iid: ([], []) for iid, _ in node_specs}
    for iid_a, orient_a, iid_b, orient_b in edge_specs:
        # forward-strand handle for iid_a has iid_b on the right
        nbr_b = MagicMock()
        nbr_b.configure_mock(name=f"handle_{iid_b}_orient{orient_b}")
        nbr_b_is_rev = orient_b
        adjacency[iid_a][0].append((iid_b, nbr_b, nbr_b_is_rev))

        # and iid_b has iid_a on the left (reverse direction)
        nbr_a = MagicMock()
        nbr_a.configure_mock(name=f"handle_{iid_a}_orient{orient_a}_from{iid_b}")
        nbr_a_is_rev = orient_a
        adjacency[iid_b][1].append((iid_a, nbr_a, nbr_a_is_rev))

    # Track extra neighbor handles so get_id/get_is_reverse work for them
    extra_handles: dict[int, list[tuple[MagicMock, bool]]] = {}
    for iid_a, orient_a, iid_b, orient_b in edge_specs:
        extra_handles.setdefault(iid_b, [])
        extra_handles.setdefault(iid_a, [])

    # Rebuild get_id to cover neighbour handles too
    nbr_id_map: dict[int, tuple[MagicMock, bool]] = {}  # id(mock) -> (int_id, is_rev)
    for iid_a, orient_a, iid_b, orient_b in edge_specs:
        right_list = adjacency[iid_a][0]
        for entry in right_list:
            if entry[0] == iid_b:
                nbr_id_map[id(entry[1])] = (iid_b, orient_b)
                break
        left_list = adjacency[iid_b][1]
        for entry in left_list:
            if entry[0] == iid_a:
                nbr_id_map[id(entry[1])] = (iid_a, orient_a)
                break

    def _get_id2(h):
        if id(h) in nbr_id_map:
            return nbr_id_map[id(h)][0]
        for iid, hh in handles.items():
            if hh is h:
                return iid
        raise KeyError(h)

    def _get_is_reverse2(h):
        if id(h) in nbr_id_map:
            return nbr_id_map[id(h)][1]
        return False

    g.get_id.side_effect = _get_id2
    g.get_is_reverse.side_effect = _get_is_reverse2

    def _follow_edges(h, go_left, fn):
        iid = _get_id2(h)
        neighbours = adjacency.get(iid, ([], []))[1 if go_left else 0]
        for _, nbr_mock, _ in neighbours:
            fn(nbr_mock)

    g.follow_edges.side_effect = _follow_edges

    # Paths
    path_handles = []
    for path_name, steps in path_specs:
        ph = MagicMock()
        ph.configure_mock(name=f"path_{path_name}")
        path_handles.append((ph, path_name, steps))

    def _for_each_path_handle(fn):
        for ph, _, _ in path_handles:
            fn(ph)

    g.for_each_path_handle.side_effect = _for_each_path_handle
    g.get_path_name.side_effect = lambda ph: next(
        name for p, name, _ in path_handles if p is ph
    )

    def _for_each_step_in_path(ph, fn):
        steps = next(s for p, _, s in path_handles if p is ph)
        for int_id, is_rev in steps:
            sh = MagicMock()
            nh = MagicMock()
            nbr_id_map[id(nh)] = (int_id, is_rev)
            sh.configure_mock(name=f"step_{int_id}")
            g.get_handle_of_step.side_effect = lambda s: s._nh
            sh._nh = nh
            fn(sh)

    g.get_handle_of_step.side_effect = lambda sh: sh._nh
    g.for_each_step_in_path.side_effect = _for_each_step_in_path

    return g


# --- Fixtures ---

NODE_SPECS = [(1, "ACGT"), (2, "TTGC")]
EDGE_SPECS = [(1, False, 2, False)]  # 1+ -> 2+
PATH_SPECS = [("p1", [(1, False), (2, False)])]


def _parser():
    return PackedGraphParser()


def test_node_count():
    g = _parser()._extract(_make_mock_graph(NODE_SPECS, EDGE_SPECS, PATH_SPECS))
    assert len(g.nodes) == 2


def test_node_ids_are_sequence_derived():
    g = _parser()._extract(_make_mock_graph(NODE_SPECS, EDGE_SPECS, PATH_SPECS))
    ids = {n.id for n in g.nodes}
    assert _node_id("ACGT") in ids
    assert _node_id("TTGC") in ids


def test_edges():
    g = _parser()._extract(_make_mock_graph(NODE_SPECS, EDGE_SPECS, PATH_SPECS))
    assert len(g.edges) == 1
    edge = g.edges[0]
    assert {edge.node_a, edge.node_b} == {_node_id("ACGT"), _node_id("TTGC")}
    assert edge.overlap is None  # PackedGraph has no overlap info


def test_paths():
    g = _parser()._extract(_make_mock_graph(NODE_SPECS, EDGE_SPECS, PATH_SPECS))
    assert len(g.paths) == 1
    path = g.paths[0]
    assert path.name == "p1"
    assert len(path.steps) == 2
    assert path.steps[0].node_id == _node_id("ACGT")
    assert path.steps[0].orient == "+"
    assert path.steps[1].node_id == _node_id("TTGC")
    assert path.steps[1].orient == "+"


def test_reverse_step_orientation():
    path_specs = [("p1", [(1, False), (2, True)])]
    g = _parser()._extract(_make_mock_graph(NODE_SPECS, EDGE_SPECS, path_specs))
    assert g.paths[0].steps[1].orient == "-"


def test_no_duplicate_edges():
    # A linear 3-node graph: each edge is visited from both ends (go_right from one,
    # go_left from the other) — deduplication should keep exactly 2 unique edges.
    node_specs = [(1, "ACGT"), (2, "TTGC"), (3, "GGAA")]
    edge_specs = [(1, False, 2, False), (2, False, 3, False)]
    g = _parser()._extract(_make_mock_graph(node_specs, edge_specs, []))
    assert len(g.edges) == 2


def test_import_error_message():
    parser = PackedGraphParser()
    with patch.dict("sys.modules", {"bdsg": None, "bdsg.bdsg": None}):
        with pytest.raises(ImportError, match="libbdsg"):
            parser.parse("fake.pg")


def test_overlap_policy_stored():
    parser = PackedGraphParser(overlap_policy="full_cigar")
    g = parser._extract(_make_mock_graph(NODE_SPECS, EDGE_SPECS, PATH_SPECS))
    assert g.overlap_policy == "full_cigar"


# ---------------------------------------------------------------------------
# Integration test with a real .pg file (skipped if libbdsg not installed)
# ---------------------------------------------------------------------------

def test_real_pg_file(tmp_path):
    """Build a small PackedGraph with bdsg, serialize to disk, parse with PackedGraphParser."""
    bdsg = pytest.importorskip("bdsg.bdsg", reason="libbdsg not installed")
    PackedGraph = bdsg.PackedGraph

    # Build a graph: node1 (ACGT) → node2 (TTGC), path p1: node1+ node2+
    g = PackedGraph()
    h1 = g.create_handle("ACGT")
    h2 = g.create_handle("TTGC")
    g.create_edge(h1, h2)
    ph = g.create_path_handle("p1")
    g.append_step(ph, h1)
    g.append_step(ph, h2)

    pg_file = str(tmp_path / "test.pg")
    g.serialize(pg_file)

    result = PackedGraphParser().parse(pg_file)

    assert len(result.nodes) == 2
    ids = {n.id for n in result.nodes}
    assert _node_id("ACGT") in ids
    assert _node_id("TTGC") in ids

    assert len(result.edges) == 1
    edge = result.edges[0]
    assert {edge.node_a, edge.node_b} == {_node_id("ACGT"), _node_id("TTGC")}
    assert edge.overlap is None

    assert len(result.paths) == 1
    path = result.paths[0]
    assert path.name == "p1"
    assert len(path.steps) == 2
    assert path.steps[0].node_id == _node_id("ACGT")
    assert path.steps[0].orient == "+"
    assert path.steps[1].node_id == _node_id("TTGC")
    assert path.steps[1].orient == "+"

