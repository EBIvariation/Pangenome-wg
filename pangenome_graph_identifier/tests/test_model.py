"""Tests for model and base parser utilities."""

import hashlib
import warnings

import pytest

from pangenome_id.parsers.base import BaseParser
from pangenome_id.parsers.gfa1 import GFA1Parser


class ConcreteParser(BaseParser):
    """Minimal concrete subclass for testing BaseParser methods."""
    def parse(self, filepath):
        raise NotImplementedError


def test_normalize_sequence_uppercase():
    p = ConcreteParser()
    assert p.normalize_sequence("acgt") == "ACGT"


def test_normalize_sequence_u_to_t():
    p = ConcreteParser()
    assert p.normalize_sequence("AUCG") == "ATCG"


def test_normalize_sequence_strips_whitespace():
    p = ConcreteParser()
    assert p.normalize_sequence("  ACGT  ") == "ACGT"


def test_node_id_from_sequence():
    p = ConcreteParser()
    node_id, normalized = p.node_id_from_sequence("ACGT", "fallback")
    expected = hashlib.sha256("ACGT".encode("ascii")).hexdigest()[:16]
    assert node_id == expected
    assert normalized == "ACGT"


def test_node_id_lowercase_normalized():
    p = ConcreteParser()
    node_id_lower, _ = p.node_id_from_sequence("acgt", "x")
    node_id_upper, _ = p.node_id_from_sequence("ACGT", "x")
    assert node_id_lower == node_id_upper


def test_node_id_u_to_t_normalized():
    p = ConcreteParser()
    node_id_u, _ = p.node_id_from_sequence("AUCG", "x")
    node_id_t, _ = p.node_id_from_sequence("ATCG", "x")
    assert node_id_u == node_id_t


def test_node_id_star_uses_fallback():
    p = ConcreteParser()
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        node_id, normalized = p.node_id_from_sequence("*", "my_node")
    assert node_id == "my_node"
    assert normalized == ""
    assert len(w) == 1
    assert "my_node" in str(w[0].message)


def test_flip_orient():
    p = ConcreteParser()
    assert p.flip_orient("+") == "-"
    assert p.flip_orient("-") == "+"


def test_canonical_edge_is_deterministic():
    p = ConcreteParser()
    # Edge and its reverse complement should produce the same canonical edge
    e1 = p.canonical_edge("aaa", "+", "bbb", "+")
    e2 = p.canonical_edge("bbb", "-", "aaa", "-")
    assert e1.node_a == e2.node_a
    assert e1.orient_a == e2.orient_a
    assert e1.node_b == e2.node_b
    assert e1.orient_b == e2.orient_b
