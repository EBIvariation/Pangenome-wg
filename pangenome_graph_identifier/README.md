# pangenome-id

A Python library and CLI tool that produces a **stable, deterministic, format-independent** hash identifier for pangenome graphs. The same graph, regardless of whether it is encoded in GFA v1 or GFA v2, always produces the same identifier.

---

## How it works

The algorithm runs in three phases:

```
GFA file  ──[parse]──►  Abstract Graph Model  ──[canonicalize]──►  bytes  ──[SHA-512]──►  identifier
```

1. **Parse** — the input file is read into a format-agnostic Abstract Graph Model (nodes, edges, paths).
2. **Canonicalize** — the model is serialized to a deterministic binary byte string (nodes and edges are sorted; path step order is preserved).
3. **Hash** — the byte string is hashed with SHA-512 to produce the identifier.

### Why format-independent?

Node identities are derived from their **sequence content**, not from the name assigned in the source file:

```
node.id = sha512(normalized_sequence)[:16]
```

Two formats that represent the same genomic sequence under different node names therefore produce the same node id, and thus the same graph hash.

---

## Installation

```bash
cd pangenome_graph_identifier
pip install -e .
```

Python ≥ 3.10 required. No external dependencies.

---

## CLI usage

```
pangenome-id [-h] [--format {gfa1,gfa2,auto}]
             [--overlap-policy {discard,length_only,full_cigar}]
             [--style {hex,ga4gh}]
             [--verbose]
             file
```

**Examples:**

```bash
# Auto-detect format, output full hex identifier
pangenome-id graph.gfa

# Compact GA4GH-style identifier
pangenome-id graph.gfa --style ga4gh

# GFA v2 file with verbose stats
pangenome-id graph.gfa2 --format gfa2 --verbose

# Include overlap lengths in the hash
pangenome-id graph.gfa --overlap-policy length_only
```

**Verbose output** (printed to stderr):

```
Format:         GFA v1
Nodes:          1042
Edges:          1387
Paths:          5
Overlap policy: discard
Canonical size: 48291 bytes
```

---

## Output styles

| Style | Example |
|---|---|
| `hex` (default) | `3a7f1b9c2e4d8f0a...` (128 hex chars) |
| `ga4gh` | `ga4gh:pg.SXaB3kqLm9vNpRt7YwZcAe` |

---

## Overlap policy

The overlap policy controls how CIGAR overlap strings on edges affect the identifier. It is embedded in the canonical bytes, so graphs hashed under different policies always produce different identifiers.

| Policy | Behaviour |
|---|---|
| `discard` (default) | Overlaps are ignored. Use for variation graphs and blunt-joined assemblies. |
| `length_only` | CIGAR is normalized to its total consumed length (`5M2I3M` → `10`). |
| `full_cigar` | CIGAR string is included verbatim. |

---

## Format support

### GFA v1
Processes `S` (segment), `L` (link), and `P` (path) lines. All other lines (`H`, `C`, `#`, …) are skipped.

### GFA v2
Processes `S` (segment), `E` (edge), and `O` (ordered group) lines. `F` (fragment) and `U` (unordered group) records are skipped as they have no equivalent in the topology model. Edge ids are discarded.

Format is auto-detected from the `VN:Z:` header tag, or from the `.gfa` file extension.

---

## Library usage

```python
from pangenome_id import identify_from_string
from pangenome_id.parsers.gfa1 import GFA1Parser
from pangenome_id.canonicalize import serialize
from pangenome_id.hasher import compute_identifier, identify_graph

# From a file
graph = GFA1Parser().parse("graph.gfa")
identifier = identify_graph(graph, style="ga4gh")

# From a string (useful in tests)
identifier = identify_from_string(gfa_text, format="gfa1", style="hex")

# Low-level: canonical bytes first
canonical = serialize(graph)
identifier = compute_identifier(canonical, style="hex")
```

---

## Project layout

```
pangenome_graph_identifier/
├── pangenome_id/
│   ├── model.py          # Abstract Graph Model dataclasses
│   ├── canonicalize.py   # AGM → deterministic canonical bytes
│   ├── hasher.py         # canonical bytes → identifier string
│   ├── parsers/
│   │   ├── base.py       # shared normalization and edge canonicalization
│   │   ├── gfa1.py       # GFA v1 parser
│   │   └── gfa2.py       # GFA v2 parser
│   └── cli.py            # CLI entry point
└── tests/
    ├── test_model.py
    ├── test_canonicalize.py
    ├── test_gfa1.py
    ├── test_gfa2.py
    └── test_cross_format.py   # cross-format identity tests
```

---

## Running the tests

```bash
pip install pytest
python -m pytest tests/ -v
```
