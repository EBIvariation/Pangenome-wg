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
node.id = sha512t24u(normalized_sequence)
```

`sha512t24u` is the GA4GH standard digest: SHA-512 of the UTF-8 sequence, first 24 bytes, URL-safe base64 encoded (32 ASCII chars).

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
             [--verbose]
             file
```

**Examples:**

```bash
# Auto-detect format
pangenome-id graph.gfa

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

## Output format

Identifiers are always in GA4GH format:

```
ga4gh:pg.SXaB3kqLm9vNpRt7YwZcAeQf
```

The suffix is the `sha512t24u` digest of the canonical bytes: first 24 bytes of SHA-512, URL-safe base64 encoded (32 chars).

---

## What contributes to the identifier

The identifier is the SHA-512 hash of a canonical byte string assembled from these components, in order:

### Delimiters

Three control bytes act as separators. They never appear in node ids (URL-safe base64), orientations (`+`/`-`), overlap strings, or path names.

| Byte | Role |
|---|---|
| `\x00` | Field separator — between fields within a record |
| `\x01` | Record separator — between individual nodes / edges / paths |
| `\x02` | Section separator — between the nodes, edges, and paths sections |

### Header (2 bytes)

| Field | Size | Notes |
|---|---|---|
| Format version | 1 B | Currently `0x02` |
| Overlap policy | 1 B | `discard=0x00`, `length_only=0x01`, `full_cigar=0x02` |

### Nodes section

Each node sorted by id, followed by the section separator:

```
id \x01  id \x01  ...  \x02
```

Node id is `sha512t24u(normalized_sequence)` — 32 URL-safe base64 chars. For `*`-sequence nodes, falls back to the source node name.

### Edges section

Each edge sorted by `(node_a, orient_a, node_b, orient_b)`, followed by the section separator:

```
node_a \x00 orient_a \x00 node_b \x00 orient_b \x00 overlap \x01  ...  \x02
```

`overlap` is empty when policy is `discard`, a length string when `length_only`, raw CIGAR when `full_cigar`.

### Paths section

Each path sorted by name (step order preserved), followed by the section separator:

```
name \x00 is_reference \x00 step_id \x00 orient [ \x00 step_id \x00 orient ... ] \x01  ...  \x02
```

`is_reference` is the ASCII character `1` or `0`.

> **What is NOT included:** source file names, original node names (they are replaced by sequence-derived ids), comment/header lines, fragment records (`F`), unordered groups (`U`), and edge ids from GFA v2.

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

### GFA v1 / v1.2

| Record | Processed as | Notes |
|---|---|---|
| `S` | Node | Sequence normalized to uppercase DNA (U→T); id = `sha512t24u(sequence)` |
| `L` | Edge | CIGAR overlap resolved per overlap policy |
| `J` | Edge | Jump distance (plain integer) resolved per overlap policy |
| `P` | Path | Comma-separated `seg+`/`seg-` segments |
| `W` | Path | Walk string `>seg<seg…`; name built as `sample_id#hap_index#seq_id` (PanSN) |

All other lines (`H`, `C`, `#`, …) are skipped.

### GFA v2
| Record | Processed as | Notes |
|---|---|---|
| `S` | Node | |
| `E` | Edge | Edge ids are discarded |
| `O` | Path | Ordered group |

`F` (fragment) and `U` (unordered group) records are skipped as they have no equivalent in the topology model.

Format is auto-detected from the `VN:Z:` header tag, or from the `.gfa` file extension.

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
