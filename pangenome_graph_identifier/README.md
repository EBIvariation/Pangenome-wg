# pangenome-id

A Python library and CLI tool that produces a **stable, deterministic, format-independent** hash identifier for pangenome graphs. The same graph, regardless of whether it is encoded in GFA v1 or GFA v2, always produces the same identifier.

---

## How it works

The algorithm runs in three phases:

```
GFA file ──[parse]──► Abstract Graph Model ──[serialize]──► RFC-8785 JSON ──[sha512t24u]──► digest
```

1. **Parse** — the input file is read into a format-agnostic Abstract Graph Model (nodes, edges, paths).
2. **Serialize** — the model is serialized to three RFC-8785 canonical JSON documents: a topology document, per-path documents, and a complete graph document.
3. **Hash** — each document is hashed with SHA-512 (sha512t24u) to produce digests.

### Three JSON documents

**Topology document** (hashed → topology digest):
```json
{
  "edges": [{"from": "<node_id>", "from_orient": "+", "overlap": null, "to": "<node_id>", "to_orient": "+"}],
  "nodes": ["<node_id>", ...],
  "overlap_policy": "discard"
}
```

**Per-path document** (hashed → path topology digest):
```json
{
  "is_reference": false,
  "steps": [{"id": "<node_id>", "orient": "+"}, ...]
}
```

**Complete graph document** (the final output):
```json
{
  "graph_topology": "<32-char-b64>",
  "names": ["path_name_a", "path_name_b"],
  "paths": ["<topo_digest_a>", "<topo_digest_b>"],
  "sequences": ["<seq_digest_a>", "<seq_digest_b>"]
}
```

### Why format-independent?

Node identities are derived from their **sequence content**, not from the name assigned in the source file:

```
node.id = sha512t24u(normalized_sequence)
```

`sha512t24u` is the GA4GH standard digest: SHA-512 of the UTF-8 sequence, first 24 bytes, URL-safe base64 encoded (32 ASCII chars). Sequences are normalized to uppercase only — no other transformation is applied.

Two formats that represent the same genomic sequence under different node names therefore produce the same node id, and thus the same graph hash.

---

## Installation

```bash
cd pangenome_graph_identifier
pip install -e .
```

To also parse vg PackedGraph (`.pg`) files:

```bash
pip install -e ".[pg]"
```

Python ≥ 3.10 required.

---

## CLI usage

```
pangenome-id [-h] [--format {gfa1,gfa2,pg,auto}]
             [--overlap-policy {discard,length_only,full_cigar}]
             [--verbose]
             file
```

**Examples:**

```bash
# Auto-detect format (GFA v1/v2 from header; .pg by extension)
pangenome-id graph.gfa

# GFA v2 file with verbose stats
pangenome-id graph.gfa2 --format gfa2 --verbose

# vg PackedGraph file (requires pip install pangenome-id[pg])
pangenome-id graph.pg

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
```

---

## Output format

The CLI writes the RFC-8785 canonical JSON document to stdout:

```json
{"graph_topology":"SXaB3kqLm9vNpRt7YwZcAeQf12345678","names":["path1","path2"],"paths":["Ab1Cd2Ef3Gh4Ij5Kl6Mn7Op8Qr9St0Uv","Wx1Yz2Ab3Cd4Ef5Gh6Ij7Kl8Mn9Op0Qr"],"sequences":["Pq2Rs3Tu4Vw5Xy6Za7Bc8De9Fg0Hi1Jk","Lm2No3Pq4Rs5Tu6Vw7Xy8Za9Bc0De1Fg"]}
```

All digest values are bare 32-character URL-safe base64 strings (sha512t24u: first 24 bytes of SHA-512).

---

## What contributes to the identifier

### Topology document fields

| Field | Description |
|---|---|
| `nodes` | Sorted list of node ids (sha512t24u of sequence) |
| `edges` | Sorted list of edges with from/to node, orientations, and overlap value |
| `overlap_policy` | `"discard"`, `"length_only"`, or `"full_cigar"` |

```json
{
  "edges": [
    {"from": "SXaB3kqLm9vNpRt7YwZcAeQf12345678", "from_orient": "+", "overlap": null, "to": "Ab1Cd2Ef3Gh4Ij5Kl6Mn7Op8Qr9St0Uv", "to_orient": "+"},
    {"from": "SXaB3kqLm9vNpRt7YwZcAeQf12345678", "from_orient": "+", "overlap": null, "to": "Wx1Yz2Ab3Cd4Ef5Gh6Ij7Kl8Mn9Op0Qr", "to_orient": "-"}
  ],
  "nodes": [
    "Ab1Cd2Ef3Gh4Ij5Kl6Mn7Op8Qr9St0Uv",
    "SXaB3kqLm9vNpRt7YwZcAeQf12345678",
    "Wx1Yz2Ab3Cd4Ef5Gh6Ij7Kl8Mn9Op0Qr"
  ],
  "overlap_policy": "discard"
}
```

### Per-path document fields

| Field | Description |
|---|---|
| `is_reference` | Whether the path is a reference path |
| `steps` | Ordered list of `{id, orient}` step objects (traversal order preserved) |

```json
{
  "is_reference": false,
  "steps": [
    {"id": "SXaB3kqLm9vNpRt7YwZcAeQf12345678", "orient": "+"},
    {"id": "Ab1Cd2Ef3Gh4Ij5Kl6Mn7Op8Qr9St0Uv", "orient": "+"},
    {"id": "Wx1Yz2Ab3Cd4Ef5Gh6Ij7Kl8Mn9Op0Qr", "orient": "-"}
  ]
}
```

The path name is **not** included in the per-path document — it appears only in the graph document's `names` array, linked by position to the corresponding path digest.

### Complete graph document fields

| Field | Description |
|---|---|
| `graph_topology` | sha512t24u digest of the topology document |
| `names` | Path names, in the order determined by `paths` |
| `paths` | sha512t24u digests of per-path topology documents, sorted lexicographically |
| `sequences` | sha512t24u digests of each path's concatenated nucleotide sequence, parallel to `paths` |

The `sequences` array is [GA4GH refget](https://samtools.github.io/hts-specs/refget.html)-compatible: each entry is `sha512t24u` of the path's full nucleotide sequence, assembled by concatenating node sequences in traversal order (applying reverse complement for `"-"` oriented steps). A `null` entry indicates the path contains a node with no sequence (`"*"` wildcard).

> **What is NOT included:** source file names, original node names (replaced by sequence-derived ids), comment/header lines, fragment records (`F`), unordered groups (`U`), and edge ids from GFA v2.

---

## Overlap policy

The overlap policy controls how CIGAR overlap strings on edges affect the identifier. It is embedded in the topology document as `overlap_policy` and reflected in edge `overlap` values, so graphs hashed under different policies always produce different topology digests.

| Policy | Behaviour |
|---|---|
| `discard` (default) | Overlaps are ignored (`overlap` is `null`). Use for variation graphs and blunt-joined assemblies. |
| `length_only` | CIGAR is normalized to its total consumed length (`5M2I3M` → `10`). |
| `full_cigar` | CIGAR string is included verbatim. |

---

## Design notes

### Edges and paths are independent

Edges and paths are parsed from separate line types and stored as independent collections in the Abstract Graph Model — there is no cross-validation between them. Edges come from `L`/`J` lines (GFA v1), `E` lines (GFA v2), or the graph's edge store (PackedGraph). Paths come from `P`/`W` lines (GFA v1) or `O` lines (GFA v2).

### Implicit topology

Consecutive steps in a path imply a traversal edge. For example, a path containing `seg1+ → seg2+` implies an edge from `seg1` (forward) to `seg2` (forward). This edge may or may not be declared as an explicit `L` line in the source file.

### Identifier impact

Because the topology digest is built only from explicit edge records, two graphs that are biologically equivalent — same paths, same topology — will produce **different identifiers** if one has explicit `L`/`E` lines and the other does not.

### Open question

It is an open design question whether parsers should infer missing edges from path steps (i.e. add an edge for every consecutive step pair not already present in the edge set). Doing so would make the identifier stable regardless of whether `L` lines were included in the source file, at the cost of making the topology digest dependent on path content.

---

## Format support

### GFA v1 / v1.2

| Record | Processed as | Notes |
|--------|--------------|-------|
| `S`    | Node  | Sequence normalized to uppercase; id = `sha512t24u(sequence)` |
| `L`    | Edge  | CIGAR overlap resolved per overlap policy |
| `J`    | Edge  | Jump distance (plain integer) resolved per overlap policy |
| `P`    | Path  | Comma-separated `seg+`/`seg-` segments |
| `W`    | Path  | Walk string `>seg<seg…`; name built as `sample_id#hap_index#seq_id` (PanSN) |

All other lines (`H`, `C`, `#`, …) are skipped.

### GFA v2

| Record | Processed as | Notes |
|--------|--------------|-------|
| `S`    | Node  | Sequence normalized to uppercase; id = `sha512t24u(sequence)` |
| `E`    | Edge  | Edge ids are discarded |
| `O`    | Path  | Ordered group |

`F` (fragment) and `U` (unordered group) records are skipped as they have no equivalent in the topology model.

### vg PackedGraph (`.pg`)

Requires `pip install pangenome-id[pg]` (`libbdsg` package).

| Element | Processed as | Notes |
|---------|--------------|-------|
| Handle  | Node  | Forward-strand sequence; id = `sha512t24u(sequence)` |
| Edge    | Edge  | No overlap information; `overlap` is always `None` |
| Path    | Path  | Step orientation preserved |

Format is auto-detected from the `.pg` file extension. GFA formats are auto-detected from the `VN:Z:` header tag or the `.gfa` extension.

---

## Project layout

```
pangenome_graph_identifier/
├── pangenome_id/
│   ├── model.py          # Abstract Graph Model dataclasses
│   ├── canonicalize.py   # AGM → RFC-8785 canonical JSON bytes
│   ├── hasher.py         # canonical JSON → identifier dict
│   ├── parsers/
│   │   ├── base.py           # shared normalization and edge canonicalization
│   │   ├── gfa1.py           # GFA v1 / v1.2 parser
│   │   ├── gfa2.py           # GFA v2 parser
│   │   └── packed_graph.py   # vg PackedGraph (.pg) parser (requires libbdsg)
│   └── cli.py            # CLI entry point
└── tests/
    ├── test_model.py
    ├── test_canonicalize.py
    ├── test_gfa1.py
    ├── test_gfa2.py
    ├── test_cross_format.py   # cross-format identity tests
    └── test_packed_graph.py   # PackedGraph parser (mock-based + optional real-file)
```

---

## Running the tests

```bash
pip install pytest
python -m pytest tests/ -v
```
