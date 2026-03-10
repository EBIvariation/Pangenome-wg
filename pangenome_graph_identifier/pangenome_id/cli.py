"""CLI entry point for pangenome-id."""

import argparse
import sys


def _detect_format(filepath: str) -> str:
    """Auto-detect format from header line or file extension."""
    if filepath.lower().endswith(".pg"):
        return "pg"
    with open(filepath, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            if "VN:Z:2" in line:
                return "gfa2"
            if "VN:Z:1" in line:
                return "gfa1"
            break
    if filepath.lower().endswith(".gfa"):
        return "gfa1"
    raise ValueError(
        "Cannot auto-detect format. Please specify --format {gfa1,gfa2}."
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="pangenome-id",
        description="Produce a stable, deterministic identifier for a pangenome graph.",
    )
    parser.add_argument("file", help="Path to the pangenome graph file")
    parser.add_argument(
        "--format",
        choices=["gfa1", "gfa2", "pg", "auto"],
        default="auto",
        help="File format (default: auto)",
    )
    parser.add_argument(
        "--overlap-policy",
        choices=["discard", "length_only", "full_cigar"],
        default="discard",
        help="How to treat overlaps (default: discard)",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print stats to stderr",
    )
    args = parser.parse_args()

    fmt = args.format
    if fmt == "auto":
        fmt = _detect_format(args.file)

    from pangenome_id.parsers.gfa1 import GFA1Parser
    from pangenome_id.parsers.gfa2 import GFA2Parser
    from pangenome_id.canonicalize import serialize

    if fmt == "gfa1":
        gfa_parser = GFA1Parser(overlap_policy=args.overlap_policy)
        format_label = "GFA v1"
    elif fmt == "pg":
        from pangenome_id.parsers.packed_graph import PackedGraphParser
        gfa_parser = PackedGraphParser(overlap_policy=args.overlap_policy)
        format_label = "vg PackedGraph"
    else:
        gfa_parser = GFA2Parser(overlap_policy=args.overlap_policy)
        format_label = "GFA v2"

    graph = gfa_parser.parse(args.file)
    canonical = serialize(graph)
    sys.stdout.buffer.write(canonical + b"\n")

    if args.verbose:
        print(f"Format:         {format_label}", file=sys.stderr)
        print(f"Nodes:          {len(graph.nodes)}", file=sys.stderr)
        print(f"Edges:          {len(graph.edges)}", file=sys.stderr)
        print(f"Paths:          {len(graph.paths)}", file=sys.stderr)
        print(f"Overlap policy: {args.overlap_policy}", file=sys.stderr)


if __name__ == "__main__":
    main()
