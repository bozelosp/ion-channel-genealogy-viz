#!/usr/bin/env python3

"""Recompute high-level network statistics from network_data.json.

The script mirrors the logic used by the visualizer: it excludes the catch-all
"Other" ion class, counts unique .mod files (unique_modelDB_mod_id), network
nodes, similarity edges restricted to the filtered node set, and reports the
redundancy/average reuse factor figures used in documentation.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, Iterable, Set


OTHER_ION_CLASS = "other"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "network_json",
        nargs="?",
        default=Path("public/data/network_data.json"),
        type=Path,
        help="Path to network_data.json (defaults to public/data/network_data.json)",
    )
    return parser.parse_args()


def load_network(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def is_other_class(node: Dict[str, Any]) -> bool:
    ion_class = (
        node.get("original_model", {}).get("ion_class")
        or node.get("ion_class")
        or ""
    )
    return isinstance(ion_class, str) and ion_class.strip().lower() == OTHER_ION_CLASS


def collect_unique_ids(node: Dict[str, Any]) -> Set[str]:
    ids: Set[str] = set()
    original = node.get("original_model", {})
    primary = original.get("unique_modelDB_mod_id")
    if isinstance(primary, str) and primary:
        ids.add(primary)

    for entry in node.get("identical_models") or []:
        if not isinstance(entry, dict):
            continue
        candidate = entry.get("unique_modelDB_mod_id")
        if isinstance(candidate, str) and candidate:
            ids.add(candidate)

    return ids


def coerce_node_id(value: Any) -> str | None:
    if isinstance(value, (str, int)):
        return str(value)
    if isinstance(value, dict) and "id" in value:
        inner = value.get("id")
        if isinstance(inner, (str, int)):
            return str(inner)
    return None


def main() -> None:
    args = parse_args()
    data = load_network(args.network_json)

    nodes: Iterable[Dict[str, Any]] = data.get("nodes", [])
    links: Iterable[Dict[str, Any]] = data.get("links", [])

    filtered_nodes = []
    unique_files: Set[str] = set()

    for node in nodes:
        ids = collect_unique_ids(node)
        unique_files.update(ids)
        if not is_other_class(node):
            filtered_nodes.append({"node": node, "ids": ids})

    filtered_node_ids = {
        str(entry["node"].get("id"))
        for entry in filtered_nodes
        if entry["node"].get("id") is not None
    }

    filtered_unique_files: Set[str] = set()
    for entry in filtered_nodes:
        filtered_unique_files.update(entry["ids"])

    filtered_links = 0
    for link in links:
        source_id = coerce_node_id(link.get("source"))
        target_id = coerce_node_id(link.get("target"))
        if source_id in filtered_node_ids and target_id in filtered_node_ids:
            filtered_links += 1

    total_files = len(filtered_unique_files)
    unique_variants = len(filtered_nodes)
    redundant_files = total_files - unique_variants
    redundancy_pct = (redundant_files / total_files * 100) if total_files else 0.0
    reuse_factor = (total_files / unique_variants) if unique_variants else 0.0

    print(f"Network data source: {args.network_json}")
    print("-- after excluding \"Other\" ion class --")
    print(f"Unique .mod files       : {total_files:,}")
    print(f"Unique network nodes    : {unique_variants:,}")
    print(f"Similarity edges        : {filtered_links:,}")
    print(f"Redundant instances     : {redundant_files:,} ({redundancy_pct:.2f}%)")
    print(f"Average reuse factor    : {reuse_factor:.2f}× per unique implementation")


if __name__ == "__main__":
    main()

