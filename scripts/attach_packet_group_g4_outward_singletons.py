#!/usr/bin/env python3
"""Attach cached outward-hardened singleton witnesses to a g=4 mesh."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from decimal import Decimal
from pathlib import Path

import numpy as np

import certify_packet_group_g4_anchor_mesh as g4
import certify_packet_group_triangle_ledger as g2
from outward_log2 import Interval, log2_int
from probe_packet_group_g4_stellar_cover import optimize_values


def load_cache(directory: Path) -> dict[str, dict]:
    result = {}
    for path in directory.glob("*.json"):
        row = json.loads(path.read_text(encoding="utf-8"))
        if row.get("schema") != "g2-outward-hardened-witness-cache-v1":
            continue
        hardened = g2._deserialize_hardened(row["hardened"])
        name = hardened["name"]
        if name in result:
            raise ValueError(f"duplicate hardened witness {name!r}")
        result[name] = hardened
    return result


def logsumexp2(values: list[float]) -> float:
    maximum = max(values)
    return maximum + math.log2(math.fsum(2.0 ** (value - maximum) for value in values))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--atlas-source", type=Path, action="append", default=[])
    args = parser.parse_args()

    ledger = json.loads(args.ledger.read_text(encoding="utf-8"))
    source_by_name = {
        Path(str(row["path"]).replace("\\", "/")).name: row
        for row in ledger["sources"]
    }
    for path in args.atlas_source:
        source_by_name[path.name] = {
            "path": str(path.resolve()),
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        }
    ledger["sources"] = [source_by_name[name] for name in sorted(source_by_name)]
    hardened = load_cache(args.cache_dir)
    allowed = {Path(str(row["path"]).replace("\\", "/")).name for row in ledger["sources"]}
    hardened = {
        name: value for name, value in hardened.items()
        if name.split(":", 1)[0] in allowed
    }
    hardened["full_bijection"] = g4.harden_full_bijection()
    names = sorted(hardened)
    constants = np.asarray([float(hardened[name]["constant"].hi) for name in names])
    charge_lowers = np.asarray([
        [float("nan") if charge is None else float(charge.lo) for charge in hardened[name]["charges"]]
        for name in names
    ])

    all_terms = []
    changed = 0
    for stratum in ledger["support_strata"]:
        if not stratum.get("feasible"):
            continue
        support = tuple(map(int, stratum["support"]))
        eligible = np.asarray([
            all(hardened[name]["charges"][coordinate] is not None for coordinate in support)
            for name in names
        ])
        candidates = np.flatnonzero(eligible)
        if not len(candidates):
            raise ValueError(f"support {support} has no outward-eligible witnesses")
        anchors = g4.parse_anchor_profiles(
            stratum["anchor_profiles"], f"support mask {stratum['support_mask']}"
        )
        anchor_array = np.asarray([[float(value) for value in row] for row in anchors])
        normalization_lowers = np.asarray([
            float(g4.outward_fractional_normalization(profile).lo)
            for profile in anchors
        ])
        values = np.full((len(anchors), len(names)), np.inf)
        support_indices = np.asarray(support, dtype=np.int64)
        values[:, candidates] = (
            constants[candidates, None]
            - charge_lowers[np.ix_(candidates, support_indices)]
            @ anchor_array[:, support_indices].T
            - normalization_lowers[None, :]
        ).T
        rows = stratum["covered_cell_ledger"]
        for row in rows:
            indices = g4.stratum_cell_indices(row, len(anchors), len(support) - 1)
            optimized = optimize_values(
                values[np.asarray(indices)],
                names,
                float(row["target_log2"]),
                0.0,
                32,
                1e-12,
                candidate_mode="column-generation",
            )
            mixture = tuple(
                (names[index], weight)
                for index, weight in zip(optimized["indices"], optimized["weights"])
            )
            cell_maximum = None
            for vertex in indices:
                value, _components = g4.evaluate_profile(anchors[vertex], mixture, hardened)
                cell_maximum = value if cell_maximum is None else Interval(
                    max(cell_maximum.lo, value.lo), max(cell_maximum.hi, value.hi)
                )
            assert cell_maximum is not None
            count, omitted, ranges = g4.simplex_integer_profile_count_upper(indices, anchors)
            term = cell_maximum + log2_int(count)
            row["mixture"] = [
                {"witness": name, "weight_exact": str(weight)}
                for name, weight in mixture
            ]
            row["maximum_vertex_log2"] = float(cell_maximum.hi)
            row["cell_local_contribution_log2_upperish"] = float(term.hi)
            row["uniform_target_passed"] = cell_maximum.hi <= Decimal(str(row["target_log2"]))
            row["outward_mixture_attachment"] = True
            row["outward_mixture_components"] = len(mixture)
            row["integer_profile_count_dropped_coordinate"] = omitted
            row["inclusive_integer_coordinate_range_widths"] = [
                high - low + 1 for low, high in ranges
            ]
            all_terms.append(float(term.hi))
            changed += 1
    ledger["outward_mixture_attachment"] = {
        "cached_witnesses_considered": len(names),
        "cells_updated": changed,
        "selection": "column-generation minimax over cached outward vertex upper endpoints",
    }
    ledger["cell_local_security_diagnostic_log2"] = logsumexp2(all_terms)
    ledger["status"] = "OUTWARD_MIXTURES_ATTACHED_PENDING_INDEPENDENT_REPLAY"
    args.output.write_text(json.dumps(ledger, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "cells_updated": changed,
        "cached_witnesses_considered": len(names),
        "largest_term": max(all_terms),
        "union_upperish": logsumexp2(all_terms),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
