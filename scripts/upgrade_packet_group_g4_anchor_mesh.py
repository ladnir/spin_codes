#!/usr/bin/env python3
"""Upgrade a diagnostic g=4 Qhull ledger to the flat exact verifier schema.

The converter does not turn binary64 discovery into a proof.  It reconstructs
the mesh over exact ``Fraction`` anchors, audits the complete covered+failed
topology, resolves and hashes every atlas source, and serializes fixed mixture
weights as exact normalized rationals.  The outward verifier remains the sole
authority for every final inequality.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np

import certify_packet_group_g4_anchor_mesh as verifier
from packet_group_outer_profile import normalization_log2
from probe_packet_group_support_face_cover import parse_atlas, solve_fixed_mixture


SCHEMA = "packet-group-g4-flat-exact-anchor-mesh-v1"


def compact_sha256(value: Any) -> str:
    payload = json.dumps(value, separators=(",", ":"), sort_keys=True)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def split_exact_degenerate_cells(
    rows: list[dict[str, Any]],
    anchors: tuple[tuple[Fraction, ...], ...],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    retained = []
    dropped = []
    for row in rows:
        indices = verifier.simplex_indices(row, len(anchors))
        determinant = verifier.simplex_determinant(indices, anchors)
        if determinant:
            retained.append(row)
            continue
        dropped.append(
            {
                "cell_id": int(row["cell_id"]),
                "anchor_indices": list(indices),
                "canonical_anchor_indices": sorted(indices),
                "exact_determinant": "0",
            }
        )
    return retained, sorted(dropped, key=lambda row: row["cell_id"])


def resolve_sources(
    ledger_path: Path,
    ledger: dict[str, Any],
    overrides: list[Path],
) -> list[Path]:
    raw_sources = ledger.get("sources", [])
    logical = [
        Path(str(row.get("path", ""))) if isinstance(row, dict) else Path(str(row))
        for row in raw_sources
    ]
    override_by_name = {path.name: path for path in overrides}
    result = []
    for raw in logical:
        normalized = Path(str(raw).replace("\\", "/"))
        candidates = (
            override_by_name.get(normalized.name),
            normalized,
            ledger_path.parent / normalized,
            Path.cwd() / normalized,
        )
        resolved = next((path.resolve() for path in candidates if path is not None and path.exists()), None)
        if resolved is None:
            raise ValueError(
                f"cannot resolve atlas source {raw}; supply it with --atlas"
            )
        result.append(resolved)
    if overrides:
        expected = {path.name for path in logical}
        unexpected = sorted(set(override_by_name) - expected)
        if unexpected:
            raise ValueError("--atlas contains sources absent from ledger: " + ", ".join(unexpected))
    if len({path.name for path in result}) != len(result):
        raise ValueError("atlas source basenames must be unique for canonical witness references")
    return result


def source_records(paths: list[Path]) -> list[dict[str, str]]:
    return [
        {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
        for path in paths
    ]


def exact_weight(value: float, max_denominator: int) -> Fraction:
    if not math.isfinite(value) or value < 0:
        raise ValueError(f"invalid discovery mixture weight {value!r}")
    return Fraction.from_float(value).limit_denominator(max_denominator)


def rationalize_weights(
    names: list[str], values: list[float], max_denominator: int
) -> tuple[tuple[str, Fraction], ...]:
    if len(names) != len(values) or not names:
        raise ValueError("fixed mixture must contain matching nonempty names and weights")
    combined: dict[str, Fraction] = {}
    for name, value in zip(names, values):
        weight = exact_weight(float(value), max_denominator)
        combined[name] = combined.get(name, Fraction(0)) + weight
    total = sum(combined.values(), Fraction(0))
    if total <= 0:
        raise ValueError("fixed mixture has zero total weight")
    return tuple(sorted((name, weight / total) for name, weight in combined.items() if weight))


def diagnostic_values(
    names: list[str],
    witnesses: dict[str, dict[str, Any]],
    anchors: tuple[tuple[Fraction, ...], ...],
    indices: tuple[int, ...],
) -> np.ndarray:
    matrix = np.asarray(
        [[float(value) for value in anchors[index]] for index in indices],
        dtype=np.float64,
    )
    normalizations = normalization_log2(verifier.GROUP_BITS, matrix)
    rows = []
    for name in names:
        witness = witnesses.get(name)
        if witness is None:
            raise ValueError(f"active mixture references unknown atlas witness {name!r}")
        rows.append(
            float(witness["constant_log2"])
            - np.asarray(witness["charge"], dtype=np.float64) @ matrix.T
            - normalizations
        )
    return np.asarray(rows)


def convert_mixture(
    row: dict[str, Any],
    anchors: tuple[tuple[Fraction, ...], ...],
    witnesses: dict[str, dict[str, Any]] | None,
    max_denominator: int,
) -> tuple[list[dict[str, Any]], float | None]:
    raw = row.get("mixture")
    if not isinstance(raw, list) or not raw:
        raise ValueError(f"cell {row.get('cell_id')} lacks a discovery mixture")
    names = [str(item.get("witness", "")) for item in raw]
    weights = [float(item["weight"]) for item in raw]
    diagnostic_maximum = None
    if witnesses is not None:
        indices = verifier.simplex_indices(row, len(anchors))
        values = diagnostic_values(names, witnesses, anchors, indices)
        result = solve_fixed_mixture(values)
        if not result.success:
            raise ValueError(f"active-set mixture LP failed in cell {row.get('cell_id')}: {result.message}")
        weights = [float(value) for value in result.x[: len(names)]]
    mixture = rationalize_weights(names, weights, max_denominator)
    if witnesses is not None:
        values = diagnostic_values([name for name, _weight in mixture], witnesses, anchors, verifier.simplex_indices(row, len(anchors)))
        vector = np.asarray([float(weight) for _name, weight in mixture])
        diagnostic_maximum = float(np.max(vector @ values))
    rendered = [
        {
            "witness": name,
            "weight_exact": (
                str(weight.numerator)
                if weight.denominator == 1
                else f"{weight.numerator}/{weight.denominator}"
            ),
        }
        for name, weight in mixture
    ]
    return rendered, diagnostic_maximum


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--atlas", type=Path, action="append", default=[])
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-denominator", type=int, default=1 << 40)
    parser.add_argument(
        "--no-reoptimize-active-mixtures",
        action="store_true",
        help="skip the default small LP rerun over each cell's active witnesses",
    )
    parser.add_argument(
        "--drop-exact-degenerate",
        action="store_true",
        help=(
            "omit only exact zero-determinant QJ cells, then require the "
            "remaining complex to pass the full exact mesh audit"
        ),
    )
    args = parser.parse_args()
    if args.max_denominator <= 0:
        raise SystemExit("g4 mesh upgrade: --max-denominator must be positive")
    source = json.loads(args.input.read_text(encoding="utf-8"))
    if int(source.get("group_bits", -1)) != verifier.GROUP_BITS:
        raise SystemExit("g4 mesh upgrade: input must have group_bits=4")

    anchors = verifier.load_anchors(source)
    covered = source.get("covered_cell_ledger", [])
    failed = source.get("failed_cell_ledger", [])
    if not isinstance(covered, list) or not isinstance(failed, list):
        raise ValueError("covered and failed cell ledgers must be arrays")
    all_cells = list(covered) + list(failed)
    expected = int(source.get("simplices", -1))
    if len(all_cells) != expected:
        raise ValueError(f"input omits cells: covered+failed={len(all_cells)} != {expected}")
    identifiers = [int(row.get("cell_id", -1)) for row in all_cells]
    if set(identifiers) != set(range(expected)) or len(set(identifiers)) != expected:
        raise ValueError("input cell IDs must be exactly 0..simplices-1")
    audit_rows = sorted(all_cells, key=lambda row: int(row["cell_id"]))
    dropped_degenerate = []
    if args.drop_exact_degenerate:
        audit_rows, dropped_degenerate = split_exact_degenerate_cells(audit_rows, anchors)
    dropped_digest = compact_sha256(dropped_degenerate)
    try:
        geometry = verifier.verify_mesh(
            {"complete_cover": True, "simplices": len(audit_rows)}, anchors, audit_rows
        )
    except ValueError as error:
        if dropped_degenerate:
            raise ValueError(
                "remaining exact complex failed after dropping "
                f"{len(dropped_degenerate)} zero-determinant cells "
                f"(dropped_sha256={dropped_digest}): {error}"
            ) from error
        raise

    paths = resolve_sources(args.input, source, args.atlas)
    witnesses = None
    reoptimize = not args.no_reoptimize_active_mixtures
    if reoptimize:
        witnesses = {row["name"]: row for row in parse_atlas(paths, verifier.GROUP_BITS)}

    target = float(source.get("target_log2", math.inf))
    converted_covered = []
    converted_failed = []
    covered_ids = {int(item["cell_id"]) for item in covered}
    for owner_rank, original in enumerate(audit_rows):
        source_cell_id = int(original["cell_id"])
        row = {
            "cell_id": owner_rank,
            "source_cell_id": source_cell_id,
            "owner_rank": owner_rank,
            "anchor_indices": [int(value) for value in original["anchor_indices"]],
        }
        mixture, diagnostic_maximum = convert_mixture(
            original, anchors, witnesses, args.max_denominator
        )
        row["mixture"] = mixture
        if diagnostic_maximum is not None:
            row["rational_mixture_diagnostic_maximum_vertex_log2"] = diagnostic_maximum
        originally_covered = source_cell_id in covered_ids
        diagnostic_closes = diagnostic_maximum is None or diagnostic_maximum <= target
        if originally_covered and diagnostic_closes:
            converted_covered.append(row)
        else:
            row["diagnostic_failure_reason"] = (
                "producer_failed_cell" if not originally_covered else "rational_active_set_no_longer_closes"
            )
            converted_failed.append(row)

    complete = not converted_failed and len(converted_covered) == len(audit_rows)
    report = {
        "schema": SCHEMA,
        "status": (
            "EXACT_SCHEMA_READY_FOR_OUTWARD_G4_ANCHOR_VERIFICATION"
            if complete
            else "EXACT_GEOMETRY_G4_ANCHOR_SCHEMA_INCOMPLETE"
        ),
        "group_bits": verifier.GROUP_BITS,
        "complete_cover": complete,
        "profile_owner_rule": verifier.FLAT_OWNER_RULE,
        "simplices": len(audit_rows),
        "source_simplices": expected,
        "covered_cells": len(converted_covered),
        "failed_cells": len(converted_failed),
        "anchor_profiles": [
            [str(value) for value in profile] for profile in anchors
        ],
        "covered_cell_ledger": converted_covered,
        "failed_cell_ledger": converted_failed,
        "dropped_exact_degenerate_cells": dropped_degenerate,
        "dropped_exact_degenerate_cells_sha256": dropped_digest,
        "sources": source_records(paths),
        "exact_geometry": {"mode": "flat_exact_simplicial_mesh", **geometry},
        "conversion": {
            "source_ledger": str(args.input),
            "source_ledger_sha256": hashlib.sha256(args.input.read_bytes()).hexdigest(),
            "max_weight_denominator": args.max_denominator,
            "active_mixtures_reoptimized": reoptimize,
            "drop_exact_degenerate": args.drop_exact_degenerate,
            "dropped_exact_degenerate_cell_count": len(dropped_degenerate),
            "note": "exact weights and geometry are inputs; only the outward verifier certifies inequalities",
        },
    }
    rendered = json.dumps(report, indent=2, sort_keys=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered + "\n", encoding="utf-8")
    print(
        f"exact_mesh_audit=PASS simplices={len(audit_rows)} "
        f"dropped_exact_degenerate={len(dropped_degenerate)}"
    )
    print(f"covered_cells={len(converted_covered)} failed_cells={len(converted_failed)}")
    print(f"complete_cover={complete}")


if __name__ == "__main__":
    main()
