#!/usr/bin/env python3
"""Select a small, diverse batch from a profile-cover residual ledger."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def residual_score(row: dict) -> float:
    if "diagnostic_gap_bits" in row:
        return float(row["diagnostic_gap_bits"])
    return float(row.get("value_log2", 0.0) - row.get("target_log2", 0.0))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", type=Path, required=True)
    parser.add_argument("--count", type=int, default=128)
    parser.add_argument("--top-worst", type=int, default=32)
    parser.add_argument("--top-volume", type=int, default=32)
    parser.add_argument(
        "--source-cell-id",
        action="append",
        default=[],
        help="restrict selection to residuals owned by one of these source cells",
    )
    parser.add_argument(
        "--volume-bias",
        type=float,
        default=0.25,
        help="exponent applied to normalized cell volume during farthest-first selection",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if (
        args.count <= 0
        or not 0 <= args.top_worst <= args.count
        or not 0 <= args.top_volume <= args.count
        or args.top_worst + args.top_volume > args.count
        or args.volume_bias < 0.0
    ):
        raise SystemExit("residual selector: invalid batch sizes")

    payload = json.loads(args.ledger.read_text(encoding="utf-8"))
    candidates = payload.get("uncovered_integer_residuals", [])
    if args.source_cell_id:
        source_cells = set(args.source_cell_id)
        candidates = [
            row for row in candidates if row.get("source_cell_id") in source_cells
        ]
    unique: dict[tuple[int, ...], dict] = {}
    for row in candidates:
        profile = tuple(int(value) for value in row["profile"])
        score = residual_score(row)
        volume = max(0.0, float(row.get("volume_proxy", 0.0)))
        current = unique.get(profile)
        if current is None:
            unique[profile] = {**row, "volume_proxy": volume}
            continue
        combined_volume = float(current.get("volume_proxy", 0.0)) + volume
        current_score = residual_score(current)
        if score > current_score:
            unique[profile] = {**row, "volume_proxy": combined_volume}
        else:
            current["volume_proxy"] = combined_volume
    rows = list(unique.values())
    rows.sort(
        key=residual_score,
        reverse=True,
    )
    if len(rows) <= args.count:
        selected = rows
    else:
        points = np.asarray([row["profile"] for row in rows], dtype=np.float64)
        points /= np.maximum(1.0, np.sum(points, axis=1, keepdims=True))
        selected_indices = list(range(min(args.top_worst, len(rows))))
        selected_set = set(selected_indices)
        volume_order = np.argsort(
            -np.asarray([float(row.get("volume_proxy", 0.0)) for row in rows])
        )
        for raw_index in volume_order:
            if len(selected_indices) >= args.top_worst + args.top_volume:
                break
            index = int(raw_index)
            if index not in selected_set:
                selected_indices.append(index)
                selected_set.add(index)
        if selected_indices:
            distances = np.min(
                np.sum(
                    (points[:, None, :] - points[selected_indices][None, :, :]) ** 2,
                    axis=2,
                ),
                axis=1,
            )
        else:
            selected_indices = [0]
            selected_set = {0}
            distances = np.sum((points - points[0]) ** 2, axis=1)
        distances[list(selected_set)] = -1.0
        volumes = np.asarray(
            [float(row.get("volume_proxy", 0.0)) for row in rows], dtype=np.float64
        )
        maximum_volume = float(np.max(volumes))
        if maximum_volume > 0.0:
            volume_weights = np.power(volumes / maximum_volume, args.volume_bias)
        else:
            volume_weights = np.ones(len(rows), dtype=np.float64)
        while len(selected_indices) < args.count:
            index = int(np.argmax(distances * volume_weights))
            selected_indices.append(index)
            selected_set.add(index)
            new_distances = np.sum((points - points[index]) ** 2, axis=1)
            distances = np.minimum(distances, new_distances)
            distances[list(selected_set)] = -1.0
        selected = [rows[index] for index in selected_indices]

    report = {
        "status": "DIAGNOSTIC_DIVERSE_PROFILE_RESIDUAL_BATCH",
        "source": str(args.ledger),
        "input_rows": len(candidates),
        "unique_profiles": len(rows),
        "selected_profiles": len(selected),
        "top_worst": min(args.top_worst, len(selected)),
        "top_volume": min(args.top_volume, max(0, len(selected) - args.top_worst)),
        "volume_bias": args.volume_bias,
        "source_cell_ids": args.source_cell_id,
        "uncovered_integer_residuals": [
            {
                **row,
                "name": f"selected_{index:04d}_" + str(row.get("name", "residual")),
            }
            for index, row in enumerate(selected)
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"input_rows={len(candidates)}")
    print(f"unique_profiles={len(rows)}")
    print(f"selected_profiles={len(selected)}")


if __name__ == "__main__":
    main()
