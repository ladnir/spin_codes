#!/usr/bin/env python3
"""Freeze the g=8 support-local atlas and exact support census."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from check_packet_group_g8_census import (
    EXPECTED_FEASIBLE_SUPPORTS,
    EXPECTED_PROFILE_COUNT,
    support_profile_count,
)
from packet_group_outer_profile import N, atom_count


ROOT = Path(__file__).resolve().parents[1]
GROUP_BITS = 8
CLASSES = GROUP_BITS + 1
EXPECTED_ATLAS_SHA256 = "c48644fa62ad74619a4db30b07950b14aaf1419d040760f7d659f43288fca45e"


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_bytes(value) -> bytes:
    """Encode the restricted manifest JSON profile deterministically."""

    return json.dumps(
        value,
        ensure_ascii=False,
        allow_nan=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")


def object_sha256(value) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def source_record(source_id: str, path: Path, schema: str | None = None) -> dict:
    record = {
        "source_id": source_id,
        "path": str(path.relative_to(ROOT)).replace("\\", "/"),
        "sha256": file_sha256(path),
    }
    if schema is not None:
        record["schema"] = schema
    return record


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--atlas", type=Path, default=ROOT / "out" / "g8_support_seed_atlas.json"
    )
    parser.add_argument(
        "--output", type=Path, default=ROOT / "G8_SUPPORT_MANIFEST.json"
    )
    args = parser.parse_args()
    atlas_sha256 = file_sha256(args.atlas)
    if atlas_sha256 != EXPECTED_ATLAS_SHA256:
        raise SystemExit(
            f"atlas digest changed: expected {EXPECTED_ATLAS_SHA256}, got {atlas_sha256}"
        )
    atlas = json.loads(args.atlas.read_text(encoding="utf-8"))
    if atlas.get("schema") != "permute-conv.packet-group-g8-support-seed-atlas.v1":
        raise SystemExit("unexpected support atlas schema")
    rows = atlas.get("rows", [])
    if len(rows) != EXPECTED_FEASIBLE_SUPPORTS:
        raise SystemExit("support atlas row count changed")

    atlas_masks = [int(row["support_mask"]) for row in rows]
    if len(set(atlas_masks)) != EXPECTED_FEASIBLE_SUPPORTS:
        raise SystemExit("support atlas masks are not unique")
    feasible_masks = []
    expected_masks = []
    for mask in range(1, 1 << CLASSES):
        count, excluded = support_profile_count(mask)
        if count:
            expected_masks.append(mask)
            feasible_masks.append(
                {
                    "mask": f"0x{mask:03x}",
                    "active_classes": [
                        index for index in range(CLASSES) if mask & (1 << index)
                    ],
                    "dimension": mask.bit_count() - 1,
                    "count": str(count),
                    "weight_cut_exclusions": str(excluded),
                }
            )
    if set(atlas_masks) != set(expected_masks):
        raise SystemExit("support atlas masks do not match the exact census")

    row_records = []
    for index, row in enumerate(rows):
        row_records.append(
            {
                "row": index,
                "support_mask": f"0x{int(row['support_mask']):03x}",
                "row_sha256": object_sha256(row),
            }
        )

    fixed_sources = [
        source_record(
            "spectrum-01", ROOT / "out" / "ebch85_band01_split_spectrum.csv"
        ),
        source_record(
            "spectrum-01-punctured",
            ROOT / "out" / "ebch84_punctured_band01_split_spectrum.csv",
        ),
        source_record(
            "spectrum-12", ROOT / "out" / "ebch86_band12_split_spectrum.csv"
        ),
        source_record(
            "conditioned-row-discovery",
            ROOT / "scripts" / "probe_packet_group_conditioned_row_outer.py",
        ),
        source_record(
            "combined-outward-evaluator",
            ROOT / "scripts" / "certify_packet_group_triangle_ledger.py",
        ),
        source_record(
            "inner-bound-evaluator",
            ROOT / "scripts" / "packet_group_profile_bound.py",
        ),
    ]
    manifest = {
        "schema": "packet-group-g8-support-manifest-v1",
        "run_id": "g8-support-atlas-3934dae",
        "parameters": {
            "group_bits": GROUP_BITS,
            "N": str(N),
            "K": str(N // 2),
            "M": str(atom_count(GROUP_BITS)),
            "minimum_profile_weight": "21",
            "probability_target_log2": "-40/1",
        },
        "arithmetic": {
            "canonical_json": "sorted-compact-json-v1",
            "integer_split_rule": "primitive-affine-gap-v1",
            "outward_evaluator_source_id": "combined-outward-evaluator",
        },
        "census": {
            "rule": "exact-support-positive-compositions-weight-cut-v1",
            "feasible_masks": feasible_masks,
            "infeasible_masks": ["0x000", "0x001"],
            "feasible_mask_count": EXPECTED_FEASIBLE_SUPPORTS,
            "feasible_profile_count": str(EXPECTED_PROFILE_COUNT),
        },
        "witness_sources": [
            {
                "source_id": "atlas-3934dae",
                "path": str(args.atlas.relative_to(ROOT)).replace("\\", "/"),
                "sha256": atlas_sha256,
                "schema": atlas["schema"],
                "rows": row_records,
            }
        ],
        "fixed_sources": fixed_sources,
    }
    args.output.write_bytes(
        (
            json.dumps(manifest, indent=2, ensure_ascii=False, allow_nan=False)
            + "\n"
        ).encode("utf-8")
    )
    print(f"manifest_sha256={file_sha256(args.output)}")
    print(f"atlas_sha256={atlas_sha256}")
    print(f"feasible_masks={len(feasible_masks)}")
    print(f"witness_rows={len(row_records)}")


if __name__ == "__main__":
    main()
