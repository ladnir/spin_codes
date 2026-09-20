#!/usr/bin/env python3
"""Run the corrected 24-node plus value-15 endpoint certificates sequentially."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXPLORATIONS = ROOT / "explorations"
MAXIMAL_SUPPORTS = (
    (0x50, (10, 20), 30),
    (0x49, (1, 9, 20), 30),
    (0x31, (1, 10, 18), 29),
    (0x32, (2, 10, 18), 30),
    (0x27, (1, 2, 4, 18), 25),
    (0x47, (1, 2, 4, 20), 27),
    (0x2B, (1, 2, 9, 18), 30),
    (0x1F, (1, 2, 4, 9, 10), 26),
)
COMPONENT_DEGREES = (1, 2, 4, 9, 10, 18, 20)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def executable(mode: str, window_nodes: int, dimension: int) -> Path:
    stem = "certify" if mode == "primary" else "verify"
    return ROOT / "scripts" / (
        f"{stem}_riffle_dp_2lap_g4_component_w{window_nodes}_d{dimension}.exe"
    )


def source(mode: str) -> Path:
    return ROOT / "scripts" / "certify_riffle_dp_2lap_g4_component_24node.cpp"


def receipt_path(
    mode: str, window_nodes: int, support_mask: int, packet_value: int
) -> Path:
    return EXPLORATIONS / (
        f"riffle_dp_2lap_g4_component_w{window_nodes}_"
        f"{mode}_s{support_mask:02x}_v{packet_value:02d}.json"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("primary", "independent"), required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    rows = []
    total_candidates = 0
    total_elapsed = 0.0
    overall_result = "PASS"
    executable_hashes = {}
    schedule = ((24, tuple(range(1, 16))), (12, (15,)))
    for support_mask, degrees, dimension in MAXIMAL_SUPPORTS:
      for window_nodes, packet_values in schedule:
        binary = executable(args.mode, window_nodes, dimension)
        executable_hashes[f"w{window_nodes}_d{dimension}"] = digest(binary)
        for packet_value in packet_values:
            path = receipt_path(args.mode, window_nodes, support_mask, packet_value)
            completed = subprocess.run(
                [str(binary), hex(support_mask), str(packet_value), str(path)],
                cwd=ROOT,
                capture_output=True,
                text=True,
                check=False,
            )
            if completed.returncode not in (0, 1):
                raise RuntimeError(
                    f"component runner: executable failed\n{completed.stdout}{completed.stderr}"
                )
            receipt = json.loads(path.read_text())
            if (
                int(receipt["component_support_mask_hex"], 16) != support_mask
                or tuple(receipt["component_degrees"]) != degrees
                or receipt["code_dimension"] != dimension
                or receipt["packet_value"] != packet_value
                or receipt["window_nodes"] != window_nodes
            ):
                raise RuntimeError("component runner: receipt identity mismatch")
            total_candidates += receipt["total_information_vectors_checked"]
            total_elapsed += receipt["elapsed_seconds"]
            rows.append(
                {
                    "component_support_mask_hex": hex(support_mask),
                    "component_degrees": list(degrees),
                    "code_dimension": dimension,
                    "window_nodes": window_nodes,
                    "packet_value": packet_value,
                    "receipt": str(path.relative_to(ROOT)).replace("\\", "/"),
                    "receipt_sha256": digest(path),
                    "information_set_count": receipt["information_set_count"],
                    "information_radius": receipt["information_radius"],
                    "information_vectors_checked": receipt[
                        "total_information_vectors_checked"
                    ],
                    "minimum_enumerated_side": min(
                        receipt["minimum_enumerated_weight"],
                        receipt["minimum_enumerated_complement_weight"],
                    ),
                    "result": receipt["result"],
                }
            )
            print(
                f"mode={args.mode} window={window_nodes} support={hex(support_mask)} value={packet_value} "
                f"result={receipt['result']} elapsed={receipt['elapsed_seconds']:.6f}"
            )
            if receipt["result"] != "PASS":
                overall_result = "COUNTEREXAMPLE"
                break
        if overall_result != "PASS":
            break
      if overall_result != "PASS":
        break

    started_masks = {
        int(row["component_support_mask_hex"], 16) for row in rows
    }
    completed_masks = {
        support_mask
        for support_mask in started_masks
        if len(
            [
                row
                for row in rows
                if int(row["component_support_mask_hex"], 16) == support_mask
                and row["result"] == "PASS"
            ]
        )
        == 16
    }
    target_masks = [
        support_mask
        for support_mask in range(1, 1 << len(COMPONENT_DEGREES))
        if sum(
            degree
            for index, degree in enumerate(COMPONENT_DEGREES)
            if (support_mask >> index) & 1
        )
        <= 30
    ]
    covered_masks = [
        support_mask
        for support_mask in target_masks
        if any(support_mask & ~completed_mask == 0 for completed_mask in completed_masks)
    ]
    if overall_result == "PASS":
        evidence_label = (
            "EXACT_EXHAUSTIVE_CERTIFICATE_RUN"
            if args.mode == "primary"
            else "EXACT_INDEPENDENT_EXHAUSTIVE_VERIFICATION_RUN"
        )
    else:
        evidence_label = "EXACT_COUNTEREXAMPLE_RUN"

    payload = {
        "schema": "riffle-dp-2lap-g4-component-endpoint-run-v2",
        "candidate": "Riffle DP-2Lap g=4",
        "mode": args.mode,
        "evidence_label": evidence_label,
        "runner_source_sha256": digest(Path(__file__).resolve()),
        "certificate_source_sha256": digest(source(args.mode)),
        "executable_sha256_by_dimension": executable_hashes,
        "planned_maximal_support_count": len(MAXIMAL_SUPPORTS),
        "started_maximal_support_count": len(started_masks),
        "fully_certified_maximal_support_count": len(completed_masks),
        "target_support_count": len(target_masks),
        "covered_support_count_by_completed_certificates": len(covered_masks),
        "certificate_cases_per_support": 16,
        "case_schedule": {
            "24_node_packet_values": list(range(1, 16)),
            "12_node_packet_values": [15],
        },
        "rows": rows,
        "total_information_vectors_checked": total_candidates,
        "total_certificate_elapsed_seconds": total_elapsed,
        "result": overall_result,
        "coverage_rule": (
            "Every nonempty component support of total dimension at most 30 is "
            "contained in at least one of the eight planned maximal supports. "
            "For each maximal support, values 1..14 require 24-node distance 97; "
            "value 15 requires 24-node distance 72 and 12-node endpoint distance "
            "36. Only fully passing maximal supports contribute to coverage."
        ),
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"output={args.output}")
    print(f"total_information_vectors_checked={total_candidates}")
    print(f"status={overall_result}")
    raise SystemExit(0 if overall_result == "PASS" else 1)


if __name__ == "__main__":
    main()
