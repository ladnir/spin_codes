#!/usr/bin/env python3
"""Run the four endpoint-aware dimension-31 frontier certificates sequentially."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXPLORATIONS = ROOT / "explorations"
SOURCE = ROOT / "scripts" / "certify_riffle_dp_2lap_g4_component_24node.cpp"
FRONTIER_SUPPORTS = (
    (0x2C, (4, 9, 18)),
    (0x33, (1, 2, 10, 18)),
    (0x4A, (2, 9, 20)),
    (0x51, (1, 10, 20)),
)
SCHEDULE = ((24, tuple(range(1, 16))), (12, (15,)))
DIMENSION = 31


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def executable(mode: str, window_nodes: int) -> Path:
    stem = "certify" if mode == "primary" else "verify"
    return ROOT / "scripts" / (
        f"{stem}_riffle_dp_2lap_g4_component_w{window_nodes}_d{DIMENSION}.exe"
    )


def receipt_path(
    mode: str,
    window_nodes: int,
    support_mask: int,
    packet_value: int,
    probe_only: bool,
) -> Path:
    stage = "probe_" if probe_only else ""
    return EXPLORATIONS / (
        f"riffle_dp_2lap_g4_component31_w{window_nodes}_"
        f"{stage}{mode}_s{support_mask:02x}_v{packet_value:02d}.json"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("primary", "independent"), required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--probe-only",
        action="store_true",
        help="Run value 15 at both windows and packet values 1 and 7 at 24 nodes.",
    )
    args = parser.parse_args()

    if args.probe_only:
        schedule = ((24, (1, 7, 15)), (12, (15,)))
    else:
        schedule = SCHEDULE

    rows = []
    total_candidates = 0
    total_elapsed = 0.0
    result = "PASS"
    executable_hashes = {
        f"w{window}_d{DIMENSION}": digest(executable(args.mode, window))
        for window in (24, 12)
    }
    for support_mask, degrees in FRONTIER_SUPPORTS:
        for window_nodes, packet_values in schedule:
            binary = executable(args.mode, window_nodes)
            for packet_value in packet_values:
                path = receipt_path(
                    args.mode,
                    window_nodes,
                    support_mask,
                    packet_value,
                    args.probe_only,
                )
                completed = subprocess.run(
                    [str(binary), hex(support_mask), str(packet_value), str(path)],
                    cwd=ROOT,
                    capture_output=True,
                    text=True,
                    check=False,
                )
                if completed.returncode not in (0, 1):
                    raise RuntimeError(
                        "dimension-31 runner: executable failed\n"
                        f"{completed.stdout}{completed.stderr}"
                    )
                receipt = json.loads(path.read_text())
                if (
                    int(receipt["component_support_mask_hex"], 16) != support_mask
                    or tuple(receipt["component_degrees"]) != degrees
                    or receipt["code_dimension"] != DIMENSION
                    or receipt["window_nodes"] != window_nodes
                    or receipt["packet_value"] != packet_value
                ):
                    raise RuntimeError("dimension-31 runner: receipt identity mismatch")
                total_candidates += receipt["total_information_vectors_checked"]
                total_elapsed += receipt["elapsed_seconds"]
                rows.append(
                    {
                        "component_support_mask_hex": hex(support_mask),
                        "component_degrees": list(degrees),
                        "code_dimension": DIMENSION,
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
                    f"mode={args.mode} window={window_nodes} "
                    f"support={hex(support_mask)} value={packet_value} "
                    f"result={receipt['result']} "
                    f"elapsed={receipt['elapsed_seconds']:.6f}",
                    flush=True,
                )
                if receipt["result"] != "PASS":
                    result = "COUNTEREXAMPLE"
                    break
            if result != "PASS":
                break
        if result != "PASS":
            break

    expected_cases = len(FRONTIER_SUPPORTS) * sum(
        len(values) for _, values in schedule
    )
    completed = result == "PASS" and len(rows) == expected_cases
    payload = {
        "schema": "riffle-dp-2lap-g4-component-dimension31-run-v1",
        "candidate": "Riffle DP-2Lap g=4",
        "mode": args.mode,
        "probe_only": args.probe_only,
        "evidence_label": (
            "EXACT_OBSTRUCTION_PROBE"
            if args.probe_only and result == "PASS"
            else (
                "EXACT_EXHAUSTIVE_CERTIFICATE_RUN"
                if args.mode == "primary" and result == "PASS"
                else (
                    "EXACT_INDEPENDENT_EXHAUSTIVE_VERIFICATION_RUN"
                    if result == "PASS"
                    else "EXACT_COUNTEREXAMPLE_RUN"
                )
            )
        ),
        "runner_source_sha256": digest(Path(__file__).resolve()),
        "certificate_source_sha256": digest(SOURCE),
        "executable_sha256_by_window_and_dimension": executable_hashes,
        "frontier_support_count": len(FRONTIER_SUPPORTS),
        "frontier_support_masks_hex": [hex(mask) for mask, _ in FRONTIER_SUPPORTS],
        "expected_case_count": expected_cases,
        "completed_case_count": len(rows),
        "rows": rows,
        "total_information_vectors_checked": total_candidates,
        "total_certificate_elapsed_seconds": total_elapsed,
        "completed": completed,
        "result": result,
        "scope_limitation": (
            "This run covers only the four dimension-31 frontier supports. "
            "The complete support-through-dimension-31 conclusion also uses "
            "the audited dimension-at-most-30 certificate."
        ),
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"output={args.output}")
    print(f"total_information_vectors_checked={total_candidates}")
    print(f"status={result}")
    raise SystemExit(0 if result == "PASS" else 1)


if __name__ == "__main__":
    main()
