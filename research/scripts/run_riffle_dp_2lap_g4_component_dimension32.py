#!/usr/bin/env python3
"""Run the four endpoint-aware dimension-32 frontier certificates sequentially."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXPLORATIONS = ROOT / "explorations"
SOURCE = ROOT / "scripts" / "certify_riffle_dp_2lap_g4_component_24node.cpp"
FRONTIER = (
    (0x2D, (1, 4, 9, 18)),
    (0x34, (4, 10, 18)),
    (0x4B, (1, 2, 9, 20)),
    (0x52, (2, 10, 20)),
)
DIMENSION = 32
FULL_SCHEDULE = ((24, tuple(range(1, 16))), (12, (15,)))
PROBE_SCHEDULE = ((24, (1, 7, 15)), (12, (15,)))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def executable(mode: str, window: int) -> Path:
    stem = "certify" if mode == "primary" else "verify"
    return ROOT / "scripts" / (
        f"{stem}_riffle_dp_2lap_g4_component_w{window}_d{DIMENSION}.exe"
    )


def receipt_path(
    mode: str, window: int, mask: int, value: int, probe_only: bool
) -> Path:
    stage = "probe_" if probe_only else ""
    return EXPLORATIONS / (
        f"riffle_dp_2lap_g4_component32_w{window}_"
        f"{stage}{mode}_s{mask:02x}_v{value:02d}.json"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("primary", "independent"), required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--probe-only", action="store_true")
    args = parser.parse_args()
    schedule = PROBE_SCHEDULE if args.probe_only else FULL_SCHEDULE

    rows = []
    total_vectors = 0
    total_elapsed = 0.0
    result = "PASS"
    executable_hashes = {
        f"w{window}_d{DIMENSION}": digest(executable(args.mode, window))
        for window in (24, 12)
    }
    for mask, degrees in FRONTIER:
        for window, values in schedule:
            binary = executable(args.mode, window)
            for value in values:
                path = receipt_path(args.mode, window, mask, value, args.probe_only)
                completed = subprocess.run(
                    [str(binary), hex(mask), str(value), str(path)],
                    cwd=ROOT,
                    capture_output=True,
                    text=True,
                    check=False,
                )
                if completed.returncode not in (0, 1):
                    raise RuntimeError(
                        "dimension-32 runner: executable failed\n"
                        f"{completed.stdout}{completed.stderr}"
                    )
                receipt = json.loads(path.read_text())
                if (
                    int(receipt["component_support_mask_hex"], 16) != mask
                    or tuple(receipt["component_degrees"]) != degrees
                    or receipt["code_dimension"] != DIMENSION
                    or receipt["window_nodes"] != window
                    or receipt["packet_value"] != value
                ):
                    raise RuntimeError("dimension-32 runner: receipt identity mismatch")
                vectors = receipt["total_information_vectors_checked"]
                total_vectors += vectors
                total_elapsed += receipt["elapsed_seconds"]
                rows.append(
                    {
                        "component_support_mask_hex": hex(mask),
                        "component_degrees": list(degrees),
                        "code_dimension": DIMENSION,
                        "window_nodes": window,
                        "packet_value": value,
                        "receipt": str(path.relative_to(ROOT)).replace("\\", "/"),
                        "receipt_sha256": digest(path),
                        "information_set_count": receipt["information_set_count"],
                        "information_radius": receipt["information_radius"],
                        "information_vectors_checked": vectors,
                        "minimum_enumerated_side": min(
                            receipt["minimum_enumerated_weight"],
                            receipt["minimum_enumerated_complement_weight"],
                        ),
                        "result": receipt["result"],
                    }
                )
                print(
                    f"mode={args.mode} window={window} support={hex(mask)} "
                    f"value={value} result={receipt['result']} "
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

    expected_cases = len(FRONTIER) * sum(len(values) for _, values in schedule)
    completed = result == "PASS" and len(rows) == expected_cases
    evidence = (
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
    )
    payload = {
        "schema": "riffle-dp-2lap-g4-component-dimension32-run-v1",
        "candidate": "Riffle DP-2Lap g=4",
        "mode": args.mode,
        "probe_only": args.probe_only,
        "evidence_label": evidence,
        "runner_source_sha256": digest(Path(__file__).resolve()),
        "certificate_source_sha256": digest(SOURCE),
        "executable_sha256_by_window_and_dimension": executable_hashes,
        "frontier_support_count": len(FRONTIER),
        "frontier_support_masks_hex": [hex(mask) for mask, _ in FRONTIER],
        "expected_case_count": expected_cases,
        "completed_case_count": len(rows),
        "rows": rows,
        "total_information_vectors_checked": total_vectors,
        "total_certificate_elapsed_seconds": total_elapsed,
        "completed": completed,
        "result": result,
        "scope_limitation": (
            "This run covers only the four dimension-32 frontier supports. "
            "The conclusion through dimension 32 also uses the audited "
            "dimension-at-most-31 certificate."
        ),
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"output={args.output}")
    print(f"total_information_vectors_checked={total_vectors}")
    print(f"status={result}")
    raise SystemExit(0 if result == "PASS" else 1)


if __name__ == "__main__":
    main()
