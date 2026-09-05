#!/usr/bin/env python3
"""Audit the Goal 06 amortized accounting and sampled zero-anchor certificates."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

from analyze_riffle_packetmul_wrapmul_2lap_g4_lifted_windows import (
    build_apply,
    lifted_step,
)
from analyze_systematic_group_kernel import systematic_state_columns


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / "constructions" / "riffle_packetmul_wrapmul_2lap_g4"
RECEIPTS = CANDIDATE / "receipts"
OUTPUT = RECEIPTS / "goal06_amortized_route_audit.json"
ANCHORS = (0, 12, 23)
ZERO_NODES = 32_737
GAPS = 34
WINDOW = 24
WINDOW_BOUND = 144
TARGET = 188_765
REJECTED_WINDOW_WEIGHT = 143
MASK64 = (1 << 64) - 1
COUNTEREXAMPLE_A = 0x89CF7ABE4FBC2365
COUNTEREXAMPLE_B = 0x87452995C56BE123


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    minimum_blocks = math.ceil((ZERO_NODES - GAPS * (WINDOW - 1)) / WINDOW)
    implied_weight = minimum_blocks * WINDOW_BOUND
    if minimum_blocks != 1332 or implied_weight != 191_808:
        raise RuntimeError("goal06 audit: amortized accounting changed")
    if implied_weight <= TARGET:
        raise RuntimeError("goal06 audit: amortized target no longer closes")

    parity = build_apply(systematic_state_columns())
    a = COUNTEREXAMPLE_A
    b = COUNTEREXAMPLE_B
    counter_weights = []
    for _ in range(WINDOW):
        a, b, output = lifted_step(a, b, parity)
        counter_weights.append(output.bit_count())
    if sum(counter_weights[:18]) != 97 or sum(counter_weights) != 278:
        raise RuntimeError("goal06 audit: Goal 05 counterexample extension changed")

    native_source = (
        ROOT
        / "scripts"
        / "certify_riffle_packetmul_wrapmul_2lap_g4_goal06_zero_anchor.cpp"
    )
    native_executable = native_source.with_suffix(".exe")
    rows = []
    for anchor in ANCHORS:
        sets_path = RECEIPTS / f"goal06_zero_anchor_node_{anchor:02d}_information_sets.json"
        input_path = RECEIPTS / f"goal06_zero_anchor_node_{anchor:02d}_c24.json"
        raw_path = RECEIPTS / f"goal06_zero_anchor_node_{anchor:02d}_c24_raw.json"
        sets = json.loads(sets_path.read_text())
        prepared = json.loads(input_path.read_text())
        raw = json.loads(raw_path.read_text())
        binary_path = RECEIPTS / prepared["binary_file"]
        if digest(binary_path) != prepared["binary_sha256"]:
            raise RuntimeError("goal06 audit: binary input hash mismatch")
        if raw["result"] != "EXHAUSTED":
            raise RuntimeError("goal06 audit: zero-anchor counterexample appeared")
        set_count = prepared["information_set_count"]
        base_weight = prepared["base_restriction_weight"]
        base_sets = prepared["fixed_base_weight_set_count"]
        expected = (
            set_count
            * sum(math.comb(64, weight) for weight in range(base_weight))
            + base_sets * math.comb(64, base_weight)
        )
        if (
            base_sets != REJECTED_WINDOW_WEIGHT % set_count + 1
            or base_weight != REJECTED_WINDOW_WEIGHT // set_count
            or expected != prepared["expected_candidate_count"]
            or expected != raw["expected_candidates"]
            or expected != raw["candidates"]
        ):
            raise RuntimeError("goal06 audit: coverage profile mismatch")
        if sets["full_information_set_count"] != set_count:
            raise RuntimeError("goal06 audit: set count mismatch")
        rows.append(
            {
                "anchor_node": anchor,
                "full_information_sets": set_count,
                "base_restriction_weight": base_weight,
                "fixed_base_weight_set_count": base_sets,
                "checked_candidates": expected,
                "native_elapsed_seconds": raw["elapsed_seconds"],
                "sets_receipt_sha256": digest(sets_path),
                "input_receipt_sha256": digest(input_path),
                "binary_sha256": digest(binary_path),
                "raw_receipt_sha256": digest(raw_path),
                "result": "NO_WEIGHT_AT_MOST_143_WORD_IN_ZERO_ANCHOR_SLICE",
            }
        )

    payload = {
        "schema": "riffle-packetmul-wrapmul-2lap-g4-goal06-route-audit-v1",
        "candidate": "Riffle PacketMul-WrapMul-2Lap g=4",
        "evidence_label": "EXACT_ACCOUNTING_AND_SAMPLED_ZERO_ANCHOR_CERTIFICATES",
        "source_sha256": digest(Path(__file__).resolve()),
        "native_source_sha256": digest(native_source),
        "native_executable_sha256": digest(native_executable),
        "amortized_target": {
            "zero_nodes": ZERO_NODES,
            "maximum_gaps": GAPS,
            "window_nodes": WINDOW,
            "proposed_window_lower_bound": WINDOW_BOUND,
            "minimum_complete_windows": minimum_blocks,
            "implied_output_weight": implied_weight,
            "required_output_weight": TARGET,
            "strict_margin": implied_weight - TARGET,
        },
        "goal05_counterexample_extension": {
            "initial_a_hex": f"0x{COUNTEREXAMPLE_A:016x}",
            "initial_b_hex": f"0x{COUNTEREXAMPLE_B:016x}",
            "node_weights": counter_weights,
            "first_18_weight": sum(counter_weights[:18]),
            "first_24_weight": sum(counter_weights),
        },
        "zero_anchor_certificates": rows,
        "total_checked_candidates": sum(row["checked_candidates"] for row in rows),
        "scope_limitation": (
            "The accounting proves that D(24) >= 144 would close the zero-gap "
            "row. The exact certificates cover only anchors 0, 12, and 23 with "
            "anchor output equal to zero. They do not prove D(24) >= 144."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"minimum_complete_windows={minimum_blocks}")
    print(f"implied_output_weight={implied_weight}")
    print(f"strict_margin={implied_weight - TARGET}")
    print(f"total_checked_candidates={payload['total_checked_candidates']}")
    print(f"output={OUTPUT}")
    print("status=EXACT_AMORTIZED_ROUTE_AUDIT")


if __name__ == "__main__":
    main()
