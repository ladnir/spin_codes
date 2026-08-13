#!/usr/bin/env python3
"""Fast structural verifier for the diagnostic identity matrix ledger manifest."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DEFAULT_MANIFEST = ROOT / "identity_matrix_ledger_manifest.json"


def log2_sum(values: list[float]) -> float:
    maximum = max(values)
    return maximum + math.log2(sum(2 ** (value - maximum) for value in values))


def matrix_command(row: dict[str, object], *, power2_punctured: bool) -> str:
    mode = str(row["mode"])
    mode_flag = "" if mode == "first-active" else f" --mode {mode}"
    construction_flag = " --power2-punctured" if power2_punctured else ""
    return (
        "python scripts\\probe_identity_matrix_ledger.py"
        f"{construction_flag}{mode_flag}"
        f" --h-min {row['h_min']} --h-max {row['h_max']}"
        f" --h-step {row['h_step']} --output-pole {row['output_pole']}"
        f" --occupancy-pole {row['occupancy_pole']}"
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--print-recompute-commands", action="store_true")
    args = parser.parse_args()

    data = json.loads(args.manifest.read_text(encoding="utf-8"))
    n = int(data["N"])
    graph_codimension = int(data["graph_codimension"])
    target = float(data["target_failure_log2"])
    low = data["low_weight"]
    episode = data["episode_weight"]
    rows = data["matrix_intervals"]
    power2_punctured = int(data.get("punctured_blocks", 0)) > 0

    if (int(low["h_min"]), int(low["h_max"])) != (1, 500):
        raise SystemExit("identity matrix ledger: invalid low-weight range")
    expected_episode = (501, 2000) if power2_punctured else (502, 2000)
    expected_step = 1 if power2_punctured else 2
    if (int(episode["h_min"]), int(episode["h_max"])) != expected_episode:
        raise SystemExit("identity matrix ledger: invalid episode range")
    if int(episode["h_step"]) != expected_step:
        raise SystemExit("identity matrix ledger: invalid episode step")
    expected = 2001 if power2_punctured else 2002
    for row in rows:
        if int(row["h_step"]) != expected_step or int(row["h_min"]) != expected:
            raise SystemExit(f"identity matrix ledger: coverage gap before {row}")
        if not power2_punctured and (int(row["h_min"]) % 2 or int(row["h_max"]) % 2):
            raise SystemExit(f"identity matrix ledger: non-even endpoint in {row}")
        expected = int(row["h_max"]) + expected_step
    if expected != n + expected_step:
        raise SystemExit(
            f"identity matrix ledger: coverage stops at {expected - expected_step}, not {n}"
        )

    ambient_rows = [
        float(low["ambient_log2_upper_approx"]),
        float(episode["ambient_log2_upper_approx"]),
        *(float(row["ambient_log2_upper_approx"]) for row in rows),
    ]
    ambient = log2_sum(ambient_rows)
    after_graph = ambient - graph_codimension
    margin = target - after_graph
    stored = data["combined"]
    if abs(ambient - float(stored["ambient_log2_upper_approx"])) > 1e-9:
        raise SystemExit("identity matrix ledger: stored ambient total is stale")
    if abs(after_graph - float(stored["after_graph_log2_upper_approx"])) > 1e-9:
        raise SystemExit("identity matrix ledger: stored graph total is stale")
    if margin <= 0:
        raise SystemExit("identity matrix ledger: diagnostic target does not pass")

    print("identity_matrix_ledger_status,DIAGNOSTIC_PASS_NOT_FORMAL")
    if power2_punctured:
        print(f"supported_weight_cover,every integer weight 1--{n}")
    else:
        print(f"supported_weight_cover,1--500 and every even weight 502--{n}")
    print(f"ambient_total_log2_upper_approx,{ambient:.12f}")
    print(f"after_graph_r24_log2_upper_approx,{after_graph:.12f}")
    print(f"margin_beyond_40_bits_approx,{margin:.12f}")
    print("remaining_hardening,outward/exact arithmetic for all floating rows")
    if args.print_recompute_commands:
        print(str(low["command"]))
        for command in episode["commands"]:
            print(str(command))
        for row in rows:
            print(matrix_command(row, power2_punctured=power2_punctured))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
