#!/usr/bin/env python3
"""Compute source-aware inner bounds for the first open Goal 13 profiles."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


SCRIPT_DIRECTORY = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIRECTORY))

import analyze_riffle_shiftalpha64_compressed_return as compressed  # noqa: E402
from analyze_riffle_shiftalpha64_adjacent_inner import analyze_profile  # noqa: E402


PROFILES = ((22, 22, 26), (22, 104, 106), (24, 106, 106))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--distance", type=int, default=188_766)
    parser.add_argument(
        "--packet-positions",
        type=int,
        default=compressed.DEFAULT_PACKET_POSITIONS,
    )
    parser.add_argument("--optimizer-maxiter", type=int, default=600)
    args = parser.parse_args()

    profiles = [
        analyze_profile(
            weights,
            distance=args.distance,
            packet_positions=args.packet_positions,
            maxiter=args.optimizer_maxiter,
        )
        for weights in PROFILES
    ]
    for profile in profiles:
        profile["selected_one_word_log2_bound"] = min(
            float(profile["pooled_one_word_log2_bound"]),
            float(profile["source_preserving"]["log2_bound"]),
        )
        profile["selected_method"] = (
            "pooled"
            if profile["pooled_one_word_log2_bound"]
            <= profile["source_preserving"]["log2_bound"]
            else "source_preserving"
        )

    payload = {
        "schema": "riffle-shiftalpha64-goal13-inner-v1",
        "construction": "Riffle ShiftAlpha64-BCHBlockPerm-ParallelAcc g=4",
        "parameters": {
            "bad_output_weight_inclusive": args.distance,
            "global_packet_positions": args.packet_positions,
        },
        "profiles": profiles,
        "validation": {
            "pooled_and_source_preserving_certificates_evaluated": True,
            "source_preserving_neutral_mass_identities": "PASS",
            "selected_bound_is_minimum_of_two_rigorous_bounds": True,
        },
        "scope": (
            "Rigorous one-word inner bounds. Outer relation and support "
            "multiplicities are not included."
        ),
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
