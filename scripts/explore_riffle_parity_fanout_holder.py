#!/usr/bin/env python3
"""Screen sparse parity-fanout shears against the all-active Holder bound."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from analyze_riffle_allone_tail_body_holder import logsumexp
from analyze_riffle_bitshuffle_splitstate_regular_bulk import (
    LOG2,
    log_choose,
    log_two_power_minus_one,
    modeled_even_floor_spectrum_logs,
)
from derive_riffle_parity_break_outer_spectrum import (
    expected_output_subset_spectrum,
)


def reference_inner_log(bulk_paths: list[Path], occupation: int) -> float:
    for path in bulk_paths:
        receipt = json.loads(path.read_text(encoding="utf-8"))
        envelopes = {
            float(row["probability"]): float(row["log2_mass"]) * LOG2
            for row in receipt["envelope"]["bernoulli_envelopes"]
        }
        for row in receipt["occupation_rows"]:
            if int(row["active_regular_outer_blocks"]) != occupation:
                continue
            probability = float(row["best_candidate_probability"])
            return (
                float(row["inner_log2_upper"]) * LOG2
                - occupation * envelopes[probability]
            )
    raise ValueError(f"occupation {occupation} not found in bulk receipts")


def all_active_margin(
    spectrum: np.ndarray,
    inner_log: float,
    occupation: int,
    p_values: np.ndarray,
    outer_bits: int,
    dimension: int,
) -> tuple[float, float]:
    log_mass = log_two_power_minus_one(dimension)
    reference_log_density = log_mass - outer_bits * LOG2
    reference: list[float] = []
    likelihood: list[float] = []
    for weight in range(1, outer_bits + 1):
        multiplicity = float(spectrum[weight])
        if not math.isfinite(multiplicity):
            continue
        shell = log_choose(outer_bits, weight)
        reference.append(shell - outer_bits * LOG2)
        likelihood.append(multiplicity - shell - reference_log_density)
    reference_array = np.asarray(reference)
    likelihood_array = np.asarray(likelihood)
    inverse_p = 1.0 / p_values
    q_values = 1.0 - inverse_p
    body_moments = np.asarray(
        [
            logsumexp(reference_array + float(p) * likelihood_array)
            for p in p_values
        ]
    )
    body_terms = log_mass + body_moments * inverse_p
    candidates = occupation * body_terms + q_values * inner_log
    best = int(np.argmin(candidates))
    return -float(candidates[best]) / LOG2, float(p_values[best])


def parse_odd_range(text: str) -> list[int]:
    start, stop = (int(value) for value in text.split(":"))
    return [value for value in range(start, stop + 1) if value & 1]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bulk", type=Path, nargs="+", required=True)
    parser.add_argument("--source-sizes", default="1:63")
    parser.add_argument("--target-sizes", default="21:63")
    parser.add_argument("--outer-bits", type=int, default=256)
    parser.add_argument("--dimension", type=int, default=128)
    parser.add_argument("--minimum-distance", type=int, default=38)
    parser.add_argument("--occupation", type=int, default=8192)
    parser.add_argument("--holder-grid-points", type=int, default=1025)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source_sizes = parse_odd_range(args.source_sizes)
    target_sizes = parse_odd_range(args.target_sizes)
    source = modeled_even_floor_spectrum_logs(
        args.outer_bits, args.minimum_distance
    )
    inner_log = reference_inner_log(args.bulk, args.occupation)
    p_values = 1.0 + np.exp2(
        np.linspace(-12.0, 28.0, args.holder_grid_points)
    )

    rows = []
    for source_size in source_sizes:
        for target_size in target_sizes:
            if source_size + target_size > args.outer_bits:
                continue
            spectrum = expected_output_subset_spectrum(
                source, source_size, target_size
            )
            margin, holder_p = all_active_margin(
                spectrum,
                inner_log,
                args.occupation,
                p_values,
                args.outer_bits,
                args.dimension,
            )
            rows.append(
                {
                    "source_size": source_size,
                    "target_size": target_size,
                    "transpose_xors_per_outer_block": (
                        source_size + target_size - 1
                    ),
                    "all_active_margin_bits": margin,
                    "best_holder_p": holder_p,
                }
            )

    rows.sort(
        key=lambda row: (
            -float(row["all_active_margin_bits"]),
            int(row["transpose_xors_per_outer_block"]),
        )
    )
    payload = {
        "schema": "riffle-parity-fanout-holder-screen-v1",
        "parameters": {
            "outer_bits": args.outer_bits,
            "dimension": args.dimension,
            "minimum_distance": args.minimum_distance,
            "occupation": args.occupation,
            "holder_grid_points": args.holder_grid_points,
            "source_sizes": source_sizes,
            "target_sizes": target_sizes,
            "bulk_uniform_reference": [str(path) for path in args.bulk],
        },
        "best_rows": rows[:50],
        "all_rows": rows,
        "scope": (
            "Nearest-binary64 screen of the all-active exact-spectrum Holder "
            "bound. Each map adds the parity on a uniformly sampled source "
            "set to every coordinate in a disjoint uniformly sampled target "
            "set. Both set sizes are odd, so the map is invertible and breaks "
            "the source parity constraint. This is a screening receipt, not "
            "an end-to-end certificate."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    for row in rows[:20]:
        print(
            f"source,{row['source_size']},target,{row['target_size']},"
            f"xors,{row['transpose_xors_per_outer_block']},"
            f"margin,{row['all_active_margin_bits']:.6f},"
            f"p,{row['best_holder_p']:.6f}"
        )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
