#!/usr/bin/env python3
"""Compare one common coefficient tilt with one optimized tilt per packet type."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

from analyze_riffle_small_exact_vs_type_bound import (
    ROOT,
    coefficient_bound_objective,
    histogram_sequence_count,
    log_fraction,
    optimize_coefficient_bound,
    outer_histogram_distribution,
)


DEFAULT_RECEIPTS = (
    ROOT
    / "constructions"
    / "riffle_shiftalpha64_bchblockperm_parallelacc_g4"
    / "receipts"
)


def single_type_objective(
    parameters: np.ndarray,
    distance: int,
    packet_positions: int,
    histogram: tuple[int, ...],
    log_total_sequences: float,
) -> float:
    row = ((histogram, 0.0, log_total_sequences),)
    return coefficient_bound_objective(parameters, distance, packet_positions, row)


def exact_log_from_receipt(data_blocks: int, parity_symbols: int, distance: int) -> float:
    path = DEFAULT_RECEIPTS / f"goal19_small_exact_vs_type_bound_b{data_blocks}.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    parity = next(
        row for row in payload["parity_results"] if row["parity_symbols"] == parity_symbols
    )
    return float(parity["exact_cumulative"][distance]["log2"])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-blocks", type=int, default=4)
    parser.add_argument("--parity-symbols", type=int, choices=(0, 1, 2), default=2)
    parser.add_argument("--distance", type=int, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    outer, _signatures = outer_histogram_distribution(
        args.data_blocks, args.parity_symbols
    )
    packet_positions = 2 * (args.data_blocks + args.parity_symbols)
    rows = tuple(
        (
            histogram,
            log_fraction(count),
            math.log(histogram_sequence_count(histogram)),
        )
        for histogram, count in outer.items()
    )
    common_log, common_parameters, common_success = optimize_coefficient_bound(
        args.distance, packet_positions, rows, None
    )

    bounds = [(-12.0, 3.0)] + [(-14.0, 14.0)] * 4
    contributions = []
    failures = 0
    for histogram, log_outer, log_total in rows:
        result = minimize(
            single_type_objective,
            common_parameters,
            args=(
                args.distance,
                packet_positions,
                histogram,
                log_total,
            ),
            method="Nelder-Mead",
            bounds=bounds,
            options={"maxiter": 2500, "xatol": 1e-9, "fatol": 1e-9},
        )
        if not result.success:
            fallback = minimize(
                single_type_objective,
                result.x,
                args=(
                    args.distance,
                    packet_positions,
                    histogram,
                    log_total,
                ),
                method="Powell",
                bounds=bounds,
                options={"maxiter": 2000, "xtol": 1e-9, "ftol": 1e-9},
            )
            if float(fallback.fun) <= float(result.fun):
                result = fallback
        failures += int(not result.success)
        contributions.append(log_outer + float(result.fun))
    typewise_log = float(logsumexp(contributions))
    exact_log2 = exact_log_from_receipt(
        args.data_blocks, args.parity_symbols, args.distance
    )
    common_log2 = common_log / math.log(2.0)
    typewise_log2 = typewise_log / math.log(2.0)
    if typewise_log2 + 1e-8 < exact_log2:
        raise RuntimeError("the typewise numerical bound fell below the exact tail")

    payload = {
        "schema": "riffle-small-typewise-bound-probe-v1",
        "evidence_label": "NUMERICAL_COEFFICIENT_UPPER_BOUND_COMPARISON",
        "data_blocks": args.data_blocks,
        "parity_symbols": args.parity_symbols,
        "packet_positions": packet_positions,
        "binary_output_length": 4 * packet_positions,
        "distance": args.distance,
        "packet_types": len(rows),
        "exact_cumulative_log2": exact_log2,
        "common_tilt_bound_log2": common_log2,
        "common_tilt_loss_bits": common_log2 - exact_log2,
        "typewise_tilt_bound_log2": typewise_log2,
        "typewise_tilt_loss_bits": typewise_log2 - exact_log2,
        "improvement_from_typewise_tilts_bits": common_log2 - typewise_log2,
        "optimizer": {
            "common_success": common_success,
            "typewise_failures": failures,
        },
        "scope": (
            "Each packet histogram receives its own numerical coefficient tilt. "
            "The result measures avoidable common-tilt loss at a small exact "
            "instance. It is not a scalable algorithm or a proof certificate."
        ),
    }
    output = args.output or (
        DEFAULT_RECEIPTS
        / (
            f"goal19_typewise_b{args.data_blocks}_p{args.parity_symbols}"
            f"_d{args.distance}.json"
        )
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({"output": str(output), **payload}, indent=2))


if __name__ == "__main__":
    main()
