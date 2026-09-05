#!/usr/bin/env python3
"""Sweep the inner step and state sizes for the BCH-like spectrum model.

FieldMulStepConv with coefficients sampled from the complete field and a
uniform random Toeplitz step have the same transition law on one fixed,
nonzero input/state vector.  This script therefore evaluates both designs at
once.  The calculation covers exactly one active outer block.  It is a
floating-point diagnostic, not a complete distance certificate.
"""

from __future__ import annotations

import argparse
import contextlib
import io
import json
import math
from pathlib import Path

from analyze_riffle_spectrumperm_bitshuffle_oneblock import evaluate, self_test


DEFAULT_OUTPUT = Path(
    "constructions/riffle_bchperm_transpose_bitshuffle_structuredstepconv/"
    "receipts/one_active_t_s_sweep.json"
)

# The first four states for each step width straddle the observed 40-bit
# boundary.  Additional rows expose the s ~= t amortization regime.
DEFAULT_CASES = {
    4: (12, 16, 17, 18),
    8: (8, 14, 15, 16, 17),
    16: (12, 14, 15, 16, 17),
    32: (12, 13, 14, 15, 16, 24, 32),
}


def run_case(
    step_bits: int,
    state_bits: int,
    args: argparse.Namespace,
    *,
    relative_distance: float | None = None,
    grid_step: float | None = None,
) -> dict[str, object]:
    # The imported evaluator reports its numerical grid progress.  Suppress
    # that per-case stream so the sweep has one compact progress line.
    with contextlib.redirect_stdout(io.StringIO()):
        case = evaluate(
            message_bits=args.message_bits,
            outer_bits=args.outer_bits,
            packet_bits=step_bits,
            sigma=state_bits,
            relative_distance=(
                args.relative_distance if relative_distance is None else relative_distance
            ),
            grid_min=args.grid_min,
            grid_max=args.grid_max,
            grid_step=args.grid_step if grid_step is None else grid_step,
            spectrum_path=None,
            modeled_minimum_distance=args.modeled_minimum_distance,
        )
    output_bits = int(case["output_bits"])
    width = step_bits + state_bits
    result = {
        "step_bits": step_bits,
        "state_bits": state_bits,
        "map_width_bits": width,
        "relative_distance": case["distance"] / case["output_bits"],
        "distance_bits": case["distance"],
        "inner_steps": output_bits // step_bits,
        "coefficient_bits": (output_bits // step_bits) * width,
        "coefficient_bits_per_output_bit": width / step_bits,
        "quasilinear_proxy_per_output_bit": width * math.log2(width) / step_bits,
        "one_active_lambda_bits": case["one_active_lambda_bits_lower_float"],
        "dominant_outer_weight": case["dominant_weight"]["outer_weight"],
        "dominant_pointwise_log2": case["dominant_weight"]["pointwise_log2_upper"],
        "dominant_log_output_surprisal": case["dominant_weight"][
            "log_output_surprisal"
        ],
    }
    print(
        f"t={step_bits},s={state_bits},d={width},"
        f"lambda={float(result['one_active_lambda_bits']):.6f}",
        flush=True,
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=1 << 20)
    parser.add_argument("--outer-bits", type=int, default=256)
    parser.add_argument("--modeled-minimum-distance", type=int, default=38)
    parser.add_argument("--relative-distance", type=float, default=0.09)
    parser.add_argument("--grid-min", type=float, default=-10.0)
    parser.add_argument("--grid-max", type=float, default=1.0)
    parser.add_argument("--grid-step", type=float, default=0.5)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    checks = self_test()
    rows = [
        run_case(step_bits, state_bits, args)
        for step_bits, states in DEFAULT_CASES.items()
        for state_bits in states
    ]
    frontier = [
        run_case(
            32,
            32,
            args,
            relative_distance=relative_distance,
            grid_step=0.25,
        )
        for relative_distance in (0.13, 0.14, 0.1405, 0.141, 0.15)
    ]
    payload = {
        "schema": "riffle-bchperm-bitshuffle-structuredstep-one-active-sweep-v1",
        "construction_families": [
            "BCHPerm-TransposeBitShuffle-FieldMulStepConv(t,s)",
            "BCHPerm-TransposeBitShuffle-ToeplitzStepConv(t,s)",
        ],
        "probability_space": {
            "outer_spectrum": (
                "modeled real-valued complement-symmetric even "
                f"[{args.outer_bits},{args.outer_bits // 2},"
                f"{args.modeled_minimum_distance}]-shaped spectrum"
            ),
            "outer_coordinate_permutations": "independent by outer block",
            "region_bit_permutations": "independent by transposed region",
            "inner_maps": "independent by inner step",
        },
        "parameters": {
            "message_bits": args.message_bits,
            "outer_bits": args.outer_bits,
            "relative_distance": args.relative_distance,
            "grid": [args.grid_min, args.grid_max, args.grid_step],
        },
        "self_test": checks,
        "rows": rows,
        "preferred_point_distance_frontier": frontier,
        "scope": (
            "One active outer block only. The spectrum is a real-valued "
            "model, and the Chernoff grid uses ordinary floating-point "
            "arithmetic. FieldMul includes the zero coefficient so that its "
            "one-vector transition equals the Toeplitz transition."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
