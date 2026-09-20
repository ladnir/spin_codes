#!/usr/bin/env python3
"""Outward-rounded full-spectrum certificate for the bit-transpose variant.

For ``a`` active outer blocks, one transposed row contains a uniform
``a``-subset of the ``L`` bit positions.  The coefficient of degree ``a`` in

    [sum_r C(g,r) x^r M_r]^(L/g)

is the unnormalised tilted row transition.  This checker computes every
coefficient with upward-rounded scaled binary64 arithmetic and divides by the
exact subset count C(L,a).  Fixed Chernoff tilts come from an exploratory
receipt; the checker does no optimisation.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np

from certify_riffle_transpose_packetshuffle_full import (
    decimal,
    fixed_z_text,
    interval_lower,
    interval_upper,
    outer_log2_interval,
    scaled_add,
    scaled_from_float,
    scaled_matrix_multiply_batch,
    scaled_matrix_power_moment_log2_upper,
    scaled_multiply,
    scaled_normalize,
    scaled_weight,
    transition_upper,
)


DEFAULT_SOURCE = Path(
    "constructions/riffle_transpose_bitshuffle_randomstepconv/"
    "receipts/g8_b1024_sigma12_exploratory.json"
)
DEFAULT_OUTPUT = Path(
    "constructions/riffle_transpose_bitshuffle_randomstepconv/"
    "receipts/g8_b1024_sigma12_full_interval.json"
)


def scaled_divide_integer_lower(value, denominator: int):
    """Upper-bound value / denominator using an exact lower truncation."""
    if denominator <= 0:
        raise ValueError("denominator must be positive")
    bit_length = denominator.bit_length()
    shift = max(0, bit_length - 53)
    leading = denominator >> shift
    denominator_scaled = scaled_from_float(np.asarray(float(leading)))
    denominator_scaled = (
        denominator_scaled[0], denominator_scaled[1] + shift
    )
    numerator_mantissa, numerator_exponent = value
    denominator_mantissa, denominator_exponent = denominator_scaled
    quotient = np.divide(numerator_mantissa, denominator_mantissa)
    quotient = np.nextafter(quotient, np.float64(math.inf))
    quotient = np.where(numerator_mantissa == 0.0, 0.0, quotient)
    return scaled_normalize(
        quotient, numerator_exponent - denominator_exponent
    )


def normalized_bitshuffle_rows_scaled_upper(
    *,
    packet_bits: int,
    sigma: int,
    blocks: int,
    z_text: str,
    max_degree: int | None = None,
    reverse: bool = False,
):
    """Return upward bounds through ``max_degree``.

    With ``reverse=True``, degree ``d`` denotes original occupation
    ``blocks-d``.  This computes high-occupation tails from the cheaper end of
    the polynomial without changing matrix multiplication order.
    """
    if blocks % packet_bits:
        raise ValueError("packet width must divide the outer-block count")
    positions = blocks // packet_bits
    if max_degree is None:
        max_degree = blocks
    if not 0 <= max_degree <= blocks:
        raise ValueError("invalid maximum degree")
    matrices = [
        scaled_from_float(
            transition_upper(
                packet_bits,
                sigma,
                z_text,
                packet_bits - rank if reverse else rank,
            )
        )
        for rank in range(packet_bits + 1)
    ]
    multiplicities = [math.comb(packet_bits, rank) for rank in range(packet_bits + 1)]

    mantissa = np.zeros((max_degree + 1, 2, 2), dtype=np.float64)
    exponent = np.zeros((max_degree + 1, 2, 2), dtype=np.int64)
    mantissa[0] = np.eye(2, dtype=np.float64)
    coefficients = scaled_normalize(mantissa, exponent)

    for position in range(positions):
        old_count = min(position * packet_bits, max_degree) + 1
        updated = (
            np.zeros_like(coefficients[0]),
            np.zeros_like(coefficients[1]),
        )
        for rank in range(packet_bits + 1):
            if rank > max_degree:
                break
            source_count = min(old_count, max_degree - rank + 1)
            old = (
                coefficients[0][:source_count],
                coefficients[1][:source_count],
            )
            terms = scaled_matrix_multiply_batch(old, matrices[rank])
            terms = scaled_weight(terms, float(multiplicities[rank]))
            destination = (
                updated[0][rank : rank + source_count],
                updated[1][rank : rank + source_count],
            )
            combined = scaled_add(destination, terms)
            updated[0][rank : rank + source_count] = combined[0]
            updated[1][rank : rank + source_count] = combined[1]
        coefficients = updated

    normalized_mantissa = np.empty_like(coefficients[0])
    normalized_exponent = np.empty_like(coefficients[1])
    for degree in range(max_degree + 1):
        normalized = scaled_divide_integer_lower(
            (
                coefficients[0][degree],
                coefficients[1][degree],
            ),
            math.comb(blocks, degree),
        )
        normalized_mantissa[degree] = normalized[0]
        normalized_exponent[degree] = normalized[1]
    return normalized_mantissa, normalized_exponent


def self_test() -> dict[str, float]:
    from analyze_riffle_transpose_bitshuffle import numeric_subset_average

    maximum_gap = 0.0
    for packet_bits, blocks, sigma, z_text in (
        (2, 6, 3, "0.37"),
        (3, 6, 4, "0.61"),
    ):
        rows = normalized_bitshuffle_rows_scaled_upper(
            packet_bits=packet_bits,
            sigma=sigma,
            blocks=blocks,
            z_text=z_text,
            max_degree=blocks,
        )
        for active_blocks in range(blocks + 1):
            recovered = np.ldexp(
                rows[0][active_blocks], rows[1][active_blocks]
            )
            exact = numeric_subset_average(
                packet_bits=packet_bits,
                sigma=sigma,
                blocks=blocks,
                active_blocks=active_blocks,
                z=float(z_text),
            )
            # The reference uses ordinary nearest rounding.  Permit its last
            # bit while requiring the certified calculation to cover it.
            if np.any(np.nextafter(recovered, math.inf) < exact):
                raise AssertionError("scaled row bound fell below reference")
            maximum_gap = max(
                maximum_gap, float(np.max(recovered - exact))
            )
        reversed_rows = normalized_bitshuffle_rows_scaled_upper(
            packet_bits=packet_bits,
            sigma=sigma,
            blocks=blocks,
            z_text=z_text,
            max_degree=blocks,
            reverse=True,
        )
        for active_blocks in range(blocks + 1):
            degree = blocks - active_blocks
            recovered = np.ldexp(
                reversed_rows[0][degree], reversed_rows[1][degree]
            )
            exact = numeric_subset_average(
                packet_bits=packet_bits,
                sigma=sigma,
                blocks=blocks,
                active_blocks=active_blocks,
                z=float(z_text),
            )
            if np.any(np.nextafter(recovered, math.inf) < exact):
                raise AssertionError("reversed row bound fell below reference")
            maximum_gap = max(
                maximum_gap, float(np.max(recovered - exact))
            )
    return {"maximum_small_row_upward_gap": maximum_gap}


def evaluate(source: Path, dps: int, tilt_quantum_tenths: int) -> dict[str, object]:
    mp.mp.dps = dps
    mp.iv.dps = dps
    source_payload = json.loads(source.read_text(encoding="utf-8"))
    case = source_payload["case"]
    message_bits = int(source_payload["message_bits"])
    outer_bits = int(case["outer_bits"])
    outer_dimension = int(case["outer_dimension_bits"])
    blocks = int(case["outer_blocks"])
    packet_bits = int(case["packet_bits"])
    positions_per_row = int(case["packet_positions_per_row"])
    sigma = int(case["sigma"])
    distance = int(case["distance"])

    tilt_by_occupation: dict[int, int] = {}
    for row in case["occupation_rows"]:
        active_blocks = int(row["active_outer_blocks"])
        source_tenth = int(round(float(row["log_output_surprisal"]) * 10))
        fixed_tenth = int(
            round(source_tenth / tilt_quantum_tenths) * tilt_quantum_tenths
        )
        tilt_by_occupation[active_blocks] = fixed_tenth

    occupations_by_tilt: dict[int, list[int]] = {}
    for active_blocks, fixed_tenth in tilt_by_occupation.items():
        occupations_by_tilt.setdefault(fixed_tenth, []).append(active_blocks)

    rows: list[dict[str, object]] = []
    total = mp.iv.mpf(0)
    completed = 0
    for fixed_tenth in sorted(occupations_by_tilt):
        z_text = fixed_z_text(fixed_tenth)
        occupations = occupations_by_tilt[fixed_tenth]
        reverse = min(occupations) > blocks // 2
        degrees = [blocks - active if reverse else active for active in occupations]
        row_matrices = normalized_bitshuffle_rows_scaled_upper(
            packet_bits=packet_bits,
            sigma=sigma,
            blocks=blocks,
            z_text=z_text,
            max_degree=max(degrees),
            reverse=reverse,
        )
        correction = -distance * mp.iv.log(mp.iv.mpf(z_text)) / mp.iv.log(2)
        for active_blocks in occupations:
            degree = blocks - active_blocks if reverse else active_blocks
            row_matrix = (
                row_matrices[0][degree],
                row_matrices[1][degree],
            )
            log_moment = scaled_matrix_power_moment_log2_upper(
                row_matrix, outer_bits
            )
            inner_upper = min(
                mp.mpf(0), interval_upper(log_moment + correction)
            )
            outer_upper = interval_upper(
                outer_log2_interval(
                    blocks=blocks,
                    outer_dimension=outer_dimension,
                    outer_bits=outer_bits,
                    active_blocks=active_blocks,
                )
            )
            point_upper = interval_upper(
                mp.iv.mpf(outer_upper) + mp.iv.mpf(inner_upper)
            )
            total += mp.iv.power(2, mp.iv.mpf(point_upper))
            rows.append(
                {
                    "active_outer_blocks": active_blocks,
                    "fixed_log_surprisal": decimal(mp.mpf(fixed_tenth) / 10, 8),
                    "fixed_z": z_text,
                    "outer_log2_upper": decimal(outer_upper),
                    "inner_log2_upper": decimal(inner_upper),
                    "pointwise_log2_upper": decimal(point_upper),
                }
            )
        completed += len(occupations_by_tilt[fixed_tenth])
        print(
            f"tilt,{fixed_tenth / 10:.1f},occupations,{completed},{blocks}",
            flush=True,
        )

    total_upper = interval_upper(total)
    lambda_lower = interval_lower(
        -mp.iv.log(mp.iv.mpf(total_upper)) / mp.iv.log(2)
    )
    rows.sort(key=lambda row: int(row["active_outer_blocks"]))
    dominant = sorted(
        rows,
        key=lambda row: mp.mpf(str(row["pointwise_log2_upper"])),
        reverse=True,
    )[:10]
    return {
        "schema": "riffle-transpose-bitshuffle-full-interval-v1",
        "candidate": f"Riffle TransposeBitShuffle-RandomStepConv g={packet_bits}",
        "source_tilt_receipt": str(source),
        "arithmetic": {
            "binary_row_dp": "IEEE-754 binary64 mantissa plus integer binary exponent; nextafter after every nonnegative mantissa add, multiply, and divide",
            "coefficient_normalization": "divide by C(L,a) using an exact 53-leading-bit lower truncation",
            "matrix_power_scaling": "per-entry exact powers of two",
            "scalar_and_sum": f"mpmath interval arithmetic at {dps} decimal digits",
            "optimizer_in_checker": False,
            "tilt_quantum_tenths": tilt_quantum_tenths,
        },
        "parameters": {
            "message_bits": message_bits,
            "output_bits": outer_bits * positions_per_row * packet_bits,
            "outer_bits": outer_bits,
            "outer_dimension_bits": outer_dimension,
            "outer_blocks": blocks,
            "packet_bits": packet_bits,
            "packet_positions_per_row": positions_per_row,
            "sigma": sigma,
            "distance": distance,
            "relative_distance_floor": distance / (
                outer_bits * positions_per_row * packet_bits
            ),
        },
        "row_formula": "coefficient x^a of [sum_r C(g,r) x^r M_r]^(L/g), divided by C(L,a)",
        "self_test": self_test(),
        "occupation_count": len(rows),
        "log2_expected_bad_upper": decimal(-lambda_lower),
        "lambda_bits_lower": decimal(lambda_lower),
        "target_lambda_bits": 40,
        "passes_target": bool(lambda_lower > 40),
        "dominant_certified_rows": dominant,
        "occupation_rows": rows,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument(
        "--tilt-quantum-tenths",
        type=int,
        default=10,
        help="round source log-surprisal tilts to this many tenths",
    )
    parser.add_argument("--self-test-only", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    mp.mp.dps = args.dps
    mp.iv.dps = args.dps
    if args.self_test_only:
        print(json.dumps(self_test(), indent=2, sort_keys=True))
        return
    if args.tilt_quantum_tenths < 1:
        raise ValueError("tilt quantum must be positive")
    payload = evaluate(args.source, args.dps, args.tilt_quantum_tenths)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "lambda_bits_lower": payload["lambda_bits_lower"],
                "passes_target": payload["passes_target"],
                "dominant_certified_rows": payload["dominant_certified_rows"],
            },
            indent=2,
            sort_keys=True,
        )
    )
    print(f"wrote,{args.output}", flush=True)


if __name__ == "__main__":
    main()
