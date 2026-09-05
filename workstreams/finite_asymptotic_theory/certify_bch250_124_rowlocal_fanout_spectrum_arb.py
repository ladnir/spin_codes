#!/usr/bin/env python3
"""Arb certificate for a row-local fanout spectrum and optional tail event."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from flint import arb, arb_mat, ctx, fmpq, fmpq_mat


WORKSTREAM = Path(__file__).resolve().parent
B = 250
K = 124
SOURCE_SIZE = 31
TARGET_SIZE = 33
LAYERS = 256
TOTAL_MASS = (1 << K) - 1
DEFAULT_CAPS = (
    WORKSTREAM / "shortened_bch250_124_subcode_constant_weight_envelope.json"
)


def exact_transition() -> fmpq_mat:
    matrix = fmpq_mat(B + 1, B + 1)
    source_denominator = math.comb(B, SOURCE_SIZE)
    remaining_length = B - SOURCE_SIZE
    target_denominator = math.comb(remaining_length, TARGET_SIZE)
    for weight in range(B + 1):
        for overlap in range(SOURCE_SIZE + 1):
            if overlap > weight:
                continue
            if SOURCE_SIZE - overlap > B - weight:
                continue
            source_probability = fmpq(
                math.comb(weight, overlap)
                * math.comb(B - weight, SOURCE_SIZE - overlap),
                source_denominator,
            )
            if overlap % 2 == 0:
                matrix[weight, weight] += source_probability
                continue
            remaining_ones = weight - overlap
            for target_overlap in range(TARGET_SIZE + 1):
                if target_overlap > remaining_ones:
                    continue
                if (
                    TARGET_SIZE - target_overlap
                    > remaining_length - remaining_ones
                ):
                    continue
                target_probability = fmpq(
                    math.comb(remaining_ones, target_overlap)
                    * math.comb(
                        remaining_length - remaining_ones,
                        TARGET_SIZE - target_overlap,
                    ),
                    target_denominator,
                )
                output_weight = (
                    weight + TARGET_SIZE - 2 * target_overlap
                )
                matrix[weight, output_weight] += (
                    source_probability * target_probability
                )
        if sum(matrix[weight, column] for column in range(B + 1)) != 1:
            raise ArithmeticError(f"transition row {weight} is not stochastic")
    return matrix


def load_caps(path: Path) -> dict[int, int]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    caps = {}
    for row in payload["spectrum"]:
        weight = int(row["weight"])
        if weight == 0:
            continue
        value = int(row["multiplicity_upper"])
        if value > 0:
            caps[weight] = value
    if sum(caps.values()) < TOTAL_MASS:
        raise ArithmeticError("shell caps do not contain the exact code mass")
    return caps


def exact_float_arb(value: float) -> arb:
    numerator, denominator = value.as_integer_ratio()
    return arb(fmpq(numerator, denominator))


def dual_shell_upper(
    layered: arb_mat,
    output_weight: int,
    caps: dict[int, int],
) -> tuple[arb, int]:
    # A nearest-float sort proposes a dual threshold only.  The dual formula
    # is valid for every threshold and does not assume that this order is exact.
    order = sorted(
        caps,
        key=lambda source_weight: float(layered[source_weight, output_weight]),
        reverse=True,
    )
    remaining = TOTAL_MASS
    marginal_weight = order[-1]
    for source_weight in order:
        allocated = min(caps[source_weight], remaining)
        remaining -= allocated
        marginal_weight = source_weight
        if remaining == 0:
            break
    if remaining:
        raise ArithmeticError("dual threshold search did not fill code mass")
    threshold = exact_float_arb(
        float(layered[marginal_weight, output_weight])
    )
    objective = threshold * TOTAL_MASS
    for source_weight, cap in caps.items():
        excess = layered[source_weight, output_weight] - threshold
        positive_part = (excess + abs(excess)) / 2
        objective += cap * positive_part
    return objective, marginal_weight


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--caps", type=Path, default=DEFAULT_CAPS)
    parser.add_argument("--precision", type=int, default=384)
    parser.add_argument("--layers", type=int, default=LAYERS)
    parser.add_argument("--eta-numerator", type=int, default=1)
    parser.add_argument("--eta-denominator", type=int, default=1000)
    parser.add_argument("--central-minimum-weight", type=int, default=1)
    parser.add_argument("--central-maximum-weight", type=int, default=B)
    parser.add_argument(
        "--tail-union-rows",
        type=int,
        default=0,
        help="multiply the expected tail multiplicity by this row count",
    )
    parser.add_argument("--tail-margin-bits", type=int, default=0)
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            WORKSTREAM
            / "bch250_124_parityfanout31x33_l256_spectrum_outward.json"
        ),
    )
    args = parser.parse_args()
    if args.precision < 128:
        raise ValueError("precision must be at least 128 bits")
    if args.layers < 0:
        raise ValueError("layers must be nonnegative")
    if args.eta_numerator < 0 or args.eta_denominator <= 0:
        raise ValueError("invalid eta")
    if not 1 <= args.central_minimum_weight <= args.central_maximum_weight <= B:
        raise ValueError("invalid central weight interval")
    if args.tail_union_rows < 0 or args.tail_margin_bits < 0:
        raise ValueError("invalid tail-event parameters")
    ctx.prec = args.precision
    ctx.threads = 1

    caps = load_caps(args.caps)
    transition = exact_transition()
    layered = arb_mat(transition) ** args.layers
    rho = arb(fmpq((1 << K) - 1, (1 << B) - 1))
    eta = fmpq(args.eta_numerator, args.eta_denominator)
    target = arb(2) ** arb(eta)
    log_two = arb(2).log()
    rows = []
    all_below = True
    dominant = None
    tail_multiplicity = arb(0)
    for output_weight in range(1, B + 1):
        multiplicity, marginal_weight = dual_shell_upper(
            layered, output_weight, caps
        )
        per_vector = multiplicity / math.comb(B, output_weight)
        ratio = per_vector / rho
        central = (
            args.central_minimum_weight
            <= output_weight
            <= args.central_maximum_weight
        )
        accepted = (ratio < target) if central else None
        if central:
            all_below = all_below and bool(accepted)
        else:
            tail_multiplicity += multiplicity
        log2_ratio = ratio.log() / log_two
        multiplicity_upper = math.nextafter(
            float(multiplicity.upper()), math.inf
        )
        row = {
            "output_weight": output_weight,
            "dual_marginal_source_weight": marginal_weight,
            "multiplicity_upper_binary64_hex": multiplicity_upper.hex(),
            "log2_density_ratio_upper_display": float(log2_ratio.upper()),
            "density_ratio_upper_display": float(ratio.upper()),
            "in_central_weight_interval": central,
            "below_density_target": accepted,
        }
        rows.append(row)
        if central and (
            dominant is None
            or row["log2_density_ratio_upper_display"]
            > dominant["log2_density_ratio_upper_display"]
        ):
            dominant = row

    tail_union = tail_multiplicity * args.tail_union_rows
    tail_comparison = (
        tail_union < arb(fmpq(1, 1 << args.tail_margin_bits))
        if args.tail_union_rows
        else None
    )
    tail_log2_upper = (
        float((tail_union.log() / log_two).upper())
        if args.tail_union_rows and not tail_union.contains(0)
        else None
    )
    complete = all_below and (tail_comparison is not False)

    result = {
        "schema": "bch250-124-rowlocal-fanout-spectrum-arb-v1",
        "status": (
            "OUTWARD_SPECTRUM_CERTIFICATE"
            if complete
            else "FAILED_SPECTRUM_CERTIFICATE"
        ),
        "claim": {
            "fanout_layers": args.layers,
            "source_size": SOURCE_SIZE,
            "target_size": TARGET_SIZE,
            "central_weight_interval": [
                args.central_minimum_weight,
                args.central_maximum_weight,
            ],
            "pointwise_density_excess_bits_strictly_below": (
                f"{args.eta_numerator}/{args.eta_denominator}"
            ),
            "all_central_output_vectors_covered": all_below,
            "per_vector_interpretation_requires_uniform_local_coordinate_permutation": True,
            "tail_union_rows": args.tail_union_rows,
            "tail_expected_count_log2_upper_display": tail_log2_upper,
            "tail_expected_count_strictly_below_2^-bits": (
                args.tail_margin_bits if tail_comparison else None
            ),
            "tail_comparison_accepted": tail_comparison,
        },
        "parameters": {
            "outer_bits": B,
            "outer_dimension": K,
            "exact_nonzero_code_mass": TOTAL_MASS,
            "source_caps": str(args.caps),
        },
        "arithmetic": {
            "backend": "python-flint Arb",
            "precision_bits": args.precision,
            "transition_entries": "exact fmpq",
            "layered_transition": "Arb matrix exponentiation",
            "shell_optimization": (
                "Arb upper evaluation of a valid LP dual threshold. "
                "Nearest floats propose thresholds but prove no ordering."
            ),
        },
        "dominant_shell": dominant,
        "shells": rows,
        "limitations": [
            "This receipt certifies the expected spectrum of one independently sampled row-local wrapper.",
            "An independent uniform local coordinate permutation turns each shell bound into the stated per-vector bound.",
            "It does not prove concentration for one wrapper reused across rows.",
            "Implementation equivalence is outside this receipt.",
        ],
    }
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({key: value for key, value in result.items() if key != "shells"}, indent=2))
    print(f"output={args.output}")
    if not all_below:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
