#!/usr/bin/env python3
"""Outward shell-sensitive certificate for one active BCH250 row."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys


WORKSTREAM = Path(__file__).resolve().parent
sys.path.insert(0, str(WORKSTREAM))

import certify_golay_ba_rm2sub_finite_q2_64 as base  # noqa: E402


B = 250
K = 124
L = 8576
N = B * L
DISTANCE = 235_840
DEFAULT_SPECTRUM = (
    WORKSTREAM / "bch250_124_parityfanout31x33_l56_cutoff12_spectrum_outward.json"
)
DEFAULT_WITNESSES = (
    WORKSTREAM
    / "bch250_124_parityfanout31x33_l56_shell_sensitive_q1_d11_diagnostic.json"
)


base.B = B
base.L = L
base.DISTANCE = DISTANCE
base.EPOCHS_PER_REGION = L // base.STEP_BITS


def zero_scaled_matrix() -> base.ScaledMatrix:
    return (
        base.SCALED_ZERO,
        base.SCALED_ZERO,
        base.SCALED_ZERO,
        base.SCALED_ZERO,
    )


def matrix_add(
    left: base.ScaledMatrix, right: base.ScaledMatrix
) -> base.ScaledMatrix:
    return tuple(
        base.scaled_add_up(a, b) for a, b in zip(left, right, strict=True)
    )  # type: ignore[return-value]


def load_spectrum(path: Path) -> dict[int, float]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload["status"] != "OUTWARD_SPECTRUM_CERTIFICATE":
        raise ValueError("spectrum receipt is not outward certified")
    result = {
        int(row["output_weight"]): float.fromhex(
            row["multiplicity_upper_binary64_hex"]
        )
        for row in payload["shells"]
    }
    if set(result) != set(range(1, B + 1)):
        raise ValueError("spectrum receipt omits a nonzero shell")
    return result


def load_witnesses(path: Path) -> dict[int, float]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    result: dict[int, float] = {}
    for row in payload["weight_rows"]:
        weight = int(row["outer_weight"])
        log_surprisal = float(row["best_log_surprisal"])
        z = math.exp(-math.exp(log_surprisal))
        if not 0.0 < z < 1.0:
            raise ArithmeticError("invalid diagnostic witness")
        result[weight] = z
    if set(result) != set(range(1, B + 1)):
        raise ValueError("witness receipt omits a nonzero shell")
    return result


def weight_moments_for_z(
    z: float,
    activation: list[tuple[int, int]],
    live: list[tuple[int, int]],
) -> list[base.Scaled]:
    # With one active outer row, an active transposed region contains exactly
    # one input bit.  Passing the impulse matrices directly is the p=1 limit
    # of the candidate-position recurrence.
    impulses = base.impulse_matrices_upper(z, activation, live)
    regions = base.region_matrices_upper(impulses, 1)
    inactive = base.scaled_matrix_from_float(regions[0])
    active = base.scaled_matrix_from_float(regions[1])

    coefficients = [zero_scaled_matrix() for _ in range(B + 1)]
    coefficients[0] = base.SCALED_IDENTITY
    for completed in range(B):
        updated = [zero_scaled_matrix() for _ in range(B + 1)]
        for weight in range(completed + 1):
            unchanged = base.scaled_matrix_multiply_up(
                coefficients[weight], inactive
            )
            updated[weight] = matrix_add(updated[weight], unchanged)
            changed = base.scaled_matrix_multiply_up(
                coefficients[weight], active
            )
            updated[weight + 1] = matrix_add(updated[weight + 1], changed)
        coefficients = updated

    moments = [base.SCALED_ZERO] * (B + 1)
    for weight in range(B + 1):
        terminal = base.scaled_add_up(
            coefficients[weight][0], coefficients[weight][1]
        )
        moments[weight] = base.scaled_multiply_up(
            terminal,
            base.scaled_from_float(
                base.ratio_upper(1, math.comb(B, weight))
            ),
        )
    return moments


def scaled_log2_upper(value: base.Scaled) -> float:
    if value[0] == 0.0:
        return -math.inf
    return math.nextafter(math.log2(value[0]) + value[1], math.inf)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectrum", type=Path, default=DEFAULT_SPECTRUM)
    parser.add_argument("--witnesses", type=Path, default=DEFAULT_WITNESSES)
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            WORKSTREAM
            / "bch250_124_parityfanout31x33_l56_q1_outward_d11.json"
        ),
    )
    args = parser.parse_args()

    spectrum = load_spectrum(args.spectrum)
    witnesses = load_witnesses(args.witnesses)
    activation, live = base.load_activation_and_live()
    moment_cache = {
        z.hex(): weight_moments_for_z(z, activation, live)
        for z in sorted(set(witnesses.values()))
    }

    aggregate = base.SCALED_ZERO
    rows = []
    for weight in range(1, B + 1):
        z = witnesses[weight]
        moment = moment_cache[z.hex()][weight]
        chernoff = base.scaled_multiply_up(
            moment,
            base.scaled_power_up(
                base.scaled_from_float(base.div_up(1.0, z)), DISTANCE
            ),
        )
        inner = (
            chernoff
            if base.scaled_at_most_power_of_two(chernoff, 0)
            else base.SCALED_ONE
        )
        contribution = base.scaled_multiply_up(
            base.scaled_from_float(base.integer_upper(L)),
            base.scaled_multiply_up(
                base.scaled_from_float(spectrum[weight]), inner
            ),
        )
        aggregate = base.scaled_add_up(aggregate, contribution)
        rows.append(
            {
                "output_weight": weight,
                "z_exact_binary64": z.hex(),
                "inner_probability_upper_scaled": base.scaled_hex(inner),
                "contribution_upper_scaled": base.scaled_hex(contribution),
                "contribution_log2_upper_display": scaled_log2_upper(
                    contribution
                ),
            }
        )

    union_log2 = scaled_log2_upper(aggregate)
    comparison_40 = base.scaled_at_most_power_of_two(aggregate, 40)
    comparison_53 = base.scaled_at_most_power_of_two(aggregate, 53)
    result = {
        "schema": "bch250-124-shell-sensitive-q1-outward-v1",
        "status": (
            "OUTWARD_CERTIFICATE" if comparison_40 else "FAILED_OUTWARD_CERTIFICATE"
        ),
        "claim": {
            "occupation": 1,
            "bad_weight_at_most": DISTANCE,
            "expected_bad_count_upper_scaled": base.scaled_hex(aggregate),
            "union_log2_upper_display": union_log2,
            "union_margin_bits_display": -union_log2,
            "comparison_to_2^-40": comparison_40,
            "comparison_to_2^-53": comparison_53,
        },
        "parameters": {
            "outer_bits": B,
            "outer_dimension": K,
            "outer_rows": L,
            "output_bits": N,
            "distance": DISTANCE,
            "fanout_layers": 56,
            "spectrum_receipt": str(args.spectrum),
            "witness_receipt": str(args.witnesses),
        },
        "method": (
            "For each exact output weight, a matrix-coefficient recurrence "
            "averages the ordered RM2Sub region product over the uniform "
            "local coordinate permutation."
        ),
        "arithmetic": {
            "format": "IEEE-754 binary64 with one-ULP outward positive operations",
            "spectrum_shell_uppers": "outward Arb endpoints advanced one binary64 ULP",
            "transcendentals": "select exact binary64 witnesses and assert no inequality",
        },
        "weight_rows": rows,
        "limitations": [
            "Occupations 2 through 8576 use separate receipts.",
            "Implementation equivalence is outside this receipt.",
        ],
    }
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {key: value for key, value in result.items() if key != "weight_rows"},
            indent=2,
        )
    )
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
