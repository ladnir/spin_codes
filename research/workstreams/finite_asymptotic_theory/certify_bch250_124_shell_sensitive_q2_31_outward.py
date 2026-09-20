#!/usr/bin/env python3
"""Outward shell-sensitive certificate for occupations 2 through 31."""

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
MINIMUM_Q = 2
MAXIMUM_Q = 31
DEFAULT_SPECTRUM = (
    WORKSTREAM / "bch250_124_parityfanout31x33_l56_cutoff12_spectrum_outward.json"
)
DEFAULT_WITNESSES = (
    WORKSTREAM
    / "bch250_124_parityfanout31x33_l256_packing_expected_q2_31_d11_diagnostic.json"
)


base.B = B
base.L = L
base.DISTANCE = DISTANCE
base.EPOCHS_PER_REGION = L // base.STEP_BITS


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


def load_witnesses(path: Path) -> dict[int, tuple[float, float]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    result: dict[int, tuple[float, float]] = {}
    for row in payload["occupation_rows"]:
        occupation = int(row["active_regular_outer_blocks"])
        if not MINIMUM_Q <= occupation <= MAXIMUM_Q:
            continue
        p = float(row["best_candidate_probability"])
        z = math.exp(-math.exp(float(row["best_log_surprisal"])))
        if not 0.0 < p < 1.0 or not 0.0 < z < 1.0:
            raise ArithmeticError("invalid diagnostic witness")
        result[occupation] = (p, z)
    if set(result) != set(range(MINIMUM_Q, MAXIMUM_Q + 1)):
        raise ValueError("witness receipt omits an occupation")
    return result


def bernoulli_majorant_upper(
    spectrum: dict[int, float], p: float
) -> tuple[float, int]:
    q_lower = base.down(1.0 - p)
    best = 0.0
    best_weight = -1
    for weight, multiplicity in spectrum.items():
        denominator = base.mul_down(
            base.integer_lower(math.comb(B, weight)),
            base.mul_down(
                base.power_down(p, weight),
                base.power_down(q_lower, B - weight),
            ),
        )
        candidate = base.div_up(multiplicity, denominator)
        if candidate > best:
            best = candidate
            best_weight = weight
    return best, best_weight


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
            / "bch250_124_parityfanout31x33_l56_q2_31_outward_d11.json"
        ),
    )
    args = parser.parse_args()

    spectrum = load_spectrum(args.spectrum)
    witnesses = load_witnesses(args.witnesses)
    activation, live = base.load_activation_and_live()
    majorants: dict[str, tuple[float, int]] = {}
    moments: dict[str, list[base.Scaled]] = {}
    for p, z in sorted(set(witnesses.values())):
        majorants.setdefault(p.hex(), bernoulli_majorant_upper(spectrum, p))
        impulses = base.impulse_matrices_upper(z, activation, live)
        candidates = base.candidate_epoch_matrices(impulses, p)
        regions = base.region_matrices_upper(candidates, MAXIMUM_Q)
        moments[f"{p.hex()}|{z.hex()}"] = base.total_moments_upper(regions)

    aggregate = base.SCALED_ZERO
    rows = []
    for occupation in range(MINIMUM_Q, MAXIMUM_Q + 1):
        p, z = witnesses[occupation]
        majorant, maximizing_weight = majorants[p.hex()]
        moment = moments[f"{p.hex()}|{z.hex()}"][occupation]
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
            base.scaled_from_float(base.integer_upper(math.comb(L, occupation))),
            base.scaled_multiply_up(
                base.scaled_power_up(
                    base.scaled_from_float(majorant), occupation
                ),
                inner,
            ),
        )
        aggregate = base.scaled_add_up(aggregate, contribution)
        rows.append(
            {
                "occupation": occupation,
                "candidate_probability_exact_binary64": p.hex(),
                "z_exact_binary64": z.hex(),
                "majorant_upper_hex": majorant.hex(),
                "majorant_maximizing_weight": maximizing_weight,
                "contribution_upper_scaled": base.scaled_hex(contribution),
                "contribution_log2_upper_display": scaled_log2_upper(
                    contribution
                ),
            }
        )

    union_log2 = scaled_log2_upper(aggregate)
    comparison_40 = base.scaled_at_most_power_of_two(aggregate, 40)
    comparison_52 = base.scaled_at_most_power_of_two(aggregate, 52)
    result = {
        "schema": "bch250-124-shell-sensitive-q2-31-outward-v1",
        "status": (
            "OUTWARD_CERTIFICATE" if comparison_40 else "FAILED_OUTWARD_CERTIFICATE"
        ),
        "claim": {
            "occupations": [MINIMUM_Q, MAXIMUM_Q],
            "bad_weight_at_most": DISTANCE,
            "expected_bad_count_upper_scaled": base.scaled_hex(aggregate),
            "union_log2_upper_display": union_log2,
            "union_margin_bits_display": -union_log2,
            "comparison_to_2^-40": comparison_40,
            "comparison_to_2^-52": comparison_52,
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
        "arithmetic": {
            "format": "IEEE-754 binary64 with one-ULP outward positive operations",
            "spectrum_shell_uppers": "outward Arb endpoints advanced one binary64 ULP",
            "transcendentals": "select exact binary64 witnesses and assert no inequality",
        },
        "occupation_rows": rows,
        "limitations": [
            "Occupations 1 and 32 through 8576 use separate receipts.",
            "Implementation equivalence is outside this receipt.",
        ],
    }
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {key: value for key, value in result.items() if key != "occupation_rows"},
            indent=2,
        )
    )
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
