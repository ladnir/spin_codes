#!/usr/bin/env python3
"""Diagnose occupations 2 through 31 with the full wrapped-row spectrum.

For a row weight ``w``, the input JSON supplies an upper bound ``A_w`` on
the expected number of wrapped BCH words of that weight.  For each Bernoulli
value probability ``p``, this program uses the valid one-row majorant

    M(p) = max_w A_w / (C(250,w) p^w (1-p)^(250-w)).

It then combines ``M(p)^Q`` with the finite RM2Sub-S19 transfer for exactly
``Q`` candidate positions in every transposed region.  The spectrum input
and witness search use nearest binary64 arithmetic, so the result is a
diagnostic rather than an outward certificate.
"""

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
L = 8576
N = B * L
DISTANCE = 235_840
MINIMUM_Q = 2
MAXIMUM_Q = 31
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
    result: dict[int, float] = {}
    for row in payload["spectrum"]:
        weight = int(row["weight"])
        if weight == 0:
            continue
        value = row.get("expected_multiplicity")
        if value is None:
            log_value = row.get("log2_expected_multiplicity")
            if log_value is None:
                continue
            value = math.exp2(float(log_value))
        multiplicity = float(value)
        if multiplicity > 0.0:
            result[weight] = multiplicity
    if not result:
        raise ValueError("spectrum has no positive nonzero shell")
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
    required = set(range(MINIMUM_Q, MAXIMUM_Q + 1))
    if set(result) != required:
        raise ValueError("witness receipt omits an occupation")
    return result


def bernoulli_majorant(
    spectrum: dict[int, float], probability: float
) -> tuple[float, int]:
    complement = 1.0 - probability
    best = -math.inf
    best_weight = -1
    for weight, multiplicity in spectrum.items():
        log_candidate = math.log(multiplicity)
        log_candidate -= math.lgamma(B + 1)
        log_candidate += math.lgamma(weight + 1)
        log_candidate += math.lgamma(B - weight + 1)
        log_candidate -= weight * math.log(probability)
        log_candidate -= (B - weight) * math.log(complement)
        if log_candidate > best:
            best = log_candidate
            best_weight = weight
    return math.exp(best), best_weight


def scaled_log2(value: base.Scaled) -> float:
    if value[0] == 0.0:
        return -math.inf
    return math.log2(value[0]) + value[1]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectrum", type=Path, required=True)
    parser.add_argument("--witnesses", type=Path, default=DEFAULT_WITNESSES)
    parser.add_argument("--fanout-layers", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.fanout_layers < 0:
        parser.error("--fanout-layers must be nonnegative")

    spectrum = load_spectrum(args.spectrum)
    witnesses = load_witnesses(args.witnesses)
    activation, live = base.load_activation_and_live()

    majorants: dict[str, tuple[float, int]] = {}
    moments: dict[str, list[base.Scaled]] = {}
    for p, z in sorted(set(witnesses.values())):
        majorants.setdefault(p.hex(), bernoulli_majorant(spectrum, p))
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
                "candidate_probability": p,
                "z": z,
                "majorant_log2": math.log2(majorant),
                "majorant_maximizing_weight": maximizing_weight,
                "contribution_log2": scaled_log2(contribution),
                "contribution_margin_bits": -scaled_log2(contribution),
            }
        )

    union_log2 = scaled_log2(aggregate)
    result = {
        "schema": "bch250-124-shell-sensitive-q2-31-diagnostic-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "parameters": {
            "outer_bits": B,
            "outer_dimension": 124,
            "outer_rows": L,
            "output_bits": N,
            "distance": DISTANCE,
            "fanout_layers": args.fanout_layers,
            "spectrum": str(args.spectrum),
            "witnesses": str(args.witnesses),
        },
        "union_log2_upper_diagnostic": union_log2,
        "union_margin_bits_diagnostic": -union_log2,
        "comparison_to_2^-40_diagnostic": union_log2 < -40.0,
        "dominant_row": max(rows, key=lambda row: row["contribution_log2"]),
        "occupation_rows": rows,
        "limitations": [
            "The spectrum envelope was produced with nearest binary64 arithmetic.",
            "The stored witnesses were optimized for the earlier uniform-reference proof.",
            "The RM2Sub positive arithmetic is outward rounded, but the complete result is only diagnostic.",
        ],
    }
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"fanout_layers={args.fanout_layers}")
    print(f"union_margin_bits={-union_log2:.9f}")
    print(
        "dominant_occupation="
        f"{result['dominant_row']['occupation']},"
        "margin_bits="
        f"{result['dominant_row']['contribution_margin_bits']:.9f}"
    )
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
