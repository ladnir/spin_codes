#!/usr/bin/env python3
"""Evaluate a finite all-occupation cover from the certified RM2Sub boxes.

The fixed BCH250 outer map is repeated in every row.  Independent row-local
ParityFanout wrappers induce the supplied expected spectrum envelope.  A
pointwise density comparison transfers the existing random-outer finite
first-moment bound to this expected counting measure.

The 31 RM2Sub witnesses and their Collatz radii come from an outward interval
certificate.  This evaluator uses high-precision point arithmetic for the
remaining finite factors, so its output is diagnostic rather than outward.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path

import mpmath as mp


WORKSTREAM = Path(__file__).resolve().parent
B = 250
K = 125
L = 8448
N = B * L
EPOCHS = N // 128
DISTANCE = 232_320
DEFAULT_SPECTRUM = (
    WORKSTREAM / "bch250_parityfanout31x33_l256_packing_expected_envelope.json"
)
DEFAULT_BOXES = WORKSTREAM / "rm2sub_dense_compact_interval_d11.json"


def mp_fraction(value: Fraction) -> mp.mpf:
    return mp.mpf(value.numerator) / value.denominator


def log_choose(total: int, selected: int) -> mp.mpf:
    return (
        mp.loggamma(total + 1)
        - mp.loggamma(selected + 1)
        - mp.loggamma(total - selected + 1)
    )


def log_binomial_mass(total: int, selected: int, probability: mp.mpf) -> mp.mpf:
    if selected == 0:
        return total * mp.log1p(-probability)
    if selected == total:
        return total * mp.log(probability)
    return (
        log_choose(total, selected)
        + selected * mp.log(probability)
        + (total - selected) * mp.log1p(-probability)
    )


def density_excess_bits(
    spectrum_path: Path,
    dimension: int = K,
    minimum_weight: int = 1,
    maximum_weight: int = B,
) -> tuple[mp.mpf, int]:
    payload = json.loads(spectrum_path.read_text(encoding="utf-8"))
    random_density = (
        mp.log(mp.power(2, dimension) - 1) - mp.log(mp.power(2, B) - 1)
    )
    best = mp.ninf
    best_weight = -1
    for row in payload["spectrum"]:
        weight = int(row["weight"])
        if weight == 0 or not minimum_weight <= weight <= maximum_weight:
            continue
        value = row.get("log2_expected_multiplicity")
        if value is None:
            continue
        log_shell_count = log_choose(B, weight)
        log_per_vector = mp.mpf(str(value)) * mp.log(2) - log_shell_count
        excess = (log_per_vector - random_density) / mp.log(2)
        if excess > best:
            best = excess
            best_weight = weight
    return best, best_weight


def load_boxes(path: Path) -> list[dict[str, object]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload["status"] != "proved":
        raise ValueError("the compact-density receipt is not proved")
    boxes = payload["boxes"]
    if not boxes:
        raise ValueError("the compact-density receipt has no boxes")
    return boxes


def containing_box(boxes: list[dict[str, object]], alpha: Fraction) -> dict[str, object]:
    for box in boxes:
        lower = Fraction(str(box["lower"]))
        upper = Fraction(str(box["upper"]))
        if lower <= alpha <= upper:
            return box
    raise ValueError(f"no certified box contains alpha={alpha}")


def evaluate_row(
    occupation: int,
    box: dict[str, object],
    density_excess: mp.mpf,
) -> dict[str, object]:
    p = mp_fraction(Fraction(str(box["candidate_probability"])))
    z = mp_fraction(Fraction(str(box["z"])))
    radius = mp_fraction(Fraction(str(box["radius_upper"])))
    vector = [mp_fraction(Fraction(str(value))) for value in box["collatz_vector"]]
    matrix_prefactor = vector[0] * max(1 / value for value in vector)

    log_bound = log_choose(L, occupation)
    log_bound += occupation * mp.log(mp.power(2, K) - 1)
    log_bound -= occupation * mp.log1p(-mp.power(2, -B))
    log_bound -= B * log_binomial_mass(L, occupation, p)
    log_bound -= DISTANCE * mp.log(z)
    log_bound += EPOCHS * mp.log(radius)
    log_bound += mp.log(matrix_prefactor)
    log_bound += occupation * density_excess * mp.log(2)
    return {
        "occupation": occupation,
        "alpha": occupation / L,
        "log2_upper": float(log_bound / mp.log(2)),
        "margin_bits": float(-log_bound / mp.log(2)),
        "box_lower": str(box["lower"]),
        "box_upper": str(box["upper"]),
        "candidate_probability": str(box["candidate_probability"]),
        "z": str(box["z"]),
        "matrix_prefactor_log2": float(mp.log(matrix_prefactor) / mp.log(2)),
    }


def logaddexp(left: mp.mpf, right: mp.mpf) -> mp.mpf:
    if left == mp.ninf:
        return right
    maximum = max(left, right)
    return maximum + mp.log(mp.exp(left - maximum) + mp.exp(right - maximum))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectrum", type=Path, default=DEFAULT_SPECTRUM)
    parser.add_argument("--boxes", type=Path, default=DEFAULT_BOXES)
    parser.add_argument("--minimum-occupation", type=int, default=32)
    parser.add_argument("--minimum-output-weight", type=int, default=1)
    parser.add_argument("--maximum-output-weight", type=int, default=B)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not 1 <= args.minimum_occupation <= L:
        raise ValueError(f"minimum occupation must lie in [1,{L}]")

    mp.mp.dps = 100
    boxes = load_boxes(args.boxes)
    if not 1 <= args.minimum_output_weight <= args.maximum_output_weight <= B:
        raise ValueError("invalid output-weight interval")
    excess, maximizing_weight = density_excess_bits(
        args.spectrum,
        minimum_weight=args.minimum_output_weight,
        maximum_weight=args.maximum_output_weight,
    )
    rows = []
    union_log = mp.ninf
    for occupation in range(args.minimum_occupation, L + 1):
        box = containing_box(boxes, Fraction(occupation, L))
        row = evaluate_row(occupation, box, excess)
        rows.append(row)
        union_log = logaddexp(
            union_log, mp.mpf(str(row["log2_upper"])) * mp.log(2)
        )
    worst = min(rows, key=lambda row: float(row["margin_bits"]))
    payload = {
        "schema": "bch250-rowlocal-fanout-finite-cover-v1",
        "status": "HIGH_PRECISION_FINITE_DIAGNOSTIC",
        "parameters": {
            "outer_bits": B,
            "outer_dimension": K,
            "outer_rows": L,
            "output_bits": N,
            "distance": DISTANCE,
            "relative_distance": DISTANCE / N,
            "inner_epochs": EPOCHS,
            "minimum_occupation": args.minimum_occupation,
            "maximum_occupation": L,
            "output_weight_interval": [
                args.minimum_output_weight,
                args.maximum_output_weight,
            ],
            "spectrum": str(args.spectrum),
            "compact_interval_receipt": str(args.boxes),
        },
        "pointwise_comparison": {
            "reference": "uniform random [250,125] injection counting measure",
            "maximum_log2_density_ratio": float(excess),
            "maximizing_weight": maximizing_weight,
        },
        "union_log2_upper": float(union_log / mp.log(2)),
        "union_margin_bits": float(-union_log / mp.log(2)),
        "worst_occupation": worst,
        "rows": rows,
        "limitations": [
            "The compact-density radius bounds are outward certified.",
            "The finite logarithms and spectrum-density comparison use 100-digit point arithmetic, not directed rounding.",
            "The result covers only the recorded occupation interval; occupations below it require separate receipts.",
            "The probability space samples the row-local fanout wrapper independently for every row.",
        ],
    }
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"density_excess_bits={float(excess):.12f}")
    print(f"worst_occupation={worst['occupation']}")
    print(f"worst_margin_bits={worst['margin_bits']:.9f}")
    print(f"union_margin_bits={payload['union_margin_bits']:.9f}")
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
