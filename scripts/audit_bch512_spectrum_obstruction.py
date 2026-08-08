#!/usr/bin/env python3
"""Audit what generic finite bounds can prove for the BCH512 spectrum.

This is deliberately not a floating-point Delsarte certificate.  It verifies
the exact BCH zero-set facts, derives the parent dual-distance floor, and then
compares two exact generic coefficient bounds with the corrected even-weight
random-like projection targets (including the +40-bit sensitivity window):

* the strength-15 orthogonal-array Christoffel bound; and
* the constant-weight packing bound coming from distance at least 62.

The resulting gap is an obstruction report, not a claim that the desired BCH
spectrum envelope is false.  It quantifies how much BCH-specific cancellation
is still required after the standard code-parameter information is exhausted.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DEFAULT_ARTIFACT = ROOT / "bch512_spectrum_obstruction.json"
N_PRIMITIVE = 511
N_EXTENDED = 512
PARENT_DIMENSION = 259
SUBCODE_DIMENSION = 256
DESIGNED_DISTANCE = 61
DUAL_DESIGNED_DISTANCE = 16
PROJECTION_INFLATION_BITS = 40
PROJECTION_INFLATION_MAX_WEIGHT = 94
AUDIT_WEIGHTS = (62, 94, 118)


def canonical_sha256(payload: dict[str, object]) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def log2_fraction(value: Fraction) -> float:
    if value <= 0:
        raise ValueError("logarithm requires a positive fraction")
    return math.log2(value.numerator) - math.log2(value.denominator)


def primitive_zero_set() -> set[int]:
    zeros: set[int] = set()
    for exponent in range(1, DESIGNED_DISTANCE):
        value = exponent
        while value not in zeros:
            zeros.add(value)
            value = (2 * value) % N_PRIMITIVE
    return zeros


def dual_zero_set(primal_zeros: set[int]) -> set[int]:
    return {
        (-exponent) % N_PRIMITIVE
        for exponent in range(N_PRIMITIVE)
        if exponent not in primal_zeros
    }


def initial_zero_run(zeros: set[int]) -> int:
    value = 0
    while value in zeros:
        value += 1
    return value


def krawtchouk_values(n: int, x: int, degree: int) -> list[int]:
    values = [1]
    if degree == 0:
        return values
    values.append(n - 2 * x)
    for j in range(1, degree):
        numerator = (n - 2 * x) * values[j] - (n - j + 1) * values[j - 1]
        if numerator % (j + 1):
            raise AssertionError("Krawtchouk recurrence lost integrality")
        values.append(numerator // (j + 1))
    return values


def christoffel_bound(weight: int) -> Fraction:
    # Dual distance >=16 makes the code an orthogonal array of strength 15.
    # Squaring a degree-7 interpolant uses only moments through degree 14.
    kernel = sum(
        (
            Fraction(value * value, math.comb(N_EXTENDED, degree))
            for degree, value in enumerate(krawtchouk_values(N_EXTENDED, weight, 7))
        ),
        Fraction(0),
    )
    return Fraction(1 << PARENT_DIMENSION, 1) / kernel


def constant_weight_packing_bound(weight: int) -> Fraction:
    # Two weight-h supports intersect in at most h-31 coordinates, since their
    # symmetric difference is a nonzero codeword of weight at least 62.  Hence
    # an (h-30)-subset occurs in at most one such support.
    subset_size = weight - (DESIGNED_DISTANCE - 1) // 2
    return Fraction(
        math.comb(N_EXTENDED, subset_size), math.comb(weight, subset_size)
    )


def projected_subcode_baseline(weight: int) -> Fraction:
    # The candidate is even, so the ambient space has dimension n-1.
    return Fraction(math.comb(N_EXTENDED, weight), 1 << (N_EXTENDED - 1 - SUBCODE_DIMENSION))


def coefficient_row(weight: int) -> dict[str, object]:
    baseline = projected_subcode_baseline(weight)
    inflation_bits = (
        PROJECTION_INFLATION_BITS if weight <= PROJECTION_INFLATION_MAX_WEIGHT else 0
    )
    target = baseline * (1 << inflation_bits)
    christoffel = christoffel_bound(weight)
    packing = constant_weight_packing_bound(weight)
    generic = min(christoffel, packing)
    return {
        "weight": weight,
        "projection_inflation_bits": inflation_bits,
        "random_like_subcode_log2": log2_fraction(baseline),
        "projection_target_log2": log2_fraction(target),
        "christoffel_upper_log2": log2_fraction(christoffel),
        "constant_weight_upper_log2": log2_fraction(packing),
        "best_generic_upper_log2": log2_fraction(generic),
        "generic_gap_above_target_bits": log2_fraction(generic / target),
        "best_generic_bound": (
            "christoffel" if christoffel <= packing else "constant_weight_packing"
        ),
    }


def boundary_interpolation_row() -> dict[str, object]:
    # For primitive weight 61, normalization by the unique 61st-power scaling
    # reduces to monic degree-30 interpolation.  Counting 30-subsets gives the
    # standard rigorous list-size bound used in the exploration note.
    weight = 61
    random_like = Fraction(math.comb(N_PRIMITIVE, weight), 1 << (N_PRIMITIVE - PARENT_DIMENSION))
    interpolation = Fraction(
        N_PRIMITIVE * math.comb(N_PRIMITIVE, 30), math.comb(weight, 30)
    )
    random_log = log2_fraction(random_like)
    target_log = random_log + 53.75
    interpolation_log = log2_fraction(interpolation)
    return {
        "primitive_weight": weight,
        "random_like_parent_log2": random_log,
        "legacy_target_log2": target_log,
        "generic_interpolation_upper_log2": interpolation_log,
        "gap_above_legacy_target_bits": interpolation_log - target_log,
        "note": "The legacy 53.75-bit allowance is compared in the logarithmic domain.",
    }


def build_payload() -> dict[str, object]:
    primal_zeros = primitive_zero_set()
    dual_zeros = dual_zero_set(primal_zeros)
    dual_run = initial_zero_run(dual_zeros)
    if len(primal_zeros) != N_PRIMITIVE - PARENT_DIMENSION:
        raise SystemExit("unexpected BCH generator degree")
    if dual_run != 15:
        raise SystemExit(f"unexpected consecutive dual-zero run: {dual_run}")

    rows = [coefficient_row(weight) for weight in AUDIT_WEIGHTS]
    if min(float(row["generic_gap_above_target_bits"]) for row in rows) < 50.0:
        raise SystemExit("generic obstruction gap unexpectedly fell below 50 bits")

    payload: dict[str, object] = {
        "schema": 1,
        "status": "NOTEWORTHY_OBSTRUCTION",
        "interpretation": (
            "No exact spectrum was located or computed.  Parameter-only finite bounds "
            "remain more than 50 bits above the corrected projection targets, "
            "so a successful certificate must use BCH-specific cancellation or a real spectrum."
        ),
        "parent": {
            "primitive_parameters": [N_PRIMITIVE, PARENT_DIMENSION, DESIGNED_DISTANCE],
            "extended_parameters_floor": [N_EXTENDED, PARENT_DIMENSION, DESIGNED_DISTANCE + 1],
            "generator_degree": len(primal_zeros),
            "dual_initial_zero_run": dual_run,
            "extended_dual_distance_floor": DUAL_DESIGNED_DISTANCE,
            "extended_dual_distance_reason": (
                "B^perp has cyclic zeros 0..14, hence distance >=16.  For extension-dual "
                "words (u+1,1), u+1 retains zeros 1..14, hence weight >=15 before "
                "the extension coordinate."
            ),
        },
        "projection": {
            "subcode_dimension": SUBCODE_DIMENSION,
            "even_weight_random_like_exponent": "2^(k-(n-1))*binom(n,w)",
            "audited_uniform_low_weight_inflation_bits": PROJECTION_INFLATION_BITS,
            "inflation_window_max_weight": PROJECTION_INFLATION_MAX_WEIGHT,
            "weights": rows,
        },
        "primitive_boundary": boundary_interpolation_row(),
        "literature": {
            "minimum_distance_status_source": "https://doi.org/10.11591/ijece.v9i2.pp1232-1239",
            "source_finding": (
                "The 2019 ZSSMP study reports BCH(511,259) as one of the two "
                "length-511 cases whose true minimum distance it did not determine."
            ),
            "finite_spectrum_bound_source": "https://doi.org/10.1023/A:1011220817609",
        },
    }
    payload["sha256"] = canonical_sha256(payload)
    return payload


def write_artifact(path: Path, payload: dict[str, object]) -> None:
    path.write_text(json.dumps(payload, indent=2) + "\n")


def verify_artifact(path: Path, regenerated: dict[str, object]) -> None:
    stored = json.loads(path.read_text())
    claimed = stored.pop("sha256", None)
    actual = canonical_sha256(stored)
    if claimed != actual:
        raise SystemExit("BCH spectrum obstruction artifact SHA-256 mismatch")
    stored["sha256"] = claimed
    if stored != regenerated:
        raise SystemExit("BCH spectrum obstruction artifact differs from regenerated audit")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_ARTIFACT)
    parser.add_argument("--write-artifact", action="store_true")
    args = parser.parse_args()

    payload = build_payload()
    if args.write_artifact:
        write_artifact(args.artifact, payload)
        print(f"bch512_spectrum_obstruction_artifact_written,{args.artifact}")
    else:
        verify_artifact(args.artifact, payload)
        print("bch512_spectrum_obstruction_status,PASS")
    for row in payload["projection"]["weights"]:
        print(
            "weight_gap_bits,"
            f"{row['weight']},{row['generic_gap_above_target_bits']:.6f},"
            f"{row['best_generic_bound']}"
        )
    print(f"artifact_sha256,{payload['sha256']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
