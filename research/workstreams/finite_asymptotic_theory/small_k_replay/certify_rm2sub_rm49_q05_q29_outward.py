#!/usr/bin/env python3
"""Certify every three-band composition at occupations five through 29."""

from __future__ import annotations

import hashlib
import json
import math
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

from flint import arb, ctx

import certify_rm2sub_rm49_q1_outward as q1
import certify_rm2sub_rm49_q30_low_outward as q30


HERE = Path(__file__).resolve().parent
MANIFEST = HERE / "RM2SUB_RM49_Q05_Q29_INPUT_MANIFEST.json"
OUTPUT = HERE / "rm2sub_rm49_t64_s14_q05_q29_outward.json"
MIN_OCCUPATION = 5
MAX_OCCUPATION = 29
GROUPS = (
    ("low", 1, 47, Fraction(1, 4)),
    ("central", 48, 416, Fraction(1, 2)),
    ("high", 417, 512, Fraction(3, 4)),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_manifest() -> tuple[dict[str, object], dict[str, Path]]:
    payload = json.loads(MANIFEST.read_text(encoding="utf-8"))
    paths: dict[str, Path] = {}
    for row in payload["files"]:
        path = (HERE / row["path"]).resolve()
        if sha256(path) != row["sha256"]:
            raise AssertionError(f"hash mismatch for {row['role']}")
        paths[str(row["role"])] = path
    return payload, paths


def compositions(total: int):
    for low in range(total + 1):
        for central in range(total - low + 1):
            yield low, central, total - low - central


def all_compositions():
    for occupation in range(MIN_OCCUPATION, MAX_OCCUPATION + 1):
        yield from compositions(occupation)


def load_witnesses(
    path: Path,
) -> tuple[
    dict[tuple[int, int, int], int],
    dict[tuple[int, int, int], float],
    dict[int, float],
]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    observed_groups = tuple(
        (
            str(row["name"]),
            int(row["lower"]),
            int(row["upper"]),
            Fraction(str(row["reference_probability"])),
        )
        for row in payload["groups"]
    )
    if observed_groups != GROUPS:
        raise AssertionError("three-band reference partition changed")
    if payload["parameters"]["change_of_measure"] != "pointwise-density-envelope":
        raise AssertionError("witness receipt does not use the pointwise envelope")

    expected = set(all_compositions())
    witnesses: dict[tuple[int, int, int], int] = {}
    binary64: dict[tuple[int, int, int], float] = {}
    for row in payload["composition_rows"]:
        occupation = int(row["occupation"])
        if not MIN_OCCUPATION <= occupation <= MAX_OCCUPATION:
            continue
        composition = (
            int(row["low_rows"]),
            int(row["central_rows"]),
            int(row["high_rows"]),
        )
        if sum(composition) != occupation:
            raise AssertionError("composition has the wrong occupation")
        raw = float(row["log_surprisal"])
        witness = round(10.0 * raw)
        if abs(raw - witness / 10.0) >= 1e-12:
            raise AssertionError("composition witness is not on the exact tenth grid")
        witnesses[composition] = witness
        binary64[composition] = float(row["log2_upper_diagnostic"])
    if set(witnesses) != expected or len(witnesses) != 4925:
        raise AssertionError("witness table does not exactly cover Q=5,...,29")

    binary64_unions = {
        int(row["occupation"]): float(row["log2_upper_diagnostic"])
        for row in payload["occupation_rows"]
        if MIN_OCCUPATION <= int(row["occupation"]) <= MAX_OCCUPATION
    }
    if set(binary64_unions) != set(range(MIN_OCCUPATION, MAX_OCCUPATION + 1)):
        raise AssertionError("binary64 occupation table is incomplete")
    return witnesses, binary64, binary64_unions


def band_envelopes(
    outer_spectrum: dict[int, int],
) -> tuple[list[Fraction], list[int]]:
    envelopes: list[Fraction] = []
    maximizing_weights: list[int] = []
    for _, lower, upper, probability in GROUPS:
        best: Fraction | None = None
        best_weight = -1
        for weight, count in outer_spectrum.items():
            if not lower <= weight <= upper:
                continue
            reference = (
                Fraction(math.comb(q1.OUTER_BITS, weight), 1)
                * probability**weight
                * (1 - probability) ** (q1.OUTER_BITS - weight)
            )
            candidate = Fraction(count, 1) / reference
            if best is None or candidate > best:
                best = candidate
                best_weight = weight
        if best is None:
            raise AssertionError("empty spectrum band")
        envelopes.append(best)
        maximizing_weights.append(best_weight)
    if maximizing_weights != [32, 48, 512]:
        raise AssertionError(f"a density maximum changed: {maximizing_weights}")
    return envelopes, maximizing_weights


def convolve(left: list[Fraction], right: list[Fraction]) -> list[Fraction]:
    result = [Fraction(0) for _ in range(len(left) + len(right) - 1)]
    for first, first_value in enumerate(left):
        for second, second_value in enumerate(right):
            result[first + second] += first_value * second_value
    return result


def binomial_distribution(count: int, probability: Fraction) -> list[Fraction]:
    return [
        Fraction(math.comb(count, live), 1)
        * probability**live
        * (1 - probability) ** (count - live)
        for live in range(count + 1)
    ]


def reference_region(
    regions: list[q1.Matrix], composition: tuple[int, int, int]
) -> q1.Matrix:
    distribution = [Fraction(1)]
    for count, (_, _, _, probability) in zip(composition, GROUPS, strict=True):
        distribution = convolve(distribution, binomial_distribution(count, probability))
    if sum(distribution) != 1 or len(distribution) != sum(composition) + 1:
        raise AssertionError("reference live-count distribution has the wrong mass")
    result = q1.zero_matrix()
    for live, probability in enumerate(distribution):
        result = q1.matrix_add(
            result, q1.matrix_scale(regions[live], q1.rational(probability))
        )
    return result


def location_count(composition: tuple[int, int, int]) -> int:
    remaining = q1.OUTER_ROWS
    result = 1
    for count in composition:
        result *= math.comb(remaining, count)
        remaining -= count
    return result


def main() -> None:
    ctx.prec = 256
    manifest, local_paths = load_manifest()
    _, q1_paths = q1.load_manifest()
    outer_spectrum = q1.load_outer_spectrum(q1_paths)
    live_spectrum = q1.load_live_spectrum(q1_paths)
    nonactivation = q30.load_nonactivation(q1_paths)
    witnesses, binary64_rows, binary64_unions = load_witnesses(
        local_paths["binary64_composition_witnesses"]
    )
    envelopes, maximizing_weights = band_envelopes(outer_spectrum)

    by_witness: dict[int, list[tuple[int, int, int]]] = defaultdict(list)
    for composition, witness in witnesses.items():
        by_witness[witness].append(composition)

    contributions: dict[tuple[int, int, int], arb] = {}
    rows: list[dict[str, object]] = []
    maximum_pointwise_delta = 0.0
    for witness in sorted(by_witness):
        log_surprisal = q1.rational(Fraction(witness, 10))
        surprisal = log_surprisal.exp()
        z = (-surprisal).exp()
        impulses = q30.impulse_matrices(
            z, live_spectrum, nonactivation, MAX_OCCUPATION
        )
        regions = q30.region_matrices(impulses, MAX_OCCUPATION)
        correction = (arb(q1.BAD_WEIGHT) * surprisal).exp()
        for composition in by_witness[witness]:
            reference = reference_region(regions, composition)
            powered = q1.matrix_power(reference, q1.OUTER_BITS)
            moment = powered[0][0] + powered[0][1]
            raw_inner = moment * correction
            inner = raw_inner if raw_inner < arb(1) else arb(1)
            density = arb(1)
            for count, envelope in zip(composition, envelopes, strict=True):
                density *= q1.rational(envelope) ** count
            contribution = arb(location_count(composition)) * density * inner
            contributions[composition] = contribution
            log2_value = q1.log2_interval(contribution)
            log2_upper = q1.upper_float(log2_value)
            maximum_pointwise_delta = max(
                maximum_pointwise_delta,
                abs(log2_upper - binary64_rows[composition]),
            )
            rows.append(
                {
                    "low_rows": composition[0],
                    "central_rows": composition[1],
                    "high_rows": composition[2],
                    "occupation": sum(composition),
                    "log_surprisal_exact": f"{witness}/10",
                    "log2_upper_interval": str(log2_value),
                    "log2_upper": log2_upper,
                    "margin_bits_lower": -log2_upper,
                    "binary64_log2_upper": binary64_rows[composition],
                }
            )

    if set(contributions) != set(all_compositions()):
        raise AssertionError("outward calculation omitted a composition")
    if maximum_pointwise_delta >= 1e-6:
        raise AssertionError("outward composition does not match binary64 discovery")

    occupation_totals: dict[int, arb] = {}
    occupation_rows: list[dict[str, object]] = []
    maximum_union_delta = 0.0
    for occupation in range(MIN_OCCUPATION, MAX_OCCUPATION + 1):
        total = sum(
            (
                contribution
                for composition, contribution in contributions.items()
                if sum(composition) == occupation
            ),
            arb(0),
        )
        occupation_totals[occupation] = total
        log2_value = q1.log2_interval(total)
        log2_upper = q1.upper_float(log2_value)
        delta = abs(log2_upper - binary64_unions[occupation])
        maximum_union_delta = max(maximum_union_delta, delta)
        if delta >= 1e-6:
            raise AssertionError(f"outward Q={occupation} union changed")
        if not total < arb(1) / (1 << q1.TARGET_MARGIN_BITS):
            raise AssertionError(f"Q={occupation} does not clear 40 bits")
        occupation_rows.append(
            {
                "occupation": occupation,
                "composition_count": math.comb(occupation + 2, 2),
                "log2_failure_probability_interval": str(log2_value),
                "margin_bits_lower": -log2_upper,
                "binary64_margin_bits": -binary64_unions[occupation],
            }
        )

    combined = sum(occupation_totals.values(), arb(0))
    combined_log2 = q1.log2_interval(combined)
    combined_upper = q1.upper_float(combined_log2)
    if not combined < arb(1) / (1 << q1.TARGET_MARGIN_BITS):
        raise AssertionError("Q=5,...,29 union does not clear 40 bits")

    rows.sort(
        key=lambda row: (
            row["occupation"],
            row["low_rows"],
            row["central_rows"],
            row["high_rows"],
        )
    )
    dominant = min(rows, key=lambda row: float(row["margin_bits_lower"]))
    payload = {
        "schema": "rm2sub-rm49-q05-q29-outward-certificate-v1",
        "status": "OUTWARD_OCCUPATION_CERTIFICATE",
        "claim": {
            "event": "some occupation-five through occupation-29 message encodes to output weight at most 13107",
            "failure_probability_upper_interval": str(combined),
            "log2_failure_probability_interval": str(combined_log2),
            "margin_bits_lower": -combined_upper,
            "clears_40_bits": True,
        },
        "coverage": {
            "occupation_interval": [MIN_OCCUPATION, MAX_OCCUPATION],
            "composition_count": len(rows),
            "distinct_witness_count": len(by_witness),
            "complete": True,
        },
        "groups": [
            {
                "name": group[0],
                "weights": [group[1], group[2]],
                "reference_probability_exact": f"{group[3].numerator}/{group[3].denominator}",
                "density_envelope_definition": "A_w/(C(512,w)*p^w*(1-p)^(512-w)) at the displayed maximizing weight",
                "density_envelope_log2_interval": str(
                    q1.log2_interval(q1.rational(envelope))
                ),
                "maximizing_weight": maximizing_weight,
            }
            for group, envelope, maximizing_weight in zip(
                GROUPS, envelopes, maximizing_weights, strict=True
            )
        ],
        "occupation_rows": occupation_rows,
        "dominant_composition": dominant,
        "composition_rows": rows,
        "arithmetic": {
            "library": "python-flint Arb",
            "precision_bits": ctx.prec,
            "maximum_pointwise_log2_delta_from_binary64": maximum_pointwise_delta,
            "maximum_occupation_union_log2_delta_from_binary64": maximum_union_delta,
        },
        "inputs": {
            "manifest": str(MANIFEST),
            "manifest_sha256": sha256(MANIFEST),
            "checker_sha256": sha256(Path(__file__).resolve()),
            "files": manifest["files"],
        },
        "scope": "This certificate covers occupations five through 29 only.",
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"status,{payload['status']}")
    print(f"margin_bits_lower,{-combined_upper:.12f}")
    print(f"composition_count,{len(rows)}")
    print(f"maximum_pointwise_delta,{maximum_pointwise_delta:.3e}")
    print(f"maximum_union_delta,{maximum_union_delta:.3e}")
    print(f"wrote,{OUTPUT}")


if __name__ == "__main__":
    main()
