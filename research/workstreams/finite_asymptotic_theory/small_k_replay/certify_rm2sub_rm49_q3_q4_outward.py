#!/usr/bin/env python3
"""Certify occupations three and four with regular/all-one mixtures."""

from __future__ import annotations

import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

from flint import arb, ctx

import certify_rm2sub_rm49_q1_outward as q1
import certify_rm2sub_rm49_q30_low_outward as q30


HERE = Path(__file__).resolve().parent
MANIFEST = HERE / "RM2SUB_RM49_Q3_Q4_INPUT_MANIFEST.json"
OUTPUT = HERE / "rm2sub_rm49_t64_s14_q3_q4_outward.json"
OCCUPATIONS = (3, 4)
REFERENCE_PROBABILITY = Fraction(1, 2)


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


def load_witnesses(path: Path) -> tuple[dict[tuple[int, int], int], dict[int, float]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    rows_by_q = {int(row["occupation"]): row for row in payload["occupation_rows"]}
    witnesses: dict[tuple[int, int], int] = {}
    diagnostic_unions: dict[int, float] = {}
    for occupation in OCCUPATIONS:
        row = rows_by_q[occupation]
        diagnostic_unions[occupation] = float(row["log2_upper_diagnostic"])
        mixtures = row["mixture_rows"]
        if len(mixtures) != occupation + 1:
            raise AssertionError(f"Q={occupation} does not contain every all-one mixture")
        for mixture in mixtures:
            regular = int(mixture["regular_rows"])
            all_one = int(mixture["all_one_rows"])
            if regular + all_one != occupation:
                raise AssertionError("mixture has the wrong occupation")
            if regular and float(mixture["candidate_probability"]) != 0.5:
                raise AssertionError("regular-row reference probability changed")
            holder = mixture["holder_order"]
            if regular and holder != "infinity":
                raise AssertionError("selected mixture does not use the pointwise envelope")
            raw = float(mixture["log_surprisal"])
            witness = round(10.0 * raw)
            if abs(raw - witness / 10.0) >= 1e-12:
                raise AssertionError("output witness is not on the exact tenth grid")
            witnesses[(occupation, all_one)] = witness
    return witnesses, diagnostic_unions


def regular_density_envelope(outer_spectrum: dict[int, int]) -> tuple[Fraction, list[int]]:
    best: Fraction | None = None
    maximizing_weights: list[int] = []
    for weight, count in outer_spectrum.items():
        if weight == q1.OUTER_BITS:
            continue
        reference = Fraction(math.comb(q1.OUTER_BITS, weight), 1 << q1.OUTER_BITS)
        candidate = Fraction(count, 1) / reference
        if best is None or candidate > best:
            best = candidate
            maximizing_weights = [weight]
        elif candidate == best:
            maximizing_weights.append(weight)
    if best is None or 32 not in maximizing_weights:
        raise AssertionError("regular RM density maximum changed")
    return best, maximizing_weights


def reference_region(
    regions: list[q1.Matrix], regular: int, forced: int
) -> q1.Matrix:
    result = q1.zero_matrix()
    for live_regular in range(regular + 1):
        probability = Fraction(math.comb(regular, live_regular), 1 << regular)
        result = q1.matrix_add(
            result,
            q1.matrix_scale(regions[forced + live_regular], q1.rational(probability)),
        )
    return result


def location_count(regular: int, all_one: int) -> int:
    return math.comb(q1.OUTER_ROWS, all_one) * math.comb(
        q1.OUTER_ROWS - all_one, regular
    )


def main() -> None:
    ctx.prec = 256
    manifest, local_paths = load_manifest()
    _, q1_paths = q1.load_manifest()
    outer_spectrum = q1.load_outer_spectrum(q1_paths)
    live_spectrum = q1.load_live_spectrum(q1_paths)
    nonactivation = q30.load_nonactivation(q1_paths)
    witnesses, diagnostic_unions = load_witnesses(
        local_paths["binary64_two_colour_witnesses"]
    )
    envelope, maximizing_weights = regular_density_envelope(outer_spectrum)

    region_cache: dict[int, list[q1.Matrix]] = {}
    mixture_rows: list[dict[str, object]] = []
    occupation_totals: dict[int, arb] = {occupation: arb(0) for occupation in OCCUPATIONS}
    for occupation in OCCUPATIONS:
        for all_one in range(occupation + 1):
            regular = occupation - all_one
            witness = witnesses[(occupation, all_one)]
            if witness not in region_cache:
                log_surprisal = q1.rational(Fraction(witness, 10))
                surprisal = log_surprisal.exp()
                z = (-surprisal).exp()
                impulses = q30.impulse_matrices(
                    z, live_spectrum, nonactivation, max(OCCUPATIONS)
                )
                region_cache[witness] = q30.region_matrices(
                    impulses, max(OCCUPATIONS)
                )
            else:
                surprisal = q1.rational(Fraction(witness, 10)).exp()
            reference = reference_region(region_cache[witness], regular, all_one)
            powered = q1.matrix_power(reference, q1.OUTER_BITS)
            moment = powered[0][0] + powered[0][1]
            raw_inner = moment * (arb(q1.BAD_WEIGHT) * surprisal).exp()
            inner = raw_inner if raw_inner < arb(1) else arb(1)
            density = q1.rational(envelope) ** regular
            contribution = arb(location_count(regular, all_one)) * density * inner
            occupation_totals[occupation] += contribution
            log2_value = q1.log2_interval(contribution)
            mixture_rows.append(
                {
                    "occupation": occupation,
                    "regular_rows": regular,
                    "all_one_rows": all_one,
                    "log_surprisal_exact": f"{witness}/10",
                    "log2_upper_interval": str(log2_value),
                    "log2_upper": q1.upper_float(log2_value),
                    "margin_bits_lower": -q1.upper_float(log2_value),
                }
            )

    occupation_rows = []
    for occupation in OCCUPATIONS:
        total = occupation_totals[occupation]
        log2_value = q1.log2_interval(total)
        log2_upper = q1.upper_float(log2_value)
        if not total < arb(1) / (1 << q1.TARGET_MARGIN_BITS):
            raise AssertionError(f"Q={occupation} does not clear 40 bits")
        occupation_rows.append(
            {
                "occupation": occupation,
                "log2_failure_probability_interval": str(log2_value),
                "margin_bits_lower": -log2_upper,
                "prior_binary64_margin_bits": -diagnostic_unions[occupation],
                "improvement_from_tighter_epoch_transfer_bits": (
                    -log2_upper + diagnostic_unions[occupation]
                ),
            }
        )
    combined = sum(occupation_totals.values(), arb(0))
    combined_log2 = q1.log2_interval(combined)
    combined_upper = q1.upper_float(combined_log2)
    if not combined < arb(1) / (1 << q1.TARGET_MARGIN_BITS):
        raise AssertionError("Q3--4 union does not clear 40 bits")

    payload = {
        "schema": "rm2sub-rm49-q3-q4-outward-certificate-v1",
        "status": "OUTWARD_OCCUPATION_CERTIFICATE",
        "claim": {
            "event": "some occupation-three or occupation-four message encodes to output weight at most 13107",
            "failure_probability_upper_interval": str(combined),
            "log2_failure_probability_interval": str(combined_log2),
            "margin_bits_lower": -combined_upper,
            "clears_40_bits": True,
        },
        "reduction": "Separate the unique all-one RM word. Dominate every other routed RM word pointwise by the uniform Bernoulli-half reference. Sum every regular/all-one mixture with the exact positive support-averaged RM2Sub transfer.",
        "regular_density": {
            "reference_probability": "1/2",
            "envelope_log2_interval": str(q1.log2_interval(q1.rational(envelope))),
            "maximizing_weights": maximizing_weights,
        },
        "occupation_rows": occupation_rows,
        "mixture_rows": mixture_rows,
        "arithmetic": {
            "library": "python-flint Arb",
            "precision_bits": ctx.prec,
            "proof_change_from_discovery": "exact A-spectrum support average and exact zero-state destination split replace the older factorial-moment relaxation",
        },
        "inputs": {
            "manifest": str(MANIFEST),
            "manifest_sha256": sha256(MANIFEST),
            "checker_sha256": sha256(Path(__file__).resolve()),
            "files": manifest["files"],
        },
        "scope": "This certificate covers occupations three and four only.",
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"status,{payload['status']}")
    for row in occupation_rows:
        print(f"Q,{row['occupation']},margin_bits_lower,{row['margin_bits_lower']:.12f}")
    print(f"combined_margin_bits_lower,{-combined_upper:.12f}")
    print(f"wrote,{OUTPUT}")


if __name__ == "__main__":
    main()
