#!/usr/bin/env python3
"""Certify the all-low composition at occupation 30 with Arb intervals."""

from __future__ import annotations

import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

from flint import arb, ctx

import certify_rm2sub_rm49_q1_outward as q1


HERE = Path(__file__).resolve().parent
MANIFEST = HERE / "RM2SUB_RM49_Q30_LOW_INPUT_MANIFEST.json"
OUTPUT = HERE / "rm2sub_rm49_t64_s14_q30_low_outward.json"
OCCUPATION = 30
LOWER_WEIGHT = 1
UPPER_WEIGHT = 95
REFERENCE_PROBABILITY = Fraction(1, 4)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_local_manifest() -> tuple[dict[str, object], dict[str, Path]]:
    payload = json.loads(MANIFEST.read_text(encoding="utf-8"))
    paths: dict[str, Path] = {}
    for row in payload["files"]:
        path = (HERE / row["path"]).resolve()
        actual = sha256(path)
        if actual != row["sha256"]:
            raise AssertionError(f"hash mismatch for {row['role']}")
        paths[str(row["role"])] = path
    return payload, paths


def load_witness(path: Path) -> tuple[int, float]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    parameters = payload["parameters"]
    expected = {
        "message_bits": q1.MESSAGE_BITS,
        "output_bits": q1.OUTPUT_BITS,
        "outer_rows": q1.OUTER_ROWS,
        "outer_block_bits": q1.OUTER_BITS,
        "step_bits": q1.STEP_BITS,
        "state_bits": q1.STATE_BITS,
        "epochs_per_region": q1.EPOCHS_PER_REGION,
        "bad_weight": q1.BAD_WEIGHT,
    }
    for key, value in expected.items():
        if int(parameters[key]) != value:
            raise AssertionError(f"witness parameter mismatch for {key}")
    selected = [
        row
        for row in payload["composition_rows"]
        if int(row["low_rows"]) == OCCUPATION
        and int(row["central_rows"]) == 0
        and int(row["high_rows"]) == 0
    ]
    if len(selected) != 1:
        raise AssertionError("binary64 receipt does not contain the all-low Q30 row")
    row = selected[0]
    raw = float(row["log_surprisal"])
    witness_tenths = round(10.0 * raw)
    if abs(raw - witness_tenths / 10.0) >= 1e-12:
        raise AssertionError("composition witness is not on the exact tenth grid")
    return witness_tenths, float(row["log2_upper_diagnostic"])


def load_nonactivation(paths: dict[str, Path]) -> list[Fraction]:
    payload = json.loads(paths["rm2sub_b_kernel_spectrum"].read_text(encoding="utf-8"))
    rows = {int(row["total_weight"]): row for row in payload["by_total_weight"]}
    if set(rows) != set(range(q1.STEP_BITS + 1)):
        raise AssertionError("B-kernel receipt does not cover every epoch weight")
    values: list[Fraction] = []
    for weight in range(q1.STEP_BITS + 1):
        row = rows[weight]
        kernel_words = int(row["kernel_words"])
        shell_size = int(row["shell_size"])
        if shell_size != math.comb(q1.STEP_BITS, weight):
            raise AssertionError(f"wrong shell size at weight {weight}")
        values.append(Fraction(kernel_words, shell_size))
    return values


def support_averaged_moments(
    z: arb, live_spectrum: dict[int, int], maximum_input_weight: int
) -> list[arb]:
    state_denominator = arb((1 << q1.STATE_BITS) - 1)
    moments: list[arb] = []
    for input_weight in range(maximum_input_weight + 1):
        shell_denominator = math.comb(q1.STEP_BITS, input_weight)
        total = arb(0)
        for codeword_weight, count in live_spectrum.items():
            kernel = arb(0)
            first = max(0, input_weight - (q1.STEP_BITS - codeword_weight))
            last = min(input_weight, codeword_weight)
            for intersection in range(first, last + 1):
                multiplicity = (
                    math.comb(codeword_weight, intersection)
                    * math.comb(
                        q1.STEP_BITS - codeword_weight,
                        input_weight - intersection,
                    )
                )
                exponent = input_weight + codeword_weight - 2 * intersection
                kernel += arb(multiplicity) * (z**exponent)
            total += arb(count) * kernel / shell_denominator
        moments.append(total / state_denominator)
    return moments


def impulse_matrices(
    z: arb,
    live_spectrum: dict[int, int],
    nonactivation: list[Fraction],
    maximum_input_weight: int,
) -> list[q1.Matrix]:
    moments = support_averaged_moments(z, live_spectrum, maximum_input_weight)
    state_denominator = arb((1 << q1.STATE_BITS) - 1)
    punctured = state_denominator / arb((1 << q1.STATE_BITS) - 2)
    matrices: list[q1.Matrix] = [
        (
            (arb(1), arb(0)),
            (arb(0), punctured * moments[0]),
        )
    ]
    for weight in range(1, maximum_input_weight + 1):
        output_factor = z**weight
        nonactivation_mass = q1.rational(nonactivation[weight])
        averaged = punctured * moments[weight]
        matrices.append(
            (
                (
                    nonactivation_mass * output_factor,
                    (arb(1) - nonactivation_mass) * output_factor,
                ),
                (averaged / state_denominator, averaged),
            )
        )
    return matrices


def region_matrices(
    impulses: list[q1.Matrix], maximum_weight: int
) -> list[q1.Matrix]:
    current = impulses[: maximum_weight + 1]
    current_bits = q1.STEP_BITS
    for _ in range(1, q1.EPOCHS_PER_REGION):
        next_bits = current_bits + q1.STEP_BITS
        updated = [q1.zero_matrix() for _ in range(maximum_weight + 1)]
        for total in range(maximum_weight + 1):
            denominator = arb(math.comb(next_bits, total))
            for next_weight in range(total + 1):
                previous_weight = total - next_weight
                if previous_weight > min(current_bits, maximum_weight):
                    continue
                if next_weight > q1.STEP_BITS:
                    continue
                coefficient = arb(
                    math.comb(current_bits, previous_weight)
                    * math.comb(q1.STEP_BITS, next_weight)
                ) / denominator
                term = q1.matrix_multiply(
                    current[previous_weight], impulses[next_weight]
                )
                updated[total] = q1.matrix_add(
                    updated[total], q1.matrix_scale(term, coefficient)
                )
        current = updated
        current_bits = next_bits
    return current


def reference_region(regions: list[q1.Matrix]) -> q1.Matrix:
    result = q1.zero_matrix()
    p = REFERENCE_PROBABILITY
    for live in range(OCCUPATION + 1):
        probability = (
            Fraction(math.comb(OCCUPATION, live), 1)
            * p**live
            * (1 - p) ** (OCCUPATION - live)
        )
        result = q1.matrix_add(
            result, q1.matrix_scale(regions[live], q1.rational(probability))
        )
    return result


def low_band_envelope(outer_spectrum: dict[int, int]) -> tuple[Fraction, int]:
    best: Fraction | None = None
    best_weight = -1
    p = REFERENCE_PROBABILITY
    for weight, count in outer_spectrum.items():
        if not LOWER_WEIGHT <= weight <= UPPER_WEIGHT:
            continue
        reference = (
            Fraction(math.comb(q1.OUTER_BITS, weight), 1)
            * p**weight
            * (1 - p) ** (q1.OUTER_BITS - weight)
        )
        candidate = Fraction(count, 1) / reference
        if best is None or candidate > best:
            best = candidate
            best_weight = weight
    if best is None or best_weight != 32:
        raise AssertionError("low-band density maximum changed")
    return best, best_weight


def main() -> None:
    ctx.prec = 256
    local_manifest, local_paths = load_local_manifest()
    _, q1_paths = q1.load_manifest()
    outer_spectrum = q1.load_outer_spectrum(q1_paths)
    live_spectrum = q1.load_live_spectrum(q1_paths)
    nonactivation = load_nonactivation(q1_paths)
    witness_tenths, binary64_log2 = load_witness(
        local_paths["binary64_composition_witness"]
    )

    log_surprisal = q1.rational(Fraction(witness_tenths, 10))
    surprisal = log_surprisal.exp()
    z = (-surprisal).exp()
    impulses = impulse_matrices(z, live_spectrum, nonactivation, OCCUPATION)
    regions = region_matrices(impulses, OCCUPATION)
    reference = reference_region(regions)
    powered = q1.matrix_power(reference, q1.OUTER_BITS)
    reference_moment = powered[0][0] + powered[0][1]
    chernoff = reference_moment * (arb(q1.BAD_WEIGHT) * surprisal).exp()
    inner = chernoff if chernoff < arb(1) else arb(1)

    envelope, maximizing_weight = low_band_envelope(outer_spectrum)
    location_count = math.comb(q1.OUTER_ROWS, OCCUPATION)
    total = arb(location_count) * (q1.rational(envelope) ** OCCUPATION) * inner
    total_log2 = q1.log2_interval(total)
    log2_upper = q1.upper_float(total_log2)
    delta = abs(log2_upper - binary64_log2)
    if delta >= 1e-6:
        raise AssertionError("outward composition does not match binary64 discovery")
    if not total < arb(1) / (1 << 1000):
        raise AssertionError("representative composition did not retain ample margin")

    payload = {
        "schema": "rm2sub-rm49-q30-low-outward-certificate-v1",
        "status": "OUTWARD_REPRESENTATIVE_COMPOSITION_CERTIFICATE",
        "claim": {
            "event": "some message with exactly 30 low-band active rows and no central/high rows encodes to weight at most 13107",
            "failure_probability_upper_interval": str(total),
            "log2_failure_probability_interval": str(total_log2),
            "margin_bits_lower": -log2_upper,
        },
        "composition": {
            "occupation": OCCUPATION,
            "low_rows": OCCUPATION,
            "central_rows": 0,
            "high_rows": 0,
            "row_location_count": str(location_count),
        },
        "change_of_measure": {
            "low_band": [LOWER_WEIGHT, UPPER_WEIGHT],
            "reference_probability": "1/4",
            "density_envelope_exact": f"{envelope.numerator}/{envelope.denominator}",
            "density_envelope_log2_interval": str(q1.log2_interval(q1.rational(envelope))),
            "maximizing_weight": maximizing_weight,
        },
        "witness": {
            "log_surprisal_exact": f"{witness_tenths}/10",
            "binary64_log2_upper": binary64_log2,
            "absolute_log2_delta": delta,
        },
        "arithmetic": {
            "library": "python-flint Arb",
            "precision_bits": ctx.prec,
            "all_certificate_operations": "positive Arb interval arithmetic after exact rational input parsing",
        },
        "inputs": {
            "manifest": str(MANIFEST),
            "manifest_sha256": sha256(MANIFEST),
            "checker_sha256": sha256(Path(__file__).resolve()),
            "files": local_manifest["files"],
        },
        "scope": "This receipt validates the refined-band transfer for one dominant composition. It is not an occupation-30 union or an end-to-end certificate.",
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"status,{payload['status']}")
    print(f"margin_bits_lower,{-log2_upper:.12f}")
    print(f"binary64_delta,{delta:.3e}")
    print(f"wrote,{OUTPUT}")


if __name__ == "__main__":
    main()
