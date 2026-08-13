#!/usr/bin/env python3
"""Exact rational Collatz certificates for the systematic PA1 live kernel.

For every rational output pole z, a deterministic floating power iteration is
used only to propose a positive 64-entry integer vector.  The claimed
inequality is then checked from scratch with Fraction arithmetic against the
exact accumulator weight chain, committed exact systematic split rows, and
the rigorous row/column/ordinary-spectrum caps for every unknown split row.
Thus floating point affects witness quality but not certificate validity.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np

from analyze_nosinger_bch_kernel import accumulator_weight_entries
from analyze_systematic_group_kernel import (
    B,
    EXACT_SLICES,
    SPECTRUM,
    maximize_fraction,
)
from certificate_spectra import load_ebch128_spectrum
from probe_systematic_group_multipole import maximize_row


ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "systematic_pa1_geometric_certificates.json"


def load_rows() -> dict[int, Counter[int]]:
    rows: dict[int, Counter[int]] = {}
    with EXACT_SLICES.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            rows.setdefault(int(row["input_weight"]), Counter())[
                int(row["state_weight"])
            ] += int(row["count"])
    return rows


def propose_vector(
    pole: Fraction,
    spectrum: dict[int, int],
    rows: dict[int, Counter[int]],
    scale_bits: int,
) -> tuple[int, ...]:
    z = float(pole)
    powers = np.array([z**t for t in range(B + 1)])
    mixer = np.zeros((B + 1, B + 1))
    for source in range(1, B + 1):
        denominator = math.comb(B, source)
        for target, count in accumulator_weight_entries(source):
            mixer[source, target] = count / denominator

    def apply(values: np.ndarray) -> np.ndarray:
        split = np.zeros(B + 1)
        for t in range(1, B + 1):
            split[t] = powers[t] * maximize_row(t, values, spectrum, rows)
        return mixer @ split

    values = np.ones(B + 1)
    values[0] = 0.0
    for _ in range(192):
        values = apply(values)
        values /= np.max(values)
    scale = 1 << scale_bits
    return tuple(
        0 if state == 0 else max(1, math.ceil(values[state] * scale))
        for state in range(B + 1)
    )


def exact_image(
    vector: tuple[int, ...],
    pole: Fraction,
    spectrum: dict[int, int],
    rows: dict[int, Counter[int]],
) -> tuple[Fraction, ...]:
    split = [Fraction(0) for _ in range(B + 1)]
    for t in range(1, B + 1):
        scores = tuple(pole**t * vector[q] for q in range(B + 1))
        split[t] = maximize_fraction(
            t=t,
            scores=scores,
            spectrum=spectrum,
            exact=rows,
        )
    image = [Fraction(0) for _ in range(B + 1)]
    for source in range(1, B + 1):
        denominator = math.comb(B, source)
        image[source] = sum(
            (count * split[target] for target, count in accumulator_weight_entries(source)),
            Fraction(0),
        ) / denominator
    return tuple(image)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--poles",
        default="19/20,49/50,99/100,997/1000",
        help="comma-separated rational output poles",
    )
    parser.add_argument("--scale-bits", type=int, default=56)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    poles = [Fraction(value) for value in args.poles.split(",")]
    if any(not 0 < pole < 1 for pole in poles):
        raise SystemExit("PA1 certificate: every pole must lie in (0,1)")
    spectrum = {
        weight: count
        for weight, count in load_ebch128_spectrum(SPECTRUM)
        if weight and count
    }
    rows = load_rows()
    certificates = []
    for pole in poles:
        vector = propose_vector(pole, spectrum, rows, args.scale_bits)
        image = exact_image(vector, pole, spectrum, rows)
        rho = max(image[state] / vector[state] for state in range(1, B + 1))
        if not 0 < rho < 1:
            raise SystemExit(f"PA1 certificate: pole {pole} does not contract")
        prefactor = max(Fraction(1, vector[state]) for state in range(1, B + 1)) * max(image)
        vector_digest = hashlib.sha256(
            b"".join(value.to_bytes(16, "little") for value in vector[1:])
        ).hexdigest()
        record = {
            "pole_numerator": pole.numerator,
            "pole_denominator": pole.denominator,
            "rho_numerator": str(rho.numerator),
            "rho_denominator": str(rho.denominator),
            "prefactor_numerator": str(prefactor.numerator),
            "prefactor_denominator": str(prefactor.denominator),
            "rho_log2": math.log2(float(rho)),
            "prefactor_log2": math.log2(float(prefactor)),
            "vector": [str(value) for value in vector[1:]],
            "vector_sha256": vector_digest,
        }
        certificates.append(record)
        print(
            f"pole={pole} rho_log2={record['rho_log2']:.12f} "
            f"prefactor_log2={record['prefactor_log2']:.12f} "
            f"vector_sha256={vector_digest}"
        )
    payload = {
        "status": "EXACT_RATIONAL_COLLATZ_CERTIFICATE",
        "scale_bits": args.scale_bits,
        "turnoff_certificate": "certify_systematic_turnoff_cap.py",
        "certificates": certificates,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
