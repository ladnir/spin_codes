#!/usr/bin/env python3
"""Verify the Gleason-form weight enumerator for RM(4,9).

For the Type-II self-dual binary code RM(4,9), write

    Phi = x^8 + 14 x^4 y^4 + y^8,
    Psi = x^4 y^4 (x^4 - y^4)^4.

The weight enumerator has the form

    sum_j c_j Phi^(64-3j) Psi^j.

This script expands the univariate version at x=1 and checks it against the
committed Markov--Borissov spectrum table.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


GLEASON_COEFFS = [
    1,
    -896,
    366400,
    -90591488,
    15125312416,
    -1804090983680,
    158633285104000,
    -10462864561216000,
    521844917396785920,
    -19703740374727094272,
    560488125281758654464,
    -11885355149550800420864,
    184790978198497210204160,
    -2057511108351137030602752,
    15886111690015316194099200,
    -81428598375642719259197440,
    261053953183414744116101120,
    -481398485705822986116792320,
    451352262788439851856297984,
    -176727086611150405833326592,
    20139935110366512638066688,
    -279475807551445268430848,
]


def poly_mul(p: list[int], q: list[int], max_deg: int) -> list[int]:
    out = [0] * (max_deg + 1)
    for i, pi in enumerate(p):
        if pi == 0:
            continue
        for j, qj in enumerate(q):
            if qj != 0 and i + j <= max_deg:
                out[i + j] += pi * qj
    return out


def poly_pow(p: list[int], exp: int, max_deg: int) -> list[int]:
    out = [0] * (max_deg + 1)
    out[0] = 1
    base = p[:]
    while exp:
        if exp & 1:
            out = poly_mul(out, base, max_deg)
        exp >>= 1
        if exp:
            base = poly_mul(base, base, max_deg)
    return out


def gleason_univariate() -> list[int]:
    # Work in t=y^4.  At x=1,
    #   Phi = 1 + 14t + t^2,
    #   Psi = t(1-t)^4.
    max_deg = 128
    phi = [0] * (max_deg + 1)
    phi[0] = 1
    phi[1] = 14
    phi[2] = 1
    psi = [0] * (max_deg + 1)
    psi[1] = 1
    psi[2] = -4
    psi[3] = 6
    psi[4] = -4
    psi[5] = 1

    total = [0] * (max_deg + 1)
    for j, coeff in enumerate(GLEASON_COEFFS):
        term = poly_mul(poly_pow(phi, 64 - 3 * j, max_deg), poly_pow(psi, j, max_deg), max_deg)
        for i, value in enumerate(term):
            total[i] += coeff * value
    return total


def read_spectrum(path: Path) -> list[int]:
    vals = [0] * 129
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            weight = int(row["weight"])
            count = int(row["count"])
            if weight % 4 != 0:
                raise SystemExit(f"weight {weight} is not divisible by 4")
            vals[weight // 4] = count
    return vals


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectrum", type=Path, default=Path(__file__).with_name("rm512_256_spectrum.csv"))
    args = parser.parse_args()

    formula = gleason_univariate()
    spectrum = read_spectrum(args.spectrum)
    ok = formula == spectrum
    print(f"Gleason formula matches spectrum: {'PASS' if ok else 'FAIL'}")
    print(f"sum = {sum(formula)}")
    print(f"min nonzero weight = {next(4*i for i, v in enumerate(formula) if i > 0 and v)}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
