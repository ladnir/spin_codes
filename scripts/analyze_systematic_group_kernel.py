#!/usr/bin/env python3
"""Robust local kernel for a systematic-output fixed-half Riffle step.

Change basis in the committed EBCH encoder so that

    E_sys(v) = (v, P v).

Before a step, independently permute the 64 incoming state coordinates and
the 64 coordinates of the current striped input group.  Conditional on their
weights, the XOR drive is therefore uniform on its Hamming slice.  Given drive
weight t, the emitted weight is exactly t and the next-state weight q has an
unknown split enumerator C[t,q].  The rigorous cap

    C[t,q] <= A_EBCH[t+q]

and the exact row mass sum_q C[t,q]=C(64,t) define a pessimistic Bellman
operator.  Every row may choose its worst capped distribution independently,
so a contraction below one is a rigorous local envelope even without knowing
the exact split weight enumerator.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import math
from collections import Counter
from decimal import Decimal, ROUND_CEILING, localcontext
from fractions import Fraction
from pathlib import Path

from certificate_spectra import load_ebch128_spectrum


ROOT = Path(__file__).resolve().parent
SPECTRUM = ROOT / "EBCH128_64.wd"
EXACT_SLICES = ROOT / "ebch128_systematic_split_slices.csv"
B = 64
MASK = (1 << B) - 1
GENERATOR = 0xF4845518B9582A1F


def apply(columns: tuple[int, ...], value: int) -> int:
    result = 0
    while value:
        bit = value & -value
        result ^= columns[bit.bit_length() - 1]
        value ^= bit
    return result


def inverse_columns(columns: tuple[int, ...]) -> tuple[int, ...]:
    rows = []
    for output in range(B):
        left = sum(((columns[column] >> output) & 1) << column for column in range(B))
        rows.append(left | (1 << (B + output)))
    for column in range(B):
        pivot = next((row for row in range(column, B) if (rows[row] >> column) & 1), None)
        if pivot is None:
            raise SystemExit("systematic kernel: singular BCH half")
        rows[column], rows[pivot] = rows[pivot], rows[column]
        for row in range(B):
            if row != column and ((rows[row] >> column) & 1):
                rows[row] ^= rows[column]
    inverse_rows = tuple(row >> B for row in rows)
    return tuple(
        sum(((inverse_rows[input_bit] >> output_bit) & 1) << input_bit for input_bit in range(B))
        for output_bit in range(B)
    )


def systematic_state_columns() -> tuple[int, ...]:
    rows = tuple((GENERATOR << row) | (1 << 127) for row in range(B))
    left = tuple(row & MASK for row in rows)
    right = tuple((row >> B) & MASK for row in rows)
    inverse_left = inverse_columns(left)
    state = tuple(apply(right, inverse_left[column]) for column in range(B))
    if any((apply(left, inverse_left[column]) != 1 << column) for column in range(B)):
        raise SystemExit("systematic kernel: inverse verification failed")
    return state


def exact_rows(
    state_columns: tuple[int, ...], radius: int
) -> dict[int, Counter[int]]:
    result: dict[int, Counter[int]] = {}
    all_state = apply(state_columns, MASK)
    if all_state != MASK:
        raise SystemExit("systematic kernel: all-ones word does not split as all ones")
    for t in range(1, radius + 1):
        row: Counter[int] = Counter()
        for indices in itertools.combinations(range(B), t):
            state = 0
            for index in indices:
                state ^= state_columns[index]
            row[state.bit_count()] += 1
        if sum(row.values()) != math.comb(B, t):
            raise SystemExit("systematic kernel: exact row mass mismatch")
        result[t] = row
        result[B - t] = Counter({B - q: count for q, count in row.items()})
    result[B] = Counter({B: 1})
    return result


def load_exact_rows(path: Path) -> dict[int, Counter[int]]:
    result: dict[int, Counter[int]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            t = int(row["input_weight"])
            q = int(row["state_weight"])
            result.setdefault(t, Counter())[q] += int(row["count"])
    for t, histogram in result.items():
        if sum(histogram.values()) != math.comb(B, t):
            raise SystemExit(f"systematic kernel: committed row t={t} has wrong mass")
    return result


def maximize_decimal(
    *,
    t: int,
    scores: tuple[Decimal, ...],
    spectrum: dict[int, int],
    exact: dict[int, Counter[int]],
) -> Decimal:
    population = math.comb(B, t)
    if t in exact:
        return sum(Decimal(count) * scores[q] for q, count in exact[t].items()) / Decimal(
            population
        )
    remaining = population
    numerator = Decimal(0)
    candidates = sorted(
        (
            (scores[q], q, min(spectrum.get(t + q, 0), math.comb(B, q)))
            for q in range(1, B + 1)
            if spectrum.get(t + q, 0)
        ),
        reverse=True,
    )
    for score, _q, cap in candidates:
        take = min(remaining, cap)
        numerator += Decimal(take) * score
        remaining -= take
        if not remaining:
            break
    if remaining:
        raise SystemExit(f"systematic kernel: spectrum caps do not cover t={t}")
    return numerator / Decimal(population)


def maximize_fraction(
    *,
    t: int,
    scores: tuple[Fraction, ...],
    spectrum: dict[int, int],
    exact: dict[int, Counter[int]],
) -> Fraction:
    population = math.comb(B, t)
    if t in exact:
        return sum((count * scores[q] for q, count in exact[t].items()), Fraction(0)) / population
    remaining = population
    numerator = Fraction(0)
    candidates = sorted(
        (
            (scores[q], q, min(spectrum.get(t + q, 0), math.comb(B, q)))
            for q in range(1, B + 1)
            if spectrum.get(t + q, 0)
        ),
        reverse=True,
    )
    for score, _q, cap in candidates:
        take = min(remaining, cap)
        numerator += take * score
        remaining -= take
        if not remaining:
            break
    if remaining:
        raise SystemExit(f"systematic kernel: exact caps do not cover t={t}")
    return numerator / population


def log2_decimal(value: Decimal) -> Decimal:
    with localcontext() as context:
        context.prec = 50
        context.rounding = ROUND_CEILING
        return value.ln() / Decimal(2).ln()


def exact_geometric_certificate(
    *,
    approximate: tuple[Decimal, ...],
    pole: Fraction,
    spectrum: dict[int, int],
    exact: dict[int, Counter[int]],
) -> tuple[Fraction, Fraction]:
    maximum = max(approximate[1:])
    scale = 1 << 48
    vector = tuple(
        0 if t == 0 else max(1, int(approximate[t] / maximum * scale + 1))
        for t in range(B + 1)
    )
    transformed = [Fraction(0) for _ in range(B + 1)]
    for t in range(1, B + 1):
        scores = tuple(pole**t * vector[q] for q in range(B + 1))
        transformed[t] = maximize_fraction(
            t=t, scores=scores, spectrum=spectrum, exact=exact
        )
    rho = max(transformed[t] / vector[t] for t in range(1, B + 1))
    if not 0 < rho < 1:
        raise SystemExit("systematic kernel: geometric contraction is not below one")
    prefactor = max(Fraction(1, vector[t]) for t in range(1, B + 1)) * max(transformed)
    return rho, prefactor


def exact_burnin_geometric_certificate(
    *,
    pole: Fraction,
    spectrum: dict[int, int],
    exact: dict[int, Counter[int]],
    burnin: int,
) -> tuple[Fraction, Fraction]:
    """Exact episode envelope using the reached vector after a finite burn-in.

    Starting from the all-one terminal vector avoids the very loose global
    domination factor incurred by an arbitrary Collatz test vector.  If v is
    the exact vector after ``burnin`` steps and O(v)<=rho*v, monotonicity and
    homogeneity control every later step.  The explicit earlier maxima fix the
    prefactor for the complete sequence.
    """

    values = tuple(Fraction(1) for _ in range(B + 1))
    maxima: list[Fraction] = []
    for _ in range(burnin):
        transformed = [Fraction(0) for _ in range(B + 1)]
        for t in range(1, B + 1):
            scores = tuple(pole**t * values[q] for q in range(B + 1))
            transformed[t] = maximize_fraction(
                t=t, scores=scores, spectrum=spectrum, exact=exact
            )
        values = tuple(transformed)
        maxima.append(max(values))

    image = [Fraction(0) for _ in range(B + 1)]
    for t in range(1, B + 1):
        scores = tuple(pole**t * values[q] for q in range(B + 1))
        image[t] = maximize_fraction(
            t=t, scores=scores, spectrum=spectrum, exact=exact
        )
    rho = max(image[t] / values[t] for t in range(1, B + 1))
    if not 0 < rho < 1:
        raise SystemExit("systematic kernel: burn-in contraction is not below one")
    prefactor = max(
        maximum / rho**step for step, maximum in enumerate(maxima)
    )
    return rho, prefactor


def chernoff_crossing(rho: Fraction, prefactor: Fraction, pole: Fraction, distance: int) -> int:
    estimate = 1 + math.ceil(
        (
            distance * math.log2(pole.denominator / pole.numerator)
            + math.log2(prefactor.numerator / prefactor.denominator)
        )
        / math.log2(rho.denominator / rho.numerator)
    )

    def passes(live: int) -> bool:
        exponent = live - 1
        return (
            prefactor.numerator
            * pow(rho.numerator, exponent)
            * pow(pole.denominator, distance)
            <= prefactor.denominator
            * pow(rho.denominator, exponent)
            * pow(pole.numerator, distance)
        )

    while estimate > 1 and passes(estimate - 1):
        estimate -= 1
    while not passes(estimate):
        estimate += 1
    return estimate


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectrum", type=Path, default=SPECTRUM)
    parser.add_argument("--exact-slices", type=Path, default=EXACT_SLICES)
    parser.add_argument("--pole-num", type=int, default=2333)
    parser.add_argument("--pole-den", type=int, default=2373)
    parser.add_argument("--steps", type=int, default=24)
    parser.add_argument("--exact-radius", type=int, default=3)
    parser.add_argument("--total-length", type=int, default=2**21)
    parser.add_argument("--exact-burnin", type=int, default=12)
    args = parser.parse_args()

    loaded = load_ebch128_spectrum(args.spectrum)
    spectrum = {weight: count for weight, count in loaded if weight and count}
    if sum(spectrum.values()) != (1 << B) - 1:
        raise SystemExit("systematic kernel: invalid spectrum mass")

    state_columns = systematic_state_columns()
    if not 0 <= args.exact_radius <= 6:
        raise SystemExit("systematic kernel: exact radius must be between zero and six")
    exact = exact_rows(state_columns, args.exact_radius)
    for t, histogram in load_exact_rows(args.exact_slices).items():
        if t in exact and exact[t] != histogram:
            raise SystemExit(f"systematic kernel: regenerated row t={t} differs")
        exact[t] = histogram
    digest = hashlib.sha256(
        b"".join(column.to_bytes(8, "little") for column in state_columns)
    ).hexdigest()
    print("systematic-output fixed-half robust kernel")
    print(f"state_matrix_sha256={digest}")
    print(f"state_column_weight_min={min(map(int.bit_count, state_columns))}")
    print(f"state_column_weight_max={max(map(int.bit_count, state_columns))}")
    for t in sorted(exact):
        row = exact[t]
        if any(count > spectrum.get(t + q, 0) for q, count in row.items()):
            raise SystemExit(f"systematic kernel: exact row t={t} exceeds BCH spectrum")
        print(
            f"exact_t={t} count={sum(row.values())} min_q={min(row)} "
            f"max_q={max(row)} mean_q="
            f"{sum(q * count for q, count in row.items()) / sum(row.values()):.9f}"
        )

    with localcontext() as context:
        context.prec = 80
        context.rounding = ROUND_CEILING
        z = Decimal(args.pole_num) / Decimal(args.pole_den)
        values = tuple(Decimal(1) for _ in range(B + 1))
        for step in range(1, args.steps + 1):
            next_values = [Decimal(0) for _ in range(B + 1)]
            for t in range(1, B + 1):
                scores = tuple(z**t * values[q] for q in range(B + 1))
                next_values[t] = maximize_decimal(
                    t=t, scores=scores, spectrum=spectrum, exact=exact
                )
            values = tuple(next_values)
            worst_t = max(range(1, B + 1), key=values.__getitem__)
            worst = values[worst_t]
            print(
                f"episode_steps={step:2d} activation_mgf_upper={worst:.18E} "
                f"log2_upper={log2_decimal(worst):.12f} activation_t={worst_t}"
            )

    pole = Fraction(args.pole_num, args.pole_den)
    rho, prefactor = exact_geometric_certificate(
        approximate=values, pole=pole, spectrum=spectrum, exact=exact
    )
    print(
        f"exact_geometric_rho={rho.numerator}/{rho.denominator} "
        f"rho_log2={math.log2(rho.numerator / rho.denominator):.12f}"
    )
    print(
        f"exact_geometric_prefactor={prefactor.numerator}/{prefactor.denominator} "
        f"prefactor_log2={math.log2(prefactor.numerator / prefactor.denominator):.12f}"
    )
    crossing = chernoff_crossing(
        rho, prefactor, pole, (9 * args.total_length) // 100
    )
    print(f"exact_geometric_chernoff_crossing={crossing}")

    burnin_rho, burnin_prefactor = exact_burnin_geometric_certificate(
        pole=pole,
        spectrum=spectrum,
        exact=exact,
        burnin=args.exact_burnin,
    )
    print(
        f"exact_burnin_rho={burnin_rho.numerator}/{burnin_rho.denominator} "
        f"rho_log2={math.log2(burnin_rho.numerator / burnin_rho.denominator):.12f}"
    )
    print(
        f"exact_burnin_prefactor={burnin_prefactor.numerator}/{burnin_prefactor.denominator} "
        f"prefactor_log2="
        f"{math.log2(burnin_prefactor.numerator / burnin_prefactor.denominator):.12f}"
    )
    burnin_crossing = chernoff_crossing(
        burnin_rho, burnin_prefactor, pole, (9 * args.total_length) // 100
    )
    print(f"exact_burnin_chernoff_crossing={burnin_crossing}")


if __name__ == "__main__":
    main()
