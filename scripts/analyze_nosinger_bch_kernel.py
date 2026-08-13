#!/usr/bin/env python3
"""Rigorous local-kernel diagnostics for a Riffle inner without Singer.

For a state/input XOR that is uniform on the weight-t slice, define

    B[t,w] = #{v : wt(v)=t and wt(E(v))=w}.

The ordinary BCH spectrum gives the rigorous cap B[t,w] <= A[w].  This
script computes B[t,w] exactly near both ends of the input-weight range and
uses the caps elsewhere in a robust Bellman envelope.  At every bounded row
the adversary may redistribute the C(64,t) inputs among BCH output weights,
subject only to the exact ordinary spectrum caps.  This is deliberately more
pessimistic than any fixed encoder and therefore gives an upper bound on the
surviving low-output MGF through a zero-input episode.

All combinatorial counts and cap checks are integer exact.  Decimal MGF
arithmetic is rounded toward +infinity, so the printed envelope is outward.
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
EXACT_SLICES = ROOT / "ebch128_cyclic_input_slices.csv"
B = 64
N_LOCAL = 128
GENERATOR = 0xF4845518B9582A1F
GENERATOR_ROWS_SHA256 = "972cfc1c6de12e4ddc0c67680fd8e3cddedbd84d0a409c7061e09a86b818ca56"


def generator_rows() -> tuple[int, ...]:
    rows = tuple((GENERATOR << row) | (1 << 127) for row in range(B))
    digest = hashlib.sha256(b"".join(row.to_bytes(16, "little") for row in rows)).hexdigest()
    if digest != GENERATOR_ROWS_SHA256:
        raise SystemExit("no-Singer kernel: BCH generator-row fingerprint mismatch")
    return rows


def exact_slice(rows: tuple[int, ...], weight: int) -> Counter[int]:
    """Enumerate B[weight,*], using the complement when that is cheaper."""

    if not 0 <= weight <= B:
        raise ValueError("input weight out of range")
    complement = weight > B // 2
    choose = B - weight if complement else weight
    base = 0
    if complement:
        for row in rows:
            base ^= row
    result: Counter[int] = Counter()
    for indices in itertools.combinations(range(B), choose):
        word = base
        for index in indices:
            word ^= rows[index]
        result[word.bit_count()] += 1
    if sum(result.values()) != math.comb(B, weight):
        raise SystemExit("no-Singer kernel: exact slice count mismatch")
    return result


def load_exact_slices(path: Path) -> dict[int, Counter[int]]:
    result: dict[int, Counter[int]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            input_weight = int(row["input_weight"])
            output_weight = int(row["output_weight"])
            count = int(row["count"])
            result.setdefault(input_weight, Counter())[output_weight] += count
    for input_weight, histogram in result.items():
        if sum(histogram.values()) != math.comb(B, input_weight):
            raise SystemExit("no-Singer kernel: committed exact slice count mismatch")
    return result


def decimal_up(numerator: int, denominator: int = 1) -> Decimal:
    with localcontext() as context:
        context.prec = 80
        context.rounding = ROUND_CEILING
        return Decimal(numerator) / Decimal(denominator)


def split_rows(weights: tuple[int, ...]) -> dict[int, tuple[tuple[int, int, Decimal], ...]]:
    rows: dict[int, tuple[tuple[int, int, Decimal], ...]] = {}
    for weight in weights:
        denominator = math.comb(N_LOCAL, weight)
        terms: list[tuple[int, int, Decimal]] = []
        for state_weight in range(max(0, weight - B), min(B, weight) + 1):
            output_weight = weight - state_weight
            probability = decimal_up(
                math.comb(B, output_weight) * math.comb(B, state_weight), denominator
            )
            terms.append((output_weight, state_weight, probability))
        rows[weight] = tuple(terms)
    return rows


def split_score(
    terms: tuple[tuple[int, int, Decimal], ...],
    z_powers: tuple[Decimal, ...],
    next_values: tuple[Decimal, ...],
) -> Decimal:
    # State zero is excluded: this operator follows a still-live episode.
    value = Decimal(0)
    for output_weight, state_weight, probability in terms:
        if state_weight:
            value += probability * z_powers[output_weight] * next_values[state_weight]
    return value


def maximize_capped_distribution(
    scores: dict[int, Decimal], spectrum: dict[int, int], population: int
) -> Decimal:
    """Maximize sum p[w]score[w], with p[w] <= A[w]/population."""

    remaining = population
    numerator = Decimal(0)
    for weight in sorted(scores, key=scores.__getitem__, reverse=True):
        take = min(remaining, spectrum[weight])
        numerator += Decimal(take) * scores[weight]
        remaining -= take
        if not remaining:
            break
    if remaining:
        raise SystemExit("no-Singer kernel: spectrum caps do not cover an input slice")
    return numerator / Decimal(population)


def robust_row_value(
    input_weight: int,
    scores: dict[int, Decimal],
    spectrum: dict[int, int],
    exact: dict[int, Counter[int]],
) -> Decimal:
    population = math.comb(B, input_weight)
    if input_weight in exact:
        row = exact[input_weight]
        if any(count > spectrum.get(weight, 0) for weight, count in row.items()):
            raise SystemExit("no-Singer kernel: exact slice exceeds ordinary spectrum")
        return sum(Decimal(count) * scores[weight] for weight, count in row.items()) / Decimal(population)
    return maximize_capped_distribution(scores, spectrum, population)


def accumulator_weight_entries(input_weight: int) -> tuple[tuple[int, int], ...]:
    """Exact (output weight,count) row for the length-64 accumulator."""

    k = (input_weight + 1) // 2
    entries = tuple(
        (
            output_weight,
            math.comb(output_weight - 1, k - 1)
            * math.comb(B - output_weight, input_weight - k),
        )
        for output_weight in range(k, B - input_weight + k + 1)
    )
    if sum(count for _weight, count in entries) != math.comb(B, input_weight):
        raise SystemExit("no-Singer kernel: accumulator row count mismatch")
    return entries


def mixer_rounds(mixer: str) -> int:
    if mixer == "identity":
        return 0
    if mixer == "pap":
        return 1
    if mixer.startswith("pap") and mixer[3:].isdigit():
        rounds = int(mixer[3:])
        if 1 <= rounds <= 16:
            return rounds
    raise ValueError(f"unsupported mixer {mixer!r}")


def mixed_row_values(
    scores: dict[int, Decimal],
    spectrum: dict[int, int],
    exact: dict[int, Counter[int]],
    mixer: str,
) -> list[Decimal]:
    values = [Decimal(0)] + [
        robust_row_value(input_weight, scores, spectrum, exact)
        for input_weight in range(1, B + 1)
    ]
    for _ in range(mixer_rounds(mixer)):
        previous = values
        values = [Decimal(0)] * (B + 1)
        for input_weight in range(1, B + 1):
            denominator = math.comb(B, input_weight)
            values[input_weight] = sum(
                Decimal(count) * previous[weight]
                for weight, count in accumulator_weight_entries(input_weight)
            ) / Decimal(denominator)
    return values


def maximize_capped_fraction(
    scores: dict[int, Fraction], spectrum: dict[int, int], population: int
) -> Fraction:
    remaining = population
    numerator = Fraction(0)
    for weight in sorted(scores, key=scores.__getitem__, reverse=True):
        take = min(remaining, spectrum[weight])
        numerator += take * scores[weight]
        remaining -= take
        if not remaining:
            break
    if remaining:
        raise SystemExit("no-Singer kernel: exact spectrum caps do not cover a slice")
    return numerator / population


def mixed_row_fraction_values(
    scores: dict[int, Fraction],
    spectrum: dict[int, int],
    exact: dict[int, Counter[int]],
    mixer: str,
) -> list[Fraction]:
    """Exact-rational counterpart of :func:`mixed_row_values`."""

    direct = [Fraction(0) for _ in range(B + 1)]
    for input_weight in range(1, B + 1):
        population = math.comb(B, input_weight)
        if input_weight in exact:
            direct[input_weight] = sum(
                (count * scores[weight] for weight, count in exact[input_weight].items()),
                Fraction(0),
            ) / population
        else:
            direct[input_weight] = maximize_capped_fraction(
                scores, spectrum, population
            )
    transformed = direct
    for _ in range(mixer_rounds(mixer)):
        previous = transformed
        transformed = [Fraction(0) for _ in range(B + 1)]
        for input_weight in range(1, B + 1):
            denominator = math.comb(B, input_weight)
            transformed[input_weight] = sum(
                (
                    count * previous[weight]
                    for weight, count in accumulator_weight_entries(input_weight)
                ),
                Fraction(0),
            ) / denominator
    return transformed


def exact_finite_termination_bound(
    *,
    spectrum: dict[int, int],
    exact: dict[int, Counter[int]],
    mixer: str,
    remaining_blocks: int,
    remaining_ones: int,
) -> tuple[Fraction, int]:
    """Worst one-step episode-termination atom for a finite input slice.

    State weight zero terminates at the current split.  For state weight
    ``q > 0``, termination requires the next input block to equal that fixed
    state support.  With ``H`` ones uniformly placed in ``b*T`` coordinates,
    this exact-hit probability is

        C(b*T-b, H-q) / C(b*T, H).

    The robust BCH slice caps are then maximized exactly, including any PAP
    mixer selected by the caller.
    """

    if remaining_blocks < 1:
        raise ValueError("remaining_blocks must be positive")
    coordinates = B * remaining_blocks
    if not 0 <= remaining_ones <= coordinates:
        raise ValueError("remaining_ones out of range")
    if coordinates + 1 < B * (remaining_ones + 1):
        raise SystemExit(
            "no-Singer kernel: finite termination endpoint monotonicity condition failed"
        )

    denominator = math.comb(coordinates, remaining_ones)
    exact_hit = [Fraction(0) for _ in range(B + 1)]
    exact_hit[0] = Fraction(1)
    for state_weight in range(1, B + 1):
        residual = remaining_ones - state_weight
        if 0 <= residual <= coordinates - B:
            exact_hit[state_weight] = Fraction(
                math.comb(coordinates - B, residual), denominator
            )

    scores: dict[int, Fraction] = {}
    for weight in spectrum:
        split_denominator = math.comb(N_LOCAL, weight)
        scores[weight] = sum(
            (
                Fraction(
                    math.comb(B, weight - state_weight)
                    * math.comb(B, state_weight),
                    split_denominator,
                )
                * exact_hit[state_weight]
                for state_weight in range(
                    max(0, weight - B), min(B, weight) + 1
                )
            ),
            Fraction(0),
        )

    values = mixed_row_fraction_values(scores, spectrum, exact, mixer)
    worst_input_weight = max(range(1, B + 1), key=values.__getitem__)
    return values[worst_input_weight], worst_input_weight


def exact_geometric_certificate(
    *,
    approximate: tuple[Decimal, ...],
    spectrum: dict[int, int],
    exact: dict[int, Counter[int]],
    pole: Fraction,
    mixer: str,
    scale_bits: int = 48,
) -> tuple[Fraction, Fraction]:
    """Return exact (rho,prefactor) with episode_MGF(L)<=prefactor*rho^(L-1)."""

    maximum = max(approximate[1:])
    scale = 1 << scale_bits
    vector = tuple(
        0 if state == 0 else max(1, int(approximate[state] / maximum * scale + 1))
        for state in range(B + 1)
    )
    scores: dict[int, Fraction] = {}
    for weight, count in spectrum.items():
        denominator = math.comb(N_LOCAL, weight)
        score = Fraction(0)
        for state_weight in range(max(1, weight - B), min(B, weight) + 1):
            output_weight = weight - state_weight
            probability = Fraction(
                math.comb(B, output_weight) * math.comb(B, state_weight), denominator
            )
            score += probability * pole**output_weight * vector[state_weight]
        scores[weight] = score

    transformed = mixed_row_fraction_values(scores, spectrum, exact, mixer)
    rho = max(transformed[state] / vector[state] for state in range(1, B + 1))
    if not 0 < rho < 1:
        raise SystemExit("no-Singer kernel: exact geometric contraction is not below one")
    domination = max(Fraction(1, vector[state]) for state in range(1, B + 1))
    prefactor = domination * max(transformed[1:])
    return rho, prefactor


def exact_chernoff_crossing(
    rho: Fraction, prefactor: Fraction, pole: Fraction, distance: int
) -> int:
    """First live length whose geometric Chernoff upper bound is below one."""

    estimate = 1 + math.ceil(
        (
            distance * (math.log2(pole.denominator) - math.log2(pole.numerator))
            + math.log2(prefactor.numerator)
            - math.log2(prefactor.denominator)
        )
        / (math.log2(rho.denominator) - math.log2(rho.numerator))
    )

    def passes(live: int) -> bool:
        exponent = live - 1
        left = (
            prefactor.numerator
            * pow(rho.numerator, exponent)
            * pow(pole.denominator, distance)
        )
        right = (
            prefactor.denominator
            * pow(rho.denominator, exponent)
            * pow(pole.numerator, distance)
        )
        return left <= right

    while estimate > 1 and passes(estimate - 1):
        estimate -= 1
    while not passes(estimate):
        estimate += 1
    return estimate


def turnoff_score(terms: tuple[tuple[int, int, Decimal], ...]) -> Decimal:
    return sum(probability for _output, state, probability in terms if state == 0)


def log2_decimal(value: Decimal) -> Decimal:
    if value <= 0:
        return Decimal("-Infinity")
    with localcontext() as context:
        context.prec = 50
        context.rounding = ROUND_CEILING
        return value.ln() / Decimal(2).ln()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectrum", type=Path, default=SPECTRUM)
    parser.add_argument("--exact-slices", type=Path, default=EXACT_SLICES)
    parser.add_argument("--exact-radius", type=int, default=5)
    parser.add_argument("--steps", type=int, default=16)
    parser.add_argument("--pole-num", type=int, default=2333)
    parser.add_argument("--pole-den", type=int, default=2373)
    parser.add_argument("--total-length", type=int, default=2**21)
    parser.add_argument("--remaining-blocks", type=int, default=7527)
    parser.add_argument("--remaining-ones", type=int, default=499)
    parser.add_argument(
        "--mixer",
        choices=("identity", "pap", *(f"pap{rounds}" for rounds in range(2, 17))),
        default="identity",
    )
    args = parser.parse_args()
    if not 0 <= args.exact_radius <= B // 2:
        raise SystemExit("no-Singer kernel: invalid exact radius")

    rows = generator_rows()
    loaded = load_ebch128_spectrum(args.spectrum)
    spectrum = {weight: count for weight, count in loaded if weight > 0 and count > 0}
    if sum(spectrum.values()) != (1 << B) - 1:
        raise SystemExit("no-Singer kernel: invalid nonzero BCH spectrum mass")
    weights = tuple(sorted(spectrum))
    split = split_rows(weights)

    exact_weights = tuple(range(1, args.exact_radius + 1)) + tuple(
        range(B - args.exact_radius, B + 1)
    )
    exact = {weight: exact_slice(rows, weight) for weight in exact_weights}
    for weight, histogram in load_exact_slices(args.exact_slices).items():
        if weight in exact and exact[weight] != histogram:
            raise SystemExit("no-Singer kernel: regenerated and committed slices differ")
        exact[weight] = histogram
    exact_weights = tuple(sorted(exact))

    print("no-Singer BCH weight-conditioned kernel")
    print(f"exact_radius={args.exact_radius} pole={args.pole_num}/{args.pole_den}")
    for input_weight in exact_weights:
        row = exact[input_weight]
        population = math.comb(B, input_weight)
        mean_num = sum(weight * count for weight, count in row.items())
        low = sum(count for weight, count in row.items() if weight <= 30)
        print(
            f"slice t={input_weight:2d} count={population} min={min(row)} max={max(row)} "
            f"mean={mean_num / population:.9f} p_w_le_30={low / population:.12g}"
        )

    with localcontext() as context:
        context.prec = 80
        context.rounding = ROUND_CEILING
        z = Decimal(args.pole_num) / Decimal(args.pole_den)
        z_powers = tuple(z**power for power in range(N_LOCAL + 1))
        values = tuple(Decimal(1) for _ in range(B + 1))

        for step in range(1, args.steps + 1):
            scores = {
                weight: split_score(split[weight], z_powers, values) for weight in weights
            }
            next_values = mixed_row_values(scores, spectrum, exact, args.mixer)
            activation_weight = max(range(1, B + 1), key=next_values.__getitem__)
            activation = next_values[activation_weight]
            worst_state = max(range(1, B + 1), key=next_values.__getitem__)
            print(
                f"episode_steps={step:2d} activation_mgf_upper={activation:.18E} "
                f"log2_upper={log2_decimal(activation):.12f} "
                f"activation_t={activation_weight} worst_state_t={worst_state}"
            )
            values = tuple(next_values)

        turnoff_scores = {weight: turnoff_score(split[weight]) for weight in weights}
        turnoff_values = mixed_row_values(
            turnoff_scores, spectrum, exact, args.mixer
        )
        activation_turnoff = max(turnoff_values[1:])
        print(
            f"arbitrary_activation_turnoff_upper={activation_turnoff:.18E} "
            f"log2_upper={log2_decimal(activation_turnoff):.12f}"
        )
        for input_weight in range(1, min(16, B) + 1):
            turnoff = turnoff_values[input_weight]
            print(
                f"state_t={input_weight:2d} robust_turnoff_upper={turnoff:.18E} "
                f"log2_upper={log2_decimal(turnoff):.12f}"
            )

    rho, prefactor = exact_geometric_certificate(
        approximate=values,
        spectrum=spectrum,
        exact=exact,
        pole=Fraction(args.pole_num, args.pole_den),
        mixer=args.mixer,
    )
    print(
        f"exact_geometric_rho={rho.numerator}/{rho.denominator} "
        f"rho_log2={math.log2(rho.numerator) - math.log2(rho.denominator):.12f}"
    )
    print(
        f"exact_geometric_prefactor={prefactor.numerator}/{prefactor.denominator} "
        f"prefactor_log2={math.log2(prefactor.numerator) - math.log2(prefactor.denominator):.12f}"
    )
    crossing = exact_chernoff_crossing(
        rho,
        prefactor,
        Fraction(args.pole_num, args.pole_den),
        (9 * args.total_length) // 100,
    )
    print(f"exact_geometric_chernoff_crossing={crossing}")

    finite_termination, finite_worst_weight = exact_finite_termination_bound(
        spectrum=spectrum,
        exact=exact,
        mixer=args.mixer,
        remaining_blocks=args.remaining_blocks,
        remaining_ones=args.remaining_ones,
    )
    termination_log2 = (
        math.log2(finite_termination.numerator)
        - math.log2(finite_termination.denominator)
    )
    termination_cap_shift = math.floor(-termination_log2)
    while finite_termination.numerator << termination_cap_shift > finite_termination.denominator:
        termination_cap_shift -= 1
    while finite_termination.numerator << (termination_cap_shift + 1) <= finite_termination.denominator:
        termination_cap_shift += 1
    print(
        f"finite_termination_blocks={args.remaining_blocks} "
        f"remaining_ones_max={args.remaining_ones} "
        f"log2_upper={termination_log2:.12f} "
        f"dyadic_cap=2^-{termination_cap_shift} "
        f"worst_input_weight={finite_worst_weight} "
        f"numerator_bits={finite_termination.numerator.bit_length()} "
        f"denominator_bits={finite_termination.denominator.bit_length()}"
    )


if __name__ == "__main__":
    main()
