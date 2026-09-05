#!/usr/bin/env python3
"""Outward verifier for the dense one-shot random-constituent cover.

The diagnostic supplies only witnesses: a dyadic triangle partition, one
categorical reference at each accepted leaf, and point witnesses at the
remaining integer lattice points.  This program reconstructs the partition
exactly and evaluates every witness with Arb interval arithmetic.
"""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path

from flint import arb, ctx

from evaluate_single_random_constituent_highprob_renyi import spectrum_caps


WORKSTREAM = Path(__file__).resolve().parent
DIAGNOSTIC = WORKSTREAM / "single_random_constituent_B512_shared_two_band_dense_cover_s22.json"
OUTPUT = WORKSTREAM / "single_random_constituent_B512_shared_two_band_dense_cover_outward_s22.json"

B = 512
K = 256
L = 4096
N = 1 << 21
D = (109 * N + 999) // 1000
MEMORY = 22
MINIMUM_OCCUPATION = 160
POINTWISE_MARGIN = 80
DELTA = Fraction(1, 1 << 51)
P_DEFECT = Fraction(79, B)
P_CENTRAL = Fraction(1, 2)


def a(value: Fraction | int) -> arb:
    if isinstance(value, Fraction):
        return arb(value.numerator) / value.denominator
    return arb(value)


def fraction_of_float(value: float) -> Fraction:
    numerator, denominator = value.as_integer_ratio()
    return Fraction(numerator, denominator)


def upper_float(value: arb) -> float:
    if not value.is_finite():
        raise ArithmeticError("nonfinite Arb enclosure")
    return math.nextafter(float(value.upper()), math.inf)


def lower_float(value: arb) -> float:
    if not value.is_finite():
        raise ArithmeticError("nonfinite Arb enclosure")
    return math.nextafter(float(value.lower()), -math.inf)


def matrix_product(left: tuple[arb, arb, arb, arb], right: tuple[arb, arb, arb, arb]) -> tuple[arb, arb, arb, arb]:
    l00, l01, l10, l11 = left
    r00, r01, r10, r11 = right
    return (
        l00 * r00 + l01 * r10,
        l00 * r01 + l01 * r11,
        l10 * r00 + l11 * r10,
        l10 * r01 + l11 * r11,
    )


def matrix_power(matrix: tuple[arb, arb, arb, arb], exponent: int) -> tuple[arb, arb, arb, arb]:
    result = (arb(1), arb(0), arb(0), arb(1))
    base = matrix
    while exponent:
        if exponent & 1:
            result = matrix_product(result, base)
        exponent >>= 1
        if exponent:
            base = matrix_product(base, base)
    return result


def exact_caps() -> tuple[list[int], Fraction]:
    """Obtain cap witnesses, then verify their failure bound exactly."""
    approximate, _ = spectrum_caps(
        math.ldexp(1.0, -51), block_bits=B, dimension=K
    )
    caps = [int(value) for value in approximate]
    messages = (1 << K) - 1
    failure = Fraction(0)
    for weight in range(1, B + 1):
        mean = Fraction(messages * math.comb(B, weight), 1 << B)
        cap = caps[weight]
        if cap == 0:
            failure += mean
            continue
        deviation = Fraction(cap + 1) - mean
        if deviation <= 0:
            raise AssertionError("positive spectrum cap does not exceed its mean")
        failure += mean / (mean + deviation * deviation)
    failure += sum(Fraction(1 << index, 1 << B) for index in range(K))
    return caps, failure


def exact_band_majorant(caps: list[int], lower: int, upper: int, probability: Fraction) -> Fraction:
    result = Fraction(0)
    for weight in range(lower, upper + 1):
        cap = caps[weight]
        if not cap:
            continue
        shell_probability = (
            Fraction(math.comb(B, weight))
            * probability**weight
            * (1 - probability) ** (B - weight)
        )
        result = max(result, Fraction(cap) / shell_probability)
    if result <= 0:
        raise AssertionError("empty spectrum band")
    return result


def normalized_probabilities(values: list[float]) -> tuple[Fraction, Fraction, Fraction]:
    weights = [fraction_of_float(float(value)) for value in values]
    total = sum(weights, Fraction(0))
    if total <= 0:
        raise AssertionError("categorical witness has zero mass")
    result = tuple(value / total for value in weights)
    if any(value < 0 for value in result):
        raise AssertionError("categorical witness has negative mass")
    return result  # type: ignore[return-value]


def log_multinomial(counts: tuple[Fraction, Fraction, Fraction]) -> arb:
    total = sum(counts, Fraction(0))
    return (a(total) + 1).lgamma() - sum(
        ((a(count) + 1).lgamma() for count in counts), arb(0)
    )


def inner_log_moment(bit_probability: Fraction, surprisal: Fraction) -> arb:
    z = (-a(surprisal)).exp()
    b = (1 + z) / 2
    q = arb(1) / (1 << MEMORY)
    terminate = q * b
    survive = (1 - q) * b
    p = a(bit_probability)
    zero = (arb(1), arb(0), terminate, survive)
    active = (terminate, survive, terminate, survive)
    mixed = tuple((1 - p) * zero[i] + p * active[i] for i in range(4))
    powered = matrix_power(mixed, N)
    return (powered[0] + powered[1]).log()


def evaluate_witness(
    counts: tuple[Fraction, Fraction, Fraction],
    probabilities: tuple[Fraction, Fraction, Fraction],
    surprisal: Fraction,
    log_majorants: tuple[arb, arb, arb],
) -> tuple[arb, arb]:
    if sum(counts, Fraction(0)) != L:
        raise AssertionError("composition does not sum to the outer-row count")
    if any(count < 0 for count in counts):
        raise AssertionError("negative composition coordinate")
    log_type = log_multinomial(counts)
    outer = log_type + sum(
        (a(counts[index]) * log_majorants[index] for index in range(3)),
        arb(0),
    )
    log_conditioning = log_type
    for count, probability in zip(counts, probabilities):
        if count == 0:
            continue
        if probability <= 0:
            raise AssertionError("positive type count has zero reference probability")
        log_conditioning += a(count) * a(probability).log()
    bit_probability = probabilities[1] * P_DEFECT + probabilities[2] * P_CENTRAL
    raw = (
        inner_log_moment(bit_probability, surprisal)
        + D * a(surprisal)
        - B * log_conditioning
    )
    return outer + raw, raw


Point = tuple[Fraction, Fraction]
Triangle = tuple[Point, Point, Point]


def point(value: list[float] | tuple[float, float]) -> Point:
    return fraction_of_float(float(value[0])), fraction_of_float(float(value[1]))


def triangle(value: list[list[float]]) -> Triangle:
    return point(value[0]), point(value[1]), point(value[2])


def canonical(cell: Triangle) -> tuple[Point, Point, Point]:
    return tuple(sorted(cell))  # type: ignore[return-value]


def split_longest(cell: Triangle) -> tuple[Triangle, Triangle]:
    vertices = list(cell)
    edges = ((0, 1), (1, 2), (2, 0))
    first, second = max(
        edges,
        key=lambda pair: (
            (vertices[pair[0]][0] - vertices[pair[1]][0]) ** 2
            + (vertices[pair[0]][1] - vertices[pair[1]][1]) ** 2
        ),
    )
    third = 3 - first - second
    middle = (
        (vertices[first][0] + vertices[second][0]) / 2,
        (vertices[first][1] + vertices[second][1]) / 2,
    )
    return (
        (vertices[first], middle, vertices[third]),
        (middle, vertices[second], vertices[third]),
    )


def reconstruct_partition(
    accepted: set[tuple[Point, Point, Point]],
    unresolved: set[tuple[Point, Point, Point]],
) -> None:
    roots: list[Triangle] = [
        ((Fraction(0), Fraction(0)), (Fraction(0), Fraction(L)), (Fraction(L - MINIMUM_OCCUPATION), Fraction(MINIMUM_OCCUPATION))),
        ((Fraction(0), Fraction(0)), (Fraction(L - MINIMUM_OCCUPATION), Fraction(MINIMUM_OCCUPATION)), (Fraction(L - MINIMUM_OCCUPATION), Fraction(0))),
    ]
    visited: set[tuple[Point, Point, Point]] = set()
    pending = list(roots)
    while pending:
        cell = pending.pop()
        key = canonical(cell)
        if key in accepted or key in unresolved:
            if key in visited:
                raise AssertionError("duplicate leaf in dense partition")
            visited.add(key)
            continue
        left, right = split_longest(cell)
        pending.extend((right, left))
        if len(pending) + len(visited) > 100000:
            raise AssertionError("receipt does not describe the deterministic partition")
    expected = accepted | unresolved
    if visited != expected:
        raise AssertionError("dense leaves do not exactly reconstruct the domain")


def point_in_triangle(lattice_point: tuple[int, int], cell: Triangle) -> bool:
    x, y = map(Fraction, lattice_point)
    (x1, y1), (x2, y2), (x3, y3) = cell
    denominator = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)
    if denominator == 0:
        return False
    first = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / denominator
    second = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / denominator
    third = 1 - first - second
    return min(first, second, third) >= 0


def lattice_points(cell: Triangle) -> set[tuple[int, int]]:
    minimum_x = math.ceil(min(vertex[0] for vertex in cell))
    maximum_x = math.floor(max(vertex[0] for vertex in cell))
    minimum_y = math.ceil(min(vertex[1] for vertex in cell))
    maximum_y = math.floor(max(vertex[1] for vertex in cell))
    return {
        (x, y)
        for x in range(minimum_x, maximum_x + 1)
        for y in range(minimum_y, maximum_y + 1)
        if point_in_triangle((x, y), cell)
    }


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    ctx.prec = 192
    receipt = json.loads(DIAGNOSTIC.read_text(encoding="utf-8"))
    expected = {
        "outer_rows": L,
        "distance_cutoff": D,
        "memory_bits": MEMORY,
        "minimum_occupation": MINIMUM_OCCUPATION,
        "target_pointwise_margin_bits": float(POINTWISE_MARGIN),
    }
    for name, value in expected.items():
        if receipt["parameters"].get(name) != value:
            raise AssertionError(f"dense diagnostic parameter mismatch: {name}")

    caps, setup_failure = exact_caps()
    support = [weight for weight in range(1, B + 1) if caps[weight]]
    if (min(support), max(support)) != (42, 470):
        raise AssertionError("unexpected nonzero-cap support")
    low = exact_band_majorant(caps, 42, 79, P_DEFECT)
    high = exact_band_majorant(caps, 433, 470, 1 - P_DEFECT)
    central = exact_band_majorant(caps, 80, 432, P_CENTRAL)
    defect = 2 * max(low, high)
    log_majorants = (arb(0), a(defect).log(), a(central).log())

    accepted_rows = receipt["accepted_cells"]
    unresolved_rows = receipt["unresolved_cells"]
    accepted_cells = {canonical(triangle(row["vertices"])) for row in accepted_rows}
    unresolved_cells = {canonical(triangle(row["vertices"])) for row in unresolved_rows}
    if len(accepted_cells) != len(accepted_rows) or len(unresolved_cells) != len(unresolved_rows):
        raise AssertionError("duplicate dense cell")
    if accepted_cells & unresolved_cells:
        raise AssertionError("cell is both accepted and unresolved")
    reconstruct_partition(accepted_cells, unresolved_cells)

    worst_vertex = -math.inf
    for index, row in enumerate(accepted_rows):
        probabilities = normalized_probabilities(row["witness_probabilities"])
        surprisal = fraction_of_float(float(row["witness_surprisal"]))
        for vertex in triangle(row["vertices"]):
            inactive, defects = vertex
            counts = (inactive, defects, Fraction(L) - inactive - defects)
            value, raw = evaluate_witness(counts, probabilities, surprisal, log_majorants)
            if upper_float(raw) >= 0:
                raise AssertionError(f"accepted cell {index} lacks a negative raw exponent")
            shifted = value + POINTWISE_MARGIN * arb(2).log()
            if upper_float(shifted) > 0:
                raise AssertionError(f"accepted cell {index} misses the pointwise margin")
            worst_vertex = max(worst_vertex, upper_float(value / arb(2).log()))
        if (index + 1) % 100 == 0:
            print(f"accepted_cells,{index + 1},{len(accepted_rows)}", flush=True)

    exact_unresolved_points = set().union(*(lattice_points(cell) for cell in unresolved_cells))
    exception_rows = receipt["lattice_exceptions"]
    stated_points = {
        (int(row["counts"][0]), int(row["counts"][1])) for row in exception_rows
    }
    if exact_unresolved_points != stated_points:
        raise AssertionError("lattice exceptions do not exactly cover unresolved cells")
    worst_exception = -math.inf
    for index, row in enumerate(exception_rows):
        counts_int = tuple(int(value) for value in row["counts"])
        counts = tuple(Fraction(value) for value in counts_int)
        probabilities = normalized_probabilities(row["reference_type_probabilities"])
        surprisal = fraction_of_float(float(row["surprisal"]))
        value, raw = evaluate_witness(counts, probabilities, surprisal, log_majorants)
        if upper_float(raw) >= 0:
            raise AssertionError(f"lattice exception {index} lacks a negative raw exponent")
        shifted = value + POINTWISE_MARGIN * arb(2).log()
        if upper_float(shifted) > 0:
            raise AssertionError(f"lattice exception {index} misses the pointwise margin")
        worst_exception = max(worst_exception, upper_float(value / arb(2).log()))

    composition_count = sum(occupation + 1 for occupation in range(MINIMUM_OCCUPATION, L + 1))
    dense_log2_upper = upper_float(
        arb(composition_count).log() / arb(2).log() - POINTWISE_MARGIN
    )
    if dense_log2_upper >= -40:
        raise AssertionError("dense aggregate does not close at 40 bits")
    payload = {
        "schema": "single-random-constituent-shared-two-band-dense-outward-v1",
        "status": "OUTWARD_DENSE_CERTIFICATE",
        "claim": {
            "minimum_occupation": MINIMUM_OCCUPATION,
            "maximum_occupation": L,
            "pointwise_log2_upper": -POINTWISE_MARGIN,
            "composition_count": composition_count,
            "aggregate_log2_upper": dense_log2_upper,
            "aggregate_margin_bits_lower": -dense_log2_upper,
            "accepted_cells": len(accepted_rows),
            "unresolved_cells": len(unresolved_rows),
            "lattice_exceptions": len(exception_rows),
            "complete_integer_lattice": True,
        },
        "outward_checks": {
            "worst_accepted_vertex_log2_upper": worst_vertex,
            "worst_lattice_exception_log2_upper": worst_exception,
            "partition_reconstructed_exactly": True,
            "unresolved_lattice_reconstructed_exactly": True,
            "arithmetic": "192-bit Arb intervals; exact rational caps, probabilities, band majorants, and dyadic geometry",
        },
        "constituent_event": {
            "total_failure_upper": upper_float(a(setup_failure)),
            "margin_bits_lower": -upper_float(a(setup_failure).log() / arb(2).log()),
            "comparison_to_2^-40": setup_failure < Fraction(1, 1 << 40),
        },
        "parameters": {**expected, "output_bits": N, "block_bits": B, "dimension": K},
        "witness_source": {
            "file": DIAGNOSTIC.name,
            "sha256": file_sha256(DIAGNOSTIC),
        },
        "proof_rule": (
            "For a fixed categorical witness, both the unclipped reference exponent and, where that exponent is negative, the final log bound are convex in the three type counts. Vertex bounds therefore cover each accepted triangle. The listed point witnesses cover every integer composition in each unresolved triangle."
        ),
        "limitations": [
            "This receipt covers occupations 160 through 4096; sparse occupations use separate transfer certificates.",
            "The inner ensemble is RandomStepConv-M22, not RM2Sub.",
        ],
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(json.dumps(payload["constituent_event"], indent=2))
    print(f"output={OUTPUT}")


if __name__ == "__main__":
    main()
