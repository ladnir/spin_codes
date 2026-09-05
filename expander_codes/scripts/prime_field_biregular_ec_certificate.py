#!/usr/bin/env python3
"""Generate and verify the finite two-sided regular prime-field EC result.

Floating point selects fixed positive markers.  Verification reparses those
markers as decimal Arb balls and recomputes every bound with outward rounding.
The default certificate has rate one half, p = 2^127-1, length near 2^21,
two-sided degree 28, convolution memory three, and the floored p-ary GV
cutoff.  Its transfer counts every fresh zero equation before applying the
projective union bound.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from flint import arb, arb_mat, arb_poly, ctx

from ea_certificate import decimal_marker
from prime_field_biregular_ec_diagnostic import (
    constraint_balanced_ec_block_logbound,
    constraint_full_support_ec_logbound,
    candidate_parameters,
)
from prime_field_ea_diagnostic import gv_distance
from regular_ec_certificate import (
    combine_uniform_slice_transfers_arb,
    identity_matrix,
    zero_matrix,
)


SCHEMA = "prime-field-biregular-ec-constraint-trace-v2"
DEFAULT_PRIME = 2**127 - 1
DEFAULT_DEGREE = 28
DEFAULT_MEMORY = 3
DEFAULT_K, DEFAULT_N, DEFAULT_REGION_LENGTH = candidate_parameters(
    2**21, DEFAULT_DEGREE
)
DEFAULT_CUTOFF = math.floor(DEFAULT_N * gv_distance(DEFAULT_PRIME, 0.5))
DEFAULT_EXACT_BANDS = (
    (1, 8),
    (9, 16),
    (17, 32),
    (33, 64),
    (65, 96),
    (97, 128),
    (129, 160),
    (161, 192),
)
DEFAULT_EXACT_MARKERS = {
    (1, 8): {
        "structural": ("0.9991447773967677", "5.839835298989263e-22"),
        "field": ("0.9991692212810199", "2.5166797721728796e-21"),
        "diagnostic_log2": -615.5275492262797,
    },
    (9, 16): {
        "structural": ("0.9967622332862647", "2.4826055705733867e-55"),
        "field": ("0.9975785283404316", "5.8774717541116e-39"),
        "diagnostic_log2": -1340.2277122918545,
    },
    (17, 32): {
        "structural": ("0.9985604162688718", "3.6970904119057333e-13"),
        "field": ("0.998592495268847", "9.185341959987983e-13"),
        "diagnostic_log2": -335.24393054943755,
    },
    (33, 64): {
        "structural": ("0.9985527284526642", "9.420467108948146e-6"),
        "field": ("0.9985564845711785", "1.0295943023022584e-5"),
        "diagnostic_log2": -132.70375464543548,
    },
    (65, 96): {
        "structural": ("0.9979801845055196", "6.855802703817471e-5"),
        "field": ("0.9979785655536838", "6.945969394666614e-5"),
        "diagnostic_log2": -72.091,
    },
    (97, 128): {
        "structural": ("0.9972801872091132", "9.252994173875877e-5"),
        "field": ("0.997278476458764", "9.36589629959553e-5"),
        "diagnostic_log2": -95.593,
    },
    (129, 160): {
        "structural": ("0.9963851659450441", "1.229659570575988e-4"),
        "field": ("0.9963794932070447", "1.2426853116162396e-4"),
        "diagnostic_log2": -218.304,
    },
    (161, 192): {
        "structural": ("0.9954851021587344", "1.5392000341824723e-4"),
        "field": ("0.995488652002496", "1.5441616026343332e-4"),
        "diagnostic_log2": -347.506,
    },
}
DEFAULT_BLOCKS = (
    (193, 256),
    (257, 400),
    (401, 800),
    (801, 1_600),
    (1_601, 3_200),
    (3_201, 6_400),
    (6_401, 12_800),
    (12_801, 25_600),
    (25_601, 51_200),
    (51_201, 102_400),
    (102_401, 204_800),
    (204_801, 409_600),
    (409_601, 819_200),
    (819_201, 917_499),
    (917_500, 983_035),
    (983_036, 1_015_803),
    (1_015_804, 1_032_187),
    (1_032_188, 1_040_379),
    (1_040_380, 1_044_475),
    (1_044_476, 1_046_523),
    (1_046_524, 1_047_547),
    (1_047_548, 1_048_059),
    (1_048_060, 1_048_315),
    (1_048_316, 1_048_443),
    (1_048_444, 1_048_507),
    (1_048_508, 1_048_539),
    (1_048_540, 1_048_555),
    (1_048_556, 1_048_563),
    (1_048_564, 1_048_567),
    (1_048_568, 1_048_569),
    (1_048_570, 1_048_570),
    (1_048_571, 1_048_571),
)


def nonwrapping_trace_matrices_arb(
    *, memory: int, output_marker: arb, equation_marker: arb
) -> tuple[arb_mat, arb_mat]:
    """Outward-rounded transfers counting every fresh zero equation."""
    size = memory + 1
    empty = [[arb(0) for _ in range(size)] for _ in range(size)]
    occupied = [[arb(0) for _ in range(size)] for _ in range(size)]
    for state in range(memory):
        empty[state][0] = occupied[state][0] = output_marker
        empty[state][state + 1] = occupied[state][state + 1] = equation_marker
    empty[memory][memory] = arb(1)
    occupied[memory][0] = output_marker
    occupied[memory][memory] = equation_marker
    return arb_mat(empty), arb_mat(occupied)


def nonwrapping_singleton_trace_matrices_arb(
    *, memory: int, output_marker: arb, equation_marker: arb
) -> tuple[arb_mat, arb_mat, arb_mat]:
    """Outward-rounded empty, singleton, and collision transfers."""
    empty, collision = nonwrapping_trace_matrices_arb(
        memory=memory,
        output_marker=output_marker,
        equation_marker=equation_marker,
    )
    singleton = arb_mat(collision)
    singleton[memory, memory] = arb(0)
    return empty, singleton, collision


def zero_poly_matrix(size: int) -> list[list[arb_poly]]:
    return [[arb_poly() for _ in range(size)] for _ in range(size)]


def identity_poly_matrix(size: int) -> list[list[arb_poly]]:
    result = zero_poly_matrix(size)
    for index in range(size):
        result[index][index] = arb_poly([1])
    return result


def multiply_poly_matrices_truncated(
    left: list[list[arb_poly]], right: list[list[arb_poly]], order: int
) -> list[list[arb_poly]]:
    """Multiply small Arb polynomial matrices modulo X^order."""
    size = len(left)
    result = zero_poly_matrix(size)
    for row in range(size):
        for inner in range(size):
            if not left[row][inner]:
                continue
            for column in range(size):
                if right[inner][column]:
                    result[row][column] += (
                        left[row][inner] * right[inner][column]
                    ).truncate(order)
        for column in range(size):
            result[row][column] = result[row][column].truncate(order)
    return result


def power_poly_matrix_truncated(
    matrix: list[list[arb_poly]], exponent: int, order: int
) -> list[list[arb_poly]]:
    """Binary exponentiation of a polynomial matrix modulo X^order."""
    result = identity_poly_matrix(len(matrix))
    power = matrix
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = multiply_poly_matrices_truncated(result, power, order)
        remaining >>= 1
        if remaining:
            power = multiply_poly_matrices_truncated(power, power, order)
    return result


def singleton_region_polynomial_arb(
    *, region_length: int, group_size: int, max_weight: int, memory: int,
    output_marker: arb, equation_marker: arb,
) -> list[list[arb_poly]]:
    """Return the truncated singleton-refined regional numerator."""
    empty, singleton, collision = nonwrapping_singleton_trace_matrices_arb(
        memory=memory,
        output_marker=output_marker,
        equation_marker=equation_marker,
    )
    size = memory + 1
    group = zero_poly_matrix(size)
    for row in range(size):
        for column in range(size):
            coefficients = [empty[row, column]]
            coefficients.append(singleton[row, column] * group_size)
            coefficients.extend(
                collision[row, column] * math.comb(group_size, weight)
                for weight in range(2, group_size + 1)
            )
            group[row][column] = arb_poly(coefficients)
    return power_poly_matrix_truncated(
        group, region_length, max_weight + 1
    )


def singleton_exact_band_term_arb(
    *, prime: int, k: int, n: int, cutoff: int, region_count: int,
    memory: int, lo: int, hi: int, markers: dict[str, list[str]],
) -> arb:
    """Rigorous singleton-refined exact regional coefficient band."""
    return singleton_exact_band_branch_term_arb(
        prime=prime,
        k=k,
        n=n,
        cutoff=cutoff,
        region_count=region_count,
        memory=memory,
        lo=lo,
        hi=hi,
        values=markers["structural"],
        field=False,
    ) + singleton_exact_band_branch_term_arb(
        prime=prime,
        k=k,
        n=n,
        cutoff=cutoff,
        region_count=region_count,
        memory=memory,
        lo=lo,
        hi=hi,
        values=markers["field"],
        field=True,
    )


def singleton_exact_band_branch_term_arb(
    *, prime: int, k: int, n: int, cutoff: int, region_count: int,
    memory: int, lo: int, hi: int, values: list[str], field: bool,
) -> arb:
    """One projective-cap branch of an exact singleton support band."""
    region_length = n // region_count
    group_size = k // region_length
    s = arb(prime - 1)
    z, v = map(arb, values)
    if not (z > 0 and z <= 1 and v > 0 and v <= 1):
        raise ValueError(f"invalid exact-band markers on [{lo},{hi}]")
    if field and not v >= arb(1) / (prime - 1):
        raise ValueError(f"field marker below 1/(p-1) on [{lo},{hi}]")
    numerator = singleton_region_polynomial_arb(
        region_length=region_length,
        group_size=group_size,
        max_weight=hi,
        memory=memory,
        output_marker=z,
        equation_marker=v,
    )
    total = arb(0)
    output_factor = z ** (-cutoff)
    for r in range(lo, hi + 1):
        normalizer = arb(math.comb(k, r))
        region = arb_mat([
            [numerator[row][column][r] / normalizer
             for column in range(memory + 1)]
            for row in range(memory + 1)
        ])
        mgf = endpoint_mass(region, region_count, memory)
        cap = v ** (-r) / s if field else v ** (-(r - 1))
        total += normalizer * mgf * output_factor * cap
    return total


def balanced_occupancy_prefix_arb(
    *, region_length: int, group_size: int, limit: int
) -> list[list[arb]]:
    total_slots = region_length * group_size
    if not 0 <= limit <= total_slots:
        raise ValueError("invalid occupancy prefix limit")
    laws = [[arb(1)]]
    for used in range(limit):
        current = laws[-1]
        following = [arb(0) for _ in range(len(current) + 1)]
        remaining = total_slots - used
        for occupied, probability in enumerate(current):
            following[occupied] += (
                probability * (group_size * occupied - used) / remaining
            )
            following[occupied + 1] += (
                probability
                * group_size
                * (region_length - occupied)
                / remaining
            )
        laws.append(following)
    return laws


def uniform_occupancy_trace_transfers_arb(
    *, region_length: int, max_occupied: int, memory: int,
    output_marker: arb, equation_marker: arb,
) -> list[arb_mat]:
    empty, occupied = nonwrapping_trace_matrices_arb(
        memory=memory,
        output_marker=output_marker,
        equation_marker=equation_marker,
    )
    result: tuple[int, list[arb_mat]] = (
        0,
        [identity_matrix(memory + 1)],
    )
    power: tuple[int, list[arb_mat]] = (1, [empty, occupied])
    remaining = region_length
    while remaining:
        if remaining & 1:
            result = combine_uniform_slice_transfers_arb(
                result, power, max_occupied
            )
        remaining >>= 1
        if remaining:
            power = combine_uniform_slice_transfers_arb(
                power, power, max_occupied
            )
    return result[1]


def endpoint_mass(matrix: arb_mat, exponent: int, start_state: int) -> arb:
    powered = matrix**exponent
    return sum(
        (powered[start_state, column] for column in range(powered.ncols())),
        arb(0),
    )


def exact_band_term_arb(
    *, prime: int, k: int, n: int, cutoff: int, region_count: int,
    memory: int, lo: int, hi: int, markers: dict[str, list[str]],
) -> arb:
    region_length = n // region_count
    group_size = k // region_length
    laws = balanced_occupancy_prefix_arb(
        region_length=region_length,
        group_size=group_size,
        limit=hi,
    )
    s = arb(prime - 1)

    def branch(values: list[str], field: bool) -> arb:
        z, v = map(arb, values)
        if not (z > 0 and z <= 1 and v > 0 and v <= 1):
            raise ValueError(f"invalid exact-band markers on [{lo},{hi}]")
        if field and not v >= arb(1) / (prime - 1):
            raise ValueError(f"field marker below 1/(p-1) on [{lo},{hi}]")
        slices = uniform_occupancy_trace_transfers_arb(
            region_length=region_length,
            max_occupied=hi,
            memory=memory,
            output_marker=z,
            equation_marker=v,
        )
        total = arb(0)
        output_factor = z ** (-cutoff)
        for r in range(lo, hi + 1):
            region = zero_matrix(memory + 1)
            for occupied, probability in enumerate(laws[r]):
                if not probability.is_zero():
                    region += slices[occupied] * probability
            mgf = endpoint_mass(region, region_count, memory)
            cap = v ** (-r) / s if field else v ** (-(r - 1))
            total += arb(math.comb(k, r)) * mgf * output_factor * cap
        return total

    return branch(markers["structural"], False) + branch(
        markers["field"], True
    )


def block_term_arb(
    *, prime: int, k: int, n: int, cutoff: int, region_count: int,
    memory: int, lo: int, hi: int,
    markers: dict[str, list[str]],
) -> arb:
    group_size = k // (n // region_count)
    common = arb(hi - lo + 1)

    def branch(values: list[str], field: bool) -> arb:
        z, v, x = map(arb, values)
        if not (x > 0 and z > 0 and z <= 1 and v > 0 and v <= 1):
            raise ValueError(f"invalid block markers on [{lo},{hi}]")
        if field and not v >= arb(1) / (prime - 1):
            raise ValueError(f"field marker below 1/(p-1) on [{lo},{hi}]")
        empty, occupied = nonwrapping_trace_matrices_arb(
            memory=memory,
            output_marker=z,
            equation_marker=v,
        )
        matrix = empty + occupied * ((1 + x) ** group_size - 1)
        mgf = endpoint_mass(matrix, n, memory)
        b = x**region_count * v
        low = arb(math.comb(k, lo)) ** (1 - region_count) * b ** (-lo)
        if lo == hi:
            endpoint = low
        else:
            high = (
                arb(math.comb(k, hi)) ** (1 - region_count)
                * b ** (-hi)
            )
            endpoint = low + high
        result = common * endpoint * mgf * z ** (-cutoff)
        return result / (prime - 1) if field else result * v

    return branch(markers["structural"], False) + branch(
        markers["field"], True
    )


def singleton_block_term_arb(
    *, prime: int, k: int, n: int, cutoff: int, region_count: int,
    memory: int, lo: int, hi: int,
    markers: dict[str, list[str]],
) -> arb:
    """Rigorous singleton-refined positive-saddle block bound."""
    group_size = k // (n // region_count)
    common = arb(hi - lo + 1)

    def branch(values: list[str], field: bool) -> arb:
        z, v, x = map(arb, values)
        if not (x > 0 and z > 0 and z <= 1 and v > 0 and v <= 1):
            raise ValueError(f"invalid block markers on [{lo},{hi}]")
        if field and not v >= arb(1) / (prime - 1):
            raise ValueError(f"field marker below 1/(p-1) on [{lo},{hi}]")
        empty, singleton, collision = nonwrapping_singleton_trace_matrices_arb(
            memory=memory,
            output_marker=z,
            equation_marker=v,
        )
        singleton_factor = group_size * x
        collision_factor = (1 + x) ** group_size - 1 - singleton_factor
        matrix = (
            empty
            + singleton * singleton_factor
            + collision * collision_factor
        )
        mgf = endpoint_mass(matrix, n, memory)
        b = x**region_count * v
        low = arb(math.comb(k, lo)) ** (1 - region_count) * b ** (-lo)
        if lo == hi:
            endpoint = low
        else:
            high = (
                arb(math.comb(k, hi)) ** (1 - region_count)
                * b ** (-hi)
            )
            endpoint = low + high
        result = common * endpoint * mgf * z ** (-cutoff)
        return result / (prime - 1) if field else result * v

    return branch(markers["structural"], False) + branch(
        markers["field"], True
    )


def full_support_markers_float(
    *, prime: int, n: int, cutoff: int, memory: int, message_weight: int
) -> dict[str, list[str]]:
    bound = constraint_full_support_ec_logbound(
        prime=prime,
        n=n,
        cutoff=cutoff,
        memory=memory,
        message_weight=message_weight,
    )
    z_a, v_a = bound.structural_markers
    z_b, v_b = bound.field_markers
    v_b = max(v_b, 5.8774717541116e-39)
    return {
        "structural": [decimal_marker(z_a), decimal_marker(v_a)],
        "field": [decimal_marker(z_b), decimal_marker(v_b)],
    }


def full_support_term_arb(
    *, prime: int, n: int, cutoff: int, memory: int, message_weight: int,
    markers: dict[str, list[str]],
) -> arb:
    r = message_weight

    def branch(values: list[str], field: bool) -> arb:
        z, v = map(arb, values)
        if field and not v >= arb(1) / (prime - 1):
            raise ValueError("full-support field marker is below 1/(p-1)")
        _, occupied = nonwrapping_trace_matrices_arb(
            memory=memory,
            output_marker=z,
            equation_marker=v,
        )
        mgf = endpoint_mass(occupied, n, memory)
        result = mgf * z ** (-cutoff)
        if field:
            return result * v ** (-r) / (prime - 1)
        return result * v ** (-(r - 1))

    return branch(markers["structural"], False) + branch(
        markers["field"], True
    )


def generate_certificate(
    *, prime: int, k: int, n: int, cutoff: int, region_count: int,
    memory: int, target_bits: int, precision_bits: int,
) -> dict[str, Any]:
    if (prime, k, n, cutoff, region_count, memory) != (
        DEFAULT_PRIME,
        DEFAULT_K,
        DEFAULT_N,
        DEFAULT_CUTOFF,
        DEFAULT_DEGREE,
        DEFAULT_MEMORY,
    ):
        raise ValueError("the optimized support partition is for the defaults")

    exact_bands = []
    for lo, hi in DEFAULT_EXACT_BANDS:
        selected = DEFAULT_EXACT_MARKERS[(lo, hi)]
        exact_bands.append({
            "lo": lo,
            "hi": hi,
            "markers": {
                "structural": list(selected["structural"]),
                "field": list(selected["field"]),
            },
            "diagnostic_log2": selected["diagnostic_log2"],
        })
        print(f"optimized exact band [{lo},{hi}]", flush=True)

    blocks = []
    for lo, hi in DEFAULT_BLOCKS:
        bound = constraint_balanced_ec_block_logbound(
            prime=prime,
            k=k,
            n=n,
            cutoff=cutoff,
            region_count=region_count,
            memory=memory,
            support_start=lo,
            support_limit=hi,
        )
        markers = {
            "structural": [
                decimal_marker(value) for value in bound.structural_markers
            ],
            "field": [
                decimal_marker(value) for value in bound.field_markers
            ],
        }
        blocks.append({
            "lo": lo,
            "hi": hi,
            "markers": markers,
            "diagnostic_log2": bound.total_log2,
        })
        print(f"optimized block [{lo},{hi}]", flush=True)

    return {
        "schema": SCHEMA,
        "parameters": {
            "prime": str(prime),
            "k": k,
            "n": n,
            "cutoff": cutoff,
            "region_count": region_count,
            "memory": memory,
            "target_bits": target_bits,
        },
        "verification": {"precision_bits": precision_bits},
        "exact_bands": exact_bands,
        "blocks": blocks,
        "full_support": {
            "r": k,
            "markers": full_support_markers_float(
                prime=prime,
                n=n,
                cutoff=cutoff,
                memory=memory,
                message_weight=k,
            ),
        },
    }


def validate_coverage(certificate: dict[str, Any]) -> None:
    k = int(certificate["parameters"]["k"])
    expected = 1
    for section in ("exact_bands", "blocks"):
        for item in certificate[section]:
            lo, hi = int(item["lo"]), int(item["hi"])
            if lo != expected or hi < lo:
                raise ValueError(f"coverage breaks at r={expected}")
            expected = hi + 1
    if expected != k:
        raise ValueError(f"coverage stops at r={expected - 1}")
    if int(certificate["full_support"]["r"]) != k:
        raise ValueError("full-support entry is missing")


@dataclass(frozen=True)
class VerificationResult:
    success: bool
    total_bound: arb
    security_bits: arb
    exact_bound: arb
    block_bound: arb
    full_support_bound: arb
    largest_exact_band: tuple[int, int]
    largest_exact_term: arb
    largest_block: tuple[int, int]
    largest_block_term: arb


def verify_certificate(certificate: dict[str, Any]) -> VerificationResult:
    if certificate.get("schema") != SCHEMA:
        raise ValueError("unsupported certificate schema")
    validate_coverage(certificate)
    parameters = certificate["parameters"]
    prime = int(parameters["prime"])
    k, n = int(parameters["k"]), int(parameters["n"])
    cutoff = int(parameters["cutoff"])
    region_count = int(parameters["region_count"])
    memory = int(parameters["memory"])
    target_bits = int(parameters["target_bits"])
    ctx.prec = int(certificate["verification"]["precision_bits"])

    exact_total = arb(0)
    largest_exact_band, largest_exact_term = (0, 0), arb(0)
    for item in certificate["exact_bands"]:
        lo, hi = int(item["lo"]), int(item["hi"])
        term = exact_band_term_arb(
            prime=prime,
            k=k,
            n=n,
            cutoff=cutoff,
            region_count=region_count,
            memory=memory,
            lo=lo,
            hi=hi,
            markers=item["markers"],
        )
        exact_total += term
        if largest_exact_band == (0, 0) or term > largest_exact_term:
            largest_exact_band, largest_exact_term = (lo, hi), term
        print(f"verified exact band [{lo},{hi}]", flush=True)

    block_total = arb(0)
    largest_block, largest_block_term = (0, 0), arb(0)
    for item in certificate["blocks"]:
        lo, hi = int(item["lo"]), int(item["hi"])
        term = block_term_arb(
            prime=prime,
            k=k,
            n=n,
            cutoff=cutoff,
            region_count=region_count,
            memory=memory,
            lo=lo,
            hi=hi,
            markers=item["markers"],
        )
        block_total += term
        if largest_block == (0, 0) or term > largest_block_term:
            largest_block, largest_block_term = (lo, hi), term
        print(f"verified block [{lo},{hi}]", flush=True)

    full = certificate["full_support"]
    full_term = full_support_term_arb(
        prime=prime,
        n=n,
        cutoff=cutoff,
        memory=memory,
        message_weight=k,
        markers=full["markers"],
    )
    total = exact_total + block_total + full_term
    security_bits = -total.log() / arb(2).log()
    return VerificationResult(
        success=bool(total < arb(2) ** (-target_bits)),
        total_bound=total,
        security_bits=security_bits,
        exact_bound=exact_total,
        block_bound=block_total,
        full_support_bound=full_term,
        largest_exact_band=largest_exact_band,
        largest_exact_term=largest_exact_term,
        largest_block=largest_block,
        largest_block_term=largest_block_term,
    )


def result_json(result: VerificationResult) -> dict[str, Any]:
    return {
        "success": result.success,
        "total_bound": str(result.total_bound),
        "security_bits": str(result.security_bits),
        "exact_bound": str(result.exact_bound),
        "block_bound": str(result.block_bound),
        "full_support_bound": str(result.full_support_bound),
        "largest_exact_band": {
            "lo": result.largest_exact_band[0],
            "hi": result.largest_exact_band[1],
            "bound": str(result.largest_exact_term),
        },
        "largest_block": {
            "lo": result.largest_block[0],
            "hi": result.largest_block[1],
            "bound": str(result.largest_block_term),
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    generate = subparsers.add_parser("generate")
    generate.add_argument("--output", type=Path, required=True)
    generate.add_argument("--target-bits", type=int, default=20)
    generate.add_argument("--precision-bits", type=int, default=256)
    verify = subparsers.add_parser("verify")
    verify.add_argument("certificate", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.command == "generate":
        certificate = generate_certificate(
            prime=DEFAULT_PRIME,
            k=DEFAULT_K,
            n=DEFAULT_N,
            cutoff=DEFAULT_CUTOFF,
            region_count=DEFAULT_DEGREE,
            memory=DEFAULT_MEMORY,
            target_bits=args.target_bits,
            precision_bits=args.precision_bits,
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(certificate, indent=2) + "\n", encoding="utf-8"
        )
        print(f"wrote {args.output}")
        return
    certificate = json.loads(args.certificate.read_text(encoding="utf-8"))
    result = verify_certificate(certificate)
    print(json.dumps(result_json(result), indent=2))
    if not result.success:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
