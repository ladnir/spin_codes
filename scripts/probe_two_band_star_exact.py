#!/usr/bin/env python3
"""Exact-joint support moments for the dominant two-band star patterns.

The product-marginal relaxation in probe_two_band_spreading.py is dominated
by partition pairs in which all active blocks share one left tile and the
right tiles have degree at most two.  This script retains the exact joint
64/64 weight law J(a,b) of every permuted EBCH word.  Conditional on (a,b),
the two half supports are uniform and independent, so hypergeometric union
transitions give the exact support moment for a fixed star pattern.

Long-double arithmetic is used throughout.  This is a design diagnostic, not
yet an outward-rounded certificate.
"""

from __future__ import annotations

import argparse
import math

import numpy as np

from punctured_ebch_outer import full_spectrum


BLOCKS = 16384
TILES = 256
BAND_COLUMNS = 64
NONTRIVIAL_CARDINALITY = (1 << 64) - 2
GRAPH_CODIMENSION = 24


def joint_half_weight_law() -> np.ndarray:
    """J[a,b] for a uniform nonzero, non-all-ones permuted EBCH word."""

    spectrum = full_spectrum()
    joint = np.zeros((BAND_COLUMNS + 1, BAND_COLUMNS + 1), dtype=np.longdouble)
    for left in range(BAND_COLUMNS + 1):
        for right in range(BAND_COLUMNS + 1):
            total = left + right
            if not 0 < total < 2 * BAND_COLUMNS or not spectrum[total]:
                continue
            joint[left, right] = (
                np.longdouble(spectrum[total])
                / NONTRIVIAL_CARDINALITY
                * np.longdouble(math.comb(BAND_COLUMNS, left))
                * np.longdouble(math.comb(BAND_COLUMNS, right))
                / np.longdouble(math.comb(2 * BAND_COLUMNS, total))
            )
    if abs(float(np.sum(joint)) - 1.0) > 1e-15:
        raise SystemExit("two-band star: joint law does not sum to one")
    return joint


def union_transitions() -> np.ndarray:
    """T[u,w,v] = Pr[|U union W|=v | |U|=u, |W|=w]."""

    transition = np.zeros(
        (BAND_COLUMNS + 1, BAND_COLUMNS + 1, BAND_COLUMNS + 1),
        dtype=np.longdouble,
    )
    for union in range(BAND_COLUMNS + 1):
        for weight in range(BAND_COLUMNS + 1):
            denominator = math.comb(BAND_COLUMNS, weight)
            for intersection in range(
                max(0, union + weight - BAND_COLUMNS), min(union, weight) + 1
            ):
                result = union + weight - intersection
                transition[union, weight, result] += (
                    np.longdouble(math.comb(union, intersection))
                    * np.longdouble(
                        math.comb(BAND_COLUMNS - union, weight - intersection)
                    )
                    / np.longdouble(denominator)
                )
    return transition


def star_kernels(pole: float) -> tuple[np.ndarray, np.ndarray]:
    """Left-union kernels after closing a degree-one or degree-two leaf."""

    joint = joint_half_weight_law()
    transition = union_transitions()
    powers = np.array(
        [np.longdouble(pole) ** weight for weight in range(BAND_COLUMNS + 1)],
        dtype=np.longdouble,
    )

    # A singleton right tile closes immediately and contributes q^b.
    left_weight = joint @ powers
    degree_one = np.einsum("uav,a->uv", transition, left_weight, optimize=True)

    # For a degree-two right tile, H[b,c] is its exact union moment.  The
    # matrix C[a,d] then retains the correlation between both left weights.
    right_union_moment = np.einsum(
        "ubv,v->ub", transition, powers, optimize=True
    )
    pair_left_weights = joint @ right_union_moment @ joint.T

    # Contract in two stages.  This avoids materializing a five-dimensional
    # tensor for the two successive updates of the common left union.
    intermediate = np.einsum(
        "uav,ab->uvb", transition, pair_left_weights, optimize=True
    )
    degree_two = np.einsum(
        "uvb,vbw->uw", intermediate, transition, optimize=True
    )
    return degree_one, degree_two


def falling_factorial(n: int, k: int) -> int:
    result = 1
    for offset in range(k):
        result *= n - offset
    return result


def exact_star_moment(singletons: int, pairs: int, pole: float) -> float:
    degree_one, degree_two = star_kernels(pole)
    distribution = np.zeros(BAND_COLUMNS + 1, dtype=np.longdouble)
    distribution[0] = 1.0
    for _ in range(singletons):
        distribution = distribution @ degree_one
    for _ in range(pairs):
        distribution = distribution @ degree_two
    powers = np.array(
        [np.longdouble(pole) ** weight for weight in range(BAND_COLUMNS + 1)],
        dtype=np.longdouble,
    )
    return float(distribution @ powers)


def outer_weighted_log2(singletons: int, pairs: int, moment: float) -> float:
    blocks = singletons + 2 * pairs
    right_tiles = singletons + pairs
    left_assignment = TILES / TILES**blocks
    right_assignment = falling_factorial(TILES, right_tiles) / TILES**blocks
    return (
        blocks * (math.log2(BLOCKS) + math.log2(NONTRIVIAL_CARDINALITY))
        - math.lgamma(blocks + 1) / math.log(2)
        - GRAPH_CODIMENSION
        + math.log2(left_assignment)
        + math.log2(right_assignment)
        + math.log2(moment)
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--group-pole", type=float, default=0.181)
    parser.add_argument("--max-blocks", type=int, default=12)
    args = parser.parse_args()

    degree_one, degree_two = star_kernels(args.group_pole)
    powers = np.array(
        [np.longdouble(args.group_pole) ** weight for weight in range(65)],
        dtype=np.longdouble,
    )
    distributions: dict[tuple[int, int], np.ndarray] = {}
    seed = np.zeros(65, dtype=np.longdouble)
    seed[0] = 1.0
    distributions[(0, 0)] = seed

    print("exact-joint two-band star moments")
    print(f"group_pole={args.group_pole:.12f}")
    for blocks in range(1, args.max_blocks + 1):
        # The dominant capacity-two product bound uses as many pairs as
        # possible, with one singleton when the block count is odd.
        pairs = blocks // 2
        singletons = blocks - 2 * pairs
        distribution = seed
        for _ in range(singletons):
            distribution = distribution @ degree_one
        for _ in range(pairs):
            distribution = distribution @ degree_two
        moment = float(distribution @ powers)
        print(
            f"active_blocks={blocks:2d} right_singletons={singletons} "
            f"right_pairs={pairs} support_moment_log2={math.log2(moment):.12f} "
            f"outer_weighted_log2={outer_weighted_log2(singletons, pairs, moment):.12f}"
        )
    print("status=DIAGNOSTIC_NEEDS_OUTWARD_ARITHMETIC_AND_ALL_PARTITION_PAIRS")


if __name__ == "__main__":
    main()
