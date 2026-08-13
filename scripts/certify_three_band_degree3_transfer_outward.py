#!/usr/bin/env python3
"""Directed-outward certificate for the degree-three peeling transfer.

Exact rational rooted moments and incremental hypergeometric moments are first
rounded upward to IEEE binary64 numbers, verified by converting each float
back to ``Fraction``.  Positive three-factor products receive a 2^-48 relative
guard, comfortably above the standard two-operation binary64 error.  Each
greedy cap dot product has fewer than 37190 terms and receives a 2^-30 guard,
comfortably above its naive ``gamma_n`` error bound.  The 43 positive diagonal
sums are combined with ``math.fsum`` and rounded upward once more.

The cell-cap ordering uses these outward cell scores, so its greedy value is
an upper bound even if outward rounding changes ties or ordering.  Finally the
three-block boundary-budget DP multiplies exact ``Fraction.from_float`` values
and proves the total transfer is at most 2^-30.
"""

from __future__ import annotations

import math
from fractions import Fraction

import numpy as np

from probe_bch_three_band_cell_cap_motifs import cell_caps
from probe_bch_three_band_transport_lp import (
    BAND_SIZES,
    exact_rows,
    remove_trivial_codewords,
    triples,
)
from three_band_exact_moments import factor_table_upper


POLE = Fraction(1, 20)
ROOT_BITS = 160
SUPPORT_CUTOFF = 106
PRODUCT_GUARD = 2.0**-48
DOT_GUARD = 2.0**-30


def upward_float(value: Fraction) -> float:
    result = float(value)
    if Fraction.from_float(result) < value:
        result = math.nextafter(result, math.inf)
    if Fraction.from_float(result) < value:
        raise AssertionError("upward float conversion failed")
    return result


def incremental_exact(columns: int, existing: int, weight: int) -> Fraction:
    denominator = math.comb(columns, weight)
    return sum(
        (
            Fraction(
                math.comb(existing, intersection)
                * math.comb(columns - existing, weight - intersection),
                denominator,
            )
            * POLE ** (weight - intersection)
        )
        for intersection in range(
            max(0, existing + weight - columns), min(existing, weight) + 1
        )
    )


def main() -> None:
    state_list = triples()
    states = np.asarray(state_list, dtype=np.int16)
    rows = remove_trivial_codewords(state_list, exact_rows(state_list))
    caps, diagonal_mass, by_total = cell_caps(state_list, rows)
    caps_float = np.asarray([upward_float(Fraction(cap)) for cap in caps])
    active_by_total = {
        weight: np.asarray(
            [index for index in by_total[weight] if caps[index]], dtype=np.int32
        )
        for weight, mass in diagonal_mass.items()
        if mass
    }

    rooted = {
        band: np.asarray(
            [
                upward_float(value)
                for value in factor_table_upper(
                    BAND_SIZES[band], 3, POLE, root_bits=ROOT_BITS
                )
            ]
        )
        for band in range(3)
    }
    incremental = {
        (band, existing): np.asarray(
            [
                upward_float(incremental_exact(BAND_SIZES[band], existing, weight))
                for weight in range(BAND_SIZES[band] + 1)
            ]
        )
        for band in range(3)
        for existing in range(BAND_SIZES[band] + 1)
    }

    def factor_upper(factors: tuple[np.ndarray, np.ndarray, np.ndarray]) -> Fraction:
        scores = factors[0][states[:, 0]] * factors[1][states[:, 1]]
        scores *= factors[2][states[:, 2]]
        scores = np.nextafter(scores * (1.0 + PRODUCT_GUARD), math.inf)
        diagonal_values = []
        for weight, indices in active_by_total.items():
            order = indices[np.argsort(scores[indices])[::-1]]
            remaining = diagonal_mass[weight]
            selected_indices = []
            selected_caps = []
            for index in order:
                take = min(remaining, caps[int(index)])
                if take:
                    selected_indices.append(int(index))
                    selected_caps.append(upward_float(Fraction(take)))
                    remaining -= take
                if not remaining:
                    break
            if remaining:
                raise AssertionError(f"degree-three transfer misses diagonal {weight}")
            products = (
                np.asarray(selected_caps)
                * scores[np.asarray(selected_indices, dtype=np.int32)]
            )
            value = float(np.sum(products, dtype=np.float64))
            value = math.nextafter(value * (1.0 + DOT_GUARD), math.inf)
            diagonal_values.append(value)
        result = math.nextafter(math.fsum(diagonal_values), math.inf)
        return Fraction.from_float(result)

    best_by_cost: dict[int, Fraction] = {}
    witness_by_cost: dict[int, tuple[int, int, int]] = {}
    for low_band in range(3):
        external = [band for band in range(3) if band != low_band]
        for first_weight in range(BAND_SIZES[external[0]] + 1):
            for second_weight in range(BAND_SIZES[external[1]] + 1):
                cost = first_weight + second_weight
                if cost > SUPPORT_CUTOFF:
                    continue
                factors = tuple(
                    rooted[band]
                    if band == low_band
                    else incremental[
                        (
                            band,
                            first_weight if band == external[0] else second_weight,
                        )
                    ]
                    for band in range(3)
                )
                value = factor_upper(factors)  # type: ignore[arg-type]
                if value > best_by_cost.get(cost, Fraction(0)):
                    best_by_cost[cost] = value
                    witness_by_cost[cost] = (
                        low_band,
                        first_weight,
                        second_weight,
                    )

    dp = [Fraction(0)] * (SUPPORT_CUTOFF + 1)
    paths: list[list[tuple[int, int, int]] | None] = [None] * (
        SUPPORT_CUTOFF + 1
    )
    dp[0] = Fraction(1)
    paths[0] = []
    for _block in range(3):
        following = [Fraction(0)] * (SUPPORT_CUTOFF + 1)
        following_paths: list[list[tuple[int, int, int]] | None] = [None] * (
            SUPPORT_CUTOFF + 1
        )
        for used, prefix in enumerate(dp):
            if not prefix:
                continue
            for cost, value in best_by_cost.items():
                if used + cost > SUPPORT_CUTOFF:
                    continue
                candidate = prefix * value
                if candidate > following[used + cost]:
                    following[used + cost] = candidate
                    following_paths[used + cost] = paths[used] + [  # type: ignore[operator]
                        witness_by_cost[cost]
                    ]
        dp = following
        paths = following_paths

    value = max(dp)
    used = dp.index(value)
    envelope = Fraction(11, 16 * (1 << 30))
    if value > envelope:
        raise SystemExit("degree-three transfer exceeds (11/16) 2^-30")
    print("outward three-band degree-three transfer certificate")
    print(f"pole={POLE} support_cutoff={SUPPORT_CUTOFF} root_bits={ROOT_BITS}")
    print(f"transfer_log2_upper={math.log2(float(value)):.12f}")
    print(f"used_budget={used} witness={paths[used]}")
    print("transfer_le_(11/16)*2^-30=PASS")
    print("status=EXACT_MOMENTS_DIRECTED_OUTWARD_BINARY64")


if __name__ == "__main__":
    main()
