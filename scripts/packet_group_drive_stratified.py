#!/usr/bin/env python3
"""Generic exact-local drive tables for permutation atom width ``g``.

The 64-bit inner block is partitioned into ``64/g`` ordered atoms.  A direct
enumeration of current-input/state atom pairs costs 4**g.  Instead, fix the
drive atom W=U xor R:

* accumulator parity and emitted weight depend only on the concrete W; and
* the weighted count of compatible (U,R) pairs depends only on wt(W).

This factorization gives a polynomial accumulator-pattern DP and a
hypergeometric input/state polynomial.  The resulting atom transfers compose
through the same fixed-size 64-bit block DP used by packet-8.
"""

from __future__ import annotations

import itertools
import math

import numpy as np

from packet_group_native import load_point_caps


BLOCK_BITS = 64
SUPPORTED_GROUPS = (1, 2, 4, 8, 16, 32, 64)


def validate_group(group_bits: int) -> int:
    if group_bits not in SUPPORTED_GROUPS or BLOCK_BITS % group_bits:
        raise ValueError(f"group width must be one of {SUPPORTED_GROUPS}")
    return BLOCK_BITS // group_bits


def accumulator_pattern_counts(group_bits: int) -> np.ndarray:
    """Count concrete drive atoms by boundary parity, weight, and emission."""

    validate_group(group_bits)
    result = np.zeros(
        (2, 2, group_bits + 1, group_bits + 1), dtype=object
    )
    for initial_parity in range(2):
        # axes: current parity, drive weight, emitted weight
        dp = np.zeros((2, group_bits + 1, group_bits + 1), dtype=object)
        dp[initial_parity, 0, 0] = 1
        for position in range(group_bits):
            next_dp = np.zeros_like(dp)
            for parity in range(2):
                for drive_weight in range(position + 1):
                    for emitted in range(position + 1):
                        count = dp[parity, drive_weight, emitted]
                        if not count:
                            continue
                        # A zero drive bit preserves parity.
                        next_dp[parity, drive_weight, emitted + parity] += count
                        # A one drive bit toggles parity before emission.
                        toggled = parity ^ 1
                        next_dp[
                            toggled,
                            drive_weight + 1,
                            emitted + toggled,
                        ] += count
            dp = next_dp
        for final_parity in range(2):
            result[initial_parity, final_parity] = dp[final_parity]
    return result


def compatible_pair_polynomials(
    group_bits: int, fugacities: np.ndarray
) -> np.ndarray:
    """Weighted compatible-pair mass indexed by drive and state weight."""

    validate_group(group_bits)
    if fugacities.shape != (group_bits + 1,) or np.any(fugacities < 0.0):
        raise ValueError("fugacity vector has wrong shape or negative entry")
    result = np.zeros((group_bits + 1, group_bits + 1), dtype=np.float64)
    for drive_weight in range(group_bits + 1):
        for overlap in range(drive_weight + 1):
            for outside in range(group_bits - drive_weight + 1):
                input_weight = overlap + outside
                state_weight = drive_weight - overlap + outside
                multiplicity = math.comb(drive_weight, overlap) * math.comb(
                    group_bits - drive_weight, outside
                )
                result[drive_weight, state_weight] += (
                    multiplicity * float(fugacities[input_weight])
                )
    return result


def _up_add(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    raw = left + right
    return np.where(right > 0.0, np.nextafter(raw, np.inf), left)


def _up_mul(left: np.ndarray | float, right: np.ndarray | float) -> np.ndarray:
    left_array = np.asarray(left)
    right_array = np.asarray(right)
    raw = left_array * right_array
    # Advancing a positive rounded product by one ulp is outward only if the
    # rounded product stayed positive.  A positive exact dyadic product that
    # underflows to zero would otherwise be mistaken for structural zero.
    positive_underflow = (left_array > 0.0) & (right_array > 0.0) & (raw == 0.0)
    if np.any(positive_underflow):
        raise FloatingPointError(
            "outward multiplication underflowed a structurally positive product"
        )
    return np.where(raw > 0.0, np.nextafter(raw, np.inf), raw)


def compatible_pair_polynomials_outward(
    group_bits: int, fugacities: np.ndarray
) -> np.ndarray:
    """Componentwise upper enclosure, treating float inputs as exact dyadics."""

    validate_group(group_bits)
    if fugacities.shape != (group_bits + 1,) or np.any(fugacities < 0.0):
        raise ValueError("fugacity vector has wrong shape or negative entry")
    result = np.zeros((group_bits + 1, group_bits + 1), dtype=np.float64)
    for drive_weight in range(group_bits + 1):
        for overlap in range(drive_weight + 1):
            for outside in range(group_bits - drive_weight + 1):
                input_weight = overlap + outside
                state_weight = drive_weight - overlap + outside
                multiplicity = math.comb(drive_weight, overlap) * math.comb(
                    group_bits - drive_weight, outside
                )
                term = float(_up_mul(float(multiplicity), fugacities[input_weight]))
                result[drive_weight, state_weight] = _up_add(
                    result[drive_weight, state_weight], np.asarray(term)
                )
    return result


def atom_transfer(group_bits: int, fugacities: np.ndarray) -> np.ndarray:
    """Return transfer[p_in,p_out,state,drive,emitted]."""

    patterns = accumulator_pattern_counts(group_bits)
    pairs = compatible_pair_polynomials(group_bits, fugacities)
    transfer = np.zeros(
        (2, 2, group_bits + 1, group_bits + 1, group_bits + 1),
        dtype=np.float64,
    )
    for incoming in range(2):
        for outgoing in range(2):
            for drive in range(group_bits + 1):
                for emitted in range(group_bits + 1):
                    count = int(patterns[incoming, outgoing, drive, emitted])
                    if count:
                        transfer[incoming, outgoing, :, drive, emitted] = (
                            count * pairs[drive]
                        )
    return transfer


def atom_transfer_outward(group_bits: int, fugacities: np.ndarray) -> np.ndarray:
    patterns = accumulator_pattern_counts(group_bits)
    pairs = compatible_pair_polynomials_outward(group_bits, fugacities)
    transfer = np.zeros(
        (2, 2, group_bits + 1, group_bits + 1, group_bits + 1),
        dtype=np.float64,
    )
    for incoming in range(2):
        for outgoing in range(2):
            for drive in range(group_bits + 1):
                for emitted in range(group_bits + 1):
                    count = int(patterns[incoming, outgoing, drive, emitted])
                    if count:
                        transfer[incoming, outgoing, :, drive, emitted] = _up_mul(
                            float(count), pairs[drive]
                        )
    return transfer


def block_histograms(
    group_bits: int, fugacities: np.ndarray, *, check_mass: bool = True
) -> np.ndarray:
    """Return H[q,d,y] before division by C(64,q)."""

    atoms = validate_group(group_bits)
    transfer = atom_transfer(group_bits, fugacities)
    entries = [
        (
            incoming,
            outgoing,
            selected,
            drive,
            emitted,
            transfer[incoming, outgoing, selected, drive, emitted],
        )
        for incoming in range(2)
        for outgoing in range(2)
        for selected in range(group_bits + 1)
        for drive in range(group_bits + 1)
        for emitted in range(group_bits + 1)
        if transfer[incoming, outgoing, selected, drive, emitted]
    ]
    dp = np.zeros(
        (2, BLOCK_BITS + 1, BLOCK_BITS + 1, BLOCK_BITS + 1),
        dtype=np.float64,
    )
    dp[0, 0, 0, 0] = 1.0
    maximum = 0
    for _ in range(atoms):
        next_dp = np.zeros_like(dp)
        for incoming, outgoing, selected, drive, emitted, coefficient in entries:
            source = dp[incoming, : maximum + 1, : maximum + 1, : maximum + 1]
            next_dp[
                outgoing,
                selected : selected + maximum + 1,
                drive : drive + maximum + 1,
                emitted : emitted + maximum + 1,
            ] += coefficient * source
        dp = next_dp
        maximum += group_bits
    result = dp.sum(axis=0)
    if check_mass:
        concrete_mass = sum(
            math.comb(group_bits, weight) * float(fugacity)
            for weight, fugacity in enumerate(fugacities)
        )
        for state_weight in range(BLOCK_BITS + 1):
            expected = (
                math.comb(BLOCK_BITS, state_weight) * concrete_mass**atoms
            )
            actual = float(np.sum(result[state_weight]))
            if abs(actual - expected) > 5e-11 * max(1.0, expected):
                raise SystemExit(
                    "generic drive table: mass mismatch "
                    f"g={group_bits} q={state_weight} actual={actual} expected={expected}"
                )
    return result


def block_histograms_outward(group_bits: int, fugacities: np.ndarray) -> np.ndarray:
    atoms = validate_group(group_bits)
    transfer = atom_transfer_outward(group_bits, fugacities)
    entries = [
        (incoming, outgoing, selected, drive, emitted,
         transfer[incoming, outgoing, selected, drive, emitted])
        for incoming in range(2)
        for outgoing in range(2)
        for selected in range(group_bits + 1)
        for drive in range(group_bits + 1)
        for emitted in range(group_bits + 1)
        if transfer[incoming, outgoing, selected, drive, emitted]
    ]
    dp = np.zeros(
        (2, BLOCK_BITS + 1, BLOCK_BITS + 1, BLOCK_BITS + 1),
        dtype=np.float64,
    )
    dp[0, 0, 0, 0] = 1.0
    maximum = 0
    for _ in range(atoms):
        next_dp = np.zeros_like(dp)
        for incoming, outgoing, selected, drive, emitted, coefficient in entries:
            source = dp[incoming, : maximum + 1, : maximum + 1, : maximum + 1]
            contribution = _up_mul(coefficient, source)
            destination = next_dp[
                outgoing,
                selected : selected + maximum + 1,
                drive : drive + maximum + 1,
                emitted : emitted + maximum + 1,
            ]
            destination[...] = _up_add(destination, contribution)
        dp = next_dp
        maximum += group_bits
    return _up_add(dp[0], dp[1])


def point_caps(group_bits: int, fugacities: np.ndarray) -> np.ndarray:
    """Exact finite max c_g[q,d], evaluated with binary64 fugacities."""

    atoms = validate_group(group_bits)
    local = compatible_pair_polynomials(group_bits, fugacities)
    native = load_point_caps()
    if native is not None:
        local_flat = np.ascontiguousarray(local.reshape(-1), dtype=np.float64)
        maxima = np.zeros((BLOCK_BITS + 1, BLOCK_BITS + 1), dtype=np.float64)
        status = native(group_bits, local_flat, maxima.reshape(-1))
        if status != 0:
            raise RuntimeError(f"native point-cap builder failed with status {status}")
        for state_weight in range(BLOCK_BITS + 1):
            maxima[state_weight] /= math.comb(BLOCK_BITS, state_weight)
        return maxima
    maxima = np.zeros((BLOCK_BITS + 1, BLOCK_BITS + 1), dtype=np.float64)
    for weights in itertools.combinations_with_replacement(
        range(group_bits + 1), atoms
    ):
        polynomial = np.array([1.0])
        for weight in weights:
            polynomial = np.convolve(polynomial, local[weight])
        maxima[:, sum(weights)] = np.maximum(maxima[:, sum(weights)], polynomial)
    for state_weight in range(BLOCK_BITS + 1):
        maxima[state_weight] /= math.comb(BLOCK_BITS, state_weight)
    return maxima


def _convolve_outward(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    result = np.zeros(len(left) + len(right) - 1, dtype=np.float64)
    for left_index, left_value in enumerate(left):
        if left_value == 0.0:
            continue
        contribution = _up_mul(left_value, right)
        target = result[left_index : left_index + len(right)]
        target[...] = _up_add(target, contribution)
    return result


def point_caps_outward(group_bits: int, fugacities: np.ndarray) -> np.ndarray:
    atoms = validate_group(group_bits)
    local = compatible_pair_polynomials_outward(group_bits, fugacities)
    maxima = np.zeros((BLOCK_BITS + 1, BLOCK_BITS + 1), dtype=np.float64)
    for weights in itertools.combinations_with_replacement(
        range(group_bits + 1), atoms
    ):
        polynomial = np.array([1.0])
        for weight in weights:
            polynomial = _convolve_outward(polynomial, local[weight])
        total = sum(weights)
        maxima[:, total] = np.maximum(maxima[:, total], polynomial)
    for state_weight in range(BLOCK_BITS + 1):
        positive = maxima[state_weight] > 0.0
        quotient = maxima[state_weight, positive] / math.comb(BLOCK_BITS, state_weight)
        maxima[state_weight, positive] = np.nextafter(quotient, np.inf)
    return maxima


def profile_classes(group_bits: int) -> tuple[int, ...]:
    validate_group(group_bits)
    return tuple(math.comb(group_bits, weight) for weight in range(group_bits + 1))


def profile_count(group_bits: int, total_bits: int) -> int:
    """Count every nonnegative packet profile, including low physical weight."""

    atoms = total_bits // group_bits
    return math.comb(atoms + group_bits, group_bits)


def feasible_profile_count(
    group_bits: int, total_bits: int, minimum_physical_weight: int
) -> int:
    """Count profiles whose total physical weight reaches the given threshold."""

    atoms = total_bits // group_bits
    if minimum_physical_weight <= 0:
        return profile_count(group_bits, total_bits)
    if minimum_physical_weight > atoms * group_bits:
        return 0
    # For excluded profiles, a_0 is determined by the positive-weight counts.
    # The threshold is tiny relative to every active construction's atom count.
    maximum_excluded_weight = minimum_physical_weight - 1
    if atoms < maximum_excluded_weight:
        raise ValueError("feasible profile count requires a general mass-aware DP")
    partitions = [0] * (maximum_excluded_weight + 1)
    partitions[0] = 1
    for weight in range(1, group_bits + 1):
        for total_weight in range(weight, maximum_excluded_weight + 1):
            partitions[total_weight] += partitions[total_weight - weight]
    return profile_count(group_bits, total_bits) - sum(partitions)
