#!/usr/bin/env python3
"""Verify the packet-weight accumulator enumerator.

For packet width g, let h_k count input packets of Hamming weight k.  The
quotient enumerator tracks the Hamming weight of the accumulator state, rather
than the state itself.  This script compares that quotient against exhaustive
enumeration and against the classical scalar accumulator formula.
"""

from __future__ import annotations

import itertools
import json
import math
from collections import Counter, defaultdict
from collections.abc import Iterator


def transition_multiplicity(g: int, old_weight: int, new_weight: int, packet_weight: int) -> int:
    """Count weight-packet_weight values that realize one quotient transition."""
    numerator = old_weight + new_weight - packet_weight
    if numerator & 1:
        return 0
    overlap = numerator // 2
    if not 0 <= overlap <= old_weight:
        return 0
    added = new_weight - overlap
    if not 0 <= added <= g - old_weight:
        return 0
    return math.comb(old_weight, overlap) * math.comb(g - old_weight, added)


def quotient_enumerator(g: int, histogram: tuple[int, ...]) -> Counter[int]:
    """Return exact output-weight counts for one packet-weight histogram."""
    if len(histogram) != g + 1 or any(count < 0 for count in histogram):
        raise ValueError("histogram must contain nonnegative h_0,...,h_g")
    length = sum(histogram)
    zero_used = (0,) * (g + 1)
    current: dict[tuple[tuple[int, ...], int, int], int] = {
        (zero_used, 0, 0): 1
    }
    for _ in range(length):
        following: defaultdict[tuple[tuple[int, ...], int, int], int] = defaultdict(int)
        for (used, old_weight, output_weight), ways in current.items():
            for packet_weight in range(g + 1):
                if used[packet_weight] == histogram[packet_weight]:
                    continue
                next_used = list(used)
                next_used[packet_weight] += 1
                next_used_tuple = tuple(next_used)
                for new_weight in range(g + 1):
                    multiplicity = transition_multiplicity(
                        g, old_weight, new_weight, packet_weight
                    )
                    if multiplicity:
                        following[
                            (
                                next_used_tuple,
                                new_weight,
                                output_weight + new_weight,
                            )
                        ] += ways * multiplicity
        current = following
    result: Counter[int] = Counter()
    for (used, _state_weight, output_weight), ways in current.items():
        if used != histogram:
            raise RuntimeError("quotient enumerator ended at a partial histogram")
        result[output_weight] += ways
    return result


def histogram_sequence_count(g: int, histogram: tuple[int, ...]) -> int:
    """Count input-value sequences with the prescribed weight histogram."""
    length = sum(histogram)
    result = math.factorial(length)
    for count in histogram:
        result //= math.factorial(count)
    for packet_weight, count in enumerate(histogram):
        result *= math.comb(g, packet_weight) ** count
    return result


def scalar_accumulator_count(length: int, input_weight: int, output_weight: int) -> int:
    """Return the classical zero-start accumulator input-output enumerator."""
    if input_weight == 0:
        return int(output_weight == 0)
    if output_weight == 0:
        return 0
    return math.comb(length - output_weight, input_weight // 2) * math.comb(
        output_weight - 1, (input_weight + 1) // 2 - 1
    )


def exhaustive_enumerators(g: int, length: int) -> dict[tuple[int, ...], Counter[int]]:
    """Enumerate every input-value sequence independently of the quotient DP."""
    result: dict[tuple[int, ...], Counter[int]] = defaultdict(Counter)
    values = range(1 << g)
    value_weights = tuple(value.bit_count() for value in values)
    for sequence in itertools.product(values, repeat=length):
        histogram = [0] * (g + 1)
        state = 0
        output_weight = 0
        for value in sequence:
            histogram[value_weights[value]] += 1
            state ^= value
            output_weight += state.bit_count()
        result[tuple(histogram)][output_weight] += 1
    return result


def audit_transition_quotient(g: int) -> int:
    """Compare every quotient multiplicity with direct state enumeration."""
    comparisons = 0
    for old_weight in range(g + 1):
        state = (1 << old_weight) - 1
        for packet_weight in range(g + 1):
            direct = Counter(
                (state ^ value).bit_count()
                for value in range(1 << g)
                if value.bit_count() == packet_weight
            )
            for new_weight in range(g + 1):
                expected = transition_multiplicity(
                    g, old_weight, new_weight, packet_weight
                )
                if direct[new_weight] != expected:
                    raise RuntimeError("five-state transition quotient changed")
                comparisons += 1
    return comparisons


def weighted_partition_full(
    g: int,
    length: int,
    packet_factors: tuple[int, ...],
    state_factors: tuple[int, ...],
) -> int:
    """Evaluate the original 2^g-state transfer operator over the integers."""
    current = [0] * (1 << g)
    current[0] = 1
    for _ in range(length):
        following = [0] * (1 << g)
        for old_state, ways in enumerate(current):
            for packet in range(1 << g):
                new_state = old_state ^ packet
                following[new_state] += (
                    ways
                    * packet_factors[packet.bit_count()]
                    * state_factors[new_state.bit_count()]
                )
        current = following
    return sum(current)


def weighted_partition_quotient(
    g: int,
    length: int,
    packet_factors: tuple[int, ...],
    state_factors: tuple[int, ...],
) -> int:
    """Evaluate the (g+1)-state Hamming-weight quotient over the integers."""
    current = [0] * (g + 1)
    current[0] = 1
    for _ in range(length):
        following = [0] * (g + 1)
        for old_weight, ways in enumerate(current):
            for new_weight in range(g + 1):
                transition = 0
                for packet_weight in range(g + 1):
                    transition += packet_factors[packet_weight] * transition_multiplicity(
                        g, old_weight, new_weight, packet_weight
                    )
                following[new_weight] += (
                    ways * transition * state_factors[new_weight]
                )
        current = following
    return sum(current)


def audit_weighted_transfer_quotient() -> list[dict[str, object]]:
    rows = []
    for probe in range(1, 9):
        length = probe + 2
        packet_factors = tuple(1 + (probe + 2) * index for index in range(5))
        state_factors = tuple(2 + (2 * probe + 1) * index for index in range(5))
        full = weighted_partition_full(
            4, length, packet_factors, state_factors
        )
        quotient = weighted_partition_quotient(
            4, length, packet_factors, state_factors
        )
        if full != quotient:
            raise RuntimeError("five-state quotient differs from 16-state transfer")
        rows.append(
            {
                "length": length,
                "packet_factors": list(packet_factors),
                "state_factors": list(state_factors),
                "partition_function": full,
            }
        )
    return rows


def audit_scalar_formula(maximum_length: int) -> int:
    comparisons = 0
    for length in range(1, maximum_length + 1):
        for input_weight in range(length + 1):
            quotient = quotient_enumerator(
                1, (length - input_weight, input_weight)
            )
            for output_weight in range(length + 1):
                formula = scalar_accumulator_count(
                    length, input_weight, output_weight
                )
                if quotient[output_weight] != formula:
                    raise RuntimeError("g=1 quotient differs from scalar formula")
                comparisons += 1
            if sum(quotient.values()) != math.comb(length, input_weight):
                raise RuntimeError("g=1 histogram mass changed")
    return comparisons


def audit_exhaustive_g4(length: int) -> tuple[int, int, dict[str, int]]:
    exhaustive = exhaustive_enumerators(4, length)
    output_comparisons = 0
    for histogram, direct in exhaustive.items():
        quotient = quotient_enumerator(4, histogram)
        if quotient != direct:
            raise RuntimeError("five-state quotient differs from exhaustive g=4 scan")
        if sum(quotient.values()) != histogram_sequence_count(4, histogram):
            raise RuntimeError("g=4 quotient has the wrong histogram mass")
        output_comparisons += 4 * length + 1
    sample_histogram = (1, 2, 1, 1, 0)
    if sum(sample_histogram) != length:
        raise RuntimeError("sample histogram does not match exhaustive length")
    sample = {
        str(weight): count
        for weight, count in sorted(exhaustive[sample_histogram].items())
    }
    return len(exhaustive), output_comparisons, sample


def audit_larger_masses() -> list[dict[str, object]]:
    histograms = [
        (3, 4, 3, 2, 0),
        (1, 3, 4, 3, 1),
        (4, 2, 2, 2, 2),
        (0, 1, 3, 4, 4),
    ]
    rows = []
    for histogram in histograms:
        enumerator = quotient_enumerator(4, histogram)
        expected = histogram_sequence_count(4, histogram)
        if sum(enumerator.values()) != expected:
            raise RuntimeError("larger quotient mass differs from multinomial count")
        rows.append(
            {
                "histogram_h0_through_h4": list(histogram),
                "input_sequences": expected,
                "output_weight_support": [min(enumerator), max(enumerator)],
                "distinct_output_weights": len(enumerator),
            }
        )
    return rows


def compositions(total: int, parts: int) -> Iterator[tuple[int, ...]]:
    """Yield ordered weak compositions of total into the requested parts."""
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for suffix in compositions(total - first, parts - 1):
            yield (first,) + suffix


def audit_one_block_packetization() -> dict[str, object]:
    """Check the complete four-bit packet histogram partition of 128 bits."""
    packet_count = 32
    mass_by_binary_weight = [0] * 129
    histograms = 0
    for histogram in compositions(packet_count, 5):
        binary_weight = sum(index * count for index, count in enumerate(histogram))
        mass_by_binary_weight[binary_weight] += histogram_sequence_count(4, histogram)
        histograms += 1
    for binary_weight, mass in enumerate(mass_by_binary_weight):
        if mass != math.comb(128, binary_weight):
            raise RuntimeError("packet histograms do not partition a BCH weight slice")
    return {
        "packet_count": packet_count,
        "packet_weight_histograms": histograms,
        "binary_weights_checked": len(mass_by_binary_weight),
        "weight22_slice": mass_by_binary_weight[22],
        "weight24_slice": mass_by_binary_weight[24],
    }


def main() -> None:
    transition_comparisons = audit_transition_quotient(4)
    transfer_rows = audit_weighted_transfer_quotient()
    scalar_comparisons = audit_scalar_formula(20)
    exhaustive_length = 5
    exhaustive_histograms, exhaustive_output_comparisons, sample = (
        audit_exhaustive_g4(exhaustive_length)
    )
    larger_rows = audit_larger_masses()
    block_packetization = audit_one_block_packetization()
    payload = {
        "schema": "riffle-parallelacc-packet-weight-formula-audit-v1",
        "packet_width": 4,
        "quotient_state_count": 5,
        "packet_histogram_coordinates": ["h0", "h1", "h2", "h3", "h4"],
        "active_coordinates": ["h1", "h2", "h3", "h4"],
        "transition_comparisons": transition_comparisons,
        "weighted_16_state_transfer_comparisons": transfer_rows,
        "scalar_formula": {
            "maximum_length": 20,
            "coefficient_comparisons": scalar_comparisons,
        },
        "exhaustive_g4": {
            "length": exhaustive_length,
            "input_sequences": (1 << 4) ** exhaustive_length,
            "packet_weight_histograms": exhaustive_histograms,
            "output_coefficient_comparisons": exhaustive_output_comparisons,
            "sample_histogram_h0_through_h4": [1, 2, 1, 1, 0],
            "sample_output_enumerator": sample,
        },
        "larger_mass_checks": larger_rows,
        "one_bch_block_packetization": block_packetization,
        "validation": {
            "binomial_transition_multiplicities_match_direct_states": True,
            "g1_reduces_to_classical_accumulator_formula": True,
            "g4_quotient_matches_exhaustive_input_sequences": True,
            "five_state_quotient_matches_weighted_16_state_transfer": True,
            "all_enumerator_masses_match_multinomial_orientation_counts": True,
            "packet_histograms_partition_every_128_bit_weight_slice": True,
        },
        "scope": (
            "Exact conditional accumulator enumerator for a fixed packet-weight "
            "histogram. This receipt does not yet average the histogram over "
            "randomized BCH blocks or sum over the outer code."
        ),
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
