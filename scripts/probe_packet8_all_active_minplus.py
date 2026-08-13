#!/usr/bin/env python3
"""Min-plus lower bound for the repeated-0xff inner trajectory.

For every incoming state weight, the exact accumulator histogram identifies
which emitted weights are possible under the random state permutation when
all eight current packets are 0xff.  The rigorous systematic EBCH split caps
identify possible next-state weights.  Join these supports to form a 65-state
directed graph whose edge cost is the smallest compatible emitted weight.

The graph deliberately forgets correlations and multiplicities, so its
shortest paths are deterministic lower bounds for every actual trajectory.
It reports the reachable minimum-mean cycle and the exact min-plus cost after
32768 all-active blocks.  It then constructs an exact integer potential and
checks every one of the 256 packet masks.  The resulting inequality charges
each zero-packet defect directly while telescoping the state potential.
"""

from __future__ import annotations

import argparse
import math

import numpy as np

from analyze_systematic_group_kernel import EXACT_SLICES, SPECTRUM
from certificate_spectra import load_ebch128_spectrum
from probe_packet8_constant_value_scalar import BITS, D, INNER_BLOCKS
from probe_packet8_repeated_value_state_transfer import build_histograms
from probe_packet8_weight_transfer import build_split_caps, load_exact


INF = 10**18


def edge_costs() -> np.ndarray:
    spectrum = {
        weight: count
        for weight, count in load_ebch128_spectrum(SPECTRUM)
        if weight and count
    }
    caps = build_split_caps(spectrum, load_exact(EXACT_SLICES))
    histograms = build_histograms(0xFF)
    edges = np.full((BITS + 1, BITS + 1), INF, dtype=np.int64)
    for incoming in range(BITS + 1):
        for emitted in range(BITS + 1):
            if not histograms[incoming, 0xFF, emitted]:
                continue
            for outgoing in range(BITS + 1):
                if caps[emitted][outgoing]:
                    edges[incoming, outgoing] = min(
                        int(edges[incoming, outgoing]), emitted
                    )
    return edges


def minimum_mean_cycle(edges: np.ndarray) -> float:
    states = edges.shape[0]
    # Karp's formula, with a virtual zero-cost source reaching every state.
    dp = np.full((states + 1, states), INF, dtype=np.int64)
    dp[0] = 0
    for length in range(1, states + 1):
        dp[length] = np.min(dp[length - 1][:, None] + edges, axis=0)
    result = math.inf
    for state in range(states):
        if dp[states, state] >= INF:
            continue
        maximum = -math.inf
        for length in range(states):
            if dp[length, state] < INF:
                maximum = max(
                    maximum,
                    (int(dp[states, state]) - int(dp[length, state]))
                    / (states - length),
                )
        result = min(result, maximum)
    return result


def shortest_cost(edges: np.ndarray, blocks: int) -> int:
    values = np.full(BITS + 1, INF, dtype=np.int64)
    values[0] = 0
    for _ in range(blocks):
        values = np.min(values[:, None] + edges, axis=0)
    return int(np.min(values))


def all_active_potential(edges: np.ndarray) -> np.ndarray:
    """Shortest-path potential for edge cost minus six."""

    adjusted = edges - 6
    potential = np.zeros(BITS + 1, dtype=np.int64)
    for _ in range(BITS + 1):
        next_potential = np.minimum(
            potential, np.min(potential[:, None] + adjusted, axis=0)
        )
        if np.array_equal(next_potential, potential):
            break
        potential = next_potential
    else:
        raise SystemExit("all-active minplus: negative adjusted cycle")
    finite = edges < INF
    incoming, outgoing = np.where(finite)
    slack = (
        edges[incoming, outgoing]
        - 6
        - potential[outgoing]
        + potential[incoming]
    )
    if int(np.min(slack)) < 0:
        raise SystemExit("all-active minplus: invalid potential")
    return potential


def defect_penalty(potential: np.ndarray) -> tuple[int, tuple[int, ...]]:
    spectrum = {
        weight: count
        for weight, count in load_ebch128_spectrum(SPECTRUM)
        if weight and count
    }
    caps = build_split_caps(spectrum, load_exact(EXACT_SLICES))
    histograms = build_histograms(0xFF)
    maximum_numerator = -INF
    maximum_denominator = 1
    witness_row = None
    for mask in range(1 << 8):
        defects = 8 - mask.bit_count()
        for incoming in range(BITS + 1):
            for emitted in np.flatnonzero(histograms[incoming, mask]):
                for outgoing in range(BITS + 1):
                    if not caps[int(emitted)][outgoing]:
                        continue
                    needed = (
                        6
                        + int(potential[outgoing])
                        - int(potential[incoming])
                        - int(emitted)
                    )
                    if defects == 0:
                        if needed > 0:
                            raise SystemExit(
                                "all-active minplus: active-mask potential failed"
                            )
                        continue
                    if needed * maximum_denominator > maximum_numerator * defects:
                        maximum_numerator = needed
                        maximum_denominator = defects
                        witness_row = (
                            mask,
                            defects,
                            incoming,
                            int(emitted),
                            outgoing,
                            needed,
                        )
    assert witness_row is not None
    penalty = max(0, (maximum_numerator + maximum_denominator - 1) // maximum_denominator)
    # Verify the rounded integer penalty against every supported transition.
    for mask in range(1 << 8):
        defects = 8 - mask.bit_count()
        for incoming in range(BITS + 1):
            for emitted in np.flatnonzero(histograms[incoming, mask]):
                for outgoing in range(BITS + 1):
                    if caps[int(emitted)][outgoing] and int(emitted) < (
                        6
                        + int(potential[outgoing])
                        - int(potential[incoming])
                        - penalty * defects
                    ):
                        raise SystemExit("all-active minplus: defect inequality failed")
    return penalty, witness_row


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--defect-blocks", type=int, default=207)
    args = parser.parse_args()
    if not 0 <= args.defect_blocks <= INNER_BLOCKS:
        raise SystemExit("all-active minplus: invalid defect count")
    edges = edge_costs()
    if np.any(np.min(edges, axis=1) >= INF):
        raise SystemExit("all-active minplus: dead incoming state")
    mean = minimum_mean_cycle(edges)
    exact_cost = shortest_cost(edges, INNER_BLOCKS)
    potential = all_active_potential(edges)
    penalty, penalty_witness = defect_penalty(potential)
    deterministic_lower = (
        6 * INNER_BLOCKS
        - penalty * args.defect_blocks
        + int(np.min(potential))
        - int(potential[0])
    )
    maximum_safe_defects = (6 * INNER_BLOCKS - D - 1) // penalty
    print("packet-8 all-active min-plus trajectory probe")
    print(f"minimum_mean_cycle_emitted_weight={mean:.12f}")
    print(f"all_active_blocks={INNER_BLOCKS} shortest_emitted_weight={exact_cost}")
    print(f"potential={','.join(map(str, potential.tolist()))}")
    print(f"potential_range={int(np.max(potential)-np.min(potential))}")
    print(f"defect_packet_penalty={penalty}")
    print(f"defect_penalty_witness={penalty_witness}")
    print(
        f"zero_packet_defects={args.defect_blocks} "
        f"deterministic_emitted_lower_bound={deterministic_lower} target={D}"
    )
    print(f"maximum_deterministically_safe_zero_packets={maximum_safe_defects}")
    print(f"all_active_exceeds_target={exact_cost > D}")
    print(f"defect_trajectory_exceeds_target={deterministic_lower > D}")
    print("status=EXACT_INTEGER_SUPPORT_GRAPH_POTENTIAL_CERTIFICATE")


if __name__ == "__main__":
    main()
