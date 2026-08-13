#!/usr/bin/env python3
"""Analyze proof-relevant abstractions of tilted packet-8 inner transitions.

The input is a trace produced by
``probe_packet8_hard_profile_thermodynamic.py``.  All reported quantities are
empirical diagnostics.  In particular, conditional entropies and effective
cell counts help identify which local variables predict the next systematic
EBCH state, but they are not probability bounds.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter, defaultdict
from pathlib import Path


def entropy(counts) -> float:
    total = sum(counts)
    if not total:
        return 0.0
    return -sum(
        (count / total) * math.log2(count / total)
        for count in counts
        if count
    )


def conditional_entropy(
    counter: Counter[tuple[int, ...]], feature_indices: tuple[int, ...], outcome_index: int
) -> float:
    groups: dict[tuple[int, ...], Counter[int]] = defaultdict(Counter)
    for key, count in counter.items():
        feature = tuple(key[index] for index in feature_indices)
        groups[feature][key[outcome_index]] += count
    total = sum(counter.values())
    return sum(
        (sum(outcomes.values()) / total) * entropy(outcomes.values())
        for outcomes in groups.values()
    )


def effective_cells(counter: Counter[tuple[int, ...]]) -> float:
    total = sum(counter.values())
    denominator = sum(count * count for count in counter.values())
    return total * total / denominator if denominator else 0.0


def rows_to_counter(rows: list[dict], names: tuple[str, ...]) -> Counter[tuple[int, ...]]:
    return Counter(
        {
            tuple(int(row[name]) for name in names): int(row["count"])
            for row in rows
        }
    )


def outcome_multiplicity(
    counter: Counter[tuple[int, ...]], feature_indices: tuple[int, ...], outcome_index: int
) -> tuple[float, int, int]:
    groups: dict[tuple[int, ...], set[int]] = defaultdict(set)
    masses: Counter[tuple[int, ...]] = Counter()
    for key, count in counter.items():
        feature = tuple(key[index] for index in feature_indices)
        groups[feature].add(key[outcome_index])
        masses[feature] += count
    total = sum(masses.values())
    weighted = sum(masses[key] * len(groups[key]) for key in groups) / total
    return weighted, min(map(len, groups.values())), max(map(len, groups.values()))


def analyze_trace(trace: dict) -> dict[str, object]:
    state_names = ("state_weight", "emitted_weight", "next_state_weight")
    drive_names = (
        "state_weight",
        "input_weight",
        "overlap",
        "drive_weight",
        "emitted_weight",
        "next_state_weight",
    )
    packet_names = (
        "state_weight",
        "packet_weight_code",
        "emitted_weight",
        "next_state_weight",
    )
    concrete_names = ("state_weight", "input", "emitted_weight", "next_state_weight")
    state = rows_to_counter(trace["state_transition"], state_names)
    drive = rows_to_counter(trace["drive_transition"], drive_names)
    packet = rows_to_counter(trace["packet_transition_top"], packet_names)
    concrete = rows_to_counter(trace["concrete_transition_top"], concrete_names)
    observations = int(trace["observations"])

    next_histogram: Counter[int] = Counter()
    state_histogram: Counter[int] = Counter()
    emitted_histogram: Counter[int] = Counter()
    for (q, emitted, next_q), count in state.items():
        state_histogram[q] += count
        emitted_histogram[emitted] += count
        next_histogram[next_q] += count

    h_next = entropy(next_histogram.values())
    h_next_given_q = conditional_entropy(state, (0,), 2)
    h_next_given_q_y = conditional_entropy(state, (0, 1), 2)
    h_next_given_q_t_y = conditional_entropy(drive, (0, 1, 4), 5)
    h_next_given_q_drive_y = conditional_entropy(drive, (0, 3, 4), 5)
    h_next_given_overlap = conditional_entropy(drive, (0, 1, 2, 4), 5)
    packet_coverage = sum(packet.values()) / observations
    concrete_coverage = sum(concrete.values()) / observations
    h_next_given_packet = conditional_entropy(packet, (0, 1, 2), 3)
    h_next_given_concrete = conditional_entropy(concrete, (0, 1, 2), 3)
    weighted_outcomes, minimum_outcomes, maximum_outcomes = outcome_multiplicity(
        state, (0, 1), 2
    )

    return {
        "beta_index": trace["beta_index"],
        "beta": trace["beta"],
        "observations": observations,
        "mean_state_weight": sum(q * count for q, count in state_histogram.items())
        / observations,
        "mean_emitted_weight": sum(y * count for y, count in emitted_histogram.items())
        / observations,
        "mean_next_state_weight": sum(q * count for q, count in next_histogram.items())
        / observations,
        "state_weight_range": [min(state_histogram), max(state_histogram)],
        "emitted_weight_range": [min(emitted_histogram), max(emitted_histogram)],
        "next_state_weight_range": [min(next_histogram), max(next_histogram)],
        "entropy_next_state_weight": h_next,
        "conditional_entropy": {
            "given_state_weight": h_next_given_q,
            "given_state_and_emitted_weight": h_next_given_q_y,
            "given_state_input_and_emitted_weight": h_next_given_q_t_y,
            "given_state_drive_and_emitted_weight": h_next_given_q_drive_y,
            "given_state_input_overlap_and_emitted_weight": h_next_given_overlap,
            "given_state_ordered_packet_weights_and_emitted_weight": h_next_given_packet,
            "given_state_concrete_input_and_emitted_weight": h_next_given_concrete,
        },
        "conditional_information_beyond_state_emitted": {
            "input_weight": h_next_given_q_y - h_next_given_q_t_y,
            "drive_weight": h_next_given_q_y - h_next_given_q_drive_y,
            "input_weight_and_overlap": h_next_given_q_y - h_next_given_overlap,
            "ordered_packet_weights": h_next_given_q_y - h_next_given_packet,
            "concrete_input": h_next_given_q_y - h_next_given_concrete,
        },
        "coverage": {
            "packet_top": packet_coverage,
            "concrete_top": concrete_coverage,
        },
        "effective_cells": {
            "state_transition": effective_cells(state),
            "drive_transition": effective_cells(drive),
            "packet_transition_top": effective_cells(packet),
            "concrete_transition_top": effective_cells(concrete),
        },
        "next_state_outcomes_per_state_emitted_cell": {
            "mass_weighted_mean": weighted_outcomes,
            "minimum": minimum_outcomes,
            "maximum": maximum_outcomes,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trace", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    data = json.loads(args.trace.read_text(encoding="utf-8"))
    traces = data.get("transition_traces", [])
    if not traces:
        raise SystemExit("tilted transition analysis: input contains no traces")
    rows = [analyze_trace(trace) for trace in traces]

    print("packet-8 tilted transition abstraction analysis")
    print(f"source={args.trace}")
    for row in rows:
        conditional = row["conditional_entropy"]
        information = row["conditional_information_beyond_state_emitted"]
        coverage = row["coverage"]
        effective = row["effective_cells"]
        outcomes = row["next_state_outcomes_per_state_emitted_cell"]
        print(
            f"beta={row['beta']:.9f} observations={row['observations']} "
            f"mean_q={row['mean_state_weight']:.6f} "
            f"mean_y={row['mean_emitted_weight']:.6f} "
            f"H(qnext|q,y)={conditional['given_state_and_emitted_weight']:.6f} "
            f"info_drive={information['drive_weight']:.6f} "
            f"info_overlap={information['input_weight_and_overlap']:.6f} "
            f"info_packet={information['ordered_packet_weights']:.6f} "
            f"info_concrete={information['concrete_input']:.6f} "
            f"packet_coverage={coverage['packet_top']:.6f} "
            f"concrete_coverage={coverage['concrete_top']:.6f} "
            f"effective_state_cells={effective['state_transition']:.3f} "
            f"qnext_outcomes_mean={outcomes['mass_weighted_mean']:.3f}"
        )
    print("status=EMPIRICAL_ABSTRACTION_DIAGNOSTIC_NOT_A_BOUND")

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps({"source": str(args.trace), "rows": rows}, indent=2) + "\n",
            encoding="utf-8",
        )
        print(f"output={args.output}")


if __name__ == "__main__":
    main()
