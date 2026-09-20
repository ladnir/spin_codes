#!/usr/bin/env python3
"""Search BCH [256,128] XOR circuits by transpose cost and live pressure."""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass
from pathlib import Path

import build_bch256_128_eq3_circuit as build
import probe_bch_forward_xor_circuit as paar


@dataclass(frozen=True)
class Candidate:
    seed: int
    transpose_xors: int
    forward_xors: int
    descending_peak_live: int
    descending_live_area: int
    greedy_peak_live: int
    greedy_live_area: int
    greedy_schedule_seed: int


def dominates(left: Candidate, right: Candidate) -> bool:
    left_values = (
        left.transpose_xors,
        left.descending_peak_live,
        left.descending_live_area,
    )
    right_values = (
        right.transpose_xors,
        right.descending_peak_live,
        right.descending_live_area,
    )
    return all(a <= b for a, b in zip(left_values, right_values)) and any(
        a < b for a, b in zip(left_values, right_values)
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--trials", type=int, default=128)
    parser.add_argument("--schedule-trials", type=int, default=8)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    rows = build.build_rows()
    target = build.columns(rows)
    paar.DIMENSION = build.DIMENSION
    paar.LENGTH = build.LENGTH
    candidates: list[Candidate] = []
    for seed in range(args.seed, args.seed + args.trials):
        circuit = paar.synthesize(seed, 2, target)
        descending = paar.descending_transpose_schedule(circuit)
        descending_peak, descending_area = paar.transpose_live_stats(
            circuit, descending
        )
        greedy_seed, greedy_peak, greedy_area = min(
            (
                (
                    schedule_seed,
                    *paar.transpose_live_stats(
                        circuit,
                        paar.greedy_transpose_schedule(circuit, schedule_seed),
                    ),
                )
                for schedule_seed in range(args.schedule_trials)
            ),
            key=lambda item: (item[1], item[2]),
        )
        candidates.append(
            Candidate(
                seed=seed,
                transpose_xors=paar.transpose_xor_count(circuit),
                forward_xors=circuit.xor_count,
                descending_peak_live=descending_peak,
                descending_live_area=descending_area,
                greedy_peak_live=greedy_peak,
                greedy_live_area=greedy_area,
                greedy_schedule_seed=greedy_seed,
            )
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(output, fieldnames=list(asdict(candidates[0])))
        writer.writeheader()
        writer.writerows(asdict(candidate) for candidate in candidates)

    pareto = [
        candidate
        for candidate in candidates
        if not any(
            dominates(other, candidate)
            for other in candidates
            if other is not candidate
        )
    ]
    print(f"trials,{len(candidates)}")
    print(f"pareto_candidates,{len(pareto)}")
    print(
        "seed,transpose_xors,descending_peak_live,descending_live_area,"
        "greedy_peak_live,greedy_live_area"
    )
    for candidate in sorted(
        pareto,
        key=lambda item: (
            item.transpose_xors,
            item.descending_peak_live,
            item.descending_live_area,
        ),
    ):
        print(
            f"{candidate.seed},{candidate.transpose_xors},"
            f"{candidate.descending_peak_live},"
            f"{candidate.descending_live_area},"
            f"{candidate.greedy_peak_live},{candidate.greedy_live_area}"
        )
    print(f"output,{args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
