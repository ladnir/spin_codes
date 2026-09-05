#!/usr/bin/env python3
"""Find low-weight paper-EC rows caused by expander draw collisions.

This is an ensemble-level experiment for the archived single-pass construction.
It does not reproduce the libOTe AES-based PRNG bit for bit.  It samples the
same idealized distribution: nominal expander edges are independent uniform
draws, and every non-fixed convolution tap is an independent uniform bit.
"Wrapping" fixes the oldest feedback tap to one; block boundaries remain
non-cyclic, as in ExConvCodeOld.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np


EXPANDER_STREAM_TAG = 0x455850414E444552
CONVOLUTION_STREAM_TAG = 0x434F4E564F4C5554


@dataclass(frozen=True)
class CollisionCandidate:
    row: int
    nominal_draws: tuple[int, ...]
    effective_support: tuple[int, ...]

    @property
    def effective_weight(self) -> int:
        return len(self.effective_support)


@dataclass(frozen=True)
class CollisionWitness:
    row: int
    nominal_draws: tuple[int, ...]
    effective_support: tuple[int, ...]
    parity_weight: int
    codeword_weight: int
    relative_weight: float


def effective_support(draws: Iterable[int]) -> tuple[int, ...]:
    """Return the coordinates with odd multiplicity in one binary row."""
    parity: set[int] = set()
    for value in draws:
        coordinate = int(value)
        if coordinate in parity:
            parity.remove(coordinate)
        else:
            parity.add(coordinate)
    return tuple(sorted(parity))


def weight_three_probability_five_draws(domain_size: int) -> float:
    """Return Pr[effective weight 3] for five draws with replacement."""
    if domain_size < 1:
        raise ValueError("domain_size must be positive")
    q = domain_size
    count = 60 * math.comb(q, 3) + 240 * math.comb(q, 4)
    return count / q**5


def weight_one_probability_five_draws(domain_size: int) -> float:
    """Return Pr[effective weight 1] for five draws with replacement."""
    if domain_size < 1:
        raise ValueError("domain_size must be positive")
    q = domain_size
    count = q + 15 * q * (q - 1) + 30 * q * math.comb(q - 1, 2)
    return count / q**5


def _rng(seed: int, stream_tag: int) -> np.random.Generator:
    sequence = np.random.SeedSequence([seed, stream_tag])
    return np.random.Generator(np.random.PCG64(sequence))


def find_collision_candidates(
    *,
    message_rows: int,
    domain_size: int,
    nominal_weight: int,
    seed: int,
    chunk_rows: int = 1 << 18,
) -> list[CollisionCandidate]:
    """Stream all nominal rows and retain rows whose draws cancel."""
    if message_rows < 1 or domain_size < 1 or nominal_weight < 1:
        raise ValueError("row count, domain size, and nominal weight must be positive")
    if chunk_rows < 1:
        raise ValueError("chunk_rows must be positive")

    generator = _rng(seed, EXPANDER_STREAM_TAG)
    dtype = np.uint32 if domain_size <= np.iinfo(np.uint32).max else np.uint64
    candidates: list[CollisionCandidate] = []

    for start in range(0, message_rows, chunk_rows):
        count = min(chunk_rows, message_rows - start)
        draws = generator.integers(
            0,
            domain_size,
            size=(count, nominal_weight),
            dtype=dtype,
        )
        draws.sort(axis=1)
        duplicate_rows = np.flatnonzero(np.any(draws[:, 1:] == draws[:, :-1], axis=1))
        for local_row in duplicate_rows:
            nominal = tuple(int(value) for value in draws[local_row])
            support = effective_support(nominal)
            if len(support) < nominal_weight:
                candidates.append(
                    CollisionCandidate(
                        row=start + int(local_row),
                        nominal_draws=nominal,
                        effective_support=support,
                    )
                )
    return candidates


def generate_wrapping_taps(*, length: int, memory: int, seed: int) -> list[int]:
    """Sample source-indexed taps with the oldest tap fixed to one."""
    if length < 1:
        raise ValueError("length must be positive")
    if not (1 <= memory <= 63):
        raise ValueError("memory must lie in [1,63]")
    generator = _rng(seed, CONVOLUTION_STREAM_TAG)
    random_tap_count = memory - 1
    if random_tap_count:
        taps = generator.integers(
            0,
            1 << random_tap_count,
            size=length,
            dtype=np.uint64,
        )
    else:
        taps = np.zeros(length, dtype=np.uint64)
    taps |= np.uint64(1 << (memory - 1))
    return [int(value) for value in taps]


def simulate_wrapping_candidates(
    supports: Sequence[Sequence[int]],
    taps: Sequence[int],
    *,
    memory: int,
) -> list[int]:
    """Return output weights for up to 63 inputs under one shared convolution.

    The computation follows the archived scatter implementation.  At source
    coordinate i, every active candidate is XORed into the selected future
    coordinates.  The fixed oldest tap is included in ``taps[i]``.  Updates
    past the block boundary are discarded rather than wrapped to its start.
    """
    length = len(taps)
    if not (1 <= len(supports) <= 63):
        raise ValueError("simulate between one and 63 candidates at a time")
    if not (1 <= memory <= 63):
        raise ValueError("memory must lie in [1,63]")

    words = [0] * length
    for candidate, support in enumerate(supports):
        bit = 1 << candidate
        for coordinate in support:
            position = int(coordinate)
            if not (0 <= position < length):
                raise ValueError("support coordinate outside convolution domain")
            words[position] ^= bit

    weights = [0] * len(supports)
    for source, tap_mask in enumerate(taps):
        active = words[source]
        if not active:
            continue

        pending_active = active
        while pending_active:
            bit = pending_active & -pending_active
            weights[bit.bit_length() - 1] += 1
            pending_active ^= bit

        pending_taps = int(tap_mask)
        while pending_taps:
            tap_bit = pending_taps & -pending_taps
            offset = tap_bit.bit_length()
            target = source + offset
            if target < length:
                words[target] ^= active
            pending_taps ^= tap_bit
    return weights


def evaluate_candidates(
    candidates: Sequence[CollisionCandidate],
    taps: Sequence[int],
    *,
    memory: int,
    systematic_length: int,
) -> list[CollisionWitness]:
    """Evaluate collision candidates in bit-sliced batches."""
    witnesses: list[CollisionWitness] = []
    for start in range(0, len(candidates), 63):
        batch = candidates[start : start + 63]
        parity_weights = simulate_wrapping_candidates(
            [candidate.effective_support for candidate in batch],
            taps,
            memory=memory,
        )
        for candidate, parity_weight in zip(batch, parity_weights, strict=True):
            codeword_weight = 1 + parity_weight
            witnesses.append(
                CollisionWitness(
                    row=candidate.row,
                    nominal_draws=candidate.nominal_draws,
                    effective_support=candidate.effective_support,
                    parity_weight=parity_weight,
                    codeword_weight=codeword_weight,
                    relative_weight=codeword_weight / systematic_length,
                )
            )
    return witnesses


def run_trial(
    *,
    k: int,
    nominal_weight: int,
    memory: int,
    seed: int,
    chunk_rows: int,
) -> dict[str, object]:
    candidates = find_collision_candidates(
        message_rows=k,
        domain_size=k,
        nominal_weight=nominal_weight,
        seed=seed,
        chunk_rows=chunk_rows,
    )
    taps = generate_wrapping_taps(length=k, memory=memory, seed=seed)
    witnesses = evaluate_candidates(
        candidates,
        taps,
        memory=memory,
        systematic_length=2 * k,
    ) if candidates else []
    witnesses.sort(key=lambda witness: (witness.codeword_weight, witness.row))

    return {
        "seed": seed,
        "k": k,
        "n": 2 * k,
        "nominal_weight": nominal_weight,
        "memory": memory,
        "boundary": "zero",
        "wrapping": "oldest feedback tap fixed to one",
        "collision_candidates": len(candidates),
        "effective_weight_histogram": {
            str(weight): sum(candidate.effective_weight == weight for candidate in candidates)
            for weight in sorted({candidate.effective_weight for candidate in candidates})
        },
        "best_witness": asdict(witnesses[0]) if witnesses else None,
        "witnesses": [asdict(witness) for witness in witnesses],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--k-log2", type=int, default=20)
    parser.add_argument("--nominal-weight", type=int, default=5)
    parser.add_argument("--memory", type=int, default=25)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--trials", type=int, default=1)
    parser.add_argument("--chunk-rows", type=int, default=1 << 18)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    if not (1 <= args.k_log2 < 63):
        raise SystemExit("--k-log2 must lie in [1,62]")
    if args.trials < 1:
        raise SystemExit("--trials must be positive")
    k = 1 << args.k_log2

    trials = [
        run_trial(
            k=k,
            nominal_weight=args.nominal_weight,
            memory=args.memory,
            seed=args.seed + trial,
            chunk_rows=args.chunk_rows,
        )
        for trial in range(args.trials)
    ]
    payload = {
        "experiment": "paper-ec-collision-finder",
        "distribution": {
            "expander": "independent uniform draws with replacement",
            "convolution": "independent random taps with fixed oldest tap",
            "prng": "NumPy PCG64; ensemble-equivalent, not libOTe-replayable",
        },
        "trials": trials,
    }
    encoded = json.dumps(payload, indent=2, sort_keys=True)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded + "\n", encoding="utf-8")
    print(encoded)


if __name__ == "__main__":
    main()
