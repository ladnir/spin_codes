#!/usr/bin/env python3
"""Attack the paper-like exact-row Expand--Convolute ensemble.

The systematic code is ``[I | B C]`` at rate one half.  Every row of ``B`` is
a uniform subset of fixed size, sampled without replacement.  ``C`` is one
wrapping convolution pass: its oldest feedback tap is fixed to one and its
other taps are independent uniform bits.

The attack partitions the parity block into contiguous regions and groups
expander rows by their region multiset.  A packed reverse pass through ``C``
maps each expander coordinate to every region-boundary state.  Gaussian
elimination then finds combinations whose boundary states vanish.  Their
parity outputs are confined to the selected regions.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

from ec_collision_finder import _rng, generate_wrapping_taps


EXACT_ROW_STREAM_TAG = 0x4558414354524F57
KERNEL_STREAM_TAG = 0x4B45524E454C5354


@dataclass(frozen=True)
class Bucket:
    key: int
    regions: tuple[int, ...]
    row_count: int
    constraint_count: int
    allowed_coordinates: int
    start: int


@dataclass(frozen=True)
class RegionWitness:
    method: str
    bucket_key: int
    regions: tuple[int, ...]
    message_rows: tuple[int, ...]
    expander_support: tuple[int, ...]
    message_weight: int
    parity_weight: int
    codeword_weight: int
    relative_weight: float
    outside_region_weight: int
    codeword_sha256: str


def sample_exact_rows(
    *,
    row_count: int,
    domain_size: int,
    row_weight: int,
    seed: int,
) -> np.ndarray:
    """Sample independent uniform row subsets without replacement."""
    if row_count < 1 or domain_size < 1 or row_weight < 1:
        raise ValueError("row count, domain size, and row weight must be positive")
    if row_weight > domain_size:
        raise ValueError("row weight exceeds the expander domain")

    generator = _rng(seed, EXACT_ROW_STREAM_TAG)
    dtype = np.uint32 if domain_size <= np.iinfo(np.uint32).max else np.uint64
    rows = generator.integers(
        0,
        domain_size,
        size=(row_count, row_weight),
        dtype=dtype,
    )

    # Reject an entire ordered tuple when it contains a duplicate.  Conditioned
    # on acceptance, the tuple is a uniform ordered sample without replacement.
    while True:
        ordered = np.sort(rows, axis=1)
        duplicate = np.any(ordered[:, 1:] == ordered[:, :-1], axis=1)
        bad = np.flatnonzero(duplicate)
        if not len(bad):
            return ordered
        rows[bad] = generator.integers(
            0,
            domain_size,
            size=(len(bad), row_weight),
            dtype=dtype,
        )


def region_boundaries(length: int, region_count: int) -> tuple[int, ...]:
    """Return boundaries compatible with ``floor(position*regions/length)``."""
    if length < 1 or not (1 <= region_count <= length):
        raise ValueError("invalid region partition")
    return tuple(
        0 if region == 0 else (region * length + region_count - 1) // region_count
        for region in range(region_count + 1)
    )


def encode_region_keys(rows: np.ndarray, region_count: int, domain_size: int) -> np.ndarray:
    """Encode each sorted region multiset as one injective integer key."""
    if rows.ndim != 2 or rows.shape[1] < 1:
        raise ValueError("rows must be a nonempty matrix")
    row_weight = rows.shape[1]
    if region_count**row_weight > np.iinfo(np.uint64).max:
        raise ValueError("region signature does not fit in 64 bits")
    regions = (rows.astype(np.uint64) * np.uint64(region_count)) // np.uint64(domain_size)
    powers = np.array([region_count**index for index in range(row_weight)], dtype=np.uint64)
    return np.sum(regions * powers, axis=1, dtype=np.uint64)


def decode_region_key(key: int, region_count: int, row_weight: int) -> tuple[int, ...]:
    regions: list[int] = []
    value = int(key)
    for _ in range(row_weight):
        regions.append(value % region_count)
        value //= region_count
    return tuple(regions)


def sort_buckets(
    rows: np.ndarray,
    *,
    region_count: int,
    domain_size: int,
    memory: int,
) -> tuple[np.ndarray, list[Bucket]]:
    """Sort row indices by signature and return dependency-capable buckets."""
    keys = encode_region_keys(rows, region_count, domain_size)
    order = np.argsort(keys, kind="stable")
    sorted_keys = keys[order]
    changes = np.flatnonzero(sorted_keys[1:] != sorted_keys[:-1]) + 1
    starts = np.concatenate((np.array([0]), changes))
    ends = np.concatenate((changes, np.array([len(sorted_keys)])))

    buckets: list[Bucket] = []
    row_weight = rows.shape[1]
    boundaries = region_boundaries(domain_size, region_count)
    for start, end in zip(starts, ends, strict=True):
        key = int(sorted_keys[start])
        regions = decode_region_key(key, region_count, row_weight)
        selected = set(regions)
        constrained = sum(
            region < region_count - 1 and region + 1 not in selected
            for region in selected
        )
        constraint_count = constrained * memory
        row_count = int(end - start)
        if row_count > constraint_count:
            allowed_coordinates = sum(
                boundaries[region + 1] - boundaries[region]
                for region in selected
            )
            buckets.append(
                Bucket(
                    key=key,
                    regions=regions,
                    row_count=row_count,
                    constraint_count=constraint_count,
                    allowed_coordinates=allowed_coordinates,
                    start=int(start),
                )
            )
    buckets.sort(
        key=lambda bucket: (
            bucket.allowed_coordinates,
            -bucket.row_count,
            bucket.key,
        )
    )
    return order, buckets


def boundary_response_labels(
    taps: Sequence[int],
    *,
    memory: int,
    boundaries: Sequence[int],
) -> list[int]:
    """Map every input coordinate to all nonterminal boundary-state bits.

    The forward convolution performs operations ``x[target] ^= x[source]``.
    Reversing those operations propagates packed output observations back to
    the input coordinates.  Bit ``r*memory+j`` observes output coordinate
    ``boundaries[r+1]-memory+j``.
    """
    length = len(taps)
    if boundaries[0] != 0 or boundaries[-1] != length:
        raise ValueError("boundaries do not cover the convolution block")
    if any(right - left < memory for left, right in zip(boundaries, boundaries[1:])):
        raise ValueError("each region must contain at least one convolution state")

    labels = [0] * length
    for region, end in enumerate(boundaries[1:-1]):
        first = end - memory
        bit_base = region * memory
        for offset in range(memory):
            labels[first + offset] ^= 1 << (bit_base + offset)

    for source in range(length - 1, -1, -1):
        pending = int(taps[source])
        response = labels[source]
        while pending:
            tap_bit = pending & -pending
            target = source + tap_bit.bit_length()
            if target < length:
                response ^= labels[target]
            pending ^= tap_bit
        labels[source] = response
    return labels


def constraint_mask(regions: Iterable[int], *, region_count: int, memory: int) -> int:
    mask = 0
    state_mask = (1 << memory) - 1
    selected = set(int(value) for value in regions)
    for region in selected:
        # A state constraint is necessary only before an unselected gap.
        # Adjacent selected regions form one allowed output interval.
        if region < region_count - 1 and region + 1 not in selected:
            mask |= state_mask << (region * memory)
    return mask


def gf2_kernel_dependencies(columns: Sequence[int]) -> list[int]:
    """Return an independent basis of dependencies among packed columns."""
    pivots: dict[int, tuple[int, int]] = {}
    dependencies: list[int] = []
    for column_index, value in enumerate(columns):
        vector = int(value)
        combination = 1 << column_index
        while vector:
            pivot = vector.bit_length() - 1
            prior = pivots.get(pivot)
            if prior is None:
                pivots[pivot] = (vector, combination)
                break
            vector ^= prior[0]
            combination ^= prior[1]
        if not vector:
            dependencies.append(combination)
    return dependencies


def augment_kernel_candidates(
    dependencies: Sequence[int],
    *,
    sample_count: int,
    seed: int,
    bucket_key: int,
) -> list[int]:
    """Add deterministic random combinations of a kernel basis."""
    candidates = set(int(value) for value in dependencies)
    if sample_count and dependencies:
        generator = _rng(seed ^ bucket_key, KERNEL_STREAM_TAG)
        dimension = len(dependencies)
        word_count = (dimension + 63) // 64
        for _ in range(sample_count):
            selector_words = generator.integers(
                0,
                1 << 64,
                size=word_count,
                dtype=np.uint64,
            )
            combination = 0
            selected = False
            for index, dependency in enumerate(dependencies):
                if int(selector_words[index // 64]) & (1 << (index % 64)):
                    combination ^= int(dependency)
                    selected = True
            if selected:
                candidates.add(combination)
    return sorted(candidates, key=lambda value: (value.bit_count(), value))


def selected_indices(mask: int) -> tuple[int, ...]:
    indices: list[int] = []
    pending = int(mask)
    while pending:
        bit = pending & -pending
        indices.append(bit.bit_length() - 1)
        pending ^= bit
    return tuple(indices)


def xor_row_support(rows: np.ndarray, row_indices: Sequence[int]) -> tuple[int, ...]:
    support: set[int] = set()
    for row_index in row_indices:
        for coordinate in rows[int(row_index)]:
            value = int(coordinate)
            if value in support:
                support.remove(value)
            else:
                support.add(value)
    return tuple(sorted(support))


def merge_region_intervals(
    regions: Iterable[int], boundaries: Sequence[int]
) -> tuple[tuple[int, int], ...]:
    intervals: list[tuple[int, int]] = []
    for region in sorted(set(int(value) for value in regions)):
        left, right = boundaries[region], boundaries[region + 1]
        if intervals and intervals[-1][1] == left:
            intervals[-1] = (intervals[-1][0], right)
        else:
            intervals.append((left, right))
    return tuple(intervals)


def simulate_candidate_batch(
    supports: Sequence[Sequence[int]],
    taps: Sequence[int],
    *,
    allowed_intervals: Sequence[tuple[int, int]],
) -> tuple[list[int], list[int]]:
    """Return total and outside-region weights for at most 63 candidates."""
    if not (1 <= len(supports) <= 63):
        raise ValueError("simulate between one and 63 candidates")
    length = len(taps)
    words = [0] * length
    for candidate, support in enumerate(supports):
        bit = 1 << candidate
        for coordinate in support:
            words[int(coordinate)] ^= bit

    total = [0] * len(supports)
    outside = [0] * len(supports)
    interval_index = 0
    for source in range(length):
        while (
            interval_index < len(allowed_intervals)
            and source >= allowed_intervals[interval_index][1]
        ):
            interval_index += 1
        allowed = (
            interval_index < len(allowed_intervals)
            and allowed_intervals[interval_index][0] <= source
        )

        active = words[source]
        pending_active = active
        while pending_active:
            bit = pending_active & -pending_active
            index = bit.bit_length() - 1
            total[index] += 1
            if not allowed:
                outside[index] += 1
            pending_active ^= bit

        if active:
            pending_taps = int(taps[source])
            while pending_taps:
                tap_bit = pending_taps & -pending_taps
                target = source + tap_bit.bit_length()
                if target < length:
                    words[target] ^= active
                pending_taps ^= tap_bit
    return total, outside


def materialize_witness(
    *,
    message_rows: Sequence[int],
    expander_support: Sequence[int],
    taps: Sequence[int],
    allowed_intervals: Sequence[tuple[int, int]],
) -> tuple[int, int, str]:
    """Recompute one witness and hash its packed systematic codeword."""
    length = len(taps)
    values = bytearray(length)
    for coordinate in expander_support:
        values[int(coordinate)] ^= 1
    for source, tap_mask in enumerate(taps):
        if not values[source]:
            continue
        pending = int(tap_mask)
        while pending:
            tap_bit = pending & -pending
            target = source + tap_bit.bit_length()
            if target < length:
                values[target] ^= 1
            pending ^= tap_bit

    parity_weight = sum(values)
    outside_weight = 0
    interval_index = 0
    for position, value in enumerate(values):
        while (
            interval_index < len(allowed_intervals)
            and position >= allowed_intervals[interval_index][1]
        ):
            interval_index += 1
        allowed = (
            interval_index < len(allowed_intervals)
            and allowed_intervals[interval_index][0] <= position
        )
        if value and not allowed:
            outside_weight += 1

    message_packed = bytearray((length + 7) // 8)
    for row in message_rows:
        message_packed[int(row) // 8] ^= 1 << (int(row) % 8)
    parity_packed = np.packbits(
        np.frombuffer(values, dtype=np.uint8),
        bitorder="little",
    ).tobytes()
    digest = hashlib.sha256(bytes(message_packed) + parity_packed).hexdigest()
    return parity_weight, outside_weight, digest


def evaluate_late_rows(
    rows: np.ndarray,
    taps: Sequence[int],
    *,
    region_count: int,
    candidate_count: int,
) -> list[RegionWitness]:
    """Evaluate rows with the latest first expander coordinate."""
    if candidate_count < 1:
        return []
    candidate_count = min(candidate_count, len(rows))
    if candidate_count == len(rows):
        row_indices = np.arange(len(rows))
    else:
        row_indices = np.argpartition(rows[:, 0], -candidate_count)[-candidate_count:]
    row_indices = row_indices[np.argsort(rows[row_indices, 0])[::-1]]
    supports = [tuple(int(value) for value in rows[int(row)]) for row in row_indices]

    witnesses: list[RegionWitness] = []
    for start in range(0, len(supports), 63):
        batch_supports = supports[start : start + 63]
        first = min(support[0] for support in batch_supports)
        parity_weights, outside_weights = simulate_candidate_batch(
            batch_supports,
            taps,
            allowed_intervals=((first, len(taps)),),
        )
        for offset, (parity_weight, outside_weight) in enumerate(
            zip(parity_weights, outside_weights, strict=True)
        ):
            index = start + offset
            support = supports[index]
            # The batch interval begins no later than every candidate support.
            # Recheck the stronger candidate-specific confinement condition.
            if outside_weight:
                raise AssertionError("late-row output appeared before its first input")
            codeword_weight = 1 + parity_weight
            witnesses.append(
                RegionWitness(
                    method="late-row-sweep",
                    bucket_key=-1,
                    regions=tuple(
                        coordinate * region_count // len(taps)
                        for coordinate in support
                    ),
                    message_rows=(int(row_indices[index]),),
                    expander_support=support,
                    message_weight=1,
                    parity_weight=parity_weight,
                    codeword_weight=codeword_weight,
                    relative_weight=codeword_weight / (2 * len(taps)),
                    outside_region_weight=outside_weight,
                    codeword_sha256="",
                )
            )
    witnesses.sort(key=lambda witness: (witness.codeword_weight, witness.message_rows))
    return witnesses


def attack_bucket(
    *,
    bucket: Bucket,
    bucket_rows: np.ndarray,
    rows: np.ndarray,
    labels: Sequence[int],
    taps: Sequence[int],
    boundaries: Sequence[int],
    region_count: int,
    memory: int,
    seed: int,
    kernel_samples: int,
) -> tuple[list[RegionWitness], dict[str, object]]:
    mask = constraint_mask(bucket.regions, region_count=region_count, memory=memory)
    columns = [
        (
            int(labels[int(rows[int(row), 0])])
            ^ int(labels[int(rows[int(row), 1])])
            ^ int(labels[int(rows[int(row), 2])])
            ^ int(labels[int(rows[int(row), 3])])
            ^ int(labels[int(rows[int(row), 4])])
        )
        & mask
        for row in bucket_rows
    ]
    dependencies = gf2_kernel_dependencies(columns)
    candidates = augment_kernel_candidates(
        dependencies,
        sample_count=kernel_samples,
        seed=seed,
        bucket_key=bucket.key,
    )
    local_selections = [selected_indices(candidate) for candidate in candidates]
    message_rows = [
        tuple(int(bucket_rows[index]) for index in selection)
        for selection in local_selections
    ]
    supports = [xor_row_support(rows, selection) for selection in message_rows]
    intervals = merge_region_intervals(bucket.regions, boundaries)

    witnesses: list[RegionWitness] = []
    for start in range(0, len(candidates), 63):
        batch_supports = supports[start : start + 63]
        parity_weights, outside_weights = simulate_candidate_batch(
            batch_supports,
            taps,
            allowed_intervals=intervals,
        )
        for offset, (parity_weight, outside_weight) in enumerate(
            zip(parity_weights, outside_weights, strict=True)
        ):
            index = start + offset
            if outside_weight:
                raise AssertionError("boundary-nullspace witness escaped its regions")
            message_weight = len(message_rows[index])
            codeword_weight = message_weight + parity_weight
            witnesses.append(
                RegionWitness(
                    method="region-nullspace",
                    bucket_key=bucket.key,
                    regions=bucket.regions,
                    message_rows=message_rows[index],
                    expander_support=supports[index],
                    message_weight=message_weight,
                    parity_weight=parity_weight,
                    codeword_weight=codeword_weight,
                    relative_weight=codeword_weight / (2 * len(taps)),
                    outside_region_weight=outside_weight,
                    codeword_sha256="",
                )
            )

    witnesses.sort(key=lambda witness: (witness.codeword_weight, witness.message_rows))
    metadata = {
        "bucket": asdict(bucket),
        "kernel_dimension": len(dependencies),
        "evaluated_kernel_vectors": len(candidates),
    }
    return witnesses, metadata


def run_attack(
    *,
    k: int,
    row_weight: int,
    memory: int,
    region_count: int,
    seed: int,
    max_buckets: int,
    kernel_samples: int,
    late_row_candidates: int,
) -> dict[str, object]:
    if row_weight != 5:
        raise ValueError("the optimized boundary signature currently requires row weight five")
    boundaries = region_boundaries(k, region_count)
    if min(right - left for left, right in zip(boundaries, boundaries[1:])) < memory:
        raise ValueError("regions must be at least as long as the convolution memory")

    rows = sample_exact_rows(
        row_count=k,
        domain_size=k,
        row_weight=row_weight,
        seed=seed,
    )
    taps = generate_wrapping_taps(length=k, memory=memory, seed=seed)
    late_row_witnesses = evaluate_late_rows(
        rows,
        taps,
        region_count=region_count,
        candidate_count=late_row_candidates,
    )
    order, buckets = sort_buckets(
        rows,
        region_count=region_count,
        domain_size=k,
        memory=memory,
    )
    labels = boundary_response_labels(taps, memory=memory, boundaries=boundaries)

    region_witnesses: list[RegionWitness] = []
    bucket_receipts: list[dict[str, object]] = []
    for bucket in buckets[:max_buckets]:
        bucket_rows = order[bucket.start : bucket.start + bucket.row_count]
        witnesses, metadata = attack_bucket(
            bucket=bucket,
            bucket_rows=bucket_rows,
            rows=rows,
            labels=labels,
            taps=taps,
            boundaries=boundaries,
            region_count=region_count,
            memory=memory,
            seed=seed,
            kernel_samples=kernel_samples,
        )
        bucket_receipts.append(metadata)
        region_witnesses.extend(witnesses)

    all_witnesses = late_row_witnesses + region_witnesses
    all_witnesses.sort(key=lambda witness: (witness.codeword_weight, witness.message_rows))
    best = all_witnesses[0] if all_witnesses else None
    late_best = late_row_witnesses[0] if late_row_witnesses else None
    region_witnesses.sort(
        key=lambda witness: (witness.codeword_weight, witness.message_rows)
    )
    region_best = region_witnesses[0] if region_witnesses else None

    def finalize(witness: RegionWitness | None) -> RegionWitness | None:
        if witness is None:
            return None
        intervals = (
            ((witness.expander_support[0], k),)
            if witness.method == "late-row-sweep"
            else merge_region_intervals(witness.regions, boundaries)
        )
        parity_weight, outside_weight, digest = materialize_witness(
            message_rows=witness.message_rows,
            expander_support=witness.expander_support,
            taps=taps,
            allowed_intervals=intervals,
        )
        if parity_weight != witness.parity_weight or outside_weight:
            raise AssertionError("materialized witness does not match bit-sliced evaluation")
        return RegionWitness(
            **{
                **asdict(witness),
                "codeword_sha256": digest,
            }
        )

    best = finalize(best)
    late_best = best if best is not None and best.method == "late-row-sweep" else finalize(late_best)
    region_best = best if best is not None and best.method == "region-nullspace" else finalize(region_best)

    return {
        "seed": seed,
        "k": k,
        "n": 2 * k,
        "row_weight": row_weight,
        "memory": memory,
        "regions": region_count,
        "boundaries": boundaries,
        "qualifying_buckets": len(buckets),
        "processed_buckets": len(bucket_receipts),
        "late_row_candidates": len(late_row_witnesses),
        "late_row_best": asdict(late_best) if late_best else None,
        "region_best": asdict(region_best) if region_best else None,
        "bucket_receipts": bucket_receipts,
        "best_witness": asdict(best) if best is not None else None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--k-log2", type=int, default=20)
    parser.add_argument("--row-weight", type=int, default=5)
    parser.add_argument("--memory", type=int, default=25)
    parser.add_argument("--regions", type=int, default=16)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--max-buckets", type=int, default=16)
    parser.add_argument("--kernel-samples", type=int, default=32)
    parser.add_argument("--late-row-candidates", type=int, default=256)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    if not (1 <= args.k_log2 < 63):
        raise SystemExit("--k-log2 must lie in [1,62]")
    if args.max_buckets < 1 or args.kernel_samples < 0 or args.late_row_candidates < 0:
        raise SystemExit("attack work counts must be nonnegative and bucket count positive")

    result = run_attack(
        k=1 << args.k_log2,
        row_weight=args.row_weight,
        memory=args.memory,
        region_count=args.regions,
        seed=args.seed,
        max_buckets=args.max_buckets,
        kernel_samples=args.kernel_samples,
        late_row_candidates=args.late_row_candidates,
    )
    payload = {
        "experiment": "paper-exact-row-ec-region-nullspace",
        "ensemble": {
            "code": "systematic [I | B C] at rate one half",
            "expander": "independent uniform fixed-size row subsets",
            "convolution": "one pass, fixed oldest feedback tap, zero initial state",
            "prng": "NumPy PCG64 ensemble sampler",
        },
        "attack": {
            "partition": "equal contiguous regions",
            "bucket_key": "unordered region multiset of each expander row",
            "constraints": "zero convolution state before every unselected gap",
        },
        "result": result,
    }
    encoded = json.dumps(payload, indent=2, sort_keys=True)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded + "\n", encoding="utf-8")
    print(encoded)


if __name__ == "__main__":
    main()
