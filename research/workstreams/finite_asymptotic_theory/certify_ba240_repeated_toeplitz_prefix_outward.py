#!/usr/bin/env python3
"""Freeze and certify one repeated Golay--BA-3/Toeplitz construction.

The verifier rebuilds the complete routed outer generator from a stable
SplitMix64/Fisher--Yates specification.  It computes every prefix kernel
dimension over GF(2), groups the exact prefix first-moment identity by routed
regions, and evaluates a rigorous upper bound with python-flint Arb.

The only random object in the certified probability statement is the shared
lower-triangular Toeplitz kernel h[1],...,h[N-1].
"""

from __future__ import annotations

from array import array
import hashlib
import json
import math
import platform
import sys
from pathlib import Path

import numpy as np
from flint import arb, ctx


WORKSTREAM = Path(__file__).resolve().parent
OUTPUT = WORKSTREAM / "ba3_B240_repeated_toeplitz_prefix_outward.json"

MASK64 = (1 << 64) - 1
SETUP_SEED = 0x4241323430544F45  # ASCII bytes "BA240TOE", read big-endian.
SPLITMIX_GAMMA = 0x9E3779B97F4A7C15
SPLITMIX_MUL1 = 0xBF58476D1CE4E5B9
SPLITMIX_MUL2 = 0x94D049BB133111EB

B = 240
LOCAL_DIMENSION = 120
L = 8832
N = B * L
PARENT_DIMENSION = LOCAL_DIMENSION * L
TARGET_DIMENSION = 1 << 20
DISTANCE_CUTOFF = (11 * N) // 100
ARBITRARY_PRECISION_BITS = 256
CERTIFIED_MARGIN_BITS = 180

GOLAY_GENERATOR = sum(1 << i for i in (11, 9, 7, 6, 5, 1, 0))
GOLAY_SPECTRUM = {0: 1, 8: 759, 12: 2576, 16: 759, 24: 1}


class SplitMix64:
    """Stable unsigned-64-bit stream used only to name a fixed setup."""

    def __init__(self, seed: int):
        self.state = seed & MASK64
        self.words_consumed = 0
        self.rejected_words = 0

    def next_u64(self) -> int:
        self.state = (self.state + SPLITMIX_GAMMA) & MASK64
        z = self.state
        z = ((z ^ (z >> 30)) * SPLITMIX_MUL1) & MASK64
        z = ((z ^ (z >> 27)) * SPLITMIX_MUL2) & MASK64
        z ^= z >> 31
        self.words_consumed += 1
        return z

    def below(self, modulus: int) -> int:
        if not 1 <= modulus <= (1 << 64):
            raise ValueError("invalid rejection-sampling modulus")
        limit = (1 << 64) - ((1 << 64) % modulus)
        while True:
            value = self.next_u64()
            if value < limit:
                return value % modulus
            self.rejected_words += 1

    def permutation(self, size: int) -> list[int]:
        result = list(range(size))
        for index in range(size - 1, 0, -1):
            other = self.below(index + 1)
            result[index], result[other] = result[other], result[index]
        return result


def little_u16(values: list[int]) -> bytes:
    encoded = array("H", values)
    if encoded.itemsize != 2:
        raise RuntimeError("platform unsigned-short width is not 16 bits")
    if sys.byteorder != "little":
        encoded.byteswap()
    return encoded.tobytes()


def little_u32(value: int) -> bytes:
    return value.to_bytes(4, "little", signed=False)


def golay_rows() -> list[int]:
    rows = []
    for message_bit in range(12):
        word = GOLAY_GENERATOR << message_bit
        word |= (word.bit_count() & 1) << 23
        rows.append(word)

    spectrum: dict[int, int] = {}
    for message in range(1 << 12):
        word = 0
        for bit, row in enumerate(rows):
            if (message >> bit) & 1:
                word ^= row
        spectrum[word.bit_count()] = spectrum.get(word.bit_count(), 0) + 1
    if spectrum != GOLAY_SPECTRUM:
        raise ArithmeticError(f"unexpected extended Golay spectrum: {spectrum}")
    return rows


def golay_direct_sum_columns() -> list[int]:
    rows = golay_rows()
    columns = []
    for block in range(B // 24):
        shift = 12 * block
        for coordinate in range(24):
            column = 0
            for message_bit, row in enumerate(rows):
                column |= ((row >> coordinate) & 1) << (shift + message_bit)
            columns.append(column)
    return columns


def add_column(basis: list[int], column: int) -> bool:
    value = column
    while value:
        pivot = value.bit_length() - 1
        if basis[pivot]:
            value ^= basis[pivot]
        else:
            basis[pivot] = value
            return True
    return False


def rank(columns: list[int]) -> int:
    basis = [0] * LOCAL_DIMENSION
    return sum(add_column(basis, column) for column in columns)


def accumulate(columns: list[int], permutation: list[int]) -> list[int]:
    result = []
    state = 0
    for coordinate in permutation:
        state ^= columns[coordinate]
        result.append(state)
    return result


def binomial_arb(n: int, k: int) -> arb:
    """Outward enclosure of C(n,k), using Arb log-gamma arithmetic."""
    if k < 0 or k > n:
        return arb(0)
    k = min(k, n - k)
    if k == 0:
        return arb(1)
    return (
        arb(n + 1).lgamma()
        - arb(k + 1).lgamma()
        - arb(n - k + 1).lgamma()
    ).exp()


def definitely_less(left: arb, right: arb) -> bool:
    return bool(left.upper() < right.lower())


def rebuild_prefix_blocks() -> tuple[list[dict[str, int]], dict[str, object]]:
    rng = SplitMix64(SETUP_SEED)

    accumulator_hashes = []
    columns = golay_direct_sum_columns()
    if rank(columns) != LOCAL_DIMENSION:
        raise ArithmeticError("Golay direct sum has the wrong dimension")
    for _ in range(2):
        permutation = rng.permutation(B)
        accumulator_hashes.append(hashlib.sha256(little_u16(permutation)).hexdigest())
        columns = accumulate(columns, permutation)
    if rank(columns) != LOCAL_DIMENSION:
        raise ArithmeticError("BA-3 generator has the wrong dimension")

    local_hash = hashlib.sha256()
    increments = np.zeros((B, L), dtype=np.uint8)
    for row in range(L):
        permutation = rng.permutation(B)
        local_hash.update(little_u16(permutation))
        basis = [0] * LOCAL_DIMENSION
        row_rank = 0
        for position, coordinate in enumerate(permutation):
            if add_column(basis, columns[coordinate]):
                increments[position, row] = 1
                row_rank += 1
        if row_rank != LOCAL_DIMENSION:
            raise ArithmeticError(f"local row {row} ends at rank {row_rank}")

    route_hash = hashlib.sha256()
    prefix_hash = hashlib.sha256()
    prefix_hash.update(little_u32(PARENT_DIMENSION))
    current_dimension = PARENT_DIMENSION
    current_prefix = 0
    earliest_deficit = None
    maximum_deficit = 0
    positions_with_deficit = 0
    sum_deficit = 0
    blocks = []
    zero_prefix = None
    previous_kappa_plus_prefix = PARENT_DIMENSION

    for region in range(B):
        block_start = current_prefix + 1
        permutation = rng.permutation(L)
        route_hash.update(little_u16(permutation))
        for row in permutation:
            current_dimension -= int(increments[region, row])
            current_prefix += 1
            prefix_hash.update(little_u32(current_dimension))
            kappa_plus_prefix = current_dimension + current_prefix
            if kappa_plus_prefix < previous_kappa_plus_prefix:
                raise ArithmeticError("kappa_t+t must be monotone")
            previous_kappa_plus_prefix = kappa_plus_prefix
            ideal_dimension = max(PARENT_DIMENSION - current_prefix, 0)
            deficit = current_dimension - ideal_dimension
            if deficit:
                if earliest_deficit is None:
                    earliest_deficit = current_prefix
                positions_with_deficit += 1
                sum_deficit += deficit
                maximum_deficit = max(maximum_deficit, deficit)
            if current_dimension == 0 and zero_prefix is None:
                zero_prefix = current_prefix
        block_end = current_prefix
        blocks.append(
            {
                "region": region,
                "prefix_start": block_start,
                "prefix_end": block_end,
                "ending_kernel_dimension": current_dimension,
                # kappa_t+t is monotone, so this is the block maximum.
                "max_kappa_plus_prefix": current_dimension + current_prefix,
            }
        )

    if current_prefix != N or current_dimension != 0 or zero_prefix is None:
        raise ArithmeticError("routed generator did not reach a zero kernel")

    metadata = {
        "accumulator_permutation_sha256": accumulator_hashes,
        "local_coordinate_permutations_sha256": local_hash.hexdigest(),
        "region_row_permutations_sha256": route_hash.hexdigest(),
        "prefix_dimensions_u32le_sha256": prefix_hash.hexdigest(),
        "splitmix64_words_consumed": rng.words_consumed,
        "splitmix64_rejected_words": rng.rejected_words,
        "earliest_prefix_rank_deficit": earliest_deficit,
        "maximum_prefix_rank_deficit": maximum_deficit,
        "prefix_positions_with_rank_deficit": positions_with_deficit,
        "sum_prefix_rank_deficit": sum_deficit,
        "first_zero_kernel_prefix": zero_prefix,
    }
    return blocks, metadata


def outward_first_moment(blocks: list[dict[str, int]]) -> tuple[arb, list[dict[str, object]], arb]:
    two = arb(2)

    # For X~Bin(N-1,1/2), the lower tail through D-1 is bounded by
    # C(N-1,D-1) / (1-(D-1)/(N-D+1)).  Also Z_0 < 2^K.
    last_mass = binomial_arb(N - 1, DISTANCE_CUTOFF - 1)
    tail_ratio_bound = arb(N - DISTANCE_CUTOFF + 1) / arb(
        N - 2 * DISTANCE_CUTOFF + 2
    )
    base_bound = (
        last_mass
        * tail_ratio_bound
        * (two ** (PARENT_DIMENSION - (N - 1)))
    )

    rows = []
    aggregate = base_bound
    final_prefix = N - DISTANCE_CUTOFF
    for source in blocks:
        start = source["prefix_start"]
        end = min(source["prefix_end"], final_prefix)
        if start > end:
            break
        if source["ending_kernel_dimension"] == 0:
            # The zero may occur inside this region. Conservatively retain the
            # whole clipped region; later zero-kernel regions contribute zero.
            pass

        maximum = source["max_kappa_plus_prefix"]
        upper_binomial = binomial_arb(N - start, DISTANCE_CUTOFF)
        lower_binomial = binomial_arb(N - end - 1, DISTANCE_CUTOFF)
        binomial_sum = upper_binomial - lower_binomial
        if not definitely_less(arb(0), binomial_sum):
            raise ArithmeticError(f"binomial block difference is not positive at {start}")
        contribution = binomial_sum * (two ** (maximum - N))
        aggregate += contribution
        rows.append(
            {
                "region": source["region"],
                "prefix_interval": [start, end],
                "max_kappa_plus_prefix": maximum,
                "ending_kernel_dimension": source["ending_kernel_dimension"],
                "log2_contribution_upper": str((contribution.log() / two.log()).upper()),
            }
        )
        if source["ending_kernel_dimension"] == 0:
            break

    return aggregate, rows, base_bound


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    if B % 24 or LOCAL_DIMENSION != B // 2 or N != B * L:
        raise ArithmeticError("inconsistent frozen parameters")
    if DISTANCE_CUTOFF != 233_164:
        raise ArithmeticError("unexpected 11% cutoff")

    ctx.prec = ARBITRARY_PRECISION_BITS
    print("rebuilding fixed BA-3 generator and all routed prefix ranks", flush=True)
    blocks, setup = rebuild_prefix_blocks()
    print("evaluating outward grouped prefix bound", flush=True)
    aggregate, rows, base_bound = outward_first_moment(blocks)
    log2_aggregate = aggregate.log() / arb(2).log()
    threshold = arb(2) ** (-CERTIFIED_MARGIN_BITS)
    if not definitely_less(aggregate, threshold):
        raise ArithmeticError(
            f"grouped bound does not pass {CERTIFIED_MARGIN_BITS} bits: {log2_aggregate}"
        )

    payload = {
        "schema": "ba3-b240-repeated-toeplitz-prefix-outward-v1",
        "status": "COMPLETE_OUTWARD_ALL_MESSAGE_CERTIFICATE",
        "claim": {
            "bad_nonzero_parent_messages_expected_upper": str(aggregate.upper()),
            "bad_nonzero_parent_messages_log2_upper": str(log2_aggregate.upper()),
            "margin_bits_lower_display": -float(log2_aggregate.upper()),
            "certified_integer_margin_bits": CERTIFIED_MARGIN_BITS,
            "comparison_to_2^-180": True,
            "bad_output_weight_at_most": DISTANCE_CUTOFF,
            "minimum_distance_on_success_at_least": DISTANCE_CUTOFF + 1,
            "relative_minimum_distance_on_success_lower": (DISTANCE_CUTOFF + 1) / N,
            "all_nonzero_parent_messages_covered": True,
        },
        "parameters": {
            "outer_constituent": "one [240,120] Golay-direct-sum plus two-accumulator code, repeated in every row",
            "outer_bits": B,
            "outer_dimension": LOCAL_DIMENSION,
            "outer_rows": L,
            "output_bits": N,
            "parent_dimension": PARENT_DIMENSION,
            "target_dimension_after_zero_shortening": TARGET_DIMENSION,
            "zero_shortened_parent_inputs": PARENT_DIMENSION - TARGET_DIMENSION,
            "zero_shortening_rule": "retain the first 2^20 row-major parent input coordinates and set the remaining 11264 coordinates to zero",
            "distance_cutoff": DISTANCE_CUTOFF,
            "setup_seed_hex": f"0x{SETUP_SEED:016x}",
            "inner": "one shared invertible random lower-triangular Toeplitz convolution over GF(2), h[0]=1",
        },
        "fixed_setup": setup,
        "setup_algorithm": {
            "word_generator": "SplitMix64 with unsigned 64-bit wraparound",
            "gamma_hex": f"0x{SPLITMIX_GAMMA:016x}",
            "multiplier_1_hex": f"0x{SPLITMIX_MUL1:016x}",
            "multiplier_2_hex": f"0x{SPLITMIX_MUL2:016x}",
            "bounded_integer": "reject x >= 2^64-(2^64 mod m), then return x mod m",
            "permutation": "descending Fisher-Yates: for i=size-1,...,1 swap positions i and below(i+1)",
            "stream_order": [
                "two length-240 accumulator permutations",
                "8832 length-240 local-coordinate permutations, in row order",
                "240 length-8832 row permutations, in region order",
            ],
            "hash_encoding": "permutation entries are concatenated unsigned 16-bit little-endian integers; prefix dimensions include kappa_0 and use unsigned 32-bit little-endian integers",
        },
        "outer_validation": {
            "golay_generator_polynomial_exponents": [11, 9, 7, 6, 5, 1, 0],
            "extended_golay_weight_spectrum": GOLAY_SPECTRUM,
            "base_and_final_generator_rank": LOCAL_DIMENSION,
        },
        "bound": {
            "identity": "Z_0 F_N + sum_{t=1}^{N-D} Z_t*C(N-t-1,D-1)/2^(N-t)",
            "prefix_count": "Z_t=2^kappa_t-1",
            "base_tail_bound": "F_N <= C(N-1,D-1)*(N-D+1)/(N-2D+2)/2^(N-1)",
            "group_rule": "on [a,b], C=max(kappa_t+t), so the sum is at most 2^(C-N)*(C(N-a,D)-C(N-b-1,D))",
            "base_term_log2_upper": str((base_bound.log() / arb(2).log()).upper()),
            "regions_in_outward_sum": len(rows),
            "region_bounds": rows,
        },
        "probability_space": {
            "fixed": "the outer generator, its two accumulator permutations, all local coordinate permutations, and all region row permutations",
            "random": "h[1],...,h[N-1] are independent uniform GF(2) bits and h[0]=1",
            "failure_event": "some nonzero parent message produces output Hamming weight at most D",
            "justification": "Markov's inequality bounds this failure probability by the expected number of bad nonzero parent messages",
        },
        "arithmetic": {
            "library": "python-flint Arb",
            "precision_bits": ARBITRARY_PRECISION_BITS,
            "integer_prefix_ranks": "exact GF(2) elimination with Python integers",
            "binomial_evaluation": "outward Arb log-gamma enclosure",
        },
        "versions": {
            "python": platform.python_version(),
            "numpy": np.__version__,
        },
        "dependencies": [
            {"path": Path(__file__).name, "sha256": sha256(Path(__file__).resolve())}
        ],
        "limitations": [
            "The certificate proves existence and high conditional success probability over the Toeplitz kernel for this one fixed outer/routing setup.",
            "It does not prove that a fresh BA/routing sample passes with a stated probability.",
            "It does not supply a linear-time implementation of full-length Toeplitz convolution.",
            "The parent-code bound is conservative for any fixed zero-shortened 2^20-dimensional subcode.",
        ],
    }
    with OUTPUT.open("w", encoding="utf-8", newline="\n") as output_file:
        output_file.write(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({"output": str(OUTPUT), **payload["claim"]}, indent=2))


if __name__ == "__main__":
    main()
