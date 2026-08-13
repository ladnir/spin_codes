#!/usr/bin/env python3
"""Floating diagnostic for the systematic striped group-chain ledger.

This script starts with the exactly enumerable one-data-block, graph-zero
family.  A weight-w local EBCH word occupies w distinct physical groups in the
striped layout, and those groups have a uniform order in the B-position chain.

For ordered active positions a_1<...<a_g, write x_0 for the gap before a_1
and x_i for the gap after a_i.  The vector (x_0,...,x_g) is uniform over weak
compositions of B-g.  If cancellation occurs at e selected active positions
other than a_1, the off time is x_0 plus the e following gaps.  Summing over
that composition gives the episode union bound below.  Each cancellation is
charged 2^-38; the exact local certificate proves 2^-39 and the reached-vector
restart factor is below two.

The arithmetic here is log-domain floating point and is a diagnostic, not the
final outward/exact theorem certificate.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from punctured_ebch_outer import (
    full_spectrum,
    heterogeneous_coefficients,
    punctured_spectrum,
)
from certify_rm_outer_prefix_exact import truncated_convolution


B = 32768
N = 1 << 21
DISTANCE = 9 * N // 100
DATA_BLOCKS = 16384
PUNCTURED_SLOTS = 128
GRAPH_CODIMENSION = 24
TURNOFF_SHIFT = 39
LOG2_POLE = math.log2(2333 / 2373)
LOG2_RHO = -0.624039853324
LOG2_PREFACTOR = -0.024525793210
GRAPH_SPECTRUM = Path(__file__).resolve().parent / "ebch128_graph24_spectrum.csv"
PA1_CERTIFICATES = (
    Path(__file__).resolve().parent / "systematic_pa1_geometric_certificates.json"
)


def log2_sum(values: np.ndarray | list[float]) -> float:
    array = np.asarray(values, dtype=float)
    if not array.size:
        return -math.inf
    maximum = float(np.max(array))
    if not math.isfinite(maximum):
        return -math.inf
    return maximum + math.log2(float(np.sum(np.exp2(array - maximum))))


LOG2_FACTORIAL = np.array(
    [math.lgamma(value + 1) / math.log(2) for value in range(B + 2)],
    dtype=float,
)


def log2_binomial(n: int | np.ndarray, k: int | np.ndarray) -> np.ndarray:
    n_array = np.asarray(n, dtype=np.int64)
    k_array = np.asarray(k, dtype=np.int64)
    return (
        LOG2_FACTORIAL[n_array]
        - LOG2_FACTORIAL[k_array]
        - LOG2_FACTORIAL[n_array - k_array]
    )


LIVE = np.arange(B + 1, dtype=np.int64)


def group_episode_one_pole_log2(
    g: int, *, pole: float, rho: float, prefactor: float
) -> float:
    """Union bound for a fixed word occupying g distinct nonzero groups."""

    if not 1 <= g <= B:
        raise ValueError("group support must lie in 1..B")
    slack = B - g
    off = np.arange(slack + 1, dtype=np.int64)
    denominator = float(log2_binomial(B, g))
    log2_rho = math.log2(rho)
    log2_prefactor = math.log2(prefactor)
    survival = np.minimum(
        0.0,
        -DISTANCE * math.log2(pole)
        + log2_prefactor
        + np.maximum(LIVE - 1, 0) * log2_rho,
    )
    survival[: DISTANCE // 64 + 1] = 0.0
    # Every extra episode contributes another prefactor/rho, as well as one
    # exact-support turnoff.  This retains the actual pole-dependent restart
    # factor instead of assuming it is below two for every pole.
    restart_log2 = -TURNOFF_SHIFT + log2_prefactor - log2_rho
    episode_terms: list[float] = []
    for cancellations in range(g):
        # Select cancellation positions from a_2,...,a_g.  The sum of x_0
        # and their following gaps is `off`; the other g-e gaps sum to the
        # remaining slack.  Stars-and-bars counts both pieces exactly.
        selected = (
            log2_binomial(off + cancellations, cancellations)
            + log2_binomial(
                slack - off + g - cancellations - 1,
                g - cancellations - 1,
            )
        )
        live = B - off
        episode_terms.append(
            float(log2_binomial(g - 1, cancellations))
            + restart_log2 * cancellations
            - denominator
            + log2_sum(selected + survival[live])
        )
    return min(0.0, log2_sum(episode_terms))


def episode_parameters(accumulate: bool) -> list[tuple[float, float, float]]:
    if not accumulate:
        return [(2333 / 2373, 2**LOG2_RHO, 2**LOG2_PREFACTOR)]
    payload = json.loads(PA1_CERTIFICATES.read_text(encoding="utf-8"))
    if payload.get("status") != "EXACT_RATIONAL_COLLATZ_CERTIFICATE":
        raise SystemExit("striped ledger: invalid PA1 certificate status")
    result = []
    for row in payload["certificates"]:
        pole = row["pole_numerator"] / row["pole_denominator"]
        rho = int(row["rho_numerator"]) / int(row["rho_denominator"])
        prefactor = int(row["prefactor_numerator"]) / int(
            row["prefactor_denominator"]
        )
        result.append((pole, rho, prefactor))
    return result


def group_episode_log2(
    g: int, parameters: list[tuple[float, float, float]]
) -> float:
    return min(
        group_episode_one_pole_log2(
            g, pole=pole, rho=rho, prefactor=prefactor
        )
        for pole, rho, prefactor in parameters
    )


def one_data_block_graph_zero(
    parameters: list[tuple[float, float, float]],
) -> tuple[float, list[tuple[int, float]]]:
    full = full_spectrum()
    punctured = punctured_spectrum()
    supports = sorted(
        {
            weight
            for weight, count in enumerate(full)
            if weight and count
        }
        | {
            weight
            for weight, count in enumerate(punctured)
            if weight and count
        }
    )
    inner = {weight: group_episode_log2(weight, parameters) for weight in supports}
    full_log2 = log2_sum(
        [math.log2(count) + inner[weight] for weight, count in enumerate(full) if weight and count]
    )
    punctured_log2 = log2_sum(
        [
            math.log2(count) + inner[weight]
            for weight, count in enumerate(punctured)
            if weight and count
        ]
    )
    averaged_slot = log2_sum(
        [
            math.log2(DATA_BLOCKS - PUNCTURED_SLOTS) + full_log2,
            math.log2(PUNCTURED_SLOTS) + punctured_log2,
        ]
    ) - math.log2(DATA_BLOCKS)
    family = math.log2(DATA_BLOCKS) - GRAPH_CODIMENSION + averaged_slot
    return family, [(weight, inner[weight]) for weight in supports]


def collision_free_aggregate(
    inner_rows: list[tuple[int, float]], h_max: int
) -> tuple[float, list[tuple[int, float]]]:
    """Expected low-output count if distinct outer coordinates never collide."""

    inner = dict(inner_rows)
    if any(weight not in inner for weight in range(21, h_max + 1)):
        raise ValueError("collision-free aggregate requires every inner group weight")
    data = heterogeneous_coefficients(
        full_blocks=DATA_BLOCKS - PUNCTURED_SLOTS,
        punctured_blocks=PUNCTURED_SLOTS,
        h_max=h_max,
    )
    # The user data message must be nonzero.  For every such fixed message,
    # averaging the random graph map is exactly a uniform average over the
    # 2^24 graph words.
    data[0] = 0
    graph = [0] * (h_max + 1)
    with GRAPH_SPECTRUM.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            weight = int(row["weight"])
            if weight <= h_max:
                graph[weight] = int(row["count"])
    outer = truncated_convolution(data, graph, h_max)
    terms = [
        (weight, math.log2(outer[weight]) - GRAPH_CODIMENSION + inner[weight])
        for weight in range(21, h_max + 1)
        if outer[weight]
    ]
    return log2_sum([value for _weight, value in terms]), terms


def two_data_block_graph_zero(
    inner_rows: list[tuple[int, float]], g_max: int
) -> tuple[float, float]:
    """Exact two-data-block contribution through group support ``g_max``.

    The tile probabilities use the actual layout: 128 marked tiles have 63
    full lanes and one punctured lane; the other 128 tiles have 64 full lanes.
    Conditional on co-tiling, the intersection with a full block has the exact
    128-coordinate hypergeometric law.  Punctured/punctured pairs cannot share
    a tile.
    """

    inner = dict(inner_rows)
    full = full_spectrum()
    punctured = punctured_spectrum()
    full_nonzero = [(w, c) for w, c in enumerate(full) if w and c]
    punctured_nonzero = [(w, c) for w, c in enumerate(punctured) if w and c]
    full_slots = DATA_BLOCKS - PUNCTURED_SLOTS
    punctured_slots = PUNCTURED_SLOTS
    same_ff = (
        128 * math.comb(63, 2) + 128 * math.comb(64, 2)
    ) / math.comb(full_slots, 2)
    same_fp = 63 / full_slots

    def pair_log(
        left: list[tuple[int, int]],
        right: list[tuple[int, int]],
        pair_count: int,
        same_tile: float,
    ) -> tuple[list[float], list[float]]:
        distinct_terms: list[float] = []
        collision_terms: list[float] = []
        for w1, c1 in left:
            for w2, c2 in right:
                base = math.log2(pair_count) + math.log2(c1) + math.log2(c2)
                if w1 + w2 <= g_max:
                    distinct_terms.append(
                        base + math.log2(1 - same_tile) + inner[w1 + w2]
                    )
                if not same_tile:
                    continue
                denominator = math.comb(128, w2)
                for intersection in range(
                    max(0, w1 + w2 - 128), min(w1, w2) + 1
                ):
                    union = w1 + w2 - intersection
                    if union > g_max:
                        continue
                    ways = math.comb(w1, intersection) * math.comb(
                        128 - w1, w2 - intersection
                    )
                    if ways:
                        collision_terms.append(
                            base
                            + math.log2(same_tile)
                            + math.log2(ways / denominator)
                            + inner[union]
                        )
        return distinct_terms, collision_terms

    pieces = [
        pair_log(
            full_nonzero,
            full_nonzero,
            math.comb(full_slots, 2),
            same_ff,
        ),
        pair_log(
            full_nonzero,
            punctured_nonzero,
            full_slots * punctured_slots,
            same_fp,
        ),
        pair_log(
            punctured_nonzero,
            punctured_nonzero,
            math.comb(punctured_slots, 2),
            0.0,
        ),
    ]
    distinct = log2_sum([value for part, _ in pieces for value in part]) - GRAPH_CODIMENSION
    collided = log2_sum([value for _, part in pieces for value in part]) - GRAPH_CODIMENSION
    return distinct, collided


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--accumulate",
        action="store_true",
        help="use one random-permute-plus-accumulate round before E_sys",
    )
    parser.add_argument("--aggregate-h-max", type=int, default=106)
    args = parser.parse_args()
    parameters = episode_parameters(args.accumulate)
    family, rows = one_data_block_graph_zero(parameters)
    aggregate, aggregate_rows = collision_free_aggregate(
        rows, args.aggregate_h_max
    )
    pair_distinct, pair_collided = two_data_block_graph_zero(
        rows, args.aggregate_h_max
    )
    print("systematic striped group ledger: one-data-block graph-zero family")
    print(f"B={B} N={N} distance={DISTANCE}")
    print(f"mixer={'PA1' if args.accumulate else 'identity'}")
    print(f"turnoff_charge=2^-{TURNOFF_SHIFT}")
    for weight, inner in rows:
        print(f"group_support={weight:3d} inner_log2_upper={inner:.12f}")
    print(f"family_log2_upper_approx={family:.12f}")
    print(f"margin_beyond_40_bits={-40.0-family:.12f}")
    print(
        f"collision_free_h_21_{args.aggregate_h_max}_log2_upper_approx="
        f"{aggregate:.12f}"
    )
    dominant = sorted(aggregate_rows, key=lambda row: row[1], reverse=True)[:10]
    print(
        "collision_free_dominant_terms="
        + ",".join(f"h{weight}:{value:.6f}" for weight, value in dominant)
    )
    print(
        f"collision_free_margin_beyond_40_bits={-40.0-aggregate:.12f}"
    )
    print(
        f"two_block_graph_zero_distinct_tiles_g_le_{args.aggregate_h_max}_log2="
        f"{pair_distinct:.12f}"
    )
    print(
        f"two_block_graph_zero_same_tile_g_le_{args.aggregate_h_max}_log2="
        f"{pair_collided:.12f}"
    )
    print("status=DIAGNOSTIC_ONLY_NOT_OUTWARD_CERTIFICATE")


if __name__ == "__main__":
    main()
