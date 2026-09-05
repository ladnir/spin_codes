#!/usr/bin/env python3
"""Sparse-shell late-activation bounds for LDPCSplitState packet placement."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


DEFAULT_PACKET_AUDIT = Path(
    "constructions/riffle_ldpcsplitstate_g4_t256_s64/"
    "receipts/fixed_nested_pair_depth2_packet_kernel.json"
)
DEFAULT_OUTPUT = Path(
    "constructions/riffle_ldpcsplitstate_g4_t256_s64/"
    "receipts/sparse_prefix_w38.json"
)


def log_choose(n: int, k: int) -> float:
    if not 0 <= k <= n:
        return -math.inf
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def log_add(left: float, right: float) -> float:
    if left == -math.inf:
        return right
    if right == -math.inf:
        return left
    maximum = max(left, right)
    return maximum + math.log(math.exp(left - maximum) + math.exp(right - maximum))


def modeled_log2_multiplicity(weight: int, minimum_distance: int = 38) -> float:
    allowed = range(minimum_distance, 256 - minimum_distance + 1, 2)
    denominator = sum(math.comb(256, candidate) for candidate in allowed)
    return (
        math.log2((1 << 128) - 2)
        + math.log2(math.comb(256, weight))
        - math.log2(denominator)
    )


def four_block_prefix_log_probability(
    *, regions: int, prefix: int, weight: int, four_packet_failure: float
) -> tuple[float, list[dict[str, float | int]]]:
    suffix = regions - prefix
    denominator_log = 4 * log_choose(regions, weight)
    total_log = -math.inf
    rows = []
    for common_prefix_regions in range(
        max(0, weight - suffix), min(prefix, weight) + 1
    ):
        term = (
            log_choose(prefix, common_prefix_regions)
            + common_prefix_regions * math.log(four_packet_failure)
            + 4 * log_choose(suffix, weight - common_prefix_regions)
            - denominator_log
        )
        total_log = log_add(total_log, term)
        rows.append(
            {
                "common_prefix_regions": common_prefix_regions,
                "log2_contribution": term / math.log(2.0),
            }
        )
    return total_log, rows


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    packet_audit = json.loads(args.packet_audit.read_text(encoding="utf-8"))
    worst = packet_audit["four_singleton_packet_placement"][
        "worst_lane_composition"
    ]
    p4 = float(worst["failure_probability"])
    local_log2 = modeled_log2_multiplicity(args.outer_weight)
    rows = []
    for active_blocks in range(1, 5):
        multiplicity_log2 = (
            log_choose(args.outer_blocks, active_blocks) / math.log(2.0)
            + active_blocks * local_log2
        )
        if active_blocks <= 3:
            log_probability = active_blocks * (
                log_choose(args.regions - args.prefix_regions, args.outer_weight)
                - log_choose(args.regions, args.outer_weight)
            )
            terms = []
            method = (
                "Every occupied prefix region contains at most three active "
                "packets and therefore activates; the prefix must be empty."
            )
        else:
            log_probability, terms = four_block_prefix_log_probability(
                regions=args.regions,
                prefix=args.prefix_regions,
                weight=args.outer_weight,
                four_packet_failure=p4,
            )
            method = (
                "Each prefix region is empty or is a four-way intersection "
                "whose four singleton packets hit a fixed kernel pattern."
            )
        probability_log2 = log_probability / math.log(2.0)
        rows.append(
            {
                "active_outer_blocks": active_blocks,
                "modeled_message_multiplicity_log2": multiplicity_log2,
                "prefix_nonactivation_log2_probability": probability_log2,
                "prefix_nonactivation_bits": -probability_log2,
                "union_margin_bits_before_output_tail": (
                    -probability_log2 - multiplicity_log2
                ),
                "method": method,
                "four_block_terms": terms,
            }
        )
    return {
        "schema": "riffle-ldpcsplitstate-sparse-prefix-v1",
        "candidate": "Riffle LDPCSplitState g=4 t=256 s=64",
        "parameters": {
            "outer_blocks": args.outer_blocks,
            "regions": args.regions,
            "prefix_regions": args.prefix_regions,
            "suffix_regions": args.regions - args.prefix_regions,
            "outer_weight": args.outer_weight,
            "modeled_local_weight_multiplicity_log2": local_log2,
            "four_singleton_packet_failure_probability": p4,
            "four_singleton_packet_failure_bits": -math.log2(p4),
        },
        "rows": rows,
        "scope": (
            "Rigorous prefix-placement calculation conditioned on four "
            "modeled weight-38 outer words and the recorded fixed B. The "
            "outer multiplicity uses the real-valued BCH-like spectrum model. "
            "The four-block row uses the worst lane composition and therefore "
            "upper-bounds every grouping of four active blocks. It does not "
            "include the post-activation output-weight moment."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--packet-audit", type=Path, default=DEFAULT_PACKET_AUDIT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--outer-blocks", type=int, default=8192)
    parser.add_argument("--regions", type=int, default=256)
    parser.add_argument("--prefix-regions", type=int, default=210)
    parser.add_argument("--outer-weight", type=int, default=38)
    args = parser.parse_args()
    payload = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    for row in payload["rows"]:
        print(
            f"active_blocks,{row['active_outer_blocks']},"
            f"nonactivation_bits,{row['prefix_nonactivation_bits']:.6f},"
            f"multiplicity_bits,{row['modeled_message_multiplicity_log2']:.6f},"
            f"margin_before_tail,{row['union_margin_bits_before_output_tail']:.6f}"
        )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
