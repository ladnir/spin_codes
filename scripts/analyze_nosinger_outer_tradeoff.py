#!/usr/bin/env python3
"""Exact outer-placement tradeoff for removing the Singer inner scrambler.

The identity-inner geometric certificate moves the live-block Chernoff
crossing to 7527 at the proposed N.  This script computes the exact ambient
late-placement union bound and the effect of an independent random graph
subcode of codimension r:

    m -> (m, R m),  R <- F_2^{r x K} uniform.

Every fixed nonzero ambient message lies in the graph with probability at
most 2^-r.  The construction remains binary-linear and scalar-extends to
GF(2^128).  A four-bit transposed subset-sum circuit applies R^T using
ceil(r/4) table lookups per output symbol.
"""

from __future__ import annotations

import argparse
import math
from fractions import Fraction
from pathlib import Path

from certify_fullsplit_h500_rational import load_outer_coefficients
from outward_log2 import log2_fraction


ROOT = Path(__file__).resolve().parent
LOCAL_SPECTRUM = ROOT / "ebch128_64_spectrum.csv"


def subset_table_xors(width: int, group_bits: int) -> int:
    groups, remainder = divmod(width, group_bits)
    result = groups * ((1 << group_bits) - group_bits - 1)
    if remainder:
        result += (1 << remainder) - remainder - 1
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-size", type=int, default=2**20)
    parser.add_argument("--codim", type=int, default=24)
    parser.add_argument("--late-blocks", type=int, default=7527)
    parser.add_argument("--group-bits", type=int, default=8)
    parser.add_argument("--h-max", type=int, default=500)
    parser.add_argument("--local-spectrum", type=Path, default=LOCAL_SPECTRUM)
    args = parser.parse_args()
    if not 1 <= args.codim <= 64:
        raise SystemExit("no-Singer outer tradeoff expects 1 <= codim <= 64")
    if not 1 <= args.group_bits <= 12:
        raise SystemExit("no-Singer outer tradeoff expects 1 <= group bits <= 12")

    data_blocks = (args.message_size + args.codim + 63) // 64
    outer_blocks = data_blocks + 1  # explicit componentwise parity block
    data_capacity = 64 * data_blocks
    shortened_coordinates = data_capacity - (args.message_size + args.codim)
    if not 0 <= shortened_coordinates < 64:
        raise SystemExit("no-Singer outer tradeoff: invalid shortened dimension")
    n = 128 * outer_blocks
    outer = load_outer_coefficients(
        args.local_spectrum, blocks=outer_blocks, h_max=args.h_max
    )
    for weight in range(1, 44):
        outer[weight] = 0
    ultra_late = sum(
        (
            Fraction(
                outer[weight] * math.comb(64 * args.late_blocks, weight),
                math.comb(n, weight),
            )
            for weight in range(44, args.h_max + 1)
            if outer[weight]
        ),
        Fraction(0),
    )
    expurgated = ultra_late / (1 << args.codim)

    groups = (args.codim + args.group_bits - 1) // args.group_bits
    table_xors = subset_table_xors(args.codim, args.group_bits)
    # Start from the already-computed ambient output and XOR one selected
    # table entry per group into it.
    output_xors = args.message_size * groups
    total_xors = table_xors + output_xors
    mask_bytes = args.message_size * ((args.codim + 7) // 8)

    print("no-Singer outer graph-subcode tradeoff")
    print(
        f"K={args.message_size} r={args.codim} data_blocks={data_blocks} "
        f"outer_blocks={outer_blocks} N={n} late_blocks={args.late_blocks}"
    )
    print(
        f"data_capacity={data_capacity} shortened_coordinates={shortened_coordinates} "
        f"shortened_parity_outer_dimension={args.message_size + args.codim} "
        f"graph_dimension={args.message_size}"
    )
    print(f"ambient_ultra_late_log2_upper={log2_fraction(ultra_late).hi}")
    print(f"graph_ultra_late_log2_upper={log2_fraction(expurgated).hi}")
    print("fixed_nonzero_inclusion_probability_upper=2^-r")
    print(
        f"transpose_group_bits={args.group_bits} transpose_groups={groups} "
        f"table_xors={table_xors} "
        f"output_xors={output_xors} total_xors={total_xors} mask_bytes={mask_bytes}"
    )


if __name__ == "__main__":
    main()
