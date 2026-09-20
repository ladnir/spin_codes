#!/usr/bin/env python3
"""Sum full-split episode terms using a product-GF outer upper bound.

This is the post-prefix companion to sum_fullsplit_episode_envelope.py.  Instead
of reading exact/prefix outer coefficients from a CSV, it bounds the direct-sum
RM outer coefficients by

    A_h <= min_z (W_loc(z)^M - 1) z^{-h}.

The inner values are read from full-split episode CSVs as a function of
H=h-r, one CSV per gap bucket.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from block_outer_upgrade_probe import load_local_spectrum, outer_block_gf_bounds, z_grid
from dense_largek_eval import log2_binom, log2add
from probe_block_recursive_inner import parse_int_list
from sum_fullsplit_episode_envelope import (
    envelope_value,
    load_knot_csv,
    parse_float_list,
    precompute_gap_sums,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--local-spectrum-csv", type=Path, default=Path(__file__).with_name("rm512_256_spectrum.csv"))
    parser.add_argument("--blocks", type=int, default=4096)
    parser.add_argument("--outer-block-bits", type=int, default=256)
    parser.add_argument("--local-length", type=int, default=512)
    parser.add_argument("--local-distance", type=int, default=32)
    parser.add_argument("--N", type=int, default=2**21)
    parser.add_argument("--inner-block-bits", type=int, default=64)
    parser.add_argument("--h-values", default="501:2000")
    parser.add_argument("--first-r-values", default="1:64")
    parser.add_argument("--late-blocks", type=int, default=5949)
    parser.add_argument("--gap-start", type=int, default=1)
    parser.add_argument("--gap-stop", type=int, default=12000)
    parser.add_argument("--gap-step", type=int, default=4000)
    parser.add_argument("--inner-knot-csvs", required=True)
    parser.add_argument("--knot-h-column", default="H")
    parser.add_argument("--knot-value-column", default="total_log2")
    parser.add_argument("--inner-envelope-lift-bits-by-gap")
    parser.add_argument("--require-knot-coverage", action="store_true")
    parser.add_argument("--z-min", type=float, default=1e-8)
    parser.add_argument("--z-max", type=float, default=0.999)
    parser.add_argument("--z-count", type=int, default=220)
    parser.add_argument("--output-csv", type=Path)
    args = parser.parse_args()

    h_values = parse_int_list(args.h_values)
    r_values = parse_int_list(args.first_r_values)
    h_max = max(h_values)
    spectrum = load_local_spectrum(str(args.local_spectrum_csv))
    zs = z_grid(args.z_min, args.z_max, args.z_count)
    outer_logs, outer_z = outer_block_gf_bounds(
        blocks=args.blocks,
        block_bits=args.outer_block_bits,
        local_length=args.local_length,
        d0=args.local_distance,
        h_max=h_max,
        zs=zs,
        model="spectrum-csv",
        spectrum=spectrum,
    )

    gap_ranges: list[tuple[int, int]] = []
    for gap_min in range(args.gap_start, args.gap_stop + 1, args.gap_step):
        gap_max = min(args.gap_stop, gap_min + args.gap_step - 1)
        gap_ranges.append((gap_min, gap_max))

    knot_lists = [
        load_knot_csv(
            Path(part.strip()),
            h_column=args.knot_h_column,
            value_column=args.knot_value_column,
        )
        for part in args.inner_knot_csvs.split(";")
        if part.strip()
    ]
    if len(knot_lists) != len(gap_ranges):
        raise ValueError("--inner-knot-csvs must have one path per gap bucket")
    if args.inner_envelope_lift_bits_by_gap:
        lift_bits = parse_float_list(args.inner_envelope_lift_bits_by_gap)
        if len(lift_bits) != len(gap_ranges):
            raise ValueError("--inner-envelope-lift-bits-by-gap must have one value per gap bucket")
    else:
        lift_bits = [0.0 for _ in gap_ranges]

    max_h_after_first = max((h - 1 for h in h_values), default=0)
    gap_sums = precompute_gap_sums(
        b=args.inner_block_bits,
        late_blocks=args.late_blocks,
        gap_ranges=gap_ranges,
        max_h_after_first=max_h_after_first,
    )

    rows = []
    total = float("-inf")
    by_gap = [float("-inf") for _ in gap_ranges]
    peak = None
    by_h: dict[int, float] = {}
    for h in h_values:
        outer = outer_logs[h] if h < len(outer_logs) else float("-inf")
        if outer == float("-inf"):
            continue
        for r in r_values:
            if r < 1 or r > min(h, args.inner_block_bits):
                continue
            H = h - r
            for idx, (gap_min, gap_max) in enumerate(gap_ranges):
                gap_sum = gap_sums.get((H, idx), float("-inf"))
                if gap_sum == float("-inf"):
                    continue
                placement = (
                    log2_binom(args.inner_block_bits, r)
                    + gap_sum
                    - log2_binom(args.N, h)
                )
                inner = min(
                    0.0,
                    envelope_value(
                        knot_lists[idx],
                        H,
                        require_coverage=args.require_knot_coverage,
                    )
                    + lift_bits[idx],
                )
                term = outer + placement + inner
                total = log2add(total, term)
                by_gap[idx] = log2add(by_gap[idx], term)
                by_h[h] = log2add(by_h.get(h, float("-inf")), term)
                row = {
                    "outer_weight": h,
                    "first_r": r,
                    "remaining_ones": H,
                    "gap_min": gap_min,
                    "gap_max": gap_max,
                    "outer_log2_bound": outer,
                    "outer_z": outer_z[h],
                    "placement_log2": placement,
                    "inner_log2": inner,
                    "term_log2": term,
                }
                rows.append(row)
                if peak is None or term > peak["term_log2"]:
                    peak = row

    if args.output_csv:
        fields = [
            "outer_weight",
            "first_r",
            "remaining_ones",
            "gap_min",
            "gap_max",
            "outer_log2_bound",
            "outer_z",
            "placement_log2",
            "inner_log2",
            "term_log2",
        ]
        with args.output_csv.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)

    print("Full-split episode product-GF outer sum")
    print(f"h_values,{args.h_values}")
    print(f"z_count,{args.z_count}")
    print(f"total_log2,{total:.6f}")
    print(f"margin_bits,{-total:.6f}")
    for (gap_min, gap_max), val in zip(gap_ranges, by_gap):
        print(f"gap_{gap_min}_{gap_max}_log2,{val:.6f}")
    if peak is not None:
        print(f"peak_h,{peak['outer_weight']}")
        print(f"peak_first_r,{peak['first_r']}")
        print(f"peak_gap_min,{peak['gap_min']}")
        print(f"peak_gap_max,{peak['gap_max']}")
        print(f"peak_outer_log2_bound,{peak['outer_log2_bound']:.6f}")
        print(f"peak_inner_log2,{peak['inner_log2']:.6f}")
        print(f"peak_term_log2,{peak['term_log2']:.6f}")
    if by_h:
        h_peak, v_peak = max(by_h.items(), key=lambda kv: kv[1])
        print(f"peak_aggregate_h,{h_peak}")
        print(f"peak_aggregate_h_log2,{v_peak:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
