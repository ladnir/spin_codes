#!/usr/bin/env python3
"""Sum full-split episode terms using a piecewise-linear inner envelope.

This is the fast companion to sweep_fullsplit_episode_terms.py.  The expensive
episode ledger is evaluated at selected H knots; this script interpolates those
log2 bounds as a conservative working envelope and sums outer/placement terms
over a dense h/r prefix.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

from certify_fullsplit_inner_bucket_sum import load_outer
from dense_largek_eval import log2_binom, log2add
from probe_block_recursive_inner import parse_int_list


def parse_knots(text: str) -> list[tuple[int, float]]:
    out: list[tuple[int, float]] = []
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        h_text, v_text = part.split(":", 1)
        out.append((int(h_text), float(v_text)))
    out.sort()
    return out


def parse_float_list(text: str) -> list[float]:
    return [float(part.strip()) for part in text.split(",") if part.strip()]


def load_knot_csv(path: Path, *, h_column: str, value_column: str) -> list[tuple[int, float]]:
    out: list[tuple[int, float]] = []
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            h_text = row.get(h_column)
            value_text = row.get(value_column)
            if h_text is None or value_text is None:
                raise ValueError(
                    f"{path} must contain columns {h_column!r} and {value_column!r}"
                )
            if not h_text or not value_text:
                continue
            out.append((int(h_text), float(value_text)))
    out.sort()
    return out


def envelope_value(
    knots: list[tuple[int, float]],
    x: int,
    *,
    require_coverage: bool,
) -> float:
    if not knots:
        return 0.0
    if x <= knots[0][0]:
        if require_coverage and x < knots[0][0]:
            raise ValueError(f"envelope query H={x} is below first knot H={knots[0][0]}")
        return knots[0][1]
    for (x0, y0), (x1, y1) in zip(knots, knots[1:]):
        if x <= x1:
            t = (x - x0) / (x1 - x0)
            return y0 + t * (y1 - y0)
    if require_coverage:
        raise ValueError(f"envelope query H={x} is above last knot H={knots[-1][0]}")
    x0, y0 = knots[-2]
    x1, y1 = knots[-1]
    slope = (y1 - y0) / (x1 - x0)
    return y1 + slope * (x - x1)


def precompute_gap_sums(
    *,
    b: int,
    late_blocks: int,
    gap_ranges: list[tuple[int, int]],
    max_h_after_first: int,
    mode: str = "exact",
) -> dict[tuple[int, int], float]:
    """log2 sum_{gap in bucket} C(b*(late_blocks+gap-1), H)."""

    out: dict[tuple[int, int], float] = {}
    if mode not in {"exact", "endpoint"}:
        raise ValueError(f"unknown gap-sum mode {mode!r}")
    for H in range(max_h_after_first + 1):
        for idx, (gap_min, gap_max) in enumerate(gap_ranges):
            if mode == "endpoint":
                after_coords = b * (late_blocks + gap_max - 1)
                out[(H, idx)] = (
                    math.log2(gap_max - gap_min + 1) + log2_binom(after_coords, H)
                    if H <= after_coords
                    else float("-inf")
                )
                continue
            total = float("-inf")
            for gap in range(gap_min, gap_max + 1):
                after_coords = b * (late_blocks + gap - 1)
                if H <= after_coords:
                    total = log2add(total, log2_binom(after_coords, H))
            out[(H, idx)] = total
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outer-prefix-csv", type=Path, required=True)
    parser.add_argument("--N", type=int, default=2**21)
    parser.add_argument("--block-bits", type=int, default=64)
    parser.add_argument("--h-values", default="32:500")
    parser.add_argument("--first-r-values", default="1:64")
    parser.add_argument("--late-blocks", type=int, default=5949)
    parser.add_argument("--gap-start", type=int, default=1)
    parser.add_argument("--gap-stop", type=int, default=12000)
    parser.add_argument("--gap-step", type=int, default=4000)
    parser.add_argument(
        "--inner-knots-by-gap",
        help=(
            "Semicolon-separated knot lists, one per gap bucket. "
            "Each list is H:log2,H:log2,..."
        ),
    )
    parser.add_argument(
        "--inner-knot-csvs",
        help=(
            "Semicolon-separated CSV paths, one per gap bucket.  Each CSV is "
            "read from --knot-h-column and --knot-value-column."
        ),
    )
    parser.add_argument("--knot-h-column", default="H")
    parser.add_argument("--knot-value-column", default="total_log2")
    parser.add_argument(
        "--inner-envelope-lift-bits",
        type=float,
        default=0.0,
        help="Add this many bits to each interpolated inner log2 bound.",
    )
    parser.add_argument(
        "--inner-envelope-lift-bits-by-gap",
        help="Comma-separated per-gap additive lifts; overrides --inner-envelope-lift-bits.",
    )
    parser.add_argument(
        "--require-knot-coverage",
        action="store_true",
        help="Fail instead of extrapolating outside the supplied H knots.",
    )
    parser.add_argument("--output-csv", type=Path)
    args = parser.parse_args()

    outer = load_outer(args.outer_prefix_csv)
    gap_ranges: list[tuple[int, int]] = []
    for gap_min in range(args.gap_start, args.gap_stop + 1, args.gap_step):
        gap_max = min(args.gap_stop, gap_min + args.gap_step - 1)
        gap_ranges.append((gap_min, gap_max))
    if bool(args.inner_knots_by_gap) == bool(args.inner_knot_csvs):
        raise ValueError("supply exactly one of --inner-knots-by-gap or --inner-knot-csvs")
    if args.inner_knot_csvs:
        knot_lists = [
            load_knot_csv(
                Path(part.strip()),
                h_column=args.knot_h_column,
                value_column=args.knot_value_column,
            )
            for part in args.inner_knot_csvs.split(";")
            if part.strip()
        ]
    else:
        knot_lists = [parse_knots(part) for part in args.inner_knots_by_gap.split(";")]
    if len(knot_lists) != len(gap_ranges):
        raise ValueError("inner knots must have one semicolon-separated list/path per gap bucket")
    if args.inner_envelope_lift_bits_by_gap:
        lift_bits = parse_float_list(args.inner_envelope_lift_bits_by_gap)
        if len(lift_bits) != len(gap_ranges):
            raise ValueError("--inner-envelope-lift-bits-by-gap must have one value per gap bucket")
    else:
        lift_bits = [args.inner_envelope_lift_bits for _ in gap_ranges]
    h_values = parse_int_list(args.h_values)
    r_values = parse_int_list(args.first_r_values)
    max_h_after_first = max((h - 1 for h in h_values if h in outer), default=0)
    gap_sums = precompute_gap_sums(
        b=args.block_bits,
        late_blocks=args.late_blocks,
        gap_ranges=gap_ranges,
        max_h_after_first=max_h_after_first,
    )

    rows = []
    total = float("-inf")
    peak = None
    by_gap = [float("-inf") for _ in gap_ranges]
    by_r = {}
    for h in h_values:
        if h not in outer:
            continue
        for r in r_values:
            if r < 1 or r > min(h, args.block_bits):
                continue
            H = h - r
            for idx, (gap_min, gap_max) in enumerate(gap_ranges):
                gap_sum = gap_sums.get((H, idx), float("-inf"))
                placement = (
                    log2_binom(args.block_bits, r)
                    + gap_sum
                    - log2_binom(args.N, h)
                    if gap_sum != float("-inf")
                    else float("-inf")
                )
                if placement == float("-inf"):
                    continue
                inner = min(
                    0.0,
                    envelope_value(
                        knot_lists[idx],
                        H,
                        require_coverage=args.require_knot_coverage,
                    )
                    + lift_bits[idx],
                )
                term = outer[h] + placement + inner
                total = log2add(total, term)
                by_gap[idx] = log2add(by_gap[idx], term)
                by_r[r] = log2add(by_r.get(r, float("-inf")), term)
                row = {
                    "outer_weight": h,
                    "first_r": r,
                    "remaining_ones": H,
                    "gap_min": gap_min,
                    "gap_max": gap_max,
                    "outer_log2": outer[h],
                    "placement_log2": placement,
                    "inner_log2": inner,
                    "term_log2": term,
                }
                rows.append(row)
                if peak is None or term > peak["term_log2"]:
                    peak = row

    fields = [
        "outer_weight",
        "first_r",
        "remaining_ones",
        "gap_min",
        "gap_max",
        "outer_log2",
        "placement_log2",
        "inner_log2",
        "term_log2",
    ]
    if args.output_csv:
        with args.output_csv.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)

    print("Full-split episode envelope sum")
    print(f"N={args.N}, b={args.block_bits}, h_values={args.h_values}, r_values={args.first_r_values}")
    print(f"total_log2,{total:.6f}")
    print(f"margin_bits,{-total:.6f}")
    for (gap_min, gap_max), val in zip(gap_ranges, by_gap):
        print(f"gap_{gap_min}_{gap_max}_log2,{val:.6f}")
    for r, val in sorted(by_r.items())[:8]:
        print(f"r_{r}_log2,{val:.6f}")
    if peak is not None:
        print(f"peak_h,{peak['outer_weight']}")
        print(f"peak_first_r,{peak['first_r']}")
        print(f"peak_gap_min,{peak['gap_min']}")
        print(f"peak_gap_max,{peak['gap_max']}")
        print(f"peak_inner_log2,{peak['inner_log2']:.6f}")
        print(f"peak_term_log2,{peak['term_log2']:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
