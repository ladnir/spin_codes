#!/usr/bin/env python3
"""Probe block-outer upgrades against the fixed-tap dense inner certificate.

The goal is to test finite-n implementation upgrades without changing the
main dense+dense proof path.  The default local model is intentionally
pessimistic: a systematic [2b,b,d0] block contributes no nonzero words below
d0, and has at most 2^b-1 nonzero words in total.  For 0<z<1 this gives

    W_loc(z) <= 1 + (2^b-1) z^d0.

The global direct-sum outer has M=ceil(k/b) blocks, so low-weight coefficients
are bounded by optimizing W_loc(z)^M / z^h.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt

from certify_global_episode_cover import (
    binom_cdf_half_be_log2,
    global_episode_inner_log2,
    survivor_only_exact_prefix_log2,
)
from dense_largek_eval import log2_binom, log2add


def parse_int_list(text: str) -> list[int]:
    out: list[int] = []
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        if ":" in part:
            vals = [int(x) for x in part.split(":")]
            if len(vals) == 2:
                lo, hi = vals
                step = 1
            elif len(vals) == 3:
                lo, hi, step = vals
            else:
                raise ValueError(f"bad integer range: {part!r}")
            out.extend(range(lo, hi + 1, step))
        else:
            out.append(int(part))
    return out


def z_grid(z_min: float, z_max: float, count: int) -> list[float]:
    if count <= 1:
        return [z_max]
    lo = math.log(z_min)
    hi = math.log(z_max)
    return [math.exp(lo + (hi - lo) * i / (count - 1)) for i in range(count)]


def log2_sub_one(log_x: float) -> float:
    if log_x == float("-inf"):
        return float("-inf")
    if log_x <= 1e-10:
        val = math.expm1(log_x * math.log(2.0))
        return math.log2(val) if val > 0.0 else float("-inf")
    return log_x + math.log2(1.0 - 2.0 ** (-log_x))


def format_log2(x: float) -> str:
    if x == float("-inf"):
        return "-inf"
    return f"{x:.6f}"


def local_floor_total_log2(block_bits: int, d0: int, z: float) -> float:
    """log2(1 + (2^b-1) z^d0), computed stably."""
    nonzero = block_bits + d0 * math.log2(z)
    # The -1 in 2^b-1 is irrelevant at these b, but keep it exact enough for
    # small smoke-test block sizes.
    if block_bits <= 50:
        nonzero = math.log2((2.0**block_bits - 1.0) * (z**d0))
    if nonzero < -60.0:
        return math.log2(1.0 + 2.0**nonzero)
    return math.log2(1.0 + 2.0**nonzero)


def local_random_like_log2(block_bits: int, d0: int, z: float) -> float:
    """Expected random-linear local spectrum with a hard floor below d0.

    This is a heuristic comparison lane, not a theorem input.
    """
    n = 2 * block_bits
    rate_factor = (2.0**block_bits - 1.0) / (2.0**n - 1.0) if n <= 100 else 2.0 ** (-block_bits)
    total = 1.0
    for j in range(d0, n + 1):
        total += rate_factor * math.comb(n, j) * (z**j)
    return math.log2(total)


def load_local_spectrum(path: str | None) -> list[tuple[int, int]]:
    if path is None:
        return []
    rows: list[tuple[int, int]] = []
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            rows.append((int(row["weight"]), int(row["count"])))
    return rows


def local_spectrum_log2(spectrum: list[tuple[int, int]], z: float) -> float:
    total = 0.0
    for weight, count in spectrum:
        total += count * (z**weight)
    return math.log2(total) if total > 0.0 else float("-inf")


def outer_block_gf_bounds(
    *,
    blocks: int,
    block_bits: int,
    d0: int,
    h_max: int,
    zs: list[float],
    model: str,
    spectrum: list[tuple[int, int]],
) -> tuple[list[float], list[float]]:
    vals = [float("-inf")] * (h_max + 1)
    best_z = [float("nan")] * (h_max + 1)
    local_logs: list[tuple[float, float]] = []
    for z in zs:
        if model == "floor-total":
            local = local_floor_total_log2(block_bits, d0, z)
        elif model == "random-like":
            local = local_random_like_log2(block_bits, d0, z)
        elif model == "spectrum-csv":
            local = local_spectrum_log2(spectrum, z)
        else:
            raise ValueError(f"unknown local model: {model}")
        local_logs.append((z, log2_sub_one(blocks * local)))

    for h in range(d0, h_max + 1):
        best = float("inf")
        best_here = float("nan")
        for z, global_log in local_logs:
            if global_log == float("-inf"):
                continue
            candidate = global_log - h * math.log2(z)
            if candidate < best:
                best = candidate
                best_here = z
        vals[h] = best if best != float("inf") else float("-inf")
        best_z[h] = best_here
    return vals, best_z


def singleton_block_floor_log2(blocks: int, block_bits: int, d0: int, h: int) -> float:
    """A sharper optional bound for one active block.

    This is only valid for h within one local block and uses the ambient volume.
    It is often better than the total-count floor right near h=d0.
    """
    if h < d0 or h > 2 * block_bits:
        return float("-inf")
    return math.log2(blocks) + min(block_bits, log2_binom(2 * block_bits, h))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--k", type=int, default=2**20)
    parser.add_argument("--block-bits", type=int, default=80)
    parser.add_argument("--sigma", type=int, default=32)
    parser.add_argument("--delta", type=float, default=0.106)
    parser.add_argument("--d0s", default="8:24")
    parser.add_argument("--h-max", type=int, default=240)
    parser.add_argument("--r-max", type=int, default=12)
    parser.add_argument("--z-min", type=float, default=1e-8)
    parser.add_argument("--z-max", type=float, default=0.999)
    parser.add_argument("--z-count", type=int, default=220)
    parser.add_argument("--block-ratio", type=float, default=1.02)
    parser.add_argument("--suffix-block-ratio", type=float, default=1.1)
    parser.add_argument("--model", choices=("floor-total", "random-like", "spectrum-csv"), default="floor-total")
    parser.add_argument(
        "--local-spectrum-csv",
        default=None,
        help="CSV with columns weight,count for one local block. Required for --model spectrum-csv.",
    )
    parser.add_argument(
        "--expect-spectrum-dim",
        type=int,
        default=None,
        help="When using --model spectrum-csv, require the local counts to sum to 2^this value.",
    )
    parser.add_argument(
        "--singleton-volume",
        action="store_true",
        help="For h<=2b, also intersect with the one-active-block ambient-volume bound.",
    )
    parser.add_argument("--out-prefix", default=None)
    args = parser.parse_args()

    blocks = (args.k + args.block_bits - 1) // args.block_bits
    k_eff = blocks * args.block_bits
    n = 2 * k_eff
    d = math.floor(args.delta * n)
    d0s = parse_int_list(args.d0s)
    zs = z_grid(args.z_min, args.z_max, args.z_count)
    spectrum = load_local_spectrum(args.local_spectrum_csv)
    if args.model == "spectrum-csv":
        if not spectrum:
            raise ValueError("--model spectrum-csv requires --local-spectrum-csv")
        if args.expect_spectrum_dim is not None:
            total_words = sum(count for _, count in spectrum)
            expected_words = 1 << args.expect_spectrum_dim
            if total_words != expected_words:
                raise ValueError(
                    f"local spectrum counts sum to {total_words}, expected 2^{args.expect_spectrum_dim}"
                )
        d0s = [min(weight for weight, count in spectrum if weight > 0 and count > 0)]

    if args.out_prefix is None:
        out_base = Path(__file__).with_name(
            f"block_outer_probe_k{args.k}_b{args.block_bits}_sig{args.sigma}_d{args.delta:g}_h{args.h_max}_{args.model}"
        )
    else:
        out_base = Path(args.out_prefix)
        if not out_base.is_absolute():
            out_base = Path.cwd() / out_base

    print(f"effective k={k_eff} in {blocks} blocks; N={n}, d={d}", flush=True)
    print("precomputing inner suffix tail and episode-cover values", flush=True)
    tail_cache = [binom_cdf_half_be_log2(u, d) for u in range(n + 1)]
    exact_survivor = survivor_only_exact_prefix_log2(n=n, h_max=args.h_max, tail_cache=tail_cache)

    inner_logs = [float("-inf")] * (args.h_max + 1)
    for h in range(1, args.h_max + 1):
        inner, _, _, _ = global_episode_inner_log2(
            n=n,
            h=h,
            sigma=args.sigma,
            delta=args.delta,
            r_max=args.r_max,
            block_ratio=args.block_ratio,
            suffix_block_ratio=args.suffix_block_ratio,
            stop_gap_bits=60.0,
            tail_confirm=8,
            suffix_survivor=True,
            tail_mode="be",
            exact_survivor=True,
            tail_cache=tail_cache,
            exact_survivor_value=exact_survivor[h],
        )
        inner_logs[h] = inner
        if h == 1 or h % 50 == 0 or h == args.h_max:
            print(f"  inner h={h:4d}: {format_log2(inner)}", flush=True)

    rows: list[dict[str, str | int | float]] = []
    summaries: list[dict[str, str | int | float]] = []
    margin_points: list[tuple[int, float, int]] = []
    for d0 in d0s:
        print(f"d0={d0}: evaluating block outer model", flush=True)
        outer_logs, outer_z = outer_block_gf_bounds(
            blocks=blocks,
            block_bits=args.block_bits,
            d0=d0,
            h_max=args.h_max,
            zs=zs,
            model=args.model,
            spectrum=spectrum,
        )
        if args.singleton_volume:
            singleton_hi = min(args.h_max, 2 * args.block_bits, 2 * d0 - 1)
            for h in range(d0, singleton_hi + 1):
                one = singleton_block_floor_log2(blocks, args.block_bits, d0, h)
                if one < outer_logs[h]:
                    outer_logs[h] = one
                    outer_z[h] = float("nan")

        total = float("-inf")
        peak_h = -1
        peak_term = float("-inf")
        peak_outer = float("-inf")
        peak_inner = float("-inf")
        for h in range(1, args.h_max + 1):
            outer = outer_logs[h]
            inner = inner_logs[h]
            term = outer + inner if outer != float("-inf") and inner != float("-inf") else float("-inf")
            total = log2add(total, term)
            if term > peak_term:
                peak_h = h
                peak_term = term
                peak_outer = outer
                peak_inner = inner
            rows.append(
                {
                    "k": args.k,
                    "k_eff": k_eff,
                    "N": n,
                    "delta": f"{args.delta:.12g}",
                    "d": d,
                    "sigma": args.sigma,
                    "block_bits": args.block_bits,
                    "blocks": blocks,
                    "d0": d0,
                    "h": h,
                    "outer_log2": format_log2(outer),
                    "outer_z": "" if math.isnan(outer_z[h]) else f"{outer_z[h]:.12g}",
                    "inner_log2": format_log2(inner),
                    "term_log2": format_log2(term),
                    "cum_log2": format_log2(total),
                }
            )
        margin = -total
        margin_points.append((d0, margin, peak_h))
        summaries.append(
            {
                "k": args.k,
                "k_eff": k_eff,
                "N": n,
                "delta": f"{args.delta:.12g}",
                "d": d,
                "sigma": args.sigma,
                "block_bits": args.block_bits,
                "blocks": blocks,
                "d0": d0,
                "h_max": args.h_max,
                "model": args.model,
                "singleton_volume": int(args.singleton_volume),
                "total_log2": format_log2(total),
                "margin_bits": format_log2(margin),
                "peak_h": peak_h,
                "peak_term_log2": format_log2(peak_term),
                "peak_outer_log2": format_log2(peak_outer),
                "peak_inner_log2": format_log2(peak_inner),
            }
        )
        print(f"  total={format_log2(total)} margin={margin:.3f} peak_h={peak_h} peak={peak_term:.3f}", flush=True)

    csv_path = out_base.with_suffix(".csv")
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    summary_path = out_base.with_name(out_base.name + "_summary").with_suffix(".csv")
    with summary_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)
    print(f"wrote {csv_path}")
    print(f"wrote {summary_path}")

    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    xs = [p[0] for p in margin_points]
    ys = [p[1] for p in margin_points]
    peaks = [p[2] for p in margin_points]
    ax.plot(xs, ys, marker="o", linewidth=2)
    ax.axhline(40.0, color="black", linestyle="--", linewidth=1)
    for x, y, h in zip(xs, ys, peaks):
        ax.annotate(f"h={h}", (x, y), textcoords="offset points", xytext=(0, 7), ha="center", fontsize=8)
    ax.set_xlabel("local block minimum distance d0")
    ax.set_ylabel(f"-log2 first moment, h <= {args.h_max} (bits)")
    ax.set_title(
        f"Block outer probe, b={args.block_bits}, sigma={args.sigma}, delta={args.delta:g}, {args.model}"
    )
    ax.grid(True, alpha=0.3)
    png_path = out_base.with_suffix(".png")
    fig.savefig(png_path, dpi=180)
    print(f"wrote {png_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
