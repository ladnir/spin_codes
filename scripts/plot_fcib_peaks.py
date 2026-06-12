"""Plot per-h paired-cover first-moment contributions.

This is a lightweight visualization wrapper around the theorem-safe
fixed-tap paired-cover bound in dense_largek_eval.py. It writes a CSV with
per-h log2 contributions and, when matplotlib is available, a PNG plot.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

from dense_largek_eval import (
    fixedtap_cover_internal_block_bound_log2,
    fixedtap_cover_internal_bound_log2,
    log2add,
    outer_small_h_prefix,
)


def parse_sigmas(text: str) -> list[int]:
    out: list[int] = []
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        if ":" in part:
            bits = [int(x) for x in part.split(":")]
            if len(bits) == 2:
                lo, hi = bits
                step = 1
            elif len(bits) == 3:
                lo, hi, step = bits
            else:
                raise ValueError(f"bad sigma range: {part}")
            out.extend(range(lo, hi + 1, step))
        else:
            out.append(int(part))
    return out


def format_log2(x: float) -> str:
    if x == float("-inf"):
        return "-inf"
    return f"{x:.6f}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--k", type=int, required=True, help="Message length k.")
    parser.add_argument("--sigmas", required=True, help="Comma list/ranges, e.g. 31,32 or 30:34.")
    parser.add_argument("--delta", type=float, default=0.106)
    parser.add_argument("--h-max", type=int, default=80)
    parser.add_argument("--h-min", type=int, default=1)
    parser.add_argument("--small-xi", type=float, default=0.1)
    parser.add_argument("--block-ratio", type=float, default=1.01)
    parser.add_argument(
        "--mode",
        choices=("block", "sampled"),
        default="block",
        help="'block' is theorem-safe; 'sampled' is the faster diagnostic.",
    )
    parser.add_argument(
        "--outer-smallh-mode",
        choices=("conv", "banded", "banded-fixedtap"),
        default="banded-fixedtap",
    )
    parser.add_argument("--parity-n", type=int, default=None)
    parser.add_argument("--out-prefix", default=None, help="Output prefix without extension.")
    parser.add_argument("--no-png", action="store_true", help="Only write CSV.")
    args = parser.parse_args()

    n = 2 * args.k
    parity_n = args.parity_n if args.parity_n is not None else args.k
    sigmas = parse_sigmas(args.sigmas)
    l_window = math.ceil((2.0 + args.small_xi) * args.delta * n)
    if args.out_prefix is None:
        sig_text = "_".join(str(s) for s in sigmas)
        out_base = Path(__file__).resolve().parent / (
            f"fcib_peaks_k{args.k}_sig{sig_text}_d{args.delta:g}_h{args.h_max}_{args.mode}"
        )
    else:
        out_base = Path(args.out_prefix)
        if not out_base.is_absolute():
            out_base = Path.cwd() / out_base

    rows: list[dict[str, str | int | float]] = []
    series: dict[int, list[tuple[int, float, float]]] = {}
    for sigma in sigmas:
        print(f"sigma={sigma}: computing outer prefix through h={args.h_max}")
        outer_logs = outer_small_h_prefix(
            args.k,
            sigma,
            args.h_max,
            outer_mode=args.outer_smallh_mode,
            parity_n=parity_n,
        )
        cumulative = float("-inf")
        points: list[tuple[int, float, float]] = []
        for h in range(args.h_min, args.h_max + 1):
            out = outer_logs[h]
            if out == float("-inf"):
                inner = float("-inf")
                term = float("-inf")
            else:
                if args.mode == "block":
                    inner = fixedtap_cover_internal_block_bound_log2(
                        n,
                        h,
                        sigma,
                        args.delta,
                        l_window,
                        block_ratio=args.block_ratio,
                    )
                else:
                    inner = fixedtap_cover_internal_bound_log2(n, h, sigma, args.delta, l_window)
                term = out + inner
                cumulative = log2add(cumulative, term)
            rows.append(
                {
                    "sigma": sigma,
                    "offset": sigma - math.ceil(math.log2(args.k)),
                    "h": h,
                    "outer_log2": format_log2(out),
                    "inner_log2": format_log2(inner),
                    "term_log2": format_log2(term),
                    "cum_log2": format_log2(cumulative),
                }
            )
            points.append((h, term, cumulative))
            if h == args.h_min or h % 8 == 0 or h == args.h_max:
                print(f"  h={h:3d} term={format_log2(term)} cum={format_log2(cumulative)}")
        series[sigma] = points
        finite = [(h, t) for h, t, _ in points if t != float("-inf")]
        if finite:
            peak_h, peak_t = max(finite, key=lambda x: x[1])
            boundary = " (RIGHT BOUNDARY: increase h-max before interpreting)" if peak_h == args.h_max else ""
            print(f"  peak h={peak_h}, log2 term={format_log2(peak_t)}{boundary}")

    csv_path = out_base.with_suffix(".csv")
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["sigma", "offset", "h", "outer_log2", "inner_log2", "term_log2", "cum_log2"],
        )
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {csv_path}")

    if args.no_png:
        return
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not available; skipped PNG")
        return

    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True, constrained_layout=True)
    for sigma, points in series.items():
        hs = [p[0] for p in points]
        terms = [p[1] for p in points]
        cums = [p[2] for p in points]
        label = f"sigma={sigma} (C={sigma - math.ceil(math.log2(args.k))})"
        axes[0].plot(hs, terms, marker=".", linewidth=1.4, label=label)
        axes[1].plot(hs, cums, marker=".", linewidth=1.4, label=label)
        finite = [(h, t) for h, t, _ in points if t != float("-inf")]
        if finite:
            peak_h, peak_t = max(finite, key=lambda x: x[1])
            axes[0].scatter([peak_h], [peak_t], s=48, zorder=4)
            note = " boundary" if peak_h == args.h_max else ""
            axes[0].annotate(
                f"h={peak_h}{note}",
                (peak_h, peak_t),
                xytext=(6, 8),
                textcoords="offset points",
                fontsize=8,
            )

    axes[0].axhline(0.0, color="black", linewidth=0.8, linestyle="--")
    axes[1].axhline(0.0, color="black", linewidth=0.8, linestyle="--")
    axes[0].set_ylabel("log2 contribution by h")
    axes[1].set_ylabel("log2 cumulative prefix")
    axes[1].set_xlabel("outer weight h")
    axes[0].set_title(
        f"Fixed-tap paired-cover {args.mode} curve, k={args.k}, delta={args.delta}, L={l_window}"
    )
    axes[0].legend()
    axes[1].legend()
    png_path = out_base.with_suffix(".png")
    png_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, dpi=180)
    print(f"Wrote {png_path}")


if __name__ == "__main__":
    main()
