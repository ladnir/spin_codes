"""Exploratory dense+dense burst envelope.

This script separates the first-moment calculation into proof lanes.  The first
implemented medium-weight lane is the one-burst (r=1) part of the fixed-tap
paired-cover theorem.

The medium outer envelope is the active fixed-tap banded generating-function
bound A_h <= W_ft-band(z) z^{-h}, optimized over a user-controlled grid of z
values. Uncontrolled proof lanes are reported as NOT_CERTIFIED instead of being
silently omitted.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

from dense_largek_eval import (
    ETA_CRIT,
    binom_cdf_half_entropy_log2,
    fixedtap_banded_outer_gf_log2,
    h2,
    linear_window_gap,
    log2_binom,
    log2add,
)


def format_log2(x: float) -> str:
    if x == float("-inf"):
        return "-inf"
    return f"{x:.6f}"


def parse_sigmas(text: str) -> list[int]:
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
                raise ValueError(f"bad sigma range {part!r}")
            out.extend(range(lo, hi + 1, step))
        else:
            out.append(int(part))
    return out


def one_burst_inner_block_log2(
    *,
    n: int,
    h: int,
    sigma: int,
    delta: float,
    l_window: int,
    block_ratio: float,
) -> float:
    """Blocked upper sum for the r=1 paired-cover inner term.

    For r=1 the paired-cover count is

        (N-T+1) * C(T+L-2, h-2) / C(N,h).

    The T block [a,b] uses N-a+1 for the decreasing placement factor,
    b+L-2 for the increasing support factor, and the binomial tail at a.
    """
    if h < 2:
        return float("-inf")
    if block_ratio <= 1.0:
        raise ValueError("block_ratio must be > 1")
    cut = math.floor(delta * n)
    denom = log2_binom(n, h)
    total = float("-inf")
    t_min = 2
    t_max = min(n, l_window + 1)
    a = t_min
    while a <= t_max:
        b = min(t_max, max(a, int(math.floor(a * block_ratio))))
        slots = b + l_window - 2
        if slots >= h - 2:
            cover = log2_binom(n - a + 1, 1) + log2_binom(slots, h - 2) - denom
            fair_len = max(0, a - sigma)
            fair_tail = binom_cdf_half_entropy_log2(fair_len, cut)
            term = math.log2(b - a + 1) + cover - (sigma - 1) + fair_tail
            total = log2add(total, term)
        a = b + 1
    return min(0.0, total)


def z_grid(z_min: float, z_max: float, count: int) -> list[float]:
    if count <= 1:
        return [z_max]
    lo = math.log(z_min)
    hi = math.log(z_max)
    return [math.exp(lo + (hi - lo) * i / (count - 1)) for i in range(count)]


def fixedtap_banded_coeff_bounds(
    *,
    k: int,
    parity_n: int,
    sigma: int,
    h_min: int,
    h_max: int,
    zs: list[float],
) -> tuple[list[float], list[float]]:
    logs = [float("-inf")] * (h_max + 1)
    best_z = [float("nan")] * (h_max + 1)
    for z in zs:
        log_w = fixedtap_banded_outer_gf_log2(k, parity_n, sigma, z)
        log_z = math.log2(z)
        for h in range(h_min, h_max + 1):
            val = log_w - h * log_z
            if logs[h] == float("-inf") or val < logs[h]:
                logs[h] = val
                best_z[h] = z
    return logs, best_z


def compute_sigma(
    *,
    k: int,
    sigma: int,
    delta: float,
    h_min: int,
    h_max: int,
    small_xi: float,
    z_min: float,
    z_max: float,
    z_count: int,
    block_ratio: float,
    theta: float,
    xi: float,
    eta_hi: float,
    gap_step: float,
) -> tuple[list[dict[str, str | int | float]], dict[str, str | int | float]]:
    n = 2 * k
    parity_n = k
    l_window = math.ceil((2.0 + small_xi) * delta * n)
    outer_logs, outer_zs = fixedtap_banded_coeff_bounds(
        k=k,
        parity_n=parity_n,
        sigma=sigma,
        h_min=h_min,
        h_max=h_max,
        zs=z_grid(z_min, z_max, z_count),
    )
    rows: list[dict[str, str | int | float]] = []
    total = float("-inf")
    peak_h = -1
    peak_term = float("-inf")
    for h in range(h_min, h_max + 1):
        outer = outer_logs[h]
        z_used = outer_zs[h]
        inner = one_burst_inner_block_log2(
            n=n,
            h=h,
            sigma=sigma,
            delta=delta,
            l_window=l_window,
            block_ratio=block_ratio,
        )
        term = outer + inner if outer != float("-inf") and inner != float("-inf") else float("-inf")
        total = log2add(total, term)
        if term > peak_term:
            peak_term = term
            peak_h = h
        rows.append(
            {
                "sigma": sigma,
                "offset": sigma - math.ceil(math.log2(k)),
                "h": h,
                "outer_gf_log2": format_log2(outer),
                "outer_z": f"{z_used:.12g}",
                "one_burst_inner_log2": format_log2(inner),
                "one_burst_term_log2": format_log2(term),
                "one_burst_cum_log2": format_log2(total),
            }
        )

    gap = linear_window_gap(delta, theta, xi, eta_hi, gap_step)
    lin_count = max(0, int(math.floor(min(eta_hi, 1.0 - theta) * n)) - int(math.ceil(ETA_CRIT * n)) + 1)
    log2_s_lin = float("-inf") if lin_count == 0 else math.log2(lin_count) + n * gap.worst_gap
    top_exp = h2(eta_hi) - 0.5
    top_count = max(0, n - int(math.ceil(eta_hi * n)) + 1)
    log2_s_top = float("-inf") if top_count == 0 else math.log2(top_count) + n * top_exp

    summary: dict[str, str | int | float] = {
        "sigma": sigma,
        "offset": sigma - math.ceil(math.log2(k)),
        "h_min": h_min,
        "h_max": h_max,
        "one_burst_total_log2": format_log2(total),
        "one_burst_peak_h": peak_h,
        "one_burst_peak_h_label": f"{peak_h}+" if peak_h == h_max else str(peak_h),
        "one_burst_peak_term_log2": format_log2(peak_term),
        "linear_top_log2": format_log2(log2add(log2_s_lin, log2_s_top)),
        "linear_gap_worst": gap.worst_gap,
        "linear_gap_eta": gap.worst_eta,
        "one_burst_outer_status": "FIXED_TAP_BANDED_GF_GRID",
        "multi_episode_medium": "NOT_CERTIFIED",
        "tail_after_h_max_before_linear": "NOT_CERTIFIED",
    }
    return rows, summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--k", type=int, required=True)
    parser.add_argument("--sigmas", required=True)
    parser.add_argument("--delta", type=float, default=0.106)
    parser.add_argument("--h-min", type=int, default=2)
    parser.add_argument("--h-max", type=int, default=512)
    parser.add_argument("--small-xi", type=float, default=0.1)
    parser.add_argument("--z-min", type=float, default=0.02)
    parser.add_argument("--z-max", type=float, default=0.999)
    parser.add_argument("--z-count", type=int, default=256)
    parser.add_argument("--block-ratio", type=float, default=1.02)
    parser.add_argument("--theta", type=float, default=0.005)
    parser.add_argument("--xi", type=float, default=8.0)
    parser.add_argument("--eta-hi", type=float, default=0.99)
    parser.add_argument("--gap-step", type=float, default=1e-5)
    parser.add_argument("--out-prefix", default=None)
    parser.add_argument("--no-png", action="store_true")
    args = parser.parse_args()

    sigmas = parse_sigmas(args.sigmas)
    if args.out_prefix is None:
        sig_text = "_".join(str(s) for s in sigmas)
        out_base = Path(__file__).resolve().parent / (
            f"burst_cert_k{args.k}_sig{sig_text}_d{args.delta:g}_h{args.h_min}_{args.h_max}"
        )
    else:
        out_base = Path(args.out_prefix)
        if not out_base.is_absolute():
            out_base = Path.cwd() / out_base

    all_rows: list[dict[str, str | int | float]] = []
    summaries: list[dict[str, str | int | float]] = []
    series: dict[int, list[tuple[int, float, float]]] = {}
    for sigma in sigmas:
        print(f"sigma={sigma}: one-burst envelope h={args.h_min}..{args.h_max}")
        rows, summary = compute_sigma(
            k=args.k,
            sigma=sigma,
            delta=args.delta,
            h_min=args.h_min,
            h_max=args.h_max,
            small_xi=args.small_xi,
            z_min=args.z_min,
            z_max=args.z_max,
            z_count=args.z_count,
            block_ratio=args.block_ratio,
            theta=args.theta,
            xi=args.xi,
            eta_hi=args.eta_hi,
            gap_step=args.gap_step,
        )
        all_rows.extend(rows)
        summaries.append(summary)
        points: list[tuple[int, float, float]] = []
        for row in rows:
            h = int(row["h"])
            term = float(row["one_burst_term_log2"]) if row["one_burst_term_log2"] != "-inf" else float("-inf")
            cum = float(row["one_burst_cum_log2"]) if row["one_burst_cum_log2"] != "-inf" else float("-inf")
            points.append((h, term, cum))
        series[sigma] = points
        print(
            "  one-burst total log2={one_burst_total_log2}, peak={one_burst_peak_h_label}, "
            "peak term={one_burst_peak_term_log2}, medium r>=2={multi_episode_medium}".format(**summary)
        )

    out_base.parent.mkdir(parents=True, exist_ok=True)
    csv_path = out_base.with_suffix(".csv")
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "sigma",
                "offset",
                "h",
                "outer_gf_log2",
                "outer_z",
                "one_burst_inner_log2",
                "one_burst_term_log2",
                "one_burst_cum_log2",
            ],
        )
        writer.writeheader()
        writer.writerows(all_rows)
    print(f"Wrote {csv_path}")

    summary_path = out_base.with_name(out_base.name + "_summary.csv")
    with summary_path.open("w", newline="") as f:
        fieldnames = list(summaries[0].keys()) if summaries else []
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(summaries)
    print(f"Wrote {summary_path}")

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
        axes[0].plot(hs, terms, linewidth=1.4, label=label)
        axes[1].plot(hs, cums, linewidth=1.4, label=label)
        finite = [(h, t) for h, t, _ in points if t != float("-inf")]
        if finite:
            peak_h, peak_t = max(finite, key=lambda p: p[1])
            axes[0].scatter([peak_h], [peak_t], s=42, zorder=4)
            note = "+" if peak_h == args.h_max else ""
            axes[0].annotate(
                f"h={peak_h}{note}",
                (peak_h, peak_t),
                xytext=(6, 8),
                textcoords="offset points",
                fontsize=8,
            )
    axes[0].axhline(0.0, color="black", linewidth=0.8, linestyle="--")
    axes[1].axhline(0.0, color="black", linewidth=0.8, linestyle="--")
    axes[0].set_ylabel("log2 one-burst contribution")
    axes[1].set_ylabel("log2 cumulative one-burst envelope")
    axes[1].set_xlabel("outer weight h")
    axes[0].set_title(
        f"One-burst envelope, k={args.k}, delta={args.delta}, fixed-tap banded GF grid"
    )
    axes[0].legend()
    axes[1].legend()
    png_path = out_base.with_suffix(".png")
    fig.savefig(png_path, dpi=180)
    print(f"Wrote {png_path}")


if __name__ == "__main__":
    main()
