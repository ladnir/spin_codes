"""Scan the corrected one-episode crossover past the exact low-h prefix.

This is a calibration tool, not a final theorem certificate.  It combines a
fixed-tap banded outer coefficient upper bound

    A_h <= min_z W(z) z^{-h}

with the corrected r=1 early-burst plus late-tail inner charge.  In the default
paired-budget mode, survival into the final suffix is charged against the shared
distance budget by a binomial tail over the combined early-plus-late live
length.  The late-tail sum is blocked: the final 2d suffix is summed exactly by
the hockey-stick identity, and the remaining short strip is upper-bounded on
geometric blocks.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

from dense_largek_eval import (
    binom_cdf_half_entropy_log2,
    fixedtap_banded_outer_gf_log2,
    log2_binom,
    log2add,
)


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
                raise ValueError(f"bad sigma range: {part}")
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


def format_log2(x: float) -> str:
    if x == float("-inf"):
        return "-inf"
    return f"{x:.6f}"


def load_exact_outer(path: str | None, sigma: int, h_max: int) -> list[float]:
    vals = [float("-inf")] * (h_max + 1)
    if path is None:
        return vals
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            if int(row["sigma"]) != sigma:
                continue
            h = int(row["h"])
            if 0 <= h <= h_max and "outer_log2" in row:
                vals[h] = float(row["outer_log2"])
    return vals


def outer_gf_bounds(
    *,
    k: int,
    sigma: int,
    h_min: int,
    h_max: int,
    zs: list[float],
) -> tuple[list[float], list[float]]:
    logs = [float("-inf")] * (h_max + 1)
    best_z = [float("nan")] * (h_max + 1)
    log_ws: list[tuple[float, float]] = []
    for z in zs:
        log_ws.append((z, fixedtap_banded_outer_gf_log2(k, k, sigma, z)))
    for h in range(h_min, h_max + 1):
        best = float("inf")
        best_here = float("nan")
        for z, log_w in log_ws:
            val = log_w - h * math.log2(z)
            if val < best:
                best = val
                best_here = z
        logs[h] = best
        best_z[h] = best_here
    return logs, best_z


def late_sum_blocked_logs(
    *,
    h_max: int,
    sigma: int,
    d: int,
    l_window: int,
    early_fair_len: int,
    block_ratio: float,
    budget_mode: str,
    late_termination: bool,
) -> list[float]:
    logs = [float("-inf")] * (h_max + 1)
    final_suffix = min(l_window, 2 * d)
    p_term = 2.0 ** (-(sigma - 1))
    early_tail_log = binom_cdf_half_entropy_log2(early_fair_len, d)
    early_tail = 1.0 if early_tail_log >= 0.0 else 2.0**early_tail_log
    for j in range(0, h_max + 1):
        if j == 0:
            logs[j] = math.log2(early_tail) if early_tail > 0.0 else float("-inf")
            continue
        if budget_mode == "separate":
            total = log2_binom(final_suffix, j)
        else:
            combined_log = binom_cdf_half_entropy_log2(early_fair_len + final_suffix, d)
            combined_prob = 1.0 if combined_log >= 0.0 else 2.0**combined_log
            term_prob = (j - 1) * p_term * early_tail if late_termination else 0.0
            final_prob = min(1.0, term_prob + combined_prob)
            total = (
                log2_binom(final_suffix, j) + math.log2(final_prob)
                if final_prob > 0.0
                else float("-inf")
            )
        a = final_suffix + 1
        while a <= l_window:
            b = min(l_window, max(a, int(math.floor(a * block_ratio))))
            if budget_mode == "separate":
                tail_log = binom_cdf_half_entropy_log2(a, d)
                tail_prob = 1.0 if tail_log >= 0.0 else 2.0**tail_log
                late_prob = min(1.0, (j - 1) * p_term + tail_prob)
            else:
                combined_log = binom_cdf_half_entropy_log2(early_fair_len + a, d)
                combined_prob = 1.0 if combined_log >= 0.0 else 2.0**combined_log
                term_prob = (j - 1) * p_term * early_tail if late_termination else 0.0
                late_prob = min(1.0, term_prob + combined_prob)
            if late_prob > 0.0:
                term = (
                    math.log2(b - a + 1)
                    + log2_binom(b - 1, j - 1)
                    + math.log2(late_prob)
                )
                total = log2add(total, term)
            a = b + 1
        logs[j] = total
    return logs


def corrected_inner_logs(
    *,
    n: int,
    sigma: int,
    delta: float,
    h_max: int,
    small_xi: float,
    t_factor: float,
    late_block_ratio: float,
    budget_mode: str,
    late_termination: bool,
) -> list[float]:
    d = math.floor(delta * n)
    t = math.ceil(t_factor * d)
    l_window = math.ceil((2.0 + small_xi) * d)
    prefix = n - l_window
    if t > prefix:
        raise ValueError("T exceeds early prefix")
    late_logs = late_sum_blocked_logs(
        h_max=h_max,
        sigma=sigma,
        d=d,
        l_window=l_window,
        early_fair_len=max(0, t - sigma),
        block_ratio=late_block_ratio,
        budget_mode=budget_mode,
        late_termination=late_termination,
    )
    early_cost = -(sigma - 1)
    inner = [float("-inf")] * (h_max + 1)
    early_start_log = log2_binom(prefix - t + 1, 1)
    for h in range(2, h_max + 1):
        denom = log2_binom(n, h)
        total = float("-inf")
        for j in range(0, h - 1):
            i = h - 2 - j
            if 0 <= i <= t - 2:
                term = early_start_log + log2_binom(t - 2, i) - denom + early_cost + late_logs[j]
                total = log2add(total, term)
        inner[h] = min(0.0, total)
    return inner


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--k", type=int, required=True)
    parser.add_argument("--sigmas", required=True)
    parser.add_argument("--delta", type=float, default=0.106)
    parser.add_argument("--h-min", type=int, default=2)
    parser.add_argument("--h-max", type=int, default=2000)
    parser.add_argument("--small-xi", type=float, default=0.1)
    parser.add_argument("--t-factor", type=float, default=2.0)
    parser.add_argument("--z-min", type=float, default=0.005)
    parser.add_argument("--z-max", type=float, default=0.999)
    parser.add_argument("--z-count", type=int, default=160)
    parser.add_argument("--late-block-ratio", type=float, default=1.002)
    parser.add_argument(
        "--budget-mode",
        choices=("paired", "separate"),
        default="paired",
        help="'paired' shares the distance budget across early+late episodes; 'separate' is the older loose diagnostic.",
    )
    parser.add_argument(
        "--include-crude-late-termination",
        action="store_true",
        help="Include the older crude late-termination union term. This is intentionally off by default.",
    )
    parser.add_argument("--exact-outer-csv", default=None)
    parser.add_argument("--exact-through", type=int, default=80)
    parser.add_argument("--out-prefix", default=None)
    parser.add_argument("--no-png", action="store_true")
    args = parser.parse_args()

    n = 2 * args.k
    sigmas = parse_sigmas(args.sigmas)
    zs = z_grid(args.z_min, args.z_max, args.z_count)
    if args.out_prefix is None:
        sig_text = "_".join(str(s) for s in sigmas)
        out_base = Path(__file__).resolve().parent / (
            f"corrected_crossover_k{args.k}_sig{sig_text}_d{args.delta:g}_h{args.h_min}_{args.h_max}"
        )
    else:
        out_base = Path(args.out_prefix)
        if not out_base.is_absolute():
            out_base = Path.cwd() / out_base

    rows: list[dict[str, str | int | float]] = []
    summaries: list[dict[str, str | int | float]] = []
    series: dict[int, list[tuple[int, float, float]]] = {}
    for sigma in sigmas:
        print(f"sigma={sigma}: outer GF grid and corrected inner h={args.h_min}..{args.h_max}", flush=True)
        gf_outer, gf_z = outer_gf_bounds(k=args.k, sigma=sigma, h_min=args.h_min, h_max=args.h_max, zs=zs)
        exact_outer = load_exact_outer(args.exact_outer_csv, sigma, args.h_max)
        inner = corrected_inner_logs(
            n=n,
            sigma=sigma,
            delta=args.delta,
            h_max=args.h_max,
            small_xi=args.small_xi,
            t_factor=args.t_factor,
            late_block_ratio=args.late_block_ratio,
            budget_mode=args.budget_mode,
            late_termination=args.include_crude_late_termination,
        )
        total = float("-inf")
        peak_h = -1
        peak_term = float("-inf")
        points: list[tuple[int, float, float]] = []
        for h in range(args.h_min, args.h_max + 1):
            use_exact = h <= args.exact_through and exact_outer[h] != float("-inf")
            outer = exact_outer[h] if use_exact else gf_outer[h]
            source = "exact" if use_exact else "gf"
            term = outer + inner[h] if outer != float("-inf") and inner[h] != float("-inf") else float("-inf")
            total = log2add(total, term)
            if term > peak_term:
                peak_term = term
                peak_h = h
            points.append((h, term, total))
            rows.append(
                {
                    "sigma": sigma,
                    "offset": sigma - math.ceil(math.log2(args.k)),
                    "h": h,
                    "outer_source": source,
                    "outer_log2": format_log2(outer),
                    "outer_z": "" if use_exact else f"{gf_z[h]:.12g}",
                    "inner_log2": format_log2(inner[h]),
                    "term_log2": format_log2(term),
                    "cum_log2": format_log2(total),
                }
            )
            if h == args.h_min or h % 100 == 0 or h == args.h_max:
                print(f"  h={h:5d} term={format_log2(term)} cum={format_log2(total)}", flush=True)
        boundary = "+" if peak_h == args.h_max else ""
        print(
            f"  total={format_log2(total)} peak_h={peak_h}{boundary} peak_term={format_log2(peak_term)}",
            flush=True,
        )
        summaries.append(
            {
                "sigma": sigma,
                "offset": sigma - math.ceil(math.log2(args.k)),
                "h_min": args.h_min,
                "h_max": args.h_max,
                "total_log2": format_log2(total),
                "peak_h": peak_h,
                "peak_h_label": f"{peak_h}+" if peak_h == args.h_max else str(peak_h),
                "peak_term_log2": format_log2(peak_term),
                "outer_after_exact": "FIXED_TAP_BANDED_GF_GRID",
                "inner_model": "CORRECTED_R1_T_FIXED_BLOCKED_LATE_TAIL",
                "budget_mode": args.budget_mode,
                "late_termination": "crude_union" if args.include_crude_late_termination else "off",
            }
        )
        series[sigma] = points

    csv_path = out_base.with_suffix(".csv")
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "sigma",
                "offset",
                "h",
                "outer_source",
                "outer_log2",
                "outer_z",
                "inner_log2",
                "term_log2",
                "cum_log2",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)
    summary_path = out_base.with_name(out_base.name + "_summary").with_suffix(".csv")
    with summary_path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "sigma",
                "offset",
                "h_min",
                "h_max",
                "total_log2",
                "peak_h",
                "peak_h_label",
                "peak_term_log2",
                "outer_after_exact",
                "inner_model",
                "budget_mode",
                "late_termination",
            ],
        )
        writer.writeheader()
        writer.writerows(summaries)
    print(f"Wrote {csv_path}")
    print(f"Wrote {summary_path}")

    if args.no_png:
        return
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is unavailable; skipped PNG")
        return

    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True, constrained_layout=True)
    for sigma, points in series.items():
        hs = [p[0] for p in points]
        terms = [p[1] for p in points]
        cums = [p[2] for p in points]
        axes[0].plot(hs, terms, linewidth=1.3, label=f"sigma={sigma}")
        axes[1].plot(hs, cums, linewidth=1.3, label=f"sigma={sigma}")
        finite = [(h, t) for h, t, _ in points if t != float("-inf")]
        if finite:
            peak_h, peak_t = max(finite, key=lambda x: x[1])
            axes[0].scatter([peak_h], [peak_t], s=36)
            axes[0].annotate(
                f"h={peak_h}{'+' if peak_h == args.h_max else ''}",
                (peak_h, peak_t),
                xytext=(6, 8),
                textcoords="offset points",
                fontsize=8,
            )
    axes[0].axhline(0.0, color="black", linewidth=0.8, linestyle="--")
    axes[1].axhline(0.0, color="black", linewidth=0.8, linestyle="--")
    axes[0].set_ylabel("log2 contribution")
    axes[1].set_ylabel("log2 cumulative")
    axes[1].set_xlabel("outer weight h")
    axes[0].set_title(
        f"Corrected r=1 crossover scan, k={args.k}, delta={args.delta}, T={args.t_factor}d"
    )
    axes[0].legend()
    axes[1].legend()
    png_path = out_base.with_suffix(".png")
    fig.savefig(png_path, dpi=180)
    print(f"Wrote {png_path}")


if __name__ == "__main__":
    main()
