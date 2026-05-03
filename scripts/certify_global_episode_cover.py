"""Global fixed-tap episode-cover first-moment certificate.

This script is a theorem-shaped finite-n upper bound for the dense+dense
medium-weight bridge.  It keeps the interleaver visible:

    sum_h A_h^out p_h^in(delta).

The inner bound covers every input-one position by either terminated ON
intervals, or by one final surviving suffix.  All live intervals share the same
distance budget through one binomial lower tail.  This avoids the two loose
diagnostic mistakes that made the medium bridge look dangerous:

* late support is not a free bucket after an early termination;
* separated live episodes do not get separate budgets of d output ones.

The interval count is intentionally pessimistic but keeps the final survivor's
geometry: r terminated intervals may occur before one final surviving suffix.
The old diagnostic treated that survivor as an ordinary interval anywhere in
[N], which was safe but too loose in the tiny-weight crossover range.
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


def binom_cdf_half_be_log2(n: int, k: int) -> float:
    """Berry-Esseen upper bound for P[Bin(n,1/2) <= k].

    This is useful near and above the mean, where the entropy bound returns 1.
    The constant is deliberately conservative.  For k below the mean we take
    the minimum of this bound and the entropy bound used elsewhere.
    """
    if k < 0:
        return float("-inf")
    if k >= n:
        return 0.0
    if n <= 0:
        return 0.0
    entropy = binom_cdf_half_entropy_log2(n, k)
    mu = 0.5 * n
    sigma = 0.5 * math.sqrt(n)
    z = (k + 0.5 - mu) / sigma
    normal = 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))
    # Conservative Berry-Esseen allowance for Bernoulli(1/2).
    be = min(1.0, max(0.0, normal + 0.8 / math.sqrt(n)))
    be_log = math.log2(be) if be > 0.0 else float("-inf")
    return min(entropy, be_log)


def tail_log2(n: int, k: int, mode: str) -> float:
    if mode == "entropy":
        return binom_cdf_half_entropy_log2(n, k)
    if mode == "be":
        return binom_cdf_half_be_log2(n, k)
    raise ValueError(f"unknown tail mode: {mode}")


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
                raise ValueError(f"bad sigma range: {part!r}")
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
    log_ws = [(z, fixedtap_banded_outer_gf_log2(k, k, sigma, z)) for z in zs]
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


def block_ranges(lo: int, hi: int, ratio: float) -> list[tuple[int, int]]:
    if lo > hi:
        return []
    if ratio <= 1.0:
        raise ValueError("block ratio must be > 1")
    out: list[tuple[int, int]] = []
    a = lo
    while a <= hi:
        b = min(hi, max(a, int(math.floor(a * ratio))))
        out.append((a, b))
        a = b + 1
    return out


def survivor_only_log2(*, n: int, h: int, d: int, block_ratio: float, tail_mode: str) -> float:
    """One final surviving suffix and no prior terminated intervals."""
    denom = log2_binom(n, h)
    total = float("-inf")
    # For suffix length U in [a,b], C(U-1,h-1) and the binomial lower-tail are
    # both upper-bounded at the most favorable endpoint for safety.
    for a, b in block_ranges(h, n, block_ratio):
        tail = tail_log2(a, d, tail_mode)
        term = math.log2(b - a + 1) + log2_binom(b - 1, h - 1) - denom + tail
        total = log2add(total, term)
    return min(0.0, total)


def survivor_only_exact_log2(*, n: int, h: int, tail_cache: list[float]) -> float:
    """Exact suffix-start count for one final survivor and no terminated intervals.

    The binomial lower tail may still be an upper bound, depending on the
    tail-cache construction, but the suffix placement sum is exact:

        sum_U C(U-1,h-1)/C(N,h) * tail(U).

    This is intended for small h, where the current certificate peak lives.
    """
    denom = log2_binom(n, h)
    total = float("-inf")
    log_count = 0.0  # C(h-1,h-1)
    for u in range(h, n + 1):
        total = log2add(total, log_count - denom + tail_cache[u])
        if u < n:
            log_count += math.log2(u / (u - h + 1))
    return min(0.0, total)


def terminated_only_log2(
    *,
    n: int,
    h: int,
    sigma: int,
    d: int,
    r: int,
    block_ratio: float,
    tail_mode: str,
) -> float:
    """All support is covered by r terminated intervals."""
    if h < 2 * r:
        return float("-inf")
    denom = log2_binom(n, h)
    total = float("-inf")
    u_lo = max(2 * r, h)
    u_hi = n
    for a, b in block_ranges(u_lo, u_hi, block_ratio):
        charged = max(0, a - r * sigma)
        tail = tail_log2(charged, d, tail_mode)
        # Decreasing placement at a, increasing length/support factors at b.
        term = (
            math.log2(b - a + 1)
            + log2_binom(b - r - 1, r - 1)
            + log2_binom(n - a + r, r)
            + log2_binom(b - 2 * r, h - 2 * r)
            - denom
            - r * (sigma - 1)
            + tail
        )
        total = log2add(total, term)
    return min(0.0, total)


def terminated_plus_survivor_log2(
    *,
    n: int,
    h: int,
    sigma: int,
    d: int,
    r: int,
    block_ratio: float,
    tail_mode: str,
) -> float:
    """r terminated intervals plus one final surviving interval.

    The survivor is counted as an ordinary ordered interval for placement,
    which is a safe overcount.  Total live length is U.
    """
    if h < 2 * r + 1:
        return float("-inf")
    denom = log2_binom(n, h)
    total = float("-inf")
    u_lo = max(2 * r + 1, h)
    u_hi = n
    for a, b in block_ranges(u_lo, u_hi, block_ratio):
        charged = max(0, a - r * sigma)
        tail = tail_log2(charged, d, tail_mode)
        term = (
            math.log2(b - a + 1)
            + log2_binom(b - r - 1, r)
            + log2_binom(n - a + r + 1, r + 1)
            + log2_binom(b - 2 * r - 1, h - 2 * r - 1)
            - denom
            - r * (sigma - 1)
            + tail
        )
        total = log2add(total, term)
    return min(0.0, total)


def terminated_plus_suffix_survivor_log2(
    *,
    n: int,
    h: int,
    sigma: int,
    d: int,
    r: int,
    block_ratio: float,
    tail_mode: str,
) -> float:
    """r terminated intervals before one final surviving suffix.

    The survivor has suffix length S.  The r terminated intervals have total
    live length T and lie in the prefix of length N-S.  Remaining input-one
    positions lie inside the terminated interval interiors or in the surviving
    suffix after its start.  The shared budget charges the survivor length S
    plus the nonterminal live outputs of the terminated intervals.
    """
    if h < 2 * r + 1:
        return float("-inf")
    denom = log2_binom(n, h)
    total = float("-inf")
    s_min = 1
    s_max = n - 2 * r
    for sa, sb in block_ranges(s_min, s_max, block_ratio):
        t_min = 2 * r
        t_max = n - sb
        if t_min > t_max:
            continue
        for ta, tb in block_ranges(t_min, t_max, block_ratio):
            if n - sa - ta + r < r:
                continue
            slots = sb + tb - 2 * r - 1
            need = h - 2 * r - 1
            if slots < need:
                continue
            charged = sa + max(0, ta - r * sigma)
            tail = tail_log2(charged, d, tail_mode)
            term = (
                math.log2(sb - sa + 1)
                + math.log2(tb - ta + 1)
                + log2_binom(tb - r - 1, r - 1)
                + log2_binom(n - sa - ta + r, r)
                + log2_binom(slots, need)
                - denom
                - r * (sigma - 1)
                + tail
            )
            total = log2add(total, term)
    return min(0.0, total)


def global_episode_inner_log2(
    *,
    n: int,
    h: int,
    sigma: int,
    delta: float,
    r_max: int | None,
    block_ratio: float,
    suffix_block_ratio: float,
    stop_gap_bits: float,
    tail_confirm: int,
    suffix_survivor: bool,
    tail_mode: str,
    exact_survivor: bool,
    tail_cache: list[float] | None,
) -> tuple[float, int, float]:
    d = math.floor(delta * n)
    if exact_survivor:
        if tail_cache is None:
            raise ValueError("exact survivor mode requires a tail cache")
        total = survivor_only_exact_log2(n=n, h=h, tail_cache=tail_cache)
    else:
        total = survivor_only_log2(n=n, h=h, d=d, block_ratio=block_ratio, tail_mode=tail_mode)
    best_r = 0
    best_piece = total
    max_r = h // 2 if r_max is None else min(r_max, h // 2)
    below_count = 0
    for r in range(1, max_r + 1):
        piece = terminated_only_log2(n=n, h=h, sigma=sigma, d=d, r=r, block_ratio=block_ratio, tail_mode=tail_mode)
        survivor_piece = (
            terminated_plus_suffix_survivor_log2(
                n=n,
                h=h,
                sigma=sigma,
                d=d,
                r=r,
                block_ratio=suffix_block_ratio,
                tail_mode=tail_mode,
            )
            if suffix_survivor
            else terminated_plus_survivor_log2(
                n=n, h=h, sigma=sigma, d=d, r=r, block_ratio=block_ratio, tail_mode=tail_mode
            )
        )
        piece = log2add(
            piece,
            survivor_piece,
        )
        total = log2add(total, piece)
        if piece > best_piece:
            best_piece = piece
            best_r = r
            below_count = 0
        elif best_piece != float("-inf") and piece <= best_piece - stop_gap_bits:
            below_count += 1
            if below_count >= tail_confirm:
                break
        else:
            below_count = 0
    return min(0.0, total), best_r, best_piece


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--k", type=int, required=True)
    parser.add_argument("--sigmas", required=True)
    parser.add_argument("--delta", type=float, default=0.106)
    parser.add_argument("--h-min", type=int, default=1)
    parser.add_argument("--h-max", type=int, default=2000)
    parser.add_argument("--z-min", type=float, default=0.005)
    parser.add_argument("--z-max", type=float, default=0.999)
    parser.add_argument("--z-count", type=int, default=160)
    parser.add_argument("--block-ratio", type=float, default=1.01)
    parser.add_argument(
        "--suffix-block-ratio",
        type=float,
        default=None,
        help="Optional coarser block ratio for the expensive suffix-survivor double sum.",
    )
    parser.add_argument("--r-max", type=int, default=None)
    parser.add_argument("--stop-gap-bits", type=float, default=60.0)
    parser.add_argument("--tail-confirm", type=int, default=8)
    parser.add_argument("--tail-mode", choices=("entropy", "be"), default="entropy")
    parser.add_argument(
        "--exact-survivor-through",
        type=int,
        default=0,
        help="Use the exact suffix-start survivor-only sum through this h.",
    )
    parser.add_argument(
        "--ordinary-survivor-interval",
        action="store_true",
        help="Use the older safe overcount that treats the final survivor as an arbitrary interval.",
    )
    parser.add_argument("--exact-outer-csv", default=None)
    parser.add_argument("--exact-through", type=int, default=80)
    parser.add_argument("--out-prefix", default=None)
    parser.add_argument("--no-png", action="store_true")
    args = parser.parse_args()

    n = 2 * args.k
    suffix_block_ratio = args.suffix_block_ratio if args.suffix_block_ratio is not None else args.block_ratio
    sigmas = parse_sigmas(args.sigmas)
    zs = z_grid(args.z_min, args.z_max, args.z_count)
    if args.out_prefix is None:
        sig_text = "_".join(str(s) for s in sigmas)
        out_base = Path(__file__).resolve().parent / (
            f"global_episode_cover_k{args.k}_sig{sig_text}_d{args.delta:g}_h{args.h_min}_{args.h_max}"
        )
    else:
        out_base = Path(args.out_prefix)
        if not out_base.is_absolute():
            out_base = Path.cwd() / out_base

    rows: list[dict[str, str | int | float]] = []
    summaries: list[dict[str, str | int | float]] = []
    series: dict[int, list[tuple[int, float, float]]] = {}
    tail_cache: list[float] | None = None
    if args.exact_survivor_through >= args.h_min:
        d = math.floor(args.delta * n)
        print(f"precomputing {args.tail_mode} tail cache for exact survivor sums", flush=True)
        tail_cache = [tail_log2(u, d, args.tail_mode) for u in range(n + 1)]
    for sigma in sigmas:
        print(f"sigma={sigma}: global episode cover h={args.h_min}..{args.h_max}", flush=True)
        gf_outer, gf_z = outer_gf_bounds(k=args.k, sigma=sigma, h_min=args.h_min, h_max=args.h_max, zs=zs)
        exact_outer = load_exact_outer(args.exact_outer_csv, sigma, args.h_max)
        total = float("-inf")
        peak_h = -1
        peak_term = float("-inf")
        points: list[tuple[int, float, float]] = []
        for h in range(args.h_min, args.h_max + 1):
            use_exact = h <= args.exact_through and exact_outer[h] != float("-inf")
            outer = exact_outer[h] if use_exact else gf_outer[h]
            source = "exact" if use_exact else "gf"
            inner, best_r, best_piece = global_episode_inner_log2(
                n=n,
                h=h,
                sigma=sigma,
                delta=args.delta,
                r_max=args.r_max,
                block_ratio=args.block_ratio,
                suffix_block_ratio=suffix_block_ratio,
                stop_gap_bits=args.stop_gap_bits,
                tail_confirm=args.tail_confirm,
                suffix_survivor=not args.ordinary_survivor_interval,
                tail_mode=args.tail_mode,
                exact_survivor=h <= args.exact_survivor_through,
                tail_cache=tail_cache,
            )
            term = outer + inner if outer != float("-inf") and inner != float("-inf") else float("-inf")
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
                    "inner_log2": format_log2(inner),
                    "term_log2": format_log2(term),
                    "cum_log2": format_log2(total),
                    "best_r": best_r,
                    "best_piece_log2": format_log2(best_piece),
                    "survivor_only_sum": "exact" if h <= args.exact_survivor_through else "blocked",
                }
            )
            if h == args.h_min or h % 25 == 0 or h == args.h_max:
                print(
                    f"  h={h:5d} term={format_log2(term)} cum={format_log2(total)} best_r={best_r}",
                    flush=True,
                )
        series[sigma] = points
        print(
            f"  total={format_log2(total)} peak_h={peak_h}{'+' if peak_h == args.h_max else ''} "
            f"peak_term={format_log2(peak_term)}",
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
                "inner_model": "GLOBAL_EPISODE_COVER_SHARED_BUDGET",
                "survivor_count": "ORDINARY_INTERVAL" if args.ordinary_survivor_interval else "FINAL_SUFFIX",
                "tail_mode": args.tail_mode,
                "exact_survivor_through": args.exact_survivor_through,
            }
        )

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
                "best_r",
                "best_piece_log2",
                "survivor_only_sum",
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
                "survivor_count",
                "tail_mode",
                "exact_survivor_through",
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
    axes[0].axhline(0.0, color="black", linewidth=0.8, linestyle="--")
    axes[1].axhline(0.0, color="black", linewidth=0.8, linestyle="--")
    axes[0].set_ylabel("log2 contribution")
    axes[1].set_ylabel("log2 cumulative")
    axes[1].set_xlabel("outer weight h")
    axes[0].set_title(f"Global episode-cover certificate, k={args.k}, delta={args.delta}")
    axes[0].legend()
    axes[1].legend()
    png_path = out_base.with_suffix(".png")
    fig.savefig(png_path, dpi=180)
    print(f"Wrote {png_path}")


if __name__ == "__main__":
    main()
