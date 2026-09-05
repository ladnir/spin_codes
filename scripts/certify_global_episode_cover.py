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

import numpy as np

from dense_largek_eval import (
    binom_cdf_half_entropy_log2,
    fixedtap_banded_outer_gf_log2,
    log2_binom,
    log2add,
    outer_small_h_banded_systematic_prefix,
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


def exact_fixedtap_outer_prefix(k_msg: int, parity_n: int, sigma: int, h_max: int) -> list[float]:
    vals = [float("-inf")] * (h_max + 1)
    if h_max <= 0:
        return vals
    prefix = outer_small_h_banded_systematic_prefix(k_msg, sigma, parity_n, h_max, fixed_tap=True)
    for h in range(1, h_max + 1):
        vals[h] = prefix[h]
    return vals


def outer_gf_bounds(
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
    log_ws = [(z, fixedtap_banded_outer_gf_log2(k, parity_n, sigma, z)) for z in zs]
    for h in range(h_min, h_max + 1):
        best = float("inf")
        best_here = float("nan")
        for z, log_w in log_ws:
            if log_w == float("-inf"):
                continue
            val = log_w - h * math.log2(z)
            if val < best:
                best = val
                best_here = z
        logs[h] = best if best != float("inf") else float("-inf")
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


def survivor_only_exact_prefix_log2(*, n: int, h_max: int, tail_cache: list[float]) -> list[float]:
    """Vectorized exact suffix-start sums for h=1..h_max.

    For fixed h, the suffix length U has law

        P[U=u] = C(u-1,h-1) / C(N,h).

    The recurrence from h to h+1 avoids the scalar log-sum loop for each h.
    """

    vals = [float("-inf")] * (h_max + 1)
    if h_max <= 0:
        return vals

    tails = np.exp2(np.asarray(tail_cache, dtype=np.float64))
    u = np.arange(n + 1, dtype=np.float64)
    probs = np.zeros(n + 1, dtype=np.float64)
    probs[1:] = 1.0 / n

    for h in range(1, h_max + 1):
        total = float(np.dot(probs[h:], tails[h:]))
        vals[h] = min(0.0, math.log2(total)) if total > 0.0 else float("-inf")
        if h < h_max:
            probs *= ((u - h) / h) * ((h + 1) / (n - h))
            probs[: h + 1] = 0.0

    return vals


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
    exact_survivor_value: float | None = None,
) -> tuple[float, int, float, list[tuple[int, float]]]:
    d = math.floor(delta * n)
    if exact_survivor:
        if exact_survivor_value is not None:
            total = exact_survivor_value
        elif tail_cache is None:
            raise ValueError("exact survivor mode requires a tail cache")
        else:
            total = survivor_only_exact_log2(n=n, h=h, tail_cache=tail_cache)
    else:
        total = survivor_only_log2(n=n, h=h, d=d, block_ratio=block_ratio, tail_mode=tail_mode)
    best_r = 0
    best_piece = total
    pieces = [(0, total)]
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
        pieces.append((r, piece))
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
    return min(0.0, total), best_r, best_piece, pieces


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--k", type=int, required=True)
    parser.add_argument(
        "--parity-extra",
        type=int,
        default=0,
        help="Extra terminated parity coordinates beyond k. The interleaver length is k+(k+parity_extra).",
    )
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
    parser.add_argument(
        "--exact-outer-mode",
        choices=("fixedtap-banded", "csv-only", "none"),
        default="fixedtap-banded",
        help="Source for exact low-weight outer counts before falling back to the GF envelope.",
    )
    parser.add_argument("--out-prefix", default=None)
    parser.add_argument(
        "--pieces-out",
        default=None,
        help="Optional long-form CSV of inner r-piece contributions for proof diagnostics.",
    )
    parser.add_argument("--no-png", action="store_true")
    args = parser.parse_args()

    parity_n = args.k + args.parity_extra
    n = args.k + parity_n
    d = math.floor(args.delta * n)
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
    piece_rows: list[dict[str, str | int | float]] = []
    summaries: list[dict[str, str | int | float]] = []
    series: dict[int, list[tuple[int, float, float]]] = {}
    tail_cache: list[float] | None = None
    exact_survivor_values: list[float] | None = None
    if args.exact_survivor_through >= args.h_min:
        print(f"precomputing {args.tail_mode} tail cache for exact survivor sums", flush=True)
        tail_cache = [tail_log2(u, d, args.tail_mode) for u in range(n + 1)]
        exact_hi = min(args.exact_survivor_through, args.h_max)
        print(f"precomputing exact survivor sums through h={exact_hi}", flush=True)
        exact_survivor_values = survivor_only_exact_prefix_log2(n=n, h_max=exact_hi, tail_cache=tail_cache)
    for sigma in sigmas:
        print(f"sigma={sigma}: global episode cover h={args.h_min}..{args.h_max}", flush=True)
        gf_outer, gf_z = outer_gf_bounds(
            k=args.k,
            parity_n=parity_n,
            sigma=sigma,
            h_min=args.h_min,
            h_max=args.h_max,
            zs=zs,
        )
        exact_outer = load_exact_outer(args.exact_outer_csv, sigma, args.h_max)
        if args.exact_outer_mode == "fixedtap-banded":
            exact_hi = min(args.exact_through, args.h_max)
            formula_outer = exact_fixedtap_outer_prefix(args.k, parity_n, sigma, exact_hi)
            for h in range(1, exact_hi + 1):
                exact_outer[h] = formula_outer[h]
        elif args.exact_outer_mode == "none":
            exact_outer = [float("-inf")] * (args.h_max + 1)
        total = float("-inf")
        peak_h = -1
        peak_term = float("-inf")
        points: list[tuple[int, float, float]] = []
        for h in range(args.h_min, args.h_max + 1):
            use_exact = h <= args.exact_through and args.exact_outer_mode == "fixedtap-banded"
            if not use_exact:
                use_exact = h <= args.exact_through and exact_outer[h] != float("-inf")
            outer = exact_outer[h] if use_exact else gf_outer[h]
            source = "exact" if use_exact else "gf"
            inner, best_r, best_piece, pieces = global_episode_inner_log2(
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
                exact_survivor_value=(
                    exact_survivor_values[h]
                    if exact_survivor_values is not None and h < len(exact_survivor_values)
                    else None
                ),
            )
            term = outer + inner if outer != float("-inf") and inner != float("-inf") else float("-inf")
            if args.pieces_out is not None:
                for r, piece in pieces:
                    piece_rows.append(
                        {
                            "sigma": sigma,
                            "h": h,
                            "r": r,
                            "outer_log2": format_log2(outer),
                            "inner_piece_log2": format_log2(piece),
                            "term_piece_log2": format_log2(
                                outer + piece if outer != float("-inf") and piece != float("-inf") else float("-inf")
                            ),
                        }
                    )
            total = log2add(total, term)
            if term > peak_term:
                peak_term = term
                peak_h = h
            points.append((h, term, total))
            rows.append(
                {
                    "k": args.k,
                    "N": n,
                    "parity_n": parity_n,
                    "parity_extra": args.parity_extra,
                    "delta": f"{args.delta:.12g}",
                    "d": d,
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
                "k": args.k,
                "N": n,
                "parity_n": parity_n,
                "parity_extra": args.parity_extra,
                "delta": f"{args.delta:.12g}",
                "d": d,
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
                "exact_outer_mode": args.exact_outer_mode,
                "exact_through": args.exact_through,
                "z_min": f"{args.z_min:.12g}",
                "z_max": f"{args.z_max:.12g}",
                "z_count": args.z_count,
                "block_ratio": f"{args.block_ratio:.12g}",
                "suffix_block_ratio": f"{suffix_block_ratio:.12g}",
                "stop_gap_bits": f"{args.stop_gap_bits:.12g}",
                "tail_confirm": args.tail_confirm,
                "r_max": "" if args.r_max is None else args.r_max,
            }
        )

    csv_path = out_base.with_suffix(".csv")
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "k",
                "N",
                "parity_n",
                "parity_extra",
                "delta",
                "d",
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
                "k",
                "N",
                "parity_n",
                "parity_extra",
                "delta",
                "d",
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
                "exact_outer_mode",
                "exact_through",
                "z_min",
                "z_max",
                "z_count",
                "block_ratio",
                "suffix_block_ratio",
                "stop_gap_bits",
                "tail_confirm",
                "r_max",
            ],
        )
        writer.writeheader()
        writer.writerows(summaries)
    print(f"Wrote {csv_path}")
    print(f"Wrote {summary_path}")
    if args.pieces_out is not None:
        pieces_path = Path(args.pieces_out)
        if not pieces_path.is_absolute():
            pieces_path = Path.cwd() / pieces_path
        pieces_path.parent.mkdir(parents=True, exist_ok=True)
        with pieces_path.open("w", newline="") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=["sigma", "h", "r", "outer_log2", "inner_piece_log2", "term_piece_log2"],
            )
            writer.writeheader()
            writer.writerows(piece_rows)
        print(f"Wrote {pieces_path}")

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
