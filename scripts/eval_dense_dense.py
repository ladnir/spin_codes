#!/usr/bin/env python3
"""Evaluate the dense+dense objective from the current LaTeX manuscript.

This script mirrors the finite-n formulas in:
  - innerDense.tex
  - outerDense.tex
  - integration.tex

It is intentionally self-contained and uses only the Python standard library.
"""

from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass
from typing import Iterable, List, Optional, Sequence, Tuple


LN2 = math.log(2.0)
LOG10 = math.log(10.0)


def logsumexp(log_values: Sequence[float]) -> float:
    if not log_values:
        return float("-inf")
    m = max(log_values)
    if math.isinf(m):
        return m
    total = 0.0
    for x in log_values:
        total += math.exp(x - m)
    return m + math.log(total)


def log_binom(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1.0) - math.lgamma(k + 1.0) - math.lgamma(n - k + 1.0)


def log_binom_pmf_half(n: int, k: int) -> float:
    return log_binom(n, k) - n * LN2


def safe_prob_from_log(logp: float) -> float:
    if logp == float("-inf"):
        return 0.0
    if logp > 0.0:
        return 1.0
    if logp < -745.0:
        return 0.0
    return math.exp(logp)


def format_from_log(logv: float) -> str:
    if logv == float("-inf"):
        return "0.000000e+00 (log10=-inf)"
    if logv > 700.0:
        return f"exp({logv:.6g}) (log10={fmt_log10(logv)})"
    if logv < -745.0:
        return f"0.000000e+00 (log10={fmt_log10(logv)})"
    return f"{math.exp(logv):.6e} (log10={fmt_log10(logv)})"


def fmt_log10(logp: float) -> str:
    if logp == float("-inf"):
        return "-inf"
    return f"{logp / LOG10:.6g}"


def parse_int_list(value: Optional[str]) -> Optional[List[int]]:
    if value is None:
        return None
    parts = [p.strip() for p in value.split(",") if p.strip()]
    if not parts:
        return None
    return [int(p) for p in parts]


def binom_cdf_half_upper(n: int, k: int, *, exact_limit: int = 2000, tail_terms: int = 128) -> float:
    """Upper bound / exact value for P[Binom(n,1/2) <= k].

    For moderate k this is computed by an exact recurrence. For large k we
    compute a short prefix of the tail exactly and then finish with a geometric
    remainder bound. This is a numerically stable upper bound for the lower tail
    and is sufficient for evaluating the paper's finite-length envelopes.
    """

    if k < 0:
        return 0.0
    if k >= n:
        return 1.0
    if k > n // 2:
        # Use symmetry to reduce to the lower tail.
        return 1.0 - binom_cdf_half_upper(n, n - k - 1, exact_limit=exact_limit, tail_terms=tail_terms)

    if k <= exact_limit:
        terms = [log_binom_pmf_half(n, i) for i in range(k + 1)]
        return safe_prob_from_log(logsumexp(terms))

    # Start from the endpoint and walk backwards for a bounded number of terms.
    logs = [log_binom_pmf_half(n, k)]
    cur_log = logs[0]
    cur_i = k
    for _ in range(max(0, tail_terms - 1)):
        if cur_i == 0:
            break
        ratio_log = math.log(cur_i) - math.log(n - cur_i + 1)
        cur_log += ratio_log
        cur_i -= 1
        logs.append(cur_log)
        if ratio_log < -30.0:
            break

    partial_log = logsumexp(logs)
    partial = safe_prob_from_log(partial_log)

    # Geometric tail bound from the last ratio observed.
    last_ratio = cur_i / (n - cur_i + 1) if cur_i > 0 else 0.0
    if last_ratio >= 1.0:
        return min(1.0, partial)
    remainder = 0.0
    if cur_i > 0 and last_ratio > 0.0:
        remainder = safe_prob_from_log(logs[-1]) * (last_ratio / (1.0 - last_ratio))
    return min(1.0, partial + remainder)


def run_tail_upper(n: int, w: int, run_target: int) -> float:
    """Upper bound for P[R(U) < run_target] under the exact run-tail identity.

    For run_target > w the event is certain.
    """

    if run_target <= 1:
        return 0.0
    if run_target > w:
        return 1.0

    denom_log = log_binom(n, w)
    logs = []
    for s in range(1, run_target):
        log_num = log_binom(w - 1, s - 1) + log_binom(n - w + 1, s)
        logs.append(log_num - denom_log)
    return min(1.0, safe_prob_from_log(logsumexp(logs)))


@dataclass
class InnerEnvelope:
    global_log: float
    run_log: float
    best_log: float
    best_run_target: int
    run_tail_log: float
    off_log: float
    ld_log: float
    global_components: Tuple[float, float, float]
    run_components: Tuple[float, float, float]


def dense_inner_envelope(
    n: int,
    w: int,
    delta: float,
    epsilon: float,
    xi: float,
    m: int,
) -> InnerEnvelope:
    rho = 1.0 - epsilon
    L = int(math.ceil((2.0 + xi) * delta * n))
    q = binom_cdf_half_upper(int(math.floor((1.0 - epsilon) * n)), int(math.floor(delta * n)))
    b = binom_cdf_half_upper(L, int(math.floor(delta * n)))

    global_terms = [
        w * math.log(rho) if w > 0 else 0.0,
        math.log(n) - m * LN2,
        math.log(q) if q > 0.0 else float("-inf"),
    ]
    global_log = logsumexp(global_terms)

    best_run_log = float("inf")
    best_run_target = 1
    best_run_tail = 0.0
    best_off = 0.0
    best_ld = 0.0
    best_components = (float("-inf"), float("-inf"), float("-inf"))
    best_run_components = best_components

    for run_target in range(1, w + 1):
        run_tail = run_tail_upper(n, w, run_target)
        off = 0.0 if run_target == 0 else min(1.0, (L * (2.0 ** (-m))) ** run_target)
        ld = min(1.0, run_target * b)
        run_log = logsumexp(
            [
                math.log(run_tail) if run_tail > 0.0 else float("-inf"),
                math.log(off) if off > 0.0 else float("-inf"),
                math.log(ld) if ld > 0.0 else float("-inf"),
            ]
        )
        if run_log < best_run_log:
            best_run_log = run_log
            best_run_target = run_target
            best_run_tail = run_tail
            best_off = off
            best_ld = ld
            best_run_components = (
                math.log(run_tail) if run_tail > 0.0 else float("-inf"),
                math.log(off) if off > 0.0 else float("-inf"),
                math.log(ld) if ld > 0.0 else float("-inf"),
            )

    best_log = min(global_log, best_run_log)
    return InnerEnvelope(
        global_log=global_log,
        run_log=best_run_log,
        best_log=best_log,
        best_run_target=best_run_target,
        run_tail_log=math.log(best_run_tail) if best_run_tail > 0.0 else float("-inf"),
        off_log=math.log(best_off) if best_off > 0.0 else float("-inf"),
        ld_log=math.log(best_ld) if best_ld > 0.0 else float("-inf"),
        global_components=tuple(global_terms),
        run_components=best_run_components,
    )


@dataclass
class SpanContribution:
    ell: int
    log_term: float


@dataclass
class OuterSpectrumResult:
    log_aw: float
    terms: List[SpanContribution]
    exact: bool
    k: int
    streams: int
    effective_n: int


def log_span_count(k: int, ell: int) -> float:
    if ell == 1:
        return math.log(k)
    return math.log(k - ell + 1) + (ell - 2) * LN2


def outer_expected_spectrum(
    n: int,
    outer_streams: int,
    M: int,
    w: int,
    *,
    exact_span_limit: int = 4000,
    span_window: int = 2048,
    adapt_rounds: int = 4,
) -> OuterSpectrumResult:
    if outer_streams <= 0:
        raise ValueError("outer_streams must be positive")
    k = n // outer_streams
    effective_n = k * outer_streams
    exact = k <= exact_span_limit

    if exact:
        ell_lo, ell_hi = 1, k
    else:
        # Approximate search around the binomial mean contribution.
        ell_center = int(round((2.0 * w / outer_streams) - M))
        ell_center = max(1, min(k, ell_center))
        ell_lo = max(1, ell_center - span_window)
        ell_hi = min(k, ell_center + span_window)

    def terms_in_range(lo: int, hi: int) -> List[SpanContribution]:
        out: List[SpanContribution] = []
        for ell in range(lo, hi + 1):
            count_log = log_span_count(k, ell)
            m_bits = outer_streams * (ell + M)
            pmf_log = log_binom_pmf_half(m_bits, w)
            out.append(SpanContribution(ell=ell, log_term=count_log + pmf_log))
        return out

    terms = terms_in_range(ell_lo, ell_hi)
    if not terms:
        return OuterSpectrumResult(float("-inf"), [], exact, k, outer_streams, effective_n)

    for _ in range(adapt_rounds):
        best = max(terms, key=lambda t: t.log_term)
        left = terms[0]
        right = terms[-1]
        need_left = (left.log_term > best.log_term - 25.0) and ell_lo > 1
        need_right = (right.log_term > best.log_term - 25.0) and ell_hi < k
        if not (need_left or need_right):
            break
        if need_left:
            new_lo = max(1, ell_lo - span_window)
            if new_lo < ell_lo:
                ell_lo = new_lo
        if need_right:
            new_hi = min(k, ell_hi + span_window)
            if new_hi > ell_hi:
                ell_hi = new_hi
        terms = terms_in_range(ell_lo, ell_hi)

    log_aw = logsumexp([t.log_term for t in terms])
    return OuterSpectrumResult(log_aw=log_aw, terms=terms, exact=exact, k=k, streams=outer_streams, effective_n=effective_n)


def choose_weights(args: argparse.Namespace, n: int) -> List[int]:
    if args.w_values is not None:
        weights = args.w_values
    elif args.w_min is not None or args.w_max is not None:
        if args.w_min is None or args.w_max is None:
            raise ValueError("Both --w-min and --w-max are required for a sweep.")
        step = args.w_step
        if step <= 0:
            raise ValueError("--w-step must be positive")
        weights = list(range(args.w_min, args.w_max + 1, step))
    else:
        weights = [max(1, min(n, int(round(args.delta * n))))]
    return [w for w in weights if 0 <= w <= n]


def parse_memory_values(args: argparse.Namespace) -> List[int]:
    if args.m_values is not None:
        return args.m_values
    if args.m_min is not None or args.m_max is not None:
        if args.m_min is None or args.m_max is None:
            raise ValueError("Both --m-min and --m-max are required for a sweep.")
        step = args.m_step
        if step <= 0:
            raise ValueError("--m-step must be positive")
        return list(range(args.m_min, args.m_max + 1, step))
    return [args.m]


def parse_outer_memory_values(args: argparse.Namespace) -> List[int]:
    if args.M_values is not None:
        return args.M_values
    if args.M_min is not None or args.M_max is not None:
        if args.M_min is None or args.M_max is None:
            raise ValueError("Both --M-min and --M-max are required for a sweep.")
        step = args.M_step
        if step <= 0:
            raise ValueError("--M-step must be positive")
        return list(range(args.M_min, args.M_max + 1, step))
    return [args.M]


def choose_outer_streams(args: argparse.Namespace) -> Tuple[int, float]:
    if args.outer_streams is not None:
        streams = args.outer_streams
        if streams <= 0:
            raise ValueError("--outer-streams must be positive")
        return streams, 1.0 / streams

    if args.outer_rate <= 0.0 or args.outer_rate >= 1.0:
        raise ValueError("--outer-rate must lie in (0,1)")
    reciprocal = 1.0 / args.outer_rate
    streams = int(round(reciprocal))
    if streams <= 0:
        raise ValueError("Could not derive a positive stream count from --outer-rate")
    if not args.allow_rounding and abs(reciprocal - streams) > 1e-9:
        raise ValueError(
            "--outer-rate should be the reciprocal of an integer stream count for the current formulas; "
            "pass --allow-rounding to use the nearest integer stream count."
        )
    return streams, 1.0 / streams


def format_term(logp: float) -> str:
    return format_from_log(logp)


def print_top_spans(result: OuterSpectrumResult, top_k: int) -> None:
    top_terms = sorted(result.terms, key=lambda t: t.log_term, reverse=True)[:top_k]
    if not top_terms:
        print("  outer span contributors: none")
        return
    print("  top outer span contributors:")
    for item in top_terms:
        print(f"    ell={item.ell:>8}  term={format_term(item.log_term)}")


def evaluate_single_setting(
    n: int,
    outer_rate: float,
    outer_streams: int,
    delta: float,
    M: int,
    m: int,
    epsilon: float,
    xi: float,
    w: int,
    args: argparse.Namespace,
) -> Tuple[float, float, float]:
    inner = dense_inner_envelope(n=n, w=w, delta=delta, epsilon=epsilon, xi=xi, m=m)
    outer = outer_expected_spectrum(
        n=n,
        outer_streams=outer_streams,
        M=M,
        w=w,
        exact_span_limit=args.exact_span_limit,
        span_window=args.span_window,
        adapt_rounds=args.adapt_rounds,
    )
    log_term = outer.log_aw + inner.best_log
    return outer, inner, log_term


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Evaluate the current dense+dense objective from the LaTeX manuscript."
    )
    parser.add_argument("n", type=int, help="Blocklength parameter used in the current formulas.")
    parser.add_argument("--outer-rate", type=float, default=0.5, help="Outer rate R_o (default: 0.5).")
    parser.add_argument(
        "--outer-streams",
        type=int,
        default=None,
        help="Override the number of outer output streams r (usually 1/outer-rate).",
    )
    parser.add_argument("--allow-rounding", action="store_true", help="Allow rounding to the nearest stream count.")
    parser.add_argument("--delta", type=float, required=True, help="Target relative distance delta.")
    parser.add_argument("--M", type=int, required=True, help="Outer memory parameter.")
    parser.add_argument("--m", type=int, required=True, help="Inner memory parameter.")
    parser.add_argument("--epsilon", type=float, default=0.1, help="Dense-inner epsilon parameter.")
    parser.add_argument("--xi", type=float, default=0.1, help="Dense-inner xi parameter.")
    parser.add_argument("--w", type=int, default=None, help="Single weight to evaluate.")
    parser.add_argument("--w-min", type=int, default=None, help="Minimum weight for a sweep.")
    parser.add_argument("--w-max", type=int, default=None, help="Maximum weight for a sweep.")
    parser.add_argument("--w-step", type=int, default=1, help="Step size for a weight sweep.")
    parser.add_argument("--w-values", type=str, default=None, help="Comma-separated explicit weights to evaluate.")
    parser.add_argument("--m-min", type=int, default=None, help="Minimum memory for a sweep.")
    parser.add_argument("--m-max", type=int, default=None, help="Maximum memory for a sweep.")
    parser.add_argument("--m-step", type=int, default=1, help="Step size for a memory sweep.")
    parser.add_argument("--m-values", type=str, default=None, help="Comma-separated explicit memory values.")
    parser.add_argument("--M-min", type=int, default=None, help="Minimum outer memory for a sweep.")
    parser.add_argument("--M-max", type=int, default=None, help="Maximum outer memory for a sweep.")
    parser.add_argument("--M-step", type=int, default=1, help="Step size for an outer-memory sweep.")
    parser.add_argument("--M-values", type=str, default=None, help="Comma-separated explicit outer-memory values.")
    parser.add_argument(
        "--span-window",
        type=int,
        default=2048,
        help="Local span window around the dominant outer span when k is large.",
    )
    parser.add_argument(
        "--exact-span-limit",
        type=int,
        default=4000,
        help="Compute the full outer span sum exactly when k is at most this large.",
    )
    parser.add_argument(
        "--adapt-rounds",
        type=int,
        default=4,
        help="Number of adaptive window expansions for the outer span sum.",
    )
    parser.add_argument("--top-k", type=int, default=5, help="How many dominant contributors to print.")
    args = parser.parse_args()

    if args.w is not None and (args.w_min is not None or args.w_max is not None or args.w_values is not None):
        parser.error("Use either --w for a single weight or a weight sweep, not both.")

    args.w_values = parse_int_list(args.w_values)
    args.m_values = parse_int_list(args.m_values)
    args.M_values = parse_int_list(args.M_values)

    if args.w is not None:
        args.w_values = [args.w]

    try:
        outer_streams, actual_rate = choose_outer_streams(args)
    except ValueError as e:
        print(f"error: {e}", file=sys.stderr)
        return 2

    weights = choose_weights(args, args.n)
    memory_values = parse_memory_values(args)
    outer_memory_values = parse_outer_memory_values(args)

    effective_n = (args.n // outer_streams) * outer_streams
    if effective_n != args.n and not args.allow_rounding:
        print(
            f"error: n={args.n} is not divisible by the derived stream count r={outer_streams}. "
            f"Use a compatible n or pass --allow-rounding.",
            file=sys.stderr,
        )
        return 2

    if effective_n != args.n:
        print(
            f"warning: rounding n down to the nearest multiple of r={outer_streams}: "
            f"{args.n} -> {effective_n}",
            file=sys.stderr,
        )

    print("Dense+dense finite-n evaluation")
    print(f"  target n          : {args.n}")
    print(f"  effective n       : {effective_n}")
    print(f"  outer rate        : {actual_rate:.12g}")
    print(f"  outer streams r   : {outer_streams}")
    print(f"  delta             : {args.delta:.12g}")
    print(f"  outer M           : {args.M}")
    print(f"  inner m           : {args.m}")
    print(f"  epsilon           : {args.epsilon:.12g}")
    print(f"  xi                : {args.xi:.12g}")
    print()

    total_best_log = float("-inf")
    best_row: Optional[Tuple[int, int, int, OuterSpectrumResult, InnerEnvelope, float]] = None
    rows: List[Tuple[int, int, int, OuterSpectrumResult, InnerEnvelope, float]] = []

    for M_val in outer_memory_values:
        for m_val in memory_values:
            for w in weights:
                outer, inner, log_term = evaluate_single_setting(
                    n=effective_n,
                    outer_rate=actual_rate,
                    outer_streams=outer_streams,
                    delta=args.delta,
                    M=M_val,
                    m=m_val,
                    epsilon=args.epsilon,
                    xi=args.xi,
                    w=w,
                    args=args,
                )
                row = (M_val, m_val, w, outer, inner, log_term)
                rows.append(row)
                if log_term > total_best_log:
                    total_best_log = log_term
                    best_row = row

    if not rows:
        print("No weights were selected for evaluation.", file=sys.stderr)
        return 2

    objective_log = logsumexp([row[5] for row in rows])
    objective = safe_prob_from_log(objective_log)

    print("Main bound components")
    print(f"  evaluated points   : {len(rows)}")
    print(f"  objective upper bd : {format_from_log(objective_log)}")
    print()

    # Print the most influential rows in descending contribution order.
    rows_sorted = sorted(rows, key=lambda row: row[5], reverse=True)[: max(1, args.top_k)]
    for idx, (M_val, m_val, w, outer, inner, log_term) in enumerate(rows_sorted, start=1):
        print(f"Contributor #{idx}")
        print(f"  outer M          : {M_val}")
        print(f"  memory m         : {m_val}")
        print(f"  weight w         : {w}")
        print(f"  outer A_w^dense  : {format_term(outer.log_aw)}")
        print(f"  inner B_w^dense  : {format_term(inner.best_log)}")
        print(f"  product term     : {format_term(log_term)}")
        print(f"  inner best run r : {inner.best_run_target}")
        print(
            f"  inner run tail   : {format_term(inner.run_tail_log)}"
            f"  | off term: {format_term(inner.off_log)}"
            f"  | LD term: {format_term(inner.ld_log)}"
        )
        print(f"  outer exact sum?  : {'yes' if outer.exact else 'no'}")
        print_top_spans(outer, args.top_k)
        print()

    if best_row is not None:
        M_val, m_val, w, outer, inner, log_term = best_row
        print("Best single point")
        print(f"  outer M          : {M_val}")
        print(f"  memory m         : {m_val}")
        print(f"  weight w         : {w}")
        print(f"  outer A_w^dense  : {format_term(outer.log_aw)}")
        print(f"  inner B_w^dense  : {format_term(inner.best_log)}")
        print(f"  product term     : {format_term(log_term)}")
        print(f"  best run target  : {inner.best_run_target}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
