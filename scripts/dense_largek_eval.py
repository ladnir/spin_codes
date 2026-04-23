#!/usr/bin/env python3
"""Concrete large-k evaluator for the dense+dense theorem regime.

This tool is meant for the regime where exact enumeration is no longer practical.
It evaluates the current dense+dense proof ingredients on a fixed message length k
and memory parameter sigma, and reports a concrete upper bound / proxy for

    Pr[d_min <= delta n]

for the current systematic-dense outer plus dense-inner line.

The output is split into:
  - rigorous finite-n pieces that come directly from explicit formulas;
  - an exponent-mode estimate on the linear-weight window, clearly labeled.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from typing import Iterable, List


ETA_CRIT = 1.0 - 2.0 ** -0.5
OUTER_LOW_SLOPE = math.log2(1.0 + math.sqrt(2.0))


def h2(x: float) -> float:
    if x <= 0.0 or x >= 1.0:
        if x == 0.0 or x == 1.0:
            return 0.0
        raise ValueError(f"x must lie in [0,1], got {x}")
    return -x * math.log2(x) - (1.0 - x) * math.log2(1.0 - x)


def log2add(x: float, y: float) -> float:
    if x == float("-inf"):
        return y
    if y == float("-inf"):
        return x
    hi = max(x, y)
    lo = min(x, y)
    return hi + math.log2(1.0 + 2.0 ** (lo - hi))


def log2_binom(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return (
        math.lgamma(n + 1.0)
        - math.lgamma(k + 1.0)
        - math.lgamma(n - k + 1.0)
    ) / math.log(2.0)


def binom_cdf_half_upper(n: int, k: int, *, exact_limit: int = 2000, tail_terms: int = 128) -> float:
    if k < 0:
        return 0.0
    if k >= n:
        return 1.0
    if k > n // 2:
        return 1.0 - binom_cdf_half_upper(n, n - k - 1, exact_limit=exact_limit, tail_terms=tail_terms)

    def log_binom(nn: int, kk: int) -> float:
        return math.lgamma(nn + 1.0) - math.lgamma(kk + 1.0) - math.lgamma(nn - kk + 1.0)

    def safe_prob(logp: float) -> float:
        if logp == float("-inf"):
            return 0.0
        if logp > 0.0:
            return 1.0
        if logp < -745.0:
            return 0.0
        return math.exp(logp)

    def logsumexp(vals: List[float]) -> float:
        m = max(vals)
        if math.isinf(m):
            return m
        s = sum(math.exp(v - m) for v in vals)
        return m + math.log(s)

    if k <= exact_limit:
        logs = [log_binom(n, i) - n * math.log(2.0) for i in range(k + 1)]
        return safe_prob(logsumexp(logs))

    logs = [log_binom(n, k) - n * math.log(2.0)]
    cur_log = logs[0]
    cur_i = k
    for _ in range(max(0, tail_terms - 1)):
        if cur_i == 0:
            break
        cur_log += math.log(cur_i) - math.log(n - cur_i + 1)
        cur_i -= 1
        logs.append(cur_log)
        if (math.log(cur_i + 1) - math.log(n - cur_i) if cur_i > 0 else -999.0) < -30.0:
            break

    partial = safe_prob(logsumexp(logs))
    last_ratio = cur_i / (n - cur_i + 1) if cur_i > 0 else 0.0
    if last_ratio >= 1.0 or last_ratio <= 0.0:
        return min(1.0, partial)
    remainder = safe_prob(logs[-1]) * (last_ratio / (1.0 - last_ratio))
    return min(1.0, partial + remainder)


def outer_generating_bound(k_msg: int, sigma: int, z: float) -> float:
    a = (1.0 + z) / 2.0
    q = ((1.0 + z) ** 2) / 2.0
    if not (0.0 < q < 1.0):
        raise ValueError("q(z) must lie in (0,1)")

    span1 = k_msg * z * (a ** (sigma + 1))

    # S = sum_{ell=2}^k (k-ell+1) q^ell
    k = k_msg
    geom = (q * q) * (1.0 - q ** (k - 1)) / (1.0 - q)
    arith_full = q * (1.0 - (k + 1) * (q ** k) + k * (q ** (k + 1))) / ((1.0 - q) ** 2)
    arith = arith_full - q
    s = (k + 1) * geom - arith
    span_ge_2 = (z * z) * (a ** sigma) * s / ((1.0 + z) ** 2)
    return span1 + span_ge_2


def outer_small_h_exact_log2(k_msg: int, sigma: int, h: int) -> float:
    total = float("-inf")

    first = math.log2(k_msg) + log2_binom(sigma + 1, h - 1) - (sigma + 1)
    total = log2add(total, first)

    for ell in range(2, k_msg + 1):
        mult = k_msg - ell + 1
        term = (
            math.log2(mult)
            + log2_binom(2 * ell + sigma - 2, h - 2)
            - (ell + sigma)
        )
        total = log2add(total, term)

    return total


def run_tail_exact_log2(n: int, w: int, r: int) -> float:
    if r <= 1:
        return float("-inf")
    if r > w:
        r = w

    denom = log2_binom(n, w)
    total = float("-inf")
    for s in range(1, r):
        term = log2_binom(w - 1, s - 1) + log2_binom(n - w + 1, s) - denom
        total = log2add(total, term)
    return total


def tiny_window_bound(n: int, sigma: int, delta: float, z: float, xi: float, kappa: float, w_out: float) -> tuple[int, float]:
    h0 = int(math.floor(kappa * math.log2(n)))
    if h0 < 1:
        return 0, 0.0
    L = int(math.ceil((2.0 + xi) * delta * n))
    b = binom_cdf_half_upper(L, int(math.floor(delta * n)))
    total = 0.0
    for h in range(1, h0 + 1):
        outer = w_out / (z ** h)
        inner = (h * (h - 1)) / (n - 1) + ((L * (2.0 ** (-sigma))) ** h) + h * b
        total += outer * inner
    return h0, total


def exact_smallw_inner_bound(n: int, w: int, sigma: int, delta: float, xi: float) -> float:
    if w <= 0:
        return 0.0

    L = int(math.ceil((2.0 + xi) * delta * n))
    b = binom_cdf_half_upper(L, int(math.floor(delta * n)))
    base = L * (2.0 ** (-sigma))

    best = 1.0
    for r in range(1, w + 1):
        log2_run = run_tail_exact_log2(n, w, r)
        run_term = 0.0 if log2_run == float("-inf") else 2.0 ** log2_run
        off_term = min(1.0, base ** r)
        bin_term = min(1.0, r * b)
        val = min(1.0, run_term + off_term + bin_term)
        if val < best:
            best = val
    return best


def small_h_rigorous_log2(
    *,
    n: int,
    k_msg: int,
    sigma: int,
    delta: float,
    xi: float,
    h_lo: int,
    h_hi: int,
) -> float:
    total = float("-inf")
    for h in range(h_lo, h_hi + 1):
        out = outer_small_h_exact_log2(k_msg, sigma, h)
        inn = exact_smallw_inner_bound(n, h, sigma, delta, xi)
        if inn <= 0.0:
            continue
        total = log2add(total, out + math.log2(inn))
    return total


def small_h_heuristic_log2(
    *,
    n: int,
    k_msg: int,
    sigma: int,
    h_lo: int,
    h_hi: int,
    theta: float,
) -> float:
    total = float("-inf")
    for h in range(h_lo, h_hi + 1):
        eta = h / n
        out = outer_small_h_exact_log2(k_msg, sigma, h)
        inn = psi_run(eta, theta)
        total = log2add(total, out - n * inn)
    return total


def geometric_window_bound(
    *,
    n: int,
    h_lo: int,
    h_hi: int,
    z: float,
    rho: float,
    w_out: float,
) -> float:
    if h_hi < h_lo:
        return 0.0
    ratio = rho / z
    if ratio >= 1.0:
        return float("inf")
    return 3.0 * w_out * (ratio ** h_lo) * (1.0 - ratio ** (h_hi - h_lo + 1)) / (1.0 - ratio)


def outer_exponent(eta: float) -> float:
    if eta <= ETA_CRIT:
        return OUTER_LOW_SLOPE * eta
    return h2(eta) - 0.5


def psi_run(eta: float, theta: float) -> float:
    alpha = theta * eta / (1.0 - eta)
    return h2(eta) - eta * h2(theta) - (1.0 - eta) * h2(alpha)


def e_bin(delta: float, xi: float) -> float:
    return (2.0 + xi) * delta * (1.0 - h2(1.0 / (2.0 + xi)))


@dataclass
class GapResult:
    worst_gap: float
    worst_eta: float
    outer_at_worst: float
    inner_at_worst: float
    run_at_worst: float
    bin_at_worst: float


def linear_window_gap(delta: float, theta: float, xi: float, eta_hi: float, step: float) -> GapResult:
    worst_gap = -float("inf")
    worst_eta = ETA_CRIT
    outer_val = inner_val = run_val = bin_val = 0.0
    eb = e_bin(delta, xi)

    eta = ETA_CRIT
    while eta <= min(eta_hi, 1.0 - theta):
        out = outer_exponent(eta)
        run = psi_run(eta, theta)
        inn = min(run, eb)
        gap = out - inn
        if gap > worst_gap:
            worst_gap = gap
            worst_eta = eta
            outer_val = out
            inner_val = inn
            run_val = run
            bin_val = eb
        eta += step

    return GapResult(worst_gap, worst_eta, outer_val, inner_val, run_val, bin_val)


def format_prob(p: float) -> str:
    if p <= 0.0:
        return "0"
    return f"{p:.6e} (2^{{{math.log2(p):.3f}}})"


def format_log2(log2p: float) -> str:
    if log2p == float("-inf"):
        return "0"
    if log2p < -900.0:
        return f"2^{{{log2p:.3f}}}"
    p = 2.0 ** log2p
    return f"{p:.6e} (2^{{{log2p:.3f}}})"


def sigma_values(args: argparse.Namespace) -> Iterable[int]:
    if args.sigma is not None:
        return [args.sigma]
    if args.sigma_min is None or args.sigma_max is None:
        raise ValueError("Need either --sigma or both --sigma-min/--sigma-max")
    return range(args.sigma_min, args.sigma_max + 1, args.sigma_step)


def total_length(k_msg: int, sigma: int, n_mode: str) -> int:
    if n_mode == "paper":
        return 2 * k_msg
    if n_mode == "extended":
        return 2 * k_msg + sigma
    raise ValueError(f"Unknown n_mode: {n_mode}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--k", type=int, required=True, help="Message length k (current dense paper uses n=2k).")
    parser.add_argument("--sigma", type=int, default=None, help="Single memory value sigma=M=m.")
    parser.add_argument("--sigma-min", type=int, default=None)
    parser.add_argument("--sigma-max", type=int, default=None)
    parser.add_argument("--sigma-step", type=int, default=1)
    parser.add_argument(
        "--n-mode",
        choices=("paper", "extended"),
        default="paper",
        help="Use n=2k ('paper') or n=2k+sigma ('extended', matching the enumerated sysBand line).",
    )
    parser.add_argument("--delta", type=float, default=0.12)
    parser.add_argument("--z", type=float, default=1.0 / 3.0)
    parser.add_argument("--rho", type=float, default=1.0 / 4.0)
    parser.add_argument("--kappa", type=float, default=0.5)
    parser.add_argument("--theta", type=float, default=0.005)
    parser.add_argument("--xi", type=float, default=8.0)
    parser.add_argument("--eta-hi", type=float, default=0.99)
    parser.add_argument("--gap-step", type=float, default=1e-5)
    parser.add_argument("--lambda-target", type=float, default=None, help="Optional target z <= 2^-lambda.")
    parser.add_argument(
        "--heuristic-h-cap",
        type=int,
        default=64,
        help="Highest weight h included in the exact-outer/asymptotic-inner heuristic splice.",
    )
    parser.add_argument(
        "--show-heuristic",
        action="store_true",
        help="Also report heuristic small-h and total values using exact small-h outer counts plus the asymptotic run exponent.",
    )
    parser.add_argument(
        "--show-exact-smallw",
        action="store_true",
        help="Also report a stronger finite-n small-weight bound using exact outer counts and the exact run-tail inner formula up to h_cap.",
    )
    args = parser.parse_args()

    print(f"Concrete dense+dense large-k evaluator")
    print(f"  message length k : {args.k}")
    print(f"  n mode           : {args.n_mode}")
    print(f"  delta            : {args.delta}")
    print(f"  z                : {args.z}")
    print(f"  rho              : {args.rho}")
    print(f"  kappa            : {args.kappa}")
    print(f"  theta            : {args.theta}")
    print(f"  xi               : {args.xi}")
    print()

    for sigma in sigma_values(args):
        n = total_length(args.k, sigma, args.n_mode)
        w_out = outer_generating_bound(args.k, sigma, args.z)
        h0, s_tiny = tiny_window_bound(n, sigma, args.delta, args.z, args.xi, args.kappa, w_out)
        h1 = h0 + 1
        h_mid_hi = int(math.floor(ETA_CRIT * n))
        s_low = geometric_window_bound(n=n, h_lo=h1, h_hi=h_mid_hi, z=args.z, rho=args.rho, w_out=w_out)

        gap = linear_window_gap(args.delta, args.theta, args.xi, args.eta_hi, args.gap_step)
        lin_count = max(0, int(math.floor(min(args.eta_hi, 1.0 - args.theta) * n)) - int(math.ceil(ETA_CRIT * n)) + 1)
        log2_s_lin = float("-inf") if lin_count == 0 else math.log2(lin_count) + n * gap.worst_gap
        s_lin = 0.0 if log2_s_lin < -900.0 else 2.0 ** log2_s_lin

        top_exp = h2(args.eta_hi) - 0.5
        top_count = max(0, n - int(math.ceil(args.eta_hi * n)) + 1)
        log2_s_top = float("-inf") if top_count == 0 else math.log2(top_count) + n * top_exp
        s_top = 0.0 if log2_s_top < -900.0 else 2.0 ** log2_s_top

        total = s_tiny + s_low + s_lin + s_top
        lg = math.log2(total) if total > 0.0 else float("-inf")

        h_cap = min(args.heuristic_h_cap, h_mid_hi)
        log2_h_small = small_h_heuristic_log2(
            n=n,
            k_msg=args.k,
            sigma=sigma,
            h_lo=1,
            h_hi=h_cap,
            theta=args.theta,
        )
        h_resid_lo = h_cap + 1
        h_resid = geometric_window_bound(
            n=n,
            h_lo=h_resid_lo,
            h_hi=h_mid_hi,
            z=args.z,
            rho=args.rho,
            w_out=w_out,
        )
        log2_h_total = log2_h_small
        if h_resid > 0.0:
            log2_h_total = log2add(log2_h_total, math.log2(h_resid))
        if log2_s_lin != float("-inf"):
            log2_h_total = log2add(log2_h_total, log2_s_lin)
        if log2_s_top != float("-inf"):
            log2_h_total = log2add(log2_h_total, log2_s_top)

        log2_r_small = small_h_rigorous_log2(
            n=n,
            k_msg=args.k,
            sigma=sigma,
            delta=args.delta,
            xi=args.xi,
            h_lo=1,
            h_hi=h_cap,
        )
        log2_r_total = log2_r_small
        if h_resid > 0.0:
            log2_r_total = log2add(log2_r_total, math.log2(h_resid))
        if log2_s_lin != float("-inf"):
            log2_r_total = log2add(log2_r_total, log2_s_lin)
        if log2_s_top != float("-inf"):
            log2_r_total = log2add(log2_r_total, log2_s_top)

        print(f"sigma = {sigma}")
        print(f"  c = sigma/log2(n)             : {sigma / math.log2(n):.6f}")
        print(f"  outer generating bound W(z)   : {format_prob(w_out)}")
        print(f"  tiny window cutoff h0         : {h0}")
        print(f"  S_tiny  (explicit finite-n)   : {format_prob(s_tiny)}")
        print(f"  S_low   (explicit geometric)  : {format_prob(s_low)}")
        print(f"  gap worst on [eta_crit,{args.eta_hi}] : {gap.worst_gap:.12f} at eta={gap.worst_eta:.6f}")
        print(f"  S_lin   (exponent-mode)       : {format_log2(log2_s_lin)}")
        print(f"  S_top   (outer-only exponent) : {format_log2(log2_s_top)}")
        print(f"  TOTAL   (mixed bound/proxy)   : {format_prob(total)}")
        if args.show_heuristic:
            print(f"  H_small (exact outer + asym inner, h<= {h_cap}) : {format_log2(log2_h_small)}")
            print(f"  H_tail  (explicit geometric residual)           : {format_prob(h_resid)}")
            print(f"  H_total (heuristic mixed)                       : {format_log2(log2_h_total)}")
        if args.show_exact_smallw:
            print(f"  R_small (exact outer + exact run-tail inner, h<= {h_cap}) : {format_log2(log2_r_small)}")
            print(f"  R_total (stronger finite-n mixed)                         : {format_log2(log2_r_total)}")
        if args.lambda_target is not None:
            target = -args.lambda_target
            ok = lg <= target
            print(f"  meets 2^-lambda target?       : {'yes' if ok else 'no'} (target log2 <= {target:.3f})")
        print()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
