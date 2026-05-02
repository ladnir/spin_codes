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
from functools import lru_cache
from typing import Iterable, List

import numpy as np


ETA_CRIT = 1.0 - 2.0 ** -0.5
OUTER_LOW_SLOPE = math.log2(1.0 + math.sqrt(2.0))
RUNMIX_GAMMA_INTERCEPT = 0.14461613
RUNMIX_GAMMA_SLOPE = -0.42472604
RUNMIX_GAMMA_SAFE_INTERCEPT = 0.1347
RUNMIX_GAMMA_SAFE_SLOPE = -0.3910
RUNMIX_LAMBDA_INTERCEPT = -0.003929694918965849
RUNMIX_LAMBDA_SLOPE = 0.024519978056113938


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


@lru_cache(maxsize=None)
def outer_small_h_exact_prefix(k_msg: int, sigma: int, h_max: int) -> tuple[float, ...]:
    """Exact low-weight outer counts for 1 <= h <= h_max.

    This computes the whole low-weight prefix at once. The first-span term is
    explicit, and the span>=2 contribution is propagated by the coefficient
    recurrence induced by multiplying by (1+x)^2 and dividing by 2 per added
    span position:

        g_{ell+1}[j] = 0.5 g_ell[j] + g_ell[j-1] + 0.5 g_ell[j-2]

    where g_ell[j] = C(2 ell + sigma - 2, j) / 2^ell.
    """
    if h_max <= 0:
        return tuple()

    vals = np.zeros(h_max + 1, dtype=np.float64)
    two_to_minus_sigma = 2.0 ** (-sigma)
    two_to_minus_sigma_plus_1 = 0.5 * two_to_minus_sigma

    for h in range(1, h_max + 1):
        j = h - 1
        if 0 <= j <= sigma + 1:
            vals[h] += k_msg * math.comb(sigma + 1, j) * two_to_minus_sigma_plus_1

    if h_max >= 2:
        max_j = h_max - 2
        g = np.zeros(max_j + 1, dtype=np.float64)
        for j in range(max_j + 1):
            if j <= sigma + 2:
                g[j] = math.comb(sigma + 2, j) / 4.0

        peak_add = 0.0
        decay_count = 0
        for ell in range(2, k_msg + 1):
            mult = (k_msg - ell + 1) * two_to_minus_sigma
            if mult != 0.0:
                add = mult * g
                vals[2:] += add

                add_peak = float(np.max(add))
                if add_peak > peak_add:
                    peak_add = add_peak
                    decay_count = 0
                elif peak_add > 0.0 and add_peak <= peak_add * (2.0 ** -80):
                    decay_count += 1
                    if k_msg > 10000 and decay_count >= 64:
                        break
                else:
                    decay_count = 0

            if ell != k_msg:
                nxt = np.zeros(max_j + 1, dtype=np.float64)
                if max_j >= 0:
                    nxt[0] = 0.5 * g[0]
                if max_j >= 1:
                    nxt[1] = 0.5 * g[1] + g[0]
                if max_j >= 2:
                    nxt[2:] = 0.5 * g[2:] + g[1:-1] + 0.5 * g[:-2]
                g = nxt

    logs = [float("-inf")] * (h_max + 1)
    for h in range(1, h_max + 1):
        if vals[h] > 0.0:
            logs[h] = math.log2(vals[h])
    return tuple(logs)


def outer_small_h_exact_log2(k_msg: int, sigma: int, h: int) -> float:
    if h <= 0:
        return float("-inf")
    return outer_small_h_exact_prefix(k_msg, sigma, h)[h]


@lru_cache(maxsize=None)
def balls_bins_cap_cached(balls: int, bins: int, cap: int) -> int:
    if balls < 0 or bins < 0:
        return 0
    if balls == 0:
        return 1
    if cap < 0 or bins == 0:
        return 0
    if balls * 2 > bins * cap:
        balls = bins * cap - balls
    if balls < bins * cap:
        d = 0
        for i in range(bins):
            bb = bins + balls - i * (cap + 1) - 1
            if bb < bins - 1 or bb < 0:
                break
            term = math.comb(bins, i) * math.comb(bb, bins - 1)
            d = d - term if (i & 1) else d + term
        return d
    if balls == bins * cap:
        return 1
    return 0


@lru_cache(maxsize=None)
def count_inputs_banded(w: int, q: int, runs: int, k: int, n: int, sigma: int) -> int:
    if w == 0 and q == 0 and runs == 0:
        return 1
    if w <= 0 or q <= 0 or runs <= 0:
        return 0
    if runs > w or sigma <= 0 or q < w or q > n or q > w * sigma:
        return 0

    mini_runs = w - runs
    assign_mini_runs = math.comb(w - 1, runs - 1)
    total = 0
    for c in range(sigma):
        tail_trim = min(sigma - 1 - c, n - k)
        q_tilde = q - tail_trim
        zeros_inside_runs = q_tilde - w - (runs - 1) * (sigma - 1) - c
        ways_mini_run_shapes = balls_bins_cap_cached(zeros_inside_runs, mini_runs, sigma - 2)
        if ways_mini_run_shapes == 0:
            continue

        outside_zeros = k - q_tilde
        boundary_bins = runs + 1 if c == sigma - 1 else runs
        if outside_zeros < 0 or boundary_bins < 0:
            continue
        if outside_zeros == 0 and boundary_bins == 0:
            ways_outside = 1
        elif boundary_bins <= 0:
            ways_outside = 0
        else:
            ways_outside = math.comb(outside_zeros + boundary_bins - 1, boundary_bins - 1)
        if ways_outside == 0:
            continue

        total += assign_mini_runs * ways_mini_run_shapes * ways_outside
    return total


@lru_cache(maxsize=None)
def parity_enum_banded_log2(k: int, n: int, sigma: int, w: int, h: int) -> float:
    if w == 0 and h == 0:
        return 0.0
    if w == 0 or h > n or sigma == 0:
        return float("-inf")
    if sigma == 1:
        val = math.comb(k, w) * math.comb(w, h) * (2.0 ** -w)
        return math.log2(val) if val > 0 else float("-inf")

    total = float("-inf")
    q_max = min(n, w * sigma)
    for q in range(w, q_max + 1):
        choose_qh = math.comb(q, h) if 0 <= h <= q else 0
        if choose_qh == 0:
            continue
        count_wq = 0
        for runs in range(1, w + 1):
            count_wq += count_inputs_banded(w, q, runs, k, n, sigma)
        if count_wq == 0:
            continue
        term = math.log2(count_wq) + math.log2(choose_qh) - q
        total = term if total == float("-inf") else log2add(total, term)
    return total


@lru_cache(maxsize=None)
def outer_small_h_banded_systematic_prefix(k_msg: int, sigma: int, parity_n: int, h_max: int) -> tuple[float, ...]:
    logs = [float("-inf")] * (h_max + 1)
    for h in range(1, h_max + 1):
        total = float("-inf")
        w_hi = min(k_msg, h)
        for w in range(1, w_hi + 1):
            hp = h - w
            term = parity_enum_banded_log2(k_msg, parity_n, sigma, w, hp)
            if term != float("-inf"):
                total = term if total == float("-inf") else log2add(total, term)
        logs[h] = total
    return tuple(logs)


def outer_small_h_banded_systematic_log2(k_msg: int, sigma: int, parity_n: int, h: int) -> float:
    if h <= 0:
        return float("-inf")
    return outer_small_h_banded_systematic_prefix(k_msg, sigma, parity_n, h)[h]


def outer_small_h_prefix(
    k_msg: int,
    sigma: int,
    h_max: int,
    *,
    outer_mode: str,
    parity_n: int | None = None,
) -> tuple[float, ...]:
    if outer_mode == "conv":
        return outer_small_h_exact_prefix(k_msg, sigma, h_max)
    if outer_mode == "banded":
        if parity_n is None:
            raise ValueError("banded outer mode requires parity_n")
        return outer_small_h_banded_systematic_prefix(k_msg, sigma, parity_n, h_max)
    raise ValueError(f"unknown outer_mode: {outer_mode}")


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
    outer_mode: str = "conv",
    parity_n: int | None = None,
) -> float:
    total = float("-inf")
    outer_logs = outer_small_h_prefix(k_msg, sigma, h_hi, outer_mode=outer_mode, parity_n=parity_n)
    for h in range(h_lo, h_hi + 1):
        out = outer_logs[h]
        inn = exact_smallw_inner_bound(n, h, sigma, delta, xi)
        if inn <= 0.0:
            continue
        total = log2add(total, out + math.log2(inn))
    return total


def adaptive_small_h_rigorous(
    *,
    n: int,
    k_msg: int,
    sigma: int,
    delta: float,
    xi: float,
    h_hi: int,
    outer_mode: str = "conv",
    parity_n: int | None = None,
    stop_gap_bits: float = 20.0,
    tail_confirm: int = 8,
    h_chunk: int = 16,
) -> tuple[float, int, int, float]:
    total = float("-inf")
    best_term = float("-inf")
    best_h = 0
    last_h = 0
    decay_count = 0

    done = False
    h_done = 0
    while h_done < h_hi and not done:
        h_cap = min(h_hi, h_done + h_chunk)
        outer_logs = outer_small_h_prefix(k_msg, sigma, h_cap, outer_mode=outer_mode, parity_n=parity_n)
        for h in range(h_done + 1, h_cap + 1):
            out = outer_logs[h]
            inn = exact_smallw_inner_bound(n, h, sigma, delta, xi)
            if inn <= 0.0 or out == float("-inf"):
                continue
            term = out + math.log2(inn)
            total = term if total == float("-inf") else log2add(total, term)
            last_h = h

            if term > best_term:
                best_term = term
                best_h = h
                decay_count = 0
            elif term <= best_term - stop_gap_bits:
                decay_count += 1
                if decay_count >= tail_confirm:
                    done = True
                    break
            else:
                decay_count = 0
        h_done = h_cap

    return total, last_h, best_h, best_term


def adaptive_small_h_localtail(
    *,
    n: int,
    k_msg: int,
    sigma: int,
    delta: float,
    xi: float,
    h_hi: int,
    outer_mode: str = "conv",
    parity_n: int | None = None,
    stop_gap_bits: float = 20.0,
    tail_confirm: int = 8,
    ratio_window: int = 4,
    h_chunk: int = 16,
) -> tuple[float, int, int, float, float]:
    """Exact tiny-weight sum plus a local ratio tail heuristic.

    This is not a proof object. It uses the actual exact term sequence and, once
    the terms have fallen well below their peak for several consecutive h,
    approximates the unseen tail by a geometric continuation whose ratio is the
    maximum observed ratio across the last few exact steps.
    """
    total = float("-inf")
    best_term = float("-inf")
    best_h = 0
    last_h = 0
    decay_count = 0
    term_logs: list[float] = []
    done = False
    h_done = 0
    while h_done < h_hi and not done:
        h_cap = min(h_hi, h_done + h_chunk)
        outer_logs = outer_small_h_prefix(k_msg, sigma, h_cap, outer_mode=outer_mode, parity_n=parity_n)
        for h in range(h_done + 1, h_cap + 1):
            out = outer_logs[h]
            inn = exact_smallw_inner_bound(n, h, sigma, delta, xi)
            if inn <= 0.0 or out == float("-inf"):
                continue
            term = out + math.log2(inn)
            term_logs.append(term)
            total = term if total == float("-inf") else log2add(total, term)
            last_h = h

            if term > best_term:
                best_term = term
                best_h = h
                decay_count = 0
            elif term <= best_term - stop_gap_bits:
                decay_count += 1
                if decay_count >= tail_confirm:
                    done = True
                    break
            else:
                decay_count = 0
        h_done = h_cap

    tail_log = float("-inf")
    if len(term_logs) >= ratio_window + 1:
        ratios = []
        for idx in range(len(term_logs) - ratio_window, len(term_logs)):
            prev = term_logs[idx - 1]
            cur = term_logs[idx]
            ratios.append(2.0 ** (cur - prev))
        q = max(ratios)
        if 0.0 < q < 1.0:
            last_term = term_logs[-1]
            tail_log = last_term + math.log2(q / (1.0 - q))

    total_with_tail = total if tail_log == float("-inf") else log2add(total, tail_log)
    return total_with_tail, last_h, best_h, best_term, tail_log


def small_h_heuristic_log2(
    *,
    n: int,
    k_msg: int,
    sigma: int,
    h_lo: int,
    h_hi: int,
    theta: float,
    outer_mode: str = "conv",
    parity_n: int | None = None,
) -> float:
    total = float("-inf")
    outer_logs = outer_small_h_prefix(k_msg, sigma, h_hi, outer_mode=outer_mode, parity_n=parity_n)
    for h in range(h_lo, h_hi + 1):
        eta = h / n
        out = outer_logs[h]
        inn = psi_run(eta, theta)
        total = log2add(total, out - n * inn)
    return total


def single_run_product_log2(
    *,
    k_msg: int,
    sigma: int,
    h_lo: int,
    h_hi: int,
    q_single: float,
    outer_mode: str = "conv",
    parity_n: int | None = None,
) -> tuple[float, int, float]:
    total = float("-inf")
    peak_h = 0
    peak_term = float("-inf")
    outer_logs = outer_small_h_prefix(k_msg, sigma, h_hi, outer_mode=outer_mode, parity_n=parity_n)
    q_log = math.log2(q_single)
    for h in range(h_lo, h_hi + 1):
        out = outer_logs[h]
        if out == float("-inf"):
            continue
        term = out + h * q_log
        total = term if total == float("-inf") else log2add(total, term)
        if term > peak_term:
            peak_term = term
            peak_h = h
    return total, peak_h, peak_term


def runmix_gamma_from_q(q_single: float) -> float:
    """Pairwise interaction penalty inferred from large exact inner-only tables.

    The model is anchored at the one-run factor q_single and uses a pairwise
    damping term exp(-gamma * C(s,2)) for s runs. On the largest exact inner
    tables currently available (n=120,128), the effective gamma extracted from
    the two-run slice is well fit by a simple linear law in q_single over the
    regime q_single ~= 0.17..0.25.
    """
    return max(0.0, RUNMIX_GAMMA_INTERCEPT + RUNMIX_GAMMA_SLOPE * q_single)


def runmix_gamma_safe_from_q(q_single: float) -> float:
    """Lower-envelope pairwise damping law for theorem-shaped experiments.

    This is a conservative affine lower envelope for the minimum effective
    pairwise exponent seen on the current exact inner-only overlap tables across
    s=2..6. It is meant to be closer to something one could plausibly prove than
    the least-squares central fit, while still retaining the same qualitative
    structure.
    """
    return max(0.0, RUNMIX_GAMMA_SAFE_INTERCEPT + RUNMIX_GAMMA_SAFE_SLOPE * q_single)


def runmix_lambda_from_q(q_single: float) -> float:
    """Small cubic correction inferred after the pairwise fit.

    After removing the pairwise damping term, the larger exact inner-only tables
    still show a mild additional suppression once the run count reaches roughly
    4--6. A tiny positive cubic coefficient captures that effect well while
    remaining essentially invisible on the smaller n=120 table. We therefore
    keep it as an optional refinement rather than part of the default central
    lane.
    """
    return max(0.0, RUNMIX_LAMBDA_INTERCEPT + RUNMIX_LAMBDA_SLOPE * q_single)


def run_count_weight(n: int, w: int, s: int) -> float:
    return math.comb(w - 1, s - 1) * math.comb(n - w + 1, s) / math.comb(n, w)


def run_mixture_model_log2(
    *,
    n: int,
    k_msg: int,
    sigma: int,
    h_lo: int,
    h_hi: int,
    q_single: float,
    gamma: float,
    outer_mode: str = "conv",
    parity_n: int | None = None,
) -> tuple[float, int, float]:
    total = float("-inf")
    peak_h = 0
    peak_term = float("-inf")
    outer_logs = outer_small_h_prefix(k_msg, sigma, h_hi, outer_mode=outer_mode, parity_n=parity_n)
    for h in range(h_lo, h_hi + 1):
        out = outer_logs[h]
        if out == float("-inf"):
            continue
        p = 0.0
        for s in range(1, h + 1):
            p += run_count_weight(n, h, s) * (q_single**s) * math.exp(-gamma * s * (s - 1) / 2.0)
        if p <= 0.0:
            continue
        term = out + math.log2(p)
        total = term if total == float("-inf") else log2add(total, term)
        if term > peak_term:
            peak_term = term
            peak_h = h
    return total, peak_h, peak_term


def run_mixture_cubic_model_log2(
    *,
    n: int,
    k_msg: int,
    sigma: int,
    h_lo: int,
    h_hi: int,
    q_single: float,
    gamma: float,
    lambd: float,
    outer_mode: str = "conv",
    parity_n: int | None = None,
) -> tuple[float, int, float]:
    total = float("-inf")
    peak_h = 0
    peak_term = float("-inf")
    outer_logs = outer_small_h_prefix(k_msg, sigma, h_hi, outer_mode=outer_mode, parity_n=parity_n)
    for h in range(h_lo, h_hi + 1):
        out = outer_logs[h]
        if out == float("-inf"):
            continue
        p = 0.0
        for s in range(1, h + 1):
            p += (
                run_count_weight(n, h, s)
                * (q_single**s)
                * math.exp(
                    -gamma * s * (s - 1) / 2.0
                    - lambd * s * (s - 1) * (s - 2) / 6.0
                )
            )
        if p <= 0.0:
            continue
        term = out + math.log2(p)
        total = term if total == float("-inf") else log2add(total, term)
        if term > peak_term:
            peak_term = term
            peak_h = h
    return total, peak_h, peak_term


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
    parser.add_argument(
        "--outer-smallh-mode",
        choices=("conv", "banded"),
        default="banded",
        help="Small-weight outer model: 'conv' uses the paper convolutional proxy, 'banded' uses the actual sysBand combinatorics.",
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
    parser.add_argument(
        "--show-adaptive-smallw",
        action="store_true",
        help="Also report an adaptive exact small-weight finite-n bound and the stopping point used before the rigorous residual tail.",
    )
    parser.add_argument(
        "--show-localtail-smallw",
        action="store_true",
        help="Also report an adaptive exact small-weight model with a local geometric tail fit from the exact term sequence.",
    )
    parser.add_argument(
        "--show-single-run-model",
        action="store_true",
        help="Also report a small-weight model using exact outer counts and a per-run factor q_single^h.",
    )
    parser.add_argument(
        "--show-run-mixture-model",
        action="store_true",
        help="Also report a small-weight model using exact outer counts and a run-count mixture with q_s = q^s exp(-gamma s(s-1)/2).",
    )
    parser.add_argument(
        "--show-run-mixture-safe",
        action="store_true",
        help="Also report a theorem-shaped safe pairwise run-mixture model using a conservative lower envelope for gamma(q).",
    )
    parser.add_argument(
        "--show-run-mixture-cubic",
        action="store_true",
        help="Also report a refined run-mixture model with an additional cubic run-count damping term.",
    )
    parser.add_argument(
        "--q-single",
        type=float,
        default=None,
        help="Override the single-run factor q used by --show-single-run-model. Default is 2*delta.",
    )
    parser.add_argument(
        "--gamma-runmix",
        type=float,
        default=None,
        help="Override the run-mixture damping gamma. Default uses the pairwise-calibrated law gamma(q).",
    )
    parser.add_argument(
        "--show-run-mixture-legacy",
        action="store_true",
        help="Also report the legacy run-mixture model with gamma = 0.2*q_single for comparison.",
    )
    parser.add_argument(
        "--show-bracket",
        action="store_true",
        help="Report the current optimistic / central / pessimistic bracket explicitly.",
    )
    parser.add_argument("--adaptive-stop-gap", type=float, default=20.0)
    parser.add_argument("--adaptive-tail-confirm", type=int, default=8)
    parser.add_argument("--localtail-ratio-window", type=int, default=4)
    parser.add_argument(
        "--bracket-h-cap",
        type=int,
        default=24,
        help="Maximum h used by the optimistic/central/pessimistic bracket models.",
    )
    args = parser.parse_args()

    print(f"Concrete dense+dense large-k evaluator")
    print(f"  message length k : {args.k}")
    print(f"  n mode           : {args.n_mode}")
    print(f"  outer small-h    : {args.outer_smallh_mode}")
    print(f"  delta            : {args.delta}")
    print(f"  z                : {args.z}")
    print(f"  rho              : {args.rho}")
    print(f"  kappa            : {args.kappa}")
    print(f"  theta            : {args.theta}")
    print(f"  xi               : {args.xi}")
    print()

    for sigma in sigma_values(args):
        n = total_length(args.k, sigma, args.n_mode)
        parity_n = n - args.k
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
        h_resid_lo = h_cap + 1
        h_resid = geometric_window_bound(
            n=n,
            h_lo=h_resid_lo,
            h_hi=h_mid_hi,
            z=args.z,
            rho=args.rho,
            w_out=w_out,
        )
        log2_h_small = log2_h_total = float("-inf")
        log2_r_small = log2_r_total = float("-inf")
        log2_a_small = log2_a_total = a_resid = float("-inf")
        h_stop = h_peak = 0
        h_peak_log2 = float("-inf")
        log2_l_total = log2_l_tail = float("-inf")
        l_stop = l_peak = 0
        l_peak_log2 = float("-inf")
        q_single = args.q_single if args.q_single is not None else min(1.0, 2.0 * args.delta)
        log2_q_total = q_peak_log2 = float("-inf")
        q_peak_h = 0
        gamma_runmix = args.gamma_runmix if args.gamma_runmix is not None else runmix_gamma_from_q(q_single)
        log2_m_total = m_peak_log2 = float("-inf")
        m_peak_h = 0
        gamma_runmix_safe = runmix_gamma_safe_from_q(q_single)
        log2_ms_total = ms_peak_log2 = float("-inf")
        ms_peak_h = 0
        lambda_runmix = runmix_lambda_from_q(q_single)
        log2_mc_total = mc_peak_log2 = float("-inf")
        mc_peak_h = 0
        gamma_runmix_legacy = 0.2 * q_single
        log2_m0_total = m0_peak_log2 = float("-inf")
        m0_peak_h = 0

        if args.show_heuristic:
            log2_h_small = small_h_heuristic_log2(
                n=n,
                k_msg=args.k,
                sigma=sigma,
                h_lo=1,
                h_hi=h_cap,
                theta=args.theta,
                outer_mode=args.outer_smallh_mode,
                parity_n=parity_n,
            )
            log2_h_total = log2_h_small
            if h_resid > 0.0:
                log2_h_total = log2add(log2_h_total, math.log2(h_resid))
            if log2_s_lin != float("-inf"):
                log2_h_total = log2add(log2_h_total, log2_s_lin)
            if log2_s_top != float("-inf"):
                log2_h_total = log2add(log2_h_total, log2_s_top)

        if args.show_exact_smallw:
            log2_r_small = small_h_rigorous_log2(
                n=n,
                k_msg=args.k,
                sigma=sigma,
                delta=args.delta,
                xi=args.xi,
                h_lo=1,
                h_hi=h_cap,
                outer_mode=args.outer_smallh_mode,
                parity_n=parity_n,
            )
            log2_r_total = log2_r_small
            if h_resid > 0.0:
                log2_r_total = log2add(log2_r_total, math.log2(h_resid))
            if log2_s_lin != float("-inf"):
                log2_r_total = log2add(log2_r_total, log2_s_lin)
            if log2_s_top != float("-inf"):
                log2_r_total = log2add(log2_r_total, log2_s_top)

        if args.show_adaptive_smallw:
            log2_a_small, h_stop, h_peak, h_peak_log2 = adaptive_small_h_rigorous(
                n=n,
                k_msg=args.k,
                sigma=sigma,
                delta=args.delta,
                xi=args.xi,
                h_hi=h_mid_hi,
                outer_mode=args.outer_smallh_mode,
                parity_n=parity_n,
                stop_gap_bits=args.adaptive_stop_gap,
                tail_confirm=args.adaptive_tail_confirm,
            )
            a_resid = geometric_window_bound(
                n=n,
                h_lo=h_stop + 1,
                h_hi=h_mid_hi,
                z=args.z,
                rho=args.rho,
                w_out=w_out,
            )
            log2_a_total = log2_a_small
            if a_resid > 0.0:
                log2_a_total = log2add(log2_a_total, math.log2(a_resid))
            if log2_s_lin != float("-inf"):
                log2_a_total = log2add(log2_a_total, log2_s_lin)
            if log2_s_top != float("-inf"):
                log2_a_total = log2add(log2_a_total, log2_s_top)

        need_bracketish = args.show_localtail_smallw or args.show_single_run_model or args.show_run_mixture_model or args.show_run_mixture_safe or args.show_run_mixture_cubic or args.show_run_mixture_legacy or args.show_bracket
        bracket_h_hi = min(args.bracket_h_cap, h_mid_hi)
        if args.show_localtail_smallw or args.show_bracket:
            log2_l_total, l_stop, l_peak, l_peak_log2, log2_l_tail = adaptive_small_h_localtail(
                n=n,
                k_msg=args.k,
                sigma=sigma,
                delta=args.delta,
                xi=args.xi,
                h_hi=bracket_h_hi if args.outer_smallh_mode == "banded" else h_mid_hi,
                outer_mode=args.outer_smallh_mode,
                parity_n=parity_n,
                stop_gap_bits=args.adaptive_stop_gap,
                tail_confirm=args.adaptive_tail_confirm,
                ratio_window=args.localtail_ratio_window,
            )
            if log2_s_lin != float("-inf"):
                log2_l_total = log2add(log2_l_total, log2_s_lin)
            if log2_s_top != float("-inf"):
                log2_l_total = log2add(log2_l_total, log2_s_top)

        if args.show_single_run_model or args.show_bracket:
            log2_q_total, q_peak_h, q_peak_log2 = single_run_product_log2(
                k_msg=args.k,
                sigma=sigma,
                h_lo=1,
                h_hi=bracket_h_hi,
                q_single=q_single,
                outer_mode=args.outer_smallh_mode,
                parity_n=parity_n,
            )
            if log2_s_lin != float("-inf"):
                log2_q_total = log2add(log2_q_total, log2_s_lin)
            if log2_s_top != float("-inf"):
                log2_q_total = log2add(log2_q_total, log2_s_top)

        if args.show_run_mixture_model or args.show_bracket:
            log2_m_total, m_peak_h, m_peak_log2 = run_mixture_model_log2(
                n=n,
                k_msg=args.k,
                sigma=sigma,
                h_lo=1,
                h_hi=bracket_h_hi,
                q_single=q_single,
                gamma=gamma_runmix,
                outer_mode=args.outer_smallh_mode,
                parity_n=parity_n,
            )
            if log2_s_lin != float("-inf"):
                log2_m_total = log2add(log2_m_total, log2_s_lin)
            if log2_s_top != float("-inf"):
                log2_m_total = log2add(log2_m_total, log2_s_top)

        if args.show_run_mixture_safe:
            log2_ms_total, ms_peak_h, ms_peak_log2 = run_mixture_model_log2(
                n=n,
                k_msg=args.k,
                sigma=sigma,
                h_lo=1,
                h_hi=bracket_h_hi,
                q_single=q_single,
                gamma=gamma_runmix_safe,
                outer_mode=args.outer_smallh_mode,
                parity_n=parity_n,
            )
            if log2_s_lin != float("-inf"):
                log2_ms_total = log2add(log2_ms_total, log2_s_lin)
            if log2_s_top != float("-inf"):
                log2_ms_total = log2add(log2_ms_total, log2_s_top)

        if args.show_run_mixture_cubic:
            log2_mc_total, mc_peak_h, mc_peak_log2 = run_mixture_cubic_model_log2(
                n=n,
                k_msg=args.k,
                sigma=sigma,
                h_lo=1,
                h_hi=bracket_h_hi,
                q_single=q_single,
                gamma=gamma_runmix,
                lambd=lambda_runmix,
                outer_mode=args.outer_smallh_mode,
                parity_n=parity_n,
            )
            if log2_s_lin != float("-inf"):
                log2_mc_total = log2add(log2_mc_total, log2_s_lin)
            if log2_s_top != float("-inf"):
                log2_mc_total = log2add(log2_mc_total, log2_s_top)

        if args.show_run_mixture_legacy:
            log2_m0_total, m0_peak_h, m0_peak_log2 = run_mixture_model_log2(
                n=n,
                k_msg=args.k,
                sigma=sigma,
                h_lo=1,
                h_hi=bracket_h_hi,
                q_single=q_single,
                gamma=gamma_runmix_legacy,
                outer_mode=args.outer_smallh_mode,
                parity_n=parity_n,
            )
            if log2_s_lin != float("-inf"):
                log2_m0_total = log2add(log2_m0_total, log2_s_lin)
            if log2_s_top != float("-inf"):
                log2_m0_total = log2add(log2_m0_total, log2_s_top)

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
        if args.show_adaptive_smallw:
            print(f"  A_small (adaptive exact small-weight sum)                : {format_log2(log2_a_small)}")
            print(f"  A_tail  (rigorous residual after adaptive stop)          : {format_prob(a_resid)}")
            print(f"  A_total (adaptive finite-n mixed)                        : {format_log2(log2_a_total)}")
            print(f"  A_peak  (dominant h, log2 contribution)                  : h={h_peak}, log2={h_peak_log2:.3f}")
            print(f"  A_stop  (last exact h included)                          : {h_stop}")
        if args.show_localtail_smallw:
            print(f"  L_total (adaptive exact + local tail model)              : {format_log2(log2_l_total)}")
            print(f"  L_tail  (local geometric continuation only)              : {format_log2(log2_l_tail)}")
            print(f"  L_peak  (dominant h, log2 contribution)                  : h={l_peak}, log2={l_peak_log2:.3f}")
            print(f"  L_stop  (last exact h included)                          : {l_stop}")
        if args.show_single_run_model:
            print(f"  Q_total (exact outer + q^h inner model)                  : {format_log2(log2_q_total)}")
            print(f"  Q_peak  (dominant h, log2 contribution)                  : h={q_peak_h}, log2={q_peak_log2:.3f}")
            print(f"  Q_q     (single-run factor)                              : {q_single:.6f}")
        if args.show_run_mixture_model:
            print(f"  M_total (exact outer + run-mixture inner model)          : {format_log2(log2_m_total)}")
            print(f"  M_peak  (dominant h, log2 contribution)                  : h={m_peak_h}, log2={m_peak_log2:.3f}")
            print(f"  M_q     (single-run factor)                              : {q_single:.6f}")
            print(f"  M_gamma (run-mixture damping)                            : {gamma_runmix:.6f}")
        if args.show_run_mixture_safe:
            print(f"  MS_total (safe pairwise run-mixture)                     : {format_log2(log2_ms_total)}")
            print(f"  MS_peak  (dominant h, log2 contribution)                 : h={ms_peak_h}, log2={ms_peak_log2:.3f}")
            print(f"  MS_q     (single-run factor)                             : {q_single:.6f}")
            print(f"  MS_gamma (safe pairwise damping)                         : {gamma_runmix_safe:.6f}")
        if args.show_run_mixture_cubic:
            print(f"  MC_total (run-mixture + cubic correction)                : {format_log2(log2_mc_total)}")
            print(f"  MC_peak  (dominant h, log2 contribution)                 : h={mc_peak_h}, log2={mc_peak_log2:.3f}")
            print(f"  MC_q     (single-run factor)                             : {q_single:.6f}")
            print(f"  MC_gamma (pairwise damping)                              : {gamma_runmix:.6f}")
            print(f"  MC_lambda (cubic damping)                                : {lambda_runmix:.6f}")
        if args.show_run_mixture_legacy:
            print(f"  M0_total (legacy run-mixture, gamma=0.2q)                : {format_log2(log2_m0_total)}")
            print(f"  M0_peak  (dominant h, log2 contribution)                 : h={m0_peak_h}, log2={m0_peak_log2:.3f}")
            print(f"  M0_gamma (legacy damping)                                : {gamma_runmix_legacy:.6f}")
        if args.show_bracket:
            print(f"  B_opt   (optimistic = local-tail)                        : {format_log2(log2_l_total)}")
            print(f"  B_ctr   (central = pairwise run-mixture)                 : {format_log2(log2_m_total)}")
            print(f"  B_pess  (pessimistic = single-run)                       : {format_log2(log2_q_total)}")
        if args.lambda_target is not None:
            target = -args.lambda_target
            ok = lg <= target
            print(f"  meets 2^-lambda target?       : {'yes' if ok else 'no'} (target log2 <= {target:.3f})")
        print()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
