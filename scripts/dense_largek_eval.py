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
RUNMIX_SAFE_SLACK_BITS = 0.5
RUNMIX_LAMBDA_INTERCEPT = -0.003929694918965849
RUNMIX_LAMBDA_SLOPE = 0.024519978056113938
LATEBASE_BETA_INTERCEPT = 0.12557873
LATEBASE_BETA_SLOPE = -0.36374342
ISOLATED_BETA_BITS_INTERCEPT = 0.67
ISOLATED_BETA_BITS_SLOPE = -1.92
ISOLATED_BETA_BITS_SAFE_MARGIN = 0.03
ISOLATED_BETA_HDEP_INTERCEPT = 0.46605946
ISOLATED_BETA_HDEP_SLOPE_Q = -1.39605011
ISOLATED_BETA_HDEP_SLOPE_H = 0.01055334
ISOLATED_BETA_HDEP_SAFE_MARGIN = 0.02
ISOLATED_POSITIONAL_PAIR_SAFE_A_BITS = -0.50296039
ISOLATED_POSITIONAL_PAIR_SAFE_B_BITS = 0.186


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


def binom_cdf_half_entropy_upper(n: int, k: int) -> float:
    if k < 0:
        return 0.0
    if k >= n or 2 * k >= n:
        return 1.0
    exponent = -n * (1.0 - h2(k / n))
    if exponent < -1074.0:
        return 0.0
    return min(1.0, 2.0**exponent)


def binom_cdf_half_entropy_log2(n: int, k: int) -> float:
    if k < 0:
        return float("-inf")
    if k >= n or 2 * k >= n:
        return 0.0
    return -n * (1.0 - h2(k / n))


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


def fixedtap_banded_outer_gf_log2(k_msg: int, parity_n: int, sigma: int, z: float) -> float:
    """Log2 generating function for the fixed-tap banded systematic outer.

    This computes

        sum_x z^{wt(x)} z^{r(x)} ((1+z)/2)^{q(x)-r(x)}

    for the active fixed-tap banded outer, where q is the number of active
    parity coordinates and r is the number of active parity clusters. The zero
    message is excluded. The transfer state is the distance since the last
    message one, capped at sigma; state sigma means the current parity window is
    inactive.
    """
    if not (0.0 < z < 1.0):
        raise ValueError("z must lie in (0,1)")
    if k_msg < 0 or parity_n < 0 or sigma <= 0:
        raise ValueError("invalid fixed-tap banded outer parameters")
    states = sigma + 1
    inactive = sigma
    a = (1.0 + z) / 2.0

    def normalize_vec(v: np.ndarray) -> tuple[np.ndarray, float]:
        scale = float(np.max(v))
        if scale <= 0.0:
            return v, float("-inf")
        return v / scale, math.log2(scale)

    def normalize_mat(mtx: np.ndarray) -> tuple[np.ndarray, float]:
        scale = float(np.max(mtx))
        if scale <= 0.0:
            return mtx, float("-inf")
        return mtx / scale, math.log2(scale)

    def apply_power(v: np.ndarray, v_log: float, matrix: np.ndarray, exp: int) -> tuple[np.ndarray, float]:
        base = matrix
        base_log = 0.0
        while exp > 0:
            if exp & 1:
                v = base @ v
                v, inc = normalize_vec(v)
                v_log += base_log + inc
            exp >>= 1
            if exp:
                base = base @ base
                base, inc = normalize_mat(base)
                base_log = 2.0 * base_log + inc
        return v, v_log

    # Message/parity positions 1..min(k_msg, parity_n): choose x_t, then emit
    # parity coordinate t from the updated window.
    choose = np.zeros((states, states), dtype=np.float64)
    for d in range(states):
        prev_active = d < sigma
        # Choose x_t=1: systematic z, active parity; inactive->active creates
        # the deterministic fixed-tap cluster start z, active continuation is fair.
        choose[0, d] += z * (a if prev_active else z)
        # Choose x_t=0.
        nd = min(sigma, d + 1)
        if nd < sigma:
            choose[nd, d] += a
        else:
            choose[nd, d] += 1.0

    # Tail parity positions after the message ends: x_t is forced to zero.
    zero = np.zeros((states, states), dtype=np.float64)
    for d in range(states):
        nd = min(sigma, d + 1)
        if nd < sigma:
            zero[nd, d] += a
        else:
            zero[nd, d] += 1.0

    # Message-only positions if the parity block is shorter: choose x_t, but no
    # parity coordinate is emitted.
    msg_only = np.zeros((states, states), dtype=np.float64)
    for d in range(states):
        msg_only[0, d] += z
        msg_only[min(sigma, d + 1), d] += 1.0

    v = np.zeros(states, dtype=np.float64)
    v[inactive] = 1.0
    v_log = 0.0
    common = min(k_msg, parity_n)
    v, v_log = apply_power(v, v_log, choose, common)
    if parity_n > common:
        v, v_log = apply_power(v, v_log, zero, parity_n - common)
    if k_msg > common:
        v, v_log = apply_power(v, v_log, msg_only, k_msg - common)

    v_sum = float(np.sum(v))
    if v_sum <= 0.0:
        return float("-inf")
    log_all = v_log + math.log2(v_sum)
    # Remove the zero message, whose generating contribution is exactly one.
    if log_all <= 1e-10:
        val = max(0.0, (2.0**log_all) - 1.0)
        return math.log2(val) if val > 0.0 else float("-inf")
    return log_all + math.log2(1.0 - 2.0 ** (-log_all))


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
def parity_enum_banded_log2(k: int, n: int, sigma: int, w: int, h: int, fixed_tap: bool = False) -> float:
    if w == 0 and h == 0:
        return 0.0
    if w == 0 or h > n or sigma == 0:
        return float("-inf")
    if sigma == 1:
        if fixed_tap:
            val = math.comb(k, w) if h == w else 0
        else:
            val = math.comb(k, w) * math.comb(w, h) * (2.0 ** -w)
        return math.log2(val) if val > 0 else float("-inf")

    total = float("-inf")
    q_max = min(n, w * sigma)
    for q in range(w, q_max + 1):
        for runs in range(1, w + 1):
            count_wqr = count_inputs_banded(w, q, runs, k, n, sigma)
            if count_wqr == 0:
                continue
            if fixed_tap:
                # Each sigma-cluster start has a deterministic fixed diagonal
                # parity one. The other q-runs active parity coordinates remain
                # fair. Thus parity weight is runs + Bin(q-runs, 1/2).
                fair = q - runs
                need = h - runs
                choose = math.comb(fair, need) if 0 <= need <= fair else 0
                denom_exp = fair
            else:
                choose = math.comb(q, h) if 0 <= h <= q else 0
                denom_exp = q
            if choose == 0:
                continue
            term = math.log2(count_wqr) + math.log2(choose) - denom_exp
            total = term if total == float("-inf") else log2add(total, term)
    return total


@lru_cache(maxsize=None)
def outer_small_h_banded_systematic_prefix(
    k_msg: int,
    sigma: int,
    parity_n: int,
    h_max: int,
    fixed_tap: bool = False,
) -> tuple[float, ...]:
    logs = [float("-inf")] * (h_max + 1)
    for h in range(1, h_max + 1):
        total = float("-inf")
        w_hi = min(k_msg, h)
        for w in range(1, w_hi + 1):
            hp = h - w
            term = parity_enum_banded_log2(k_msg, parity_n, sigma, w, hp, fixed_tap=fixed_tap)
            if term != float("-inf"):
                total = term if total == float("-inf") else log2add(total, term)
        logs[h] = total
    return tuple(logs)


def outer_small_h_banded_systematic_log2(
    k_msg: int,
    sigma: int,
    parity_n: int,
    h: int,
    fixed_tap: bool = False,
) -> float:
    if h <= 0:
        return float("-inf")
    return outer_small_h_banded_systematic_prefix(k_msg, sigma, parity_n, h, fixed_tap=fixed_tap)[h]


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
    if outer_mode in ("banded", "banded-fixedtap"):
        if parity_n is None:
            raise ValueError("banded outer mode requires parity_n")
        return outer_small_h_banded_systematic_prefix(
            k_msg, sigma, parity_n, h_max, fixed_tap=(outer_mode == "banded-fixedtap")
        )
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
        late = 0.0
        if L - h + 1 >= h and n - h + 1 >= h:
            late = math.exp2(log2_binom(L - h + 1, h) - log2_binom(n - h + 1, h))
        fixedtap_term_base = h * (2.0 ** (1 - sigma))
        inner = (h * (h - 1)) / (n - 1) + late + fixedtap_term_base + b
        total += outer * inner
    return h0, total


def exact_smallw_inner_bound(n: int, w: int, sigma: int, delta: float, xi: float) -> float:
    if w <= 0:
        return 0.0

    L = int(math.ceil((2.0 + xi) * delta * n))
    b = binom_cdf_half_upper(L, int(math.floor(delta * n)))
    late = 0.0
    if L - w + 1 >= w and n - w + 1 >= w:
        late = math.exp2(log2_binom(L - w + 1, w) - log2_binom(n - w + 1, w))
    # Fixed-tap dense inner: a genuine termination must be aligned with an input
    # one in the critical state, so the crude candidate count is w rather than L.
    base = w * (2.0 ** (1 - sigma))

    adjacency = w * (w - 1) / max(1, n - 1)
    return min(1.0, adjacency + late + base + b)


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


def late_start_base_prob(n: int, w: int, s: int, l: int) -> float:
    if s < 1 or s > w:
        return 0.0
    if l - w + 1 < s:
        return 0.0
    return math.comb(l - w + 1, s) / math.comb(n - w + 1, s)


def late_start_mixture_model_log2(
    *,
    n: int,
    k_msg: int,
    sigma: int,
    h_lo: int,
    h_hi: int,
    l: int,
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
            p += run_count_weight(n, h, s) * late_start_base_prob(n, h, s, l)
        if p <= 0.0:
            continue
        term = out + math.log2(p)
        total = term if total == float("-inf") else log2add(total, term)
        if term > peak_term:
            peak_term = term
            peak_h = h
    return total, peak_h, peak_term


def latebase_beta_from_q(q_single: float) -> float:
    """Residual pairwise inflation after factoring out exact late placement.

    Units are natural-log per run pair, so the correction factor is
    exp(beta * C(s,2)).
    """
    return max(0.0, LATEBASE_BETA_INTERCEPT + LATEBASE_BETA_SLOPE * q_single)


def isolated_beta_bits_from_q(q_single: float) -> float:
    """Safe upper envelope for the isolated-slice residual, in log2 bits/pair.

    This is calibrated from exact inner-only overlap tables after factoring out:
      1. the exact isolated-slice mass Pr[R(U)=h], and
      2. the exact late-placement base on that slice.

    It is intentionally conservative and should be read as a theorem-shaped
    working envelope rather than a best fit.
    """
    return max(0.0, ISOLATED_BETA_BITS_INTERCEPT + ISOLATED_BETA_BITS_SLOPE * q_single)


def isolated_beta_hdep_bits_from_qh(q_single: float, h: int) -> float:
    """Working h-dependent envelope for isolated-slice residual, in log2 bits/pair."""
    return max(
        0.0,
        ISOLATED_BETA_HDEP_INTERCEPT
        + ISOLATED_BETA_HDEP_SLOPE_Q * q_single
        + ISOLATED_BETA_HDEP_SLOPE_H * max(0, h - 2),
    )


def isolated_full_run_prob(n: int, h: int) -> float:
    if h < 0 or h > n:
        return 0.0
    return math.comb(n - h + 1, h) / math.comb(n, h)


def isolated_order_late_prob(n: int, h: int, r: int, l: int) -> float:
    """Exact isolated-slice probability that the r-th one lies in the last l positions.

    Conditioned on wt(U)=h and R(U)=h, map the isolated one positions
    X_1 < ... < X_h to the ordinary subset Y_i = X_i - i + 1 of [N-h+1].
    Then X_r > N-l iff Y_r > N-l-r+1, i.e. at most r-1 of the Y_i lie
    before that threshold.
    """
    if r < 1 or r > h:
        return 0.0
    total_space = n - h + 1
    cutoff = n - l - r + 1
    tail = l - h + r
    if cutoff < 0:
        return 1.0
    num = 0
    for j in range(r):
        if j > cutoff:
            break
        tail_take = h - j
        if tail_take < 0 or tail_take > tail:
            continue
        num += math.comb(cutoff, j) * math.comb(tail, tail_take)
    return num / math.comb(total_space, h)


def isolated_early_count_prob(n: int, h: int, l: int, j: int) -> float:
    """Exact isolated-slice probability that exactly j starts occur before the final l positions."""
    total_space = n - h + 1
    tail = l - h + 1
    prefix = n - l
    if j < 0 or j > h or j > prefix or h - j > tail or tail < 0:
        return 0.0
    return math.comb(prefix, j) * math.comb(tail, h - j) / math.comb(total_space, h)


def isolated_early_pair_mean(n: int, h: int, l: int) -> float:
    """Expected number of early-start pairs on the isolated slice."""
    total_space = n - h + 1
    prefix = n - l
    if h < 2 or prefix < 2 or total_space < 2:
        return 0.0
    return (h * (h - 1) / 2.0) * (prefix * (prefix - 1)) / (total_space * (total_space - 1))


def _span_count_sum(prefix: int, l: int, h: int, span_lo: int, span_hi: int) -> int:
    """Count weight-h words whose early support span lies in [span_lo, span_hi]."""
    if h < 2 or span_hi < span_lo:
        return 0
    span_lo = max(2, span_lo)
    span_hi = min(prefix, span_hi)
    if span_hi < span_lo:
        return 0
    r = h - 2
    t_lo = span_lo - 2 + l
    t_hi = span_hi - 2 + l

    def sum_binom(a: int, b: int, rr: int) -> int:
        if b < a or rr < 0:
            return 0
        return math.comb(b + 1, rr + 1) - math.comb(a, rr + 1)

    sum_c_r = sum_binom(t_lo, t_hi, r)
    sum_c_r1 = sum_binom(t_lo, t_hi, r + 1)
    # With t=s-2+l, the number of placements for a fixed span is
    # (prefix - s + 1) C(s-2+l, h-2) = (prefix+l-1-t) C(t,r).
    sum_t_c_r = r * sum_c_r + (r + 1) * sum_c_r1
    return (prefix + l - 1) * sum_c_r - sum_t_c_r


def fixedtap_multiepisode_span_bound(
    n: int,
    h: int,
    sigma: int,
    delta: float,
    l: int,
    *,
    entropy_tail_bound: bool = True,
) -> float:
    """Pessimistic fixed-tap multi-episode upper bound for a full weight-h slice.

    If the early input ones span s positions, episodes that terminate within
    l steps can cover at most l+1 positions of that span each. Unless one such
    episode pays the binomial lower-tail event, at least ceil(s/(l+1)) aligned fixed-tap
    terminations are needed. The span distribution is exact under the uniform
    weight-h law.
    """
    if h <= 0:
        return 0.0
    cut = math.floor(delta * n)
    if l < 0:
        return 1.0
    l = min(n, l)
    bin_tail = (
        binom_cdf_half_entropy_upper(l, cut)
        if entropy_tail_bound
        else binom_cdf_half_upper(l, cut)
    )

    prefix = n - l
    if prefix <= 0:
        return min(1.0, 1.0 + bin_tail)

    denom_log2 = log2_binom(n, h)
    total = 0.0

    no_early_log2 = log2_binom(l, h) - denom_log2
    if no_early_log2 > -1074.0:
        total += min(1.0, 2.0**no_early_log2)

    one_early_log2 = math.log2(prefix) + log2_binom(l, h - 1) - denom_log2
    if one_early_log2 > -1074.0:
        needed_terms = 1
        log2_episode = log2_binom(h, needed_terms) - needed_terms * (sigma - 1)
        episode_term = 1.0 if log2_episode >= 0.0 else 2.0**log2_episode
        total += min(1.0, 2.0**one_early_log2) * episode_term

    cover = l + 1
    max_terms = (prefix + cover - 1) // cover
    for needed_terms in range(1, max_terms + 1):
        span_lo = max(2, (needed_terms - 1) * cover + 1)
        span_hi = min(prefix, needed_terms * cover)
        count = _span_count_sum(prefix, l, h, span_lo, span_hi)
        if count <= 0:
            continue
        count_log2 = math.log2(count)
        prob_log2 = count_log2 - denom_log2
        if prob_log2 < -1074.0:
            continue
        log2_episode = log2_binom(h, needed_terms) - needed_terms * (sigma - 1)
        episode_term = 1.0 if log2_episode >= 0.0 else 2.0**log2_episode
        total += min(1.0, 2.0**prob_log2) * episode_term
    return min(1.0, total + bin_tail)


def fixedtap_cover_internal_bound_log2(
    n: int,
    h: int,
    sigma: int,
    delta: float,
    l: int,
    *,
    grid_steps: int = 160,
) -> float:
    """Cover-count diagnostic for fixed-tap multi-episode suppression.

    A bad early pattern with no length-l surviving episode must be coverable by
    short start-to-termination intervals. For r intervals of total live length T,
    this lane uses the theorem-safe paired cover count

        C(n-T+r,r) * C(T-r-1,r-1) * C(T+l-2r,h-2r) / C(n,h),

    which overcounts by allowing the cover intervals to sit anywhere in the
    block. The 2r distinguished positions are the starts and aligned fixed-tap
    terminations. We then pay 2^{-r(sigma-1)} for those terminations and charge
    a fair lower-tail on roughly T-r*sigma non-terminal live outputs. This is a
    diagnostic implementation of the paired-cover theorem; the T loop is still
    sampled, not an exact summation.
    """
    if h <= 0:
        return float("-inf")
    cut = math.floor(delta * n)
    l = min(n, max(1, l))
    cover = l + 1
    denom = log2_binom(n, h)
    total = float("-inf")

    no_early = log2_binom(l, h) - denom
    total = log2add(total, no_early)

    for r in range(1, h // 2 + 1):
        t_min = 2 * r
        t_max = min(n, r * cover)
        if t_min > t_max:
            continue
        raw = {t_min, t_max, min(t_max, max(t_min, r * sigma)), min(t_max, max(t_min, 2 * cut + r * sigma + 1))}
        if t_max > t_min:
            for idx in range(grid_steps):
                frac = idx / max(1, grid_steps - 1)
                raw.add(int(round(t_min * ((t_max / t_min) ** frac))))
        candidates = sorted(t for t in raw if t_min <= t <= t_max)
        for t in candidates:
            if n - t + r < r:
                continue
            remaining_slots = min(n, t + l) - 2 * r
            remaining_ones = h - 2 * r
            if remaining_slots < remaining_ones:
                continue
            cover_count = (
                log2_binom(n - t + r, r)
                + log2_binom(t - r - 1, r - 1)
                + log2_binom(remaining_slots, remaining_ones)
                - denom
            )
            fair_len = max(0, t - r * sigma)
            fair_tail = binom_cdf_half_entropy_log2(fair_len, cut)
            total = log2add(total, cover_count - r * (sigma - 1) + fair_tail)
    return min(0.0, total)


def fixedtap_cover_internal_block_bound_log2(
    n: int,
    h: int,
    sigma: int,
    delta: float,
    l: int,
    *,
    block_ratio: float = 1.01,
) -> float:
    """Theorem-safe block upper sum for the paired-cover/internal-tail bound.

    This implements the corollary in the manuscript with a monotone block
    bound over T. On a block [a,b], the decreasing placement factor is evaluated
    at a, the increasing length/support factors at b, and the fair-tail upper
    bound at a. The block length is then paid explicitly.
    """
    if h <= 0:
        return float("-inf")
    if block_ratio <= 1.0:
        raise ValueError("block_ratio must be > 1")
    cut = math.floor(delta * n)
    l = min(n, max(1, l))
    cover = l + 1
    denom = log2_binom(n, h)
    total = log2_binom(l, h) - denom

    for r in range(1, h // 2 + 1):
        t_min = 2 * r
        t_max = min(n, r * cover)
        if t_min > t_max:
            continue
        a = t_min
        while a <= t_max:
            b = min(t_max, max(a, int(math.floor(a * block_ratio))))
            block_len = b - a + 1
            remaining_slots = b + l - 2 * r
            remaining_ones = h - 2 * r
            if remaining_slots >= remaining_ones:
                cover_count = (
                    log2_binom(n - a + r, r)
                    + log2_binom(b - r - 1, r - 1)
                    + log2_binom(remaining_slots, remaining_ones)
                    - denom
                )
                fair_len = max(0, a - r * sigma)
                fair_tail = binom_cdf_half_entropy_log2(fair_len, cut)
                term = math.log2(block_len) + cover_count - r * (sigma - 1) + fair_tail
                total = log2add(total, term)
            a = b + 1
    return min(0.0, total)


def isolated_position_single_model_prob(n: int, h: int, cut: int, m: int) -> float:
    """Position-sensitive isolated-slice model using one-run penalties by remaining suffix length.

    For an early isolated start at shifted position y <= n-2*cut, the corresponding remaining suffix
    has length t = n-y+1. This model assigns that start the theorem-shaped one-run penalty

        tau_y = min(1, t * 2^{-m} + Pr[Binomial(t, 1/2) <= cut]),

    then averages the product of those penalties over the exact isolated h-subset law.
    """
    total_space = n - h + 1
    l = 2 * cut
    tail = l - h + 1
    prefix = n - l
    if tail < 0:
        return 0.0
    taus = []
    for y in range(1, prefix + 1):
        t = n - y + 1
        tau = min(1.0, t * (2.0 ** (-m)) + binom_cdf_half_upper(t, cut))
        taus.append(tau)

    elem = [0.0] * (h + 1)
    elem[0] = 1.0
    for tau in taus:
        for j in range(h, 0, -1):
            elem[j] += tau * elem[j - 1]

    total = 0.0
    denom = math.comb(total_space, h)
    for j in range(h + 1):
        if h - j > tail:
            continue
        total += elem[j] * math.comb(tail, h - j) / denom
    return total


def isolated_position_single_model_prefix(
    n: int,
    h_hi: int,
    cut: int,
    m: int,
    *,
    entropy_tail_bound: bool,
) -> list[float]:
    """Compute positional one-run probabilities for all h <= h_hi.

    The large-k evaluator uses the entropy-tail option: it replaces each
    binomial CDF by a standard entropy upper bound, making this a pessimistic
    theorem-shaped lane and avoiding millions of exact CDF evaluations.
    """
    l = 2 * cut
    prefix = n - l
    elem = [0.0] * (h_hi + 1)
    elem[0] = 1.0
    for y in range(1, prefix + 1):
        t = n - y + 1
        bin_tail = (
            binom_cdf_half_entropy_upper(t, cut)
            if entropy_tail_bound
            else binom_cdf_half_upper(t, cut)
        )
        tau = min(1.0, t * (2.0 ** (-m)) + bin_tail)
        for j in range(min(h_hi, y), 0, -1):
            elem[j] += tau * elem[j - 1]

    probs = [0.0] * (h_hi + 1)
    for h in range(1, h_hi + 1):
        total_space = n - h + 1
        tail = l - h + 1
        if tail < 0:
            continue
        denom = math.comb(total_space, h)
        total = 0.0
        for j in range(h + 1):
            if h - j > tail:
                continue
            total += elem[j] * math.comb(tail, h - j) / denom
        probs[h] = min(1.0, total)
    return probs


def isolated_slice_model_log2(
    *,
    n: int,
    k_msg: int,
    sigma: int,
    h_lo: int,
    h_hi: int,
    l: int,
    beta_bits_per_pair: float = 0.0,
    add_adjacency_remainder: bool = False,
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
        p_full = isolated_full_run_prob(n, h)
        late = late_start_base_prob(n, h, h, l)
        corr_bits = beta_bits_per_pair * h * (h - 1) / 2.0
        p = p_full * late * (2.0**corr_bits)
        if add_adjacency_remainder:
            p += 1.0 - p_full
        p = min(1.0, p)
        if p <= 0.0:
            continue
        term = out + math.log2(p)
        total = term if total == float("-inf") else log2add(total, term)
        if term > peak_term:
            peak_term = term
            peak_h = h
    return total, peak_h, peak_term


def isolated_window_model_log2(
    *,
    n: int,
    k_msg: int,
    sigma: int,
    delta: float,
    z: float,
    rho: float,
    xi: float,
    kappa: float,
    theta: float,
    eta_hi: float,
    gap_step: float,
    h_iso_hi: int,
    beta_bits_per_pair: float = 0.0,
    add_adjacency_remainder: bool = False,
    outer_mode: str = "conv",
    parity_n: int | None = None,
) -> tuple[float, float, int, float, int]:
    w_out = outer_generating_bound(k_msg, sigma, z)
    h0, _ = tiny_window_bound(n, sigma, delta, z, xi, kappa, w_out)
    h_iso_hi = min(h_iso_hi, h0)
    l = math.ceil(2.0 * delta * n)
    log2_iso, peak_h, peak_log2 = isolated_slice_model_log2(
        n=n,
        k_msg=k_msg,
        sigma=sigma,
        h_lo=1,
        h_hi=h_iso_hi,
        l=l,
        beta_bits_per_pair=beta_bits_per_pair,
        add_adjacency_remainder=add_adjacency_remainder,
        outer_mode=outer_mode,
        parity_n=parity_n,
    )
    h_lo_rest = h_iso_hi + 1
    h_mid_hi = int(math.floor(ETA_CRIT * n))
    log2_rest = float("-inf")
    if h_lo_rest <= h_mid_hi:
        s_rest = geometric_window_bound(n=n, h_lo=h_lo_rest, h_hi=h_mid_hi, z=z, rho=rho, w_out=w_out)
        if s_rest > 0.0:
            log2_rest = math.log2(s_rest)
    gap = linear_window_gap(delta, theta, xi, eta_hi, gap_step)
    lin_count = max(0, int(math.floor(min(eta_hi, 1.0 - theta) * n)) - int(math.ceil(ETA_CRIT * n)) + 1)
    log2_lin = float("-inf") if lin_count == 0 else math.log2(lin_count) + n * gap.worst_gap
    top_exp = h2(eta_hi) - 0.5
    top_count = max(0, n - int(math.ceil(eta_hi * n)) + 1)
    log2_top = float("-inf") if top_count == 0 else math.log2(top_count) + n * top_exp
    total = log2_iso
    if log2_rest != float("-inf"):
        total = log2add(total, log2_rest)
    if log2_lin != float("-inf"):
        total = log2add(total, log2_lin)
    if log2_top != float("-inf"):
        total = log2add(total, log2_top)
    return total, log2_iso, peak_h, peak_log2, h0


def isolated_runmix_split_model_log2(
    *,
    n: int,
    k_msg: int,
    sigma: int,
    delta: float,
    h_iso_hi: int,
    h_hi: int,
    beta_bits_per_pair: float,
    gamma_runmix_safe: float,
    outer_mode: str = "conv",
    parity_n: int | None = None,
) -> tuple[float, int, float]:
    total = float("-inf")
    peak_h = 0
    peak_term = float("-inf")
    outer_logs = outer_small_h_prefix(k_msg, sigma, h_hi, outer_mode=outer_mode, parity_n=parity_n)
    l = math.ceil(2.0 * delta * n)
    q_single = min(1.0, 2.0 * delta)
    for h in range(1, h_hi + 1):
        out = outer_logs[h]
        if out == float("-inf"):
            continue
        if h <= h_iso_hi:
            p_full = isolated_full_run_prob(n, h)
            p = p_full * late_start_base_prob(n, h, h, l) * (2.0 ** (beta_bits_per_pair * h * (h - 1) / 2.0))
            p += 1.0 - p_full
        else:
            p = 0.0
            for s in range(1, h + 1):
                p += run_count_weight(n, h, s) * (q_single**s) * math.exp(-gamma_runmix_safe * s * (s - 1) / 2.0)
        p = min(1.0, p)
        if p <= 0.0:
            continue
        term = out + math.log2(p)
        total = term if total == float("-inf") else log2add(total, term)
        if term > peak_term:
            peak_term = term
            peak_h = h
    return total, peak_h, peak_term


def isolated_slice_model_hdep_log2(
    *,
    n: int,
    k_msg: int,
    sigma: int,
    h_lo: int,
    h_hi: int,
    l: int,
    q_single: float,
    beta_margin_bits: float = 0.0,
    add_adjacency_remainder: bool = False,
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
        p_full = isolated_full_run_prob(n, h)
        late = late_start_base_prob(n, h, h, l)
        beta_bits_per_pair = isolated_beta_hdep_bits_from_qh(q_single, h) + beta_margin_bits
        corr_bits = beta_bits_per_pair * h * (h - 1) / 2.0
        p = p_full * late * (2.0**corr_bits)
        if add_adjacency_remainder:
            p += 1.0 - p_full
        p = min(1.0, p)
        if p <= 0.0:
            continue
        term = out + math.log2(p)
        total = term if total == float("-inf") else log2add(total, term)
        if term > peak_term:
            peak_term = term
            peak_h = h
    return total, peak_h, peak_term


def isolated_positional_pair_model_log2(
    *,
    n: int,
    k_msg: int,
    sigma: int,
    delta: float,
    adjacency_xi: float,
    h_lo: int,
    h_hi: int,
    pair_a_bits: float,
    pair_b_bits: float,
    entropy_tail_bound: bool = True,
    clamp_correction_nonnegative: bool = True,
    adjacency_mode: str = "none",
    outer_mode: str = "conv",
    parity_n: int | None = None,
) -> tuple[float, int, float]:
    """Position-sensitive isolated model with an early-pair residual inflation.

    The base term is the isolated-slice positional one-run model. The residual
    correction is measured in log2 bits and is affine in the exact expected
    number of early isolated pairs. This is a theorem-shaped diagnostic lane,
    not a proved bound.
    """
    total = float("-inf")
    peak_h = 0
    peak_term = float("-inf")
    cut = math.floor(delta * n)
    l = 2 * cut
    one_episode_l = min(n, max(l, math.ceil((2.0 + adjacency_xi) * delta * n)))
    one_episode_off = one_episode_l * (2.0 ** (-sigma))
    one_episode_bin = binom_cdf_half_entropy_upper(one_episode_l, cut)
    outer_logs = outer_small_h_prefix(k_msg, sigma, h_hi, outer_mode=outer_mode, parity_n=parity_n)
    pos_probs = isolated_position_single_model_prefix(
        n, h_hi, cut, sigma, entropy_tail_bound=entropy_tail_bound
    )
    for h in range(h_lo, h_hi + 1):
        out = outer_logs[h]
        if out == float("-inf"):
            continue
        p_full = isolated_full_run_prob(n, h)
        p_iso = pos_probs[h]
        epair = isolated_early_pair_mean(n, h, l)
        corr_bits = pair_a_bits + pair_b_bits * epair
        if clamp_correction_nonnegative:
            corr_bits = max(0.0, corr_bits)
        p = p_full * p_iso * (2.0**corr_bits)
        if adjacency_mode == "prob":
            p += 1.0 - p_full
        elif adjacency_mode == "late":
            for s in range(1, h):
                p += run_count_weight(n, h, s) * late_start_base_prob(n, h, s, l)
        elif adjacency_mode == "one-episode":
            p_noniso = 0.0
            for s in range(1, h):
                p_noniso += run_count_weight(n, h, s) * (
                    late_start_base_prob(n, h, s, one_episode_l) + one_episode_off + one_episode_bin
                )
            p += p_noniso
        elif adjacency_mode != "none":
            raise ValueError(f"unknown adjacency_mode: {adjacency_mode}")
        p = min(1.0, p)
        if p <= 0.0:
            continue
        term = out + math.log2(p)
        total = term if total == float("-inf") else log2add(total, term)
        if term > peak_term:
            peak_term = term
            peak_h = h
    return total, peak_h, peak_term


def fixedtap_multiepisode_model_log2(
    *,
    n: int,
    k_msg: int,
    sigma: int,
    delta: float,
    h_lo: int,
    h_hi: int,
    l: int,
    entropy_tail_bound: bool = True,
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
        p = fixedtap_multiepisode_span_bound(
            n,
            h,
            sigma,
            delta,
            l,
            entropy_tail_bound=entropy_tail_bound,
        )
        p = min(1.0, p)
        if p <= 0.0:
            continue
        term = out + math.log2(p)
        total = term if total == float("-inf") else log2add(total, term)
        if term > peak_term:
            peak_term = term
            peak_h = h
    return total, peak_h, peak_term


def fixedtap_cover_internal_model_log2(
    *,
    n: int,
    k_msg: int,
    sigma: int,
    delta: float,
    h_lo: int,
    h_hi: int,
    l: int,
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
        logp = fixedtap_cover_internal_bound_log2(n, h, sigma, delta, l)
        term = out + logp
        total = term if total == float("-inf") else log2add(total, term)
        if term > peak_term:
            peak_term = term
            peak_h = h
    return total, peak_h, peak_term


def fixedtap_cover_internal_block_model_log2(
    *,
    n: int,
    k_msg: int,
    sigma: int,
    delta: float,
    h_lo: int,
    h_hi: int,
    l: int,
    outer_mode: str = "conv",
    parity_n: int | None = None,
    block_ratio: float = 1.01,
) -> tuple[float, int, float]:
    total = float("-inf")
    peak_h = 0
    peak_term = float("-inf")
    outer_logs = outer_small_h_prefix(k_msg, sigma, h_hi, outer_mode=outer_mode, parity_n=parity_n)
    for h in range(h_lo, h_hi + 1):
        out = outer_logs[h]
        if out == float("-inf"):
            continue
        logp = fixedtap_cover_internal_block_bound_log2(
            n, h, sigma, delta, l, block_ratio=block_ratio
        )
        term = out + logp
        total = term if total == float("-inf") else log2add(total, term)
        if term > peak_term:
            peak_term = term
            peak_h = h
    return total, peak_h, peak_term


def latebase_pairwise_model_log2(
    *,
    n: int,
    k_msg: int,
    sigma: int,
    h_lo: int,
    h_hi: int,
    l: int,
    beta: float,
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
            p += run_count_weight(n, h, s) * late_start_base_prob(n, h, s, l) * math.exp(beta * s * (s - 1) / 2.0)
        if p <= 0.0:
            continue
        term = out + math.log2(p)
        total = term if total == float("-inf") else log2add(total, term)
        if term > peak_term:
            peak_term = term
            peak_h = h
    return total, peak_h, peak_term


def one_episode_correction_model_log2(
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
) -> tuple[float, int, float, int]:
    total = float("-inf")
    peak_h = 0
    peak_term = float("-inf")
    l = min(n, math.ceil((2.0 + xi) * delta * n))
    cut = math.floor(delta * n)
    off = l * (2.0 ** (-sigma))
    bin_tail = binom_cdf_half_upper(l, cut)
    outer_logs = outer_small_h_prefix(k_msg, sigma, h_hi, outer_mode=outer_mode, parity_n=parity_n)
    for h in range(h_lo, h_hi + 1):
        out = outer_logs[h]
        if out == float("-inf"):
            continue
        p = 0.0
        for s in range(1, h + 1):
            p += run_count_weight(n, h, s) * late_start_base_prob(n, h, s, l)
        p = min(1.0, p + off + bin_tail)
        if p <= 0.0:
            continue
        term = out + math.log2(p)
        total = term if total == float("-inf") else log2add(total, term)
        if term > peak_term:
            peak_term = term
            peak_h = h
    return total, peak_h, peak_term, l


def optimized_one_episode_model_log2(
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
    l_steps: int = 64,
) -> tuple[float, int, float, int]:
    cut = math.floor(delta * n)
    l_min = max(cut + 1, math.ceil(2.0 * delta * n))
    l_max = min(n, max(l_min, math.ceil((2.0 + xi) * delta * n)))
    if l_min > l_max:
        l_min = l_max = min(n, max(1, cut + 1))

    candidates = sorted({int(round(v)) for v in np.linspace(l_min, l_max, l_steps)} | {l_min, l_max})
    best = (float("inf"), 0, float("inf"), l_min)
    for l in candidates:
        total = float("-inf")
        peak_h = 0
        peak_term = float("-inf")
        off = l * (2.0 ** (-sigma))
        bin_tail = binom_cdf_half_upper(l, cut)
        outer_logs = outer_small_h_prefix(k_msg, sigma, h_hi, outer_mode=outer_mode, parity_n=parity_n)
        for h in range(h_lo, h_hi + 1):
            out = outer_logs[h]
            if out == float("-inf"):
                continue
            p = 0.0
            for s in range(1, h + 1):
                p += run_count_weight(n, h, s) * late_start_base_prob(n, h, s, l)
            p = min(1.0, p + off + bin_tail)
            if p <= 0.0:
                continue
            term = out + math.log2(p)
            total = term if total == float("-inf") else log2add(total, term)
            if term > peak_term:
                peak_term = term
                peak_h = h
        scalar = float("inf") if total == float("-inf") else total
        if scalar < best[0]:
            best = (scalar, peak_h, peak_term, l)
    return best


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
    # The episode length used in the forced-termination binomial tail cannot
    # exceed the block. Earlier exploratory checks used the formal expression
    # with (2+xi)delta > 1; capping exposes when that certificate is relying on
    # an impossible episode length.
    length_frac = min(1.0, (2.0 + xi) * delta)
    threshold_frac = delta / length_frac
    if threshold_frac >= 0.5:
        return 0.0
    return length_frac * (1.0 - h2(threshold_frac))


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
    if log2p == float("inf"):
        return "inf"
    if log2p > 900.0:
        return f"2^{{{log2p:.3f}}}"
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
        choices=("conv", "banded", "banded-fixedtap"),
        default="banded-fixedtap",
        help="Small-weight outer model: 'conv' uses the paper convolutional proxy, 'banded' uses the sysBand combinatorics, 'banded-fixedtap' fixes one parity tap per row.",
    )
    parser.add_argument("--delta", type=float, default=0.12)
    parser.add_argument("--z", type=float, default=1.0 / 3.0)
    parser.add_argument("--rho", type=float, default=1.0 / 4.0)
    parser.add_argument("--kappa", type=float, default=0.5)
    parser.add_argument("--theta", type=float, default=0.005)
    parser.add_argument("--xi", type=float, default=8.0)
    parser.add_argument(
        "--small-xi",
        type=float,
        default=0.1,
        help="Window slack used by small-weight fixed-tap bounds; kept separate from the linear-window xi.",
    )
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
        "--show-late-base-model",
        action="store_true",
        help="Also report the exact conditioned late-placement base mixture implied by Corollary first-run-late-qbase.",
    )
    parser.add_argument(
        "--show-one-episode-model",
        action="store_true",
        help="Also report the theorem-shaped one-episode correction model using L=(2+xi)delta n.",
    )
    parser.add_argument(
        "--show-one-episode-optimized",
        action="store_true",
        help="Also report the one-episode correction model after optimizing over L in [2delta n, (2+xi)delta n].",
    )
    parser.add_argument(
        "--show-fixedtap-multiepisode",
        action="store_true",
        help="Also report a pessimistic fixed-tap multi-episode bound using the exact full-slice early-span law.",
    )
    parser.add_argument(
        "--show-fixedtap-cover-internal",
        action="store_true",
        help="Also report the diagnostic cover-count model with an internal fair-output tail.",
    )
    parser.add_argument(
        "--show-fixedtap-cover-block",
        action="store_true",
        help="Also report a theorem-safe blocked upper sum for the paired-cover/internal-tail bound.",
    )
    parser.add_argument(
        "--fixedtap-cover-block-ratio",
        type=float,
        default=1.01,
        help="Multiplicative T-block ratio for --show-fixedtap-cover-block.",
    )
    parser.add_argument(
        "--show-latebase-pairwise",
        action="store_true",
        help="Also report exact late placement with a fitted pairwise residual inflation exp(beta*C(s,2)).",
    )
    parser.add_argument(
        "--show-isolated-model",
        action="store_true",
        help="Also report the exact isolated-slice late-placement base Pr[R=h] * late_base(h,h,h).",
    )
    parser.add_argument(
        "--show-isolated-safe",
        action="store_true",
        help="Also report a theorem-shaped isolated-slice model with a safe residual inflation on the R(U)=h slice.",
    )
    parser.add_argument(
        "--show-isolated-safe-adj",
        action="store_true",
        help="Also report the isolated-safe model after adding the full adjacency remainder 1-Pr[R=h].",
    )
    parser.add_argument(
        "--show-isolated-pos-pair",
        action="store_true",
        help="Also report the position-sensitive isolated one-run model with early-pair residual inflation.",
    )
    parser.add_argument(
        "--show-isolated-hybrid",
        action="store_true",
        help="Also report isolated-slice control only on the dominant tiny window, with the usual coarse handoff above it.",
    )
    parser.add_argument(
        "--show-split-safe",
        action="store_true",
        help="Also report a proof-oriented split lane: isolated-slice control on very tiny h, then safe run-mixture above it.",
    )
    parser.add_argument(
        "--show-split-safe-hdep",
        action="store_true",
        help="Also report the split-safe lane with an h-dependent isolated-slice residual envelope.",
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
        "--show-run-mixture-safe-slack",
        action="store_true",
        help="Also report the safe pairwise run-mixture model after adding an explicit multiplicative slack factor.",
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
        "--isolated-beta-bits",
        type=float,
        default=None,
        help="Override the isolated-slice residual envelope in log2 bits per isolated run pair.",
    )
    parser.add_argument(
        "--isolated-beta-safe-margin",
        type=float,
        default=ISOLATED_BETA_BITS_SAFE_MARGIN,
        help="Extra safety margin added to the isolated-slice residual envelope, in log2 bits per pair.",
    )
    parser.add_argument(
        "--runmix-safe-slack-bits",
        type=float,
        default=RUNMIX_SAFE_SLACK_BITS,
        help="Extra log2 slack added to the safe pairwise run-mixture lane. Default is the current empirical half-bit envelope.",
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
    parser.add_argument(
        "--isolated-h-cap",
        type=int,
        default=10,
        help="Maximum h controlled by the isolated-slice hybrid lane before handing off to the usual coarse windows.",
    )
    parser.add_argument(
        "--split-iso-h-cap",
        type=int,
        default=4,
        help="Maximum h controlled by the isolated slice in the split-safe lane.",
    )
    parser.add_argument(
        "--isolated-beta-hdep-safe-margin",
        type=float,
        default=ISOLATED_BETA_HDEP_SAFE_MARGIN,
        help="Extra safety margin added to the h-dependent isolated-slice residual envelope, in log2 bits per pair.",
    )
    parser.add_argument(
        "--pos-pair-a-bits",
        type=float,
        default=ISOLATED_POSITIONAL_PAIR_SAFE_A_BITS,
        help="Intercept for the positional early-pair residual lane, in log2 bits.",
    )
    parser.add_argument(
        "--pos-pair-b-bits",
        type=float,
        default=ISOLATED_POSITIONAL_PAIR_SAFE_B_BITS,
        help="Slope for the positional early-pair residual lane, in log2 bits per expected early pair.",
    )
    parser.add_argument(
        "--pos-pair-allow-negative-correction",
        action="store_true",
        help="Allow the positional early-pair affine correction to reduce the base model when negative.",
    )
    parser.add_argument(
        "--pos-pair-exact-binomial",
        action="store_true",
        help="Use exact binomial CDFs in the positional lane instead of entropy upper bounds. Intended for small exact-table checks.",
    )
    parser.add_argument(
        "--pos-pair-adjacency-mode",
        choices=("none", "prob", "late", "one-episode"),
        default="one-episode",
        help="How to bound the non-isolated R(U)<h contribution in the positional-pair lane.",
    )
    parser.add_argument(
        "--pos-pair-adjacency-xi",
        type=float,
        default=0.1,
        help="Window slack xi used only for the positional-pair non-isolated one-episode remainder.",
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
    print(f"  small xi         : {args.small_xi}")
    print()

    for sigma in sigma_values(args):
        n = total_length(args.k, sigma, args.n_mode)
        parity_n = n - args.k
        w_out = outer_generating_bound(args.k, sigma, args.z)
        h0, s_tiny = tiny_window_bound(n, sigma, args.delta, args.z, args.small_xi, args.kappa, w_out)
        h1 = h0 + 1
        h_mid_hi = int(math.floor(ETA_CRIT * n))
        s_low = geometric_window_bound(n=n, h_lo=h1, h_hi=h_mid_hi, z=args.z, rho=args.rho, w_out=w_out)

        gap = linear_window_gap(args.delta, args.theta, args.xi, args.eta_hi, args.gap_step)
        lin_count = max(0, int(math.floor(min(args.eta_hi, 1.0 - args.theta) * n)) - int(math.ceil(ETA_CRIT * n)) + 1)
        log2_s_lin = float("-inf") if lin_count == 0 else math.log2(lin_count) + n * gap.worst_gap
        s_lin = 0.0 if log2_s_lin < -900.0 else (float("inf") if log2_s_lin > 900.0 else 2.0 ** log2_s_lin)

        top_exp = h2(args.eta_hi) - 0.5
        top_count = max(0, n - int(math.ceil(args.eta_hi * n)) + 1)
        log2_s_top = float("-inf") if top_count == 0 else math.log2(top_count) + n * top_exp
        s_top = 0.0 if log2_s_top < -900.0 else (float("inf") if log2_s_top > 900.0 else 2.0 ** log2_s_top)

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
        log2_b_total = b_peak_log2 = float("-inf")
        b_peak_h = 0
        log2_e_total = e_peak_log2 = float("-inf")
        e_peak_h = 0
        e_l = 0
        log2_eo_total = eo_peak_log2 = float("-inf")
        eo_peak_h = 0
        eo_l = 0
        log2_lb_total = lb_peak_log2 = float("-inf")
        lb_peak_h = 0
        lb_l = 0
        beta_late = latebase_beta_from_q(q_single)
        beta_iso_bits = (
            args.isolated_beta_bits
            if args.isolated_beta_bits is not None
            else isolated_beta_bits_from_q(q_single) + args.isolated_beta_safe_margin
        )
        log2_i_total = i_peak_log2 = float("-inf")
        i_peak_h = 0
        log2_is_total = is_peak_log2 = float("-inf")
        is_peak_h = 0
        log2_isa_total = isa_peak_log2 = float("-inf")
        isa_peak_h = 0
        log2_ipp_total = ipp_peak_log2 = float("-inf")
        ipp_peak_h = 0
        log2_ih_total = log2_ih_iso = ih_peak_log2 = float("-inf")
        ih_peak_h = 0
        ih_h0 = 0
        log2_ss_total = ss_peak_log2 = float("-inf")
        ss_peak_h = 0
        log2_ssh_total = ssh_peak_log2 = float("-inf")
        ssh_peak_h = 0
        iso_l = math.ceil(2.0 * args.delta * n)
        gamma_runmix = args.gamma_runmix if args.gamma_runmix is not None else runmix_gamma_from_q(q_single)
        log2_m_total = m_peak_log2 = float("-inf")
        m_peak_h = 0
        gamma_runmix_safe = runmix_gamma_safe_from_q(q_single)
        log2_ms_total = ms_peak_log2 = float("-inf")
        ms_peak_h = 0
        log2_mss_total = mss_peak_log2 = float("-inf")
        mss_peak_h = 0
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
                xi=args.small_xi,
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
                xi=args.small_xi,
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

        bracket_h_hi = min(args.bracket_h_cap, h_mid_hi)
        if args.show_localtail_smallw or args.show_bracket:
            log2_l_total, l_stop, l_peak, l_peak_log2, log2_l_tail = adaptive_small_h_localtail(
                n=n,
                k_msg=args.k,
                sigma=sigma,
                delta=args.delta,
                xi=args.small_xi,
                h_hi=bracket_h_hi if args.outer_smallh_mode in ("banded", "banded-fixedtap") else h_mid_hi,
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

        if args.show_late_base_model:
            log2_b_total, b_peak_h, b_peak_log2 = late_start_mixture_model_log2(
                n=n,
                k_msg=args.k,
                sigma=sigma,
                h_lo=1,
                h_hi=bracket_h_hi,
                l=math.ceil(2.0 * args.delta * n),
                outer_mode=args.outer_smallh_mode,
                parity_n=parity_n,
            )
            if log2_s_lin != float("-inf"):
                log2_b_total = log2add(log2_b_total, log2_s_lin)
            if log2_s_top != float("-inf"):
                log2_b_total = log2add(log2_b_total, log2_s_top)

        if args.show_one_episode_model:
            log2_e_total, e_peak_h, e_peak_log2, e_l = one_episode_correction_model_log2(
                n=n,
                k_msg=args.k,
                sigma=sigma,
                delta=args.delta,
                xi=args.small_xi,
                h_lo=1,
                h_hi=bracket_h_hi,
                outer_mode=args.outer_smallh_mode,
                parity_n=parity_n,
            )
            if log2_s_lin != float("-inf"):
                log2_e_total = log2add(log2_e_total, log2_s_lin)
            if log2_s_top != float("-inf"):
                log2_e_total = log2add(log2_e_total, log2_s_top)

        if args.show_one_episode_optimized:
            log2_eo_total, eo_peak_h, eo_peak_log2, eo_l = optimized_one_episode_model_log2(
                n=n,
                k_msg=args.k,
                sigma=sigma,
                delta=args.delta,
                xi=args.small_xi,
                h_lo=1,
                h_hi=bracket_h_hi,
                outer_mode=args.outer_smallh_mode,
                parity_n=parity_n,
            )
            if log2_s_lin != float("-inf"):
                log2_eo_total = log2add(log2_eo_total, log2_s_lin)
            if log2_s_top != float("-inf"):
                log2_eo_total = log2add(log2_eo_total, log2_s_top)

        if args.show_fixedtap_multiepisode:
            fme_l = math.ceil((2.0 + args.small_xi) * args.delta * n)
            log2_fme_total, fme_peak_h, fme_peak_log2 = fixedtap_multiepisode_model_log2(
                n=n,
                k_msg=args.k,
                sigma=sigma,
                delta=args.delta,
                h_lo=1,
                h_hi=bracket_h_hi,
                l=fme_l,
                outer_mode=args.outer_smallh_mode,
                parity_n=parity_n,
            )
            if log2_s_lin != float("-inf"):
                log2_fme_total = log2add(log2_fme_total, log2_s_lin)
            if log2_s_top != float("-inf"):
                log2_fme_total = log2add(log2_fme_total, log2_s_top)

        if args.show_fixedtap_cover_internal:
            fci_l = math.ceil((2.0 + args.small_xi) * args.delta * n)
            log2_fci_total, fci_peak_h, fci_peak_log2 = fixedtap_cover_internal_model_log2(
                n=n,
                k_msg=args.k,
                sigma=sigma,
                delta=args.delta,
                h_lo=1,
                h_hi=bracket_h_hi,
                l=fci_l,
                outer_mode=args.outer_smallh_mode,
                parity_n=parity_n,
            )
            if log2_s_lin != float("-inf"):
                log2_fci_total = log2add(log2_fci_total, log2_s_lin)
            if log2_s_top != float("-inf"):
                log2_fci_total = log2add(log2_fci_total, log2_s_top)

        if args.show_fixedtap_cover_block:
            fcib_l = math.ceil((2.0 + args.small_xi) * args.delta * n)
            log2_fcib_total, fcib_peak_h, fcib_peak_log2 = fixedtap_cover_internal_block_model_log2(
                n=n,
                k_msg=args.k,
                sigma=sigma,
                delta=args.delta,
                h_lo=1,
                h_hi=bracket_h_hi,
                l=fcib_l,
                outer_mode=args.outer_smallh_mode,
                parity_n=parity_n,
                block_ratio=args.fixedtap_cover_block_ratio,
            )
            if log2_s_lin != float("-inf"):
                log2_fcib_total = log2add(log2_fcib_total, log2_s_lin)
            if log2_s_top != float("-inf"):
                log2_fcib_total = log2add(log2_fcib_total, log2_s_top)

        if args.show_latebase_pairwise:
            lb_l = math.ceil(2.0 * args.delta * n)
            log2_lb_total, lb_peak_h, lb_peak_log2 = latebase_pairwise_model_log2(
                n=n,
                k_msg=args.k,
                sigma=sigma,
                h_lo=1,
                h_hi=bracket_h_hi,
                l=lb_l,
                beta=beta_late,
                outer_mode=args.outer_smallh_mode,
                parity_n=parity_n,
            )
            if log2_s_lin != float("-inf"):
                log2_lb_total = log2add(log2_lb_total, log2_s_lin)
            if log2_s_top != float("-inf"):
                log2_lb_total = log2add(log2_lb_total, log2_s_top)

        if args.show_isolated_model:
            log2_i_total, i_peak_h, i_peak_log2 = isolated_slice_model_log2(
                n=n,
                k_msg=args.k,
                sigma=sigma,
                h_lo=1,
                h_hi=bracket_h_hi,
                l=iso_l,
                beta_bits_per_pair=0.0,
                add_adjacency_remainder=False,
                outer_mode=args.outer_smallh_mode,
                parity_n=parity_n,
            )
            if log2_s_lin != float("-inf"):
                log2_i_total = log2add(log2_i_total, log2_s_lin)
            if log2_s_top != float("-inf"):
                log2_i_total = log2add(log2_i_total, log2_s_top)

        if args.show_isolated_safe:
            log2_is_total, is_peak_h, is_peak_log2 = isolated_slice_model_log2(
                n=n,
                k_msg=args.k,
                sigma=sigma,
                h_lo=1,
                h_hi=bracket_h_hi,
                l=iso_l,
                beta_bits_per_pair=beta_iso_bits,
                add_adjacency_remainder=False,
                outer_mode=args.outer_smallh_mode,
                parity_n=parity_n,
            )
            if log2_s_lin != float("-inf"):
                log2_is_total = log2add(log2_is_total, log2_s_lin)
            if log2_s_top != float("-inf"):
                log2_is_total = log2add(log2_is_total, log2_s_top)

        if args.show_isolated_safe_adj:
            log2_isa_total, isa_peak_h, isa_peak_log2 = isolated_slice_model_log2(
                n=n,
                k_msg=args.k,
                sigma=sigma,
                h_lo=1,
                h_hi=bracket_h_hi,
                l=iso_l,
                beta_bits_per_pair=beta_iso_bits,
                add_adjacency_remainder=True,
                outer_mode=args.outer_smallh_mode,
                parity_n=parity_n,
            )
            if log2_s_lin != float("-inf"):
                log2_isa_total = log2add(log2_isa_total, log2_s_lin)
            if log2_s_top != float("-inf"):
                log2_isa_total = log2add(log2_isa_total, log2_s_top)

        if args.show_isolated_pos_pair:
            log2_ipp_total, ipp_peak_h, ipp_peak_log2 = isolated_positional_pair_model_log2(
                n=n,
                k_msg=args.k,
                sigma=sigma,
                delta=args.delta,
                adjacency_xi=args.pos_pair_adjacency_xi,
                h_lo=1,
                h_hi=bracket_h_hi,
                pair_a_bits=args.pos_pair_a_bits,
                pair_b_bits=args.pos_pair_b_bits,
                entropy_tail_bound=not args.pos_pair_exact_binomial,
                clamp_correction_nonnegative=not args.pos_pair_allow_negative_correction,
                adjacency_mode=args.pos_pair_adjacency_mode,
                outer_mode=args.outer_smallh_mode,
                parity_n=parity_n,
            )
            if log2_s_lin != float("-inf"):
                log2_ipp_total = log2add(log2_ipp_total, log2_s_lin)
            if log2_s_top != float("-inf"):
                log2_ipp_total = log2add(log2_ipp_total, log2_s_top)

        if args.show_isolated_hybrid:
            log2_ih_total, log2_ih_iso, ih_peak_h, ih_peak_log2, ih_h0 = isolated_window_model_log2(
                n=n,
                k_msg=args.k,
                sigma=sigma,
                delta=args.delta,
                z=args.z,
                rho=args.rho,
                xi=args.small_xi,
                kappa=args.kappa,
                theta=args.theta,
                eta_hi=args.eta_hi,
                gap_step=args.gap_step,
                h_iso_hi=args.isolated_h_cap,
                beta_bits_per_pair=beta_iso_bits,
                add_adjacency_remainder=True,
                outer_mode=args.outer_smallh_mode,
                parity_n=parity_n,
            )

        if args.show_split_safe:
            log2_ss_total, ss_peak_h, ss_peak_log2 = isolated_runmix_split_model_log2(
                n=n,
                k_msg=args.k,
                sigma=sigma,
                delta=args.delta,
                h_iso_hi=args.split_iso_h_cap,
                h_hi=bracket_h_hi,
                beta_bits_per_pair=beta_iso_bits,
                gamma_runmix_safe=gamma_runmix_safe,
                outer_mode=args.outer_smallh_mode,
                parity_n=parity_n,
            )
            if log2_s_lin != float("-inf"):
                log2_ss_total = log2add(log2_ss_total, log2_s_lin)
            if log2_s_top != float("-inf"):
                log2_ss_total = log2add(log2_ss_total, log2_s_top)

        if args.show_split_safe_hdep:
            # Isolated h-dependent control through split_iso_h_cap, then safe run-mixture above it.
            split_total = float("-inf")
            peak_h = 0
            peak_term = float("-inf")
            outer_logs = outer_small_h_prefix(
                args.k, sigma, bracket_h_hi, outer_mode=args.outer_smallh_mode, parity_n=parity_n
            )
            for h in range(1, bracket_h_hi + 1):
                out = outer_logs[h]
                if out == float("-inf"):
                    continue
                if h <= args.split_iso_h_cap:
                    p_full = isolated_full_run_prob(n, h)
                    beta_bits_per_pair = (
                        isolated_beta_hdep_bits_from_qh(q_single, h) + args.isolated_beta_hdep_safe_margin
                    )
                    p = p_full * late_start_base_prob(n, h, h, iso_l) * (
                        2.0 ** (beta_bits_per_pair * h * (h - 1) / 2.0)
                    )
                    p += 1.0 - p_full
                else:
                    p = 0.0
                    for s in range(1, h + 1):
                        p += run_count_weight(n, h, s) * (q_single**s) * math.exp(
                            -gamma_runmix_safe * s * (s - 1) / 2.0
                        )
                p = min(1.0, p)
                if p <= 0.0:
                    continue
                term = out + math.log2(p)
                split_total = term if split_total == float("-inf") else log2add(split_total, term)
                if term > peak_term:
                    peak_term = term
                    peak_h = h
            log2_ssh_total, ssh_peak_h, ssh_peak_log2 = split_total, peak_h, peak_term
            if log2_s_lin != float("-inf"):
                log2_ssh_total = log2add(log2_ssh_total, log2_s_lin)
            if log2_s_top != float("-inf"):
                log2_ssh_total = log2add(log2_ssh_total, log2_s_top)

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

        if args.show_run_mixture_safe_slack:
            log2_mss_total, mss_peak_h, mss_peak_log2 = run_mixture_model_log2(
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
                log2_mss_total = log2add(log2_mss_total, log2_s_lin)
            if log2_s_top != float("-inf"):
                log2_mss_total = log2add(log2_mss_total, log2_s_top)
            log2_mss_total += args.runmix_safe_slack_bits

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
        print(f"  offset sigma-ceil(log2 k): {sigma - math.ceil(math.log2(args.k))}")
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
        if args.show_late_base_model:
            print(f"  B_total (exact outer + exact late-placement base)        : {format_log2(log2_b_total)}")
            print(f"  B_peak  (dominant h, log2 contribution)                  : h={b_peak_h}, log2={b_peak_log2:.3f}")
        if args.show_one_episode_model:
            print(f"  E_total (late base + one-episode correction)             : {format_log2(log2_e_total)}")
            print(f"  E_peak  (dominant h, log2 contribution)                  : h={e_peak_h}, log2={e_peak_log2:.3f}")
            print(f"  E_L     ((2+xi)delta n window)                           : {e_l}")
        if args.show_one_episode_optimized:
            print(f"  EO_total (optimized one-episode correction)              : {format_log2(log2_eo_total)}")
            print(f"  EO_peak  (dominant h, log2 contribution)                 : h={eo_peak_h}, log2={eo_peak_log2:.3f}")
            print(f"  EO_L     (best window length)                            : {eo_l}")
        if args.show_fixedtap_multiepisode:
            print(f"  FME_total (fixed-tap multi-episode bound)                : {format_log2(log2_fme_total)}")
            print(f"  FME_peak  (dominant h, log2 contribution)                : h={fme_peak_h}, log2={fme_peak_log2:.3f}")
            print(f"  FME_L     (late-placement window)                        : {fme_l}")
            print(f"  FME_stat  (full-slice early span)                        : exact")
        if args.show_fixedtap_cover_internal:
            print(f"  FCI_total (cover + internal fair-tail diagnostic)        : {format_log2(log2_fci_total)}")
            print(f"  FCI_peak  (dominant h, log2 contribution)                : h={fci_peak_h}, log2={fci_peak_log2:.3f}")
            print(f"  FCI_L     (late-placement window)                        : {fci_l}")
            print(f"  FCI_stat  (cover-count envelope)                         : diagnostic")
        if args.show_fixedtap_cover_block:
            print(f"  FCIB_total (blocked cover theorem upper sum)             : {format_log2(log2_fcib_total)}")
            print(f"  FCIB_peak  (dominant h, log2 contribution)               : h={fcib_peak_h}, log2={fcib_peak_log2:.3f}")
            print(f"  FCIB_L     (late-placement window)                       : {fcib_l}")
            print(f"  FCIB_ratio (multiplicative T-block ratio)                : {args.fixedtap_cover_block_ratio}")
        if args.show_latebase_pairwise:
            print(f"  LB_total (late base + pairwise residual)                 : {format_log2(log2_lb_total)}")
            print(f"  LB_peak  (dominant h, log2 contribution)                 : h={lb_peak_h}, log2={lb_peak_log2:.3f}")
            print(f"  LB_L     (late-placement window)                         : {lb_l}")
            print(f"  LB_beta  (nat-log per run pair)                          : {beta_late:.6f}")
        if args.show_isolated_model:
            print(f"  I_total  (isolated slice late-placement base)            : {format_log2(log2_i_total)}")
            print(f"  I_peak   (dominant h, log2 contribution)                 : h={i_peak_h}, log2={i_peak_log2:.3f}")
            print(f"  I_L      (late-placement window)                         : {iso_l}")
        if args.show_isolated_safe:
            print(f"  IS_total (isolated slice + safe residual)                : {format_log2(log2_is_total)}")
            print(f"  IS_peak  (dominant h, log2 contribution)                 : h={is_peak_h}, log2={is_peak_log2:.3f}")
            print(f"  IS_beta  (safe residual bits per pair)                   : {beta_iso_bits:.6f}")
        if args.show_isolated_safe_adj:
            print(f"  ISA_total (isolated-safe + adjacency remainder)          : {format_log2(log2_isa_total)}")
            print(f"  ISA_peak  (dominant h, log2 contribution)                : h={isa_peak_h}, log2={isa_peak_log2:.3f}")
        if args.show_isolated_pos_pair:
            print(f"  IPP_total (positional + early-pair residual)             : {format_log2(log2_ipp_total)}")
            print(f"  IPP_peak  (dominant h, log2 contribution)                : h={ipp_peak_h}, log2={ipp_peak_log2:.3f}")
            print(f"  IPP_corr  (a + b E[pairs], bits)                         : a={args.pos_pair_a_bits:.6f}, b={args.pos_pair_b_bits:.6f}")
            print(f"  IPP_adj   (non-isolated remainder mode)                  : {args.pos_pair_adjacency_mode}")
            if args.pos_pair_adjacency_mode == "one-episode":
                print(f"  IPP_adj_xi                                                : {args.pos_pair_adjacency_xi:.6f}")
        if args.show_isolated_hybrid:
            print(f"  IH_total (isolated tiny slice + coarse handoff)          : {format_log2(log2_ih_total)}")
            print(f"  IH_iso   (isolated-controlled tiny slice only)           : {format_log2(log2_ih_iso)}")
            print(f"  IH_peak  (dominant isolated h, log2 contribution)        : h={ih_peak_h}, log2={ih_peak_log2:.3f}")
            print(f"  IH_h0    (paper tiny-window cutoff)                      : {ih_h0}")
        if args.show_split_safe:
            print(f"  SS_total (isolated tiny h + safe run-mixture above)      : {format_log2(log2_ss_total)}")
            print(f"  SS_peak  (dominant h, log2 contribution)                 : h={ss_peak_h}, log2={ss_peak_log2:.3f}")
            print(f"  SS_hiso  (isolated slice used through h)                 : {args.split_iso_h_cap}")
        if args.show_split_safe_hdep:
            print(f"  SSH_total (split-safe with h-dependent isolated fit)     : {format_log2(log2_ssh_total)}")
            print(f"  SSH_peak  (dominant h, log2 contribution)                : h={ssh_peak_h}, log2={ssh_peak_log2:.3f}")
            print(f"  SSH_hiso  (isolated slice used through h)                : {args.split_iso_h_cap}")
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
        if args.show_run_mixture_safe_slack:
            print(f"  MSS_total (safe pairwise + slack)                        : {format_log2(log2_mss_total)}")
            print(f"  MSS_peak  (dominant h, log2 contribution)                : h={mss_peak_h}, log2={mss_peak_log2 + args.runmix_safe_slack_bits:.3f}")
            print(f"  MSS_slack (added log2 slack)                             : {args.runmix_safe_slack_bits:.3f}")
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
