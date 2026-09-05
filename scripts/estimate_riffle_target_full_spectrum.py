#!/usr/bin/env python3
"""Heuristic target-size spectrum estimate for the riffle construction.

This intentionally works at the production parameters.  It sums the expected
low-weight count over outer field-symbol weights and uses a saddle approximation
for the parallel accumulator's random gap lengths.  The result is an engineering
estimate, not a distance certificate.
"""

from __future__ import annotations

import argparse
import math

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.special import log_ndtr, logsumexp


DATA_SYMBOLS = 16_384
FIELD_SIZE = 1 << 64
PACKETS_PER_BCH_BLOCK = 32
STATE_WEIGHT_LAW = np.asarray((1, 4, 6, 4, 1), dtype=np.float64) / 16.0


def mds_weight_log(symbol_weight: int, parity_symbols: int) -> float:
    """Natural log of the [n,k,p+1] MDS symbol-weight enumerator."""
    t = symbol_weight
    outer_symbols = DATA_SYMBOLS + parity_symbols
    distance = parity_symbols + 1
    # Factor q^(t-3) out of the alternating sum.  At q=2^64 the correction
    # differs from one by less than floating-point epsilon in our t range.
    correction = math.fsum(
        (-1.0 if j & 1 else 1.0)
        * math.comb(t - 1, j)
        * FIELD_SIZE ** (-j)
        for j in range(t - distance + 1)
    )
    return (
        math.lgamma(outer_symbols + 1)
        - math.lgamma(t + 1)
        - math.lgamma(outer_symbols - t + 1)
        + math.log(FIELD_SIZE - 1)
        + (t - distance) * math.log(FIELD_SIZE)
        + math.log(correction)
    )


def spacing_tail_log(relative_weight: float, active_packets: int) -> float:
    """Dirichlet-spacing saddle estimate for accumulator weight <= delta*N."""
    threshold = 4.0 * relative_weight
    coefficients = np.arange(5, dtype=np.float64) - threshold

    def cumulants(tilt: float) -> tuple[float, float, float]:
        denominators = 1.0 - tilt * coefficients
        moment = float(np.sum(STATE_WEIGHT_LAW / denominators))
        first = float(np.sum(STATE_WEIGHT_LAW * coefficients / denominators**2))
        second = float(
            np.sum(2.0 * STATE_WEIGHT_LAW * coefficients**2 / denominators**3)
        )
        cumulant = -math.log1p(tilt * threshold) + active_packets * math.log(moment)
        derivative = (
            -threshold / (1.0 + tilt * threshold)
            + active_packets * first / moment
        )
        curvature = (
            threshold**2 / (1.0 + tilt * threshold) ** 2
            + active_packets * (second / moment - (first / moment) ** 2)
        )
        return cumulant, derivative, curvature

    tilt = brentq(
        lambda value: cumulants(value)[1],
        -1.0 / threshold + 1e-12,
        -1e-15,
        xtol=1e-14,
        rtol=1e-14,
    )
    cumulant, _derivative, curvature = cumulants(tilt)
    w = -math.sqrt(-2.0 * cumulant)
    u = tilt * math.sqrt(curvature)
    log_normal = float(log_ndtr(w))
    mills = math.exp(-0.5 * w * w - 0.5 * math.log(2.0 * math.pi) - log_normal)
    correction = 1.0 + mills * (1.0 / w - 1.0 / u)
    return log_normal + (math.log(correction) if correction > 0.0 else 0.0)


def support_log_probability(symbol_weight: int, support: int) -> float:
    """Binomial packet-support proxy for t independently permuted BCH blocks."""
    trials = PACKETS_PER_BCH_BLOCK * symbol_weight
    return (
        math.lgamma(trials + 1)
        - math.lgamma(support + 1)
        - math.lgamma(trials - support + 1)
        + support * math.log(15.0)
        - trials * math.log(16.0)
    )


def inner_occupation_log(relative_weight: float, symbol_weight: int, window: int) -> tuple[float, int]:
    trials = PACKETS_PER_BCH_BLOCK * symbol_weight

    def term(support: int) -> float:
        return support_log_probability(symbol_weight, support) + spacing_tail_log(
            relative_weight, support
        )

    optimum = minimize_scalar(
        lambda value: -term(max(1, min(trials, int(round(value))))),
        bounds=(1, trials),
        method="bounded",
        options={"xatol": 0.25},
    )
    peak = max(1, min(trials, int(round(optimum.x))))
    low = max(1, peak - window)
    high = min(trials, peak + window)
    mixture = float(logsumexp([term(support) for support in range(low, high + 1)]))
    return mixture, peak


def total_logs(
    relative_weight: float,
    parity_counts: list[int],
    max_symbol_weight: int,
    window: int,
) -> dict[int, tuple[float, list[tuple[int, int, float]]]]:
    # The expensive accumulator calculation is common to every outer parity
    # count, so compute it once per (delta,t) and reuse it across the sweep.
    inner = {
        symbol_weight: inner_occupation_log(relative_weight, symbol_weight, window)
        for symbol_weight in range(min(parity_counts) + 1, max_symbol_weight + 1)
    }
    results = {}
    for parity_symbols in parity_counts:
        rows = []
        for symbol_weight in range(parity_symbols + 1, max_symbol_weight + 1):
            mixture, peak = inner[symbol_weight]
            value = mds_weight_log(symbol_weight, parity_symbols) + mixture
            rows.append((symbol_weight, peak, value / math.log(2.0)))
        total = float(logsumexp([row[2] * math.log(2.0) for row in rows])) / math.log(2.0)
        results[parity_symbols] = (
            total,
            sorted(rows, key=lambda row: row[2], reverse=True)[:5],
        )
    return results


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("deltas", nargs="+", type=float)
    parser.add_argument("--parity-symbols", nargs="+", type=int, default=[2])
    parser.add_argument("--max-symbol-weight", type=int, default=300)
    parser.add_argument("--support-window", type=int, default=180)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    print(f"message_bits,{DATA_SYMBOLS * 64}")
    print("evidence_label,NONRIGOROUS_FULL_SPECTRUM_MODEL")
    for delta in args.deltas:
        results = total_logs(
            delta,
            args.parity_symbols,
            args.max_symbol_weight,
            args.support_window,
        )
        for parity_symbols in args.parity_symbols:
            total, top = results[parity_symbols]
            binary_length = 128 * (DATA_SYMBOLS + parity_symbols)
            print(
                f"parity_symbols,{parity_symbols},delta,{delta:.9f},"
                f"binary_length,{binary_length},"
                f"distance,{math.floor(delta * binary_length)},"
                f"total_log2,{total:.12f},top,{top}"
            )


if __name__ == "__main__":
    main()
