#!/usr/bin/env python3
"""Coverage diagnostic for repeated EBCH128 plus RandomStepConv.

For Q active outer rows, let J be the number of BCH coordinate regions that
contain at least one input difference.  The exact BCH weight spectrum gives
the probability that all Q permuted supports lie in any fixed j-set.  A
positive union bound controls the rare event J<r.  On J>=r, delete all but r
covered regions and all but one difference in each retained region.  The
RandomStepConv Chernoff moment can only increase under these deletions.

The calculation is algebraically rigorous, but binary64 evaluation makes the
output a diagnostic rather than an outward-rounded certificate.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

import evaluate_ebch128_randomstepconv_g1 as base
from evaluate_ebch128_randomstepconv_q1_exact import uniform_coefficients


WORKSTREAM = Path(__file__).resolve().parent
OUTPUT = WORKSTREAM / "ebch128_randomstepconv_g1_s30_coverage_d11_diagnostic.json"


def log_choose_scalar(n: int, k: int) -> float:
    if k < 0 or k > n:
        return -math.inf
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def containment_logs(spectrum: dict[int, int]) -> np.ndarray:
    """Return log Pr[support(C) subset U_j] for a uniform nonzero word."""
    result = np.full(base.B + 1, -math.inf)
    log_nonzero_mass = math.log((1 << base.K) - 1)
    for j in range(base.B + 1):
        terms = []
        for weight, count in spectrum.items():
            if weight == 0 or weight > j or count == 0:
                continue
            terms.append(
                math.log(count)
                + log_choose_scalar(j, weight)
                - log_choose_scalar(base.B, weight)
            )
        if terms:
            result[j] = float(logsumexp(terms)) - log_nonzero_mass
    return result


def retained_region_moments(z: float, memory_bits: int) -> np.ndarray:
    """Moments after retaining exactly r singly-active coordinate regions."""
    zero, active = base.step_matrices(z, memory_bits)
    bit_regions = uniform_coefficients(
        base.log_entries(zero), base.log_entries(active), base.L, 1
    )
    coordinate_products = uniform_coefficients(
        bit_regions[0], bit_regions[1], base.B, base.B
    )
    return np.logaddexp(
        coordinate_products[:, 0, 0], coordinate_products[:, 0, 1]
    )


def coverage_moment_bounds(
    log_moments: np.ndarray,
    log_containment: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Bound E[z^W | Q active rows] and return optimizing coverage r."""
    q_values = np.arange(1, base.ACTIVE_ROWS + 1, dtype=np.float64)
    best = np.full(base.ACTIVE_ROWS, math.inf)
    best_r = np.zeros(base.ACTIVE_ROWS, dtype=np.int64)
    rare_prefix = np.full(base.ACTIVE_ROWS, -math.inf)

    # At threshold r, rare_prefix contains the sum over J<r.  The event
    # J>=r contributes at most the r-retained-region moment.
    for r in range(base.B + 1):
        candidate = np.logaddexp(rare_prefix, log_moments[r])
        improved = candidate < best
        best[improved] = candidate[improved]
        best_r[improved] = r

        if r < base.B and math.isfinite(log_containment[r]):
            log_event_upper = (
                log_choose_scalar(base.B, r) + q_values * log_containment[r]
            )
            rare_prefix = np.logaddexp(
                rare_prefix, log_event_upper + log_moments[r]
            )

    return np.minimum(best, 0.0), best_r


def exact_q1_moment(
    log_moments: np.ndarray, spectrum: dict[int, int]
) -> float:
    log_mass = math.log((1 << base.K) - 1)
    return float(
        logsumexp(
            [
                math.log(count) - log_mass + log_moments[weight]
                for weight, count in spectrum.items()
                if weight > 0 and count
            ]
        )
    )


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    spectrum = base.load_spectrum(args.spectrum)
    log_containment = containment_logs(spectrum)
    occupations = np.arange(1, base.ACTIVE_ROWS + 1, dtype=np.int64)
    outer_logs = base.log_choose(base.ACTIVE_ROWS, occupations) + occupations * math.log(
        (1 << base.K) - 1
    )
    best = np.full(base.ACTIVE_ROWS, math.inf)
    witnesses = np.full(base.ACTIVE_ROWS, math.nan)
    thresholds = np.zeros(base.ACTIVE_ROWS, dtype=np.int64)

    u_values = np.asarray(
        sorted(
            set(float(x) for x in np.arange(args.grid_min, args.grid_max + args.grid_step / 2, args.grid_step))
            | set(float(x) for x in np.arange(args.refine_min, args.refine_max + args.refine_step / 2, args.refine_step))
        )
    )
    for index, u in enumerate(u_values):
        surprisal = math.exp(float(u))
        z = math.exp(-surprisal)
        log_moments = retained_region_moments(z, args.memory_bits)
        conditional, selected_r = coverage_moment_bounds(log_moments, log_containment)
        conditional[0] = exact_q1_moment(log_moments, spectrum)
        selected_r[0] = -1
        inner = np.minimum(0.0, conditional + base.D * surprisal)
        values = outer_logs + inner
        improved = values < best
        best[improved] = values[improved]
        witnesses[improved] = float(u)
        thresholds[improved] = selected_r[improved]
        print(f"tilt,{index + 1},{len(u_values)},u,{u:.6f}", flush=True)

    aggregate = float(logsumexp(best))
    dominant = int(np.argmax(best))
    top = np.argsort(best)[-20:][::-1]
    return {
        "schema": "ebch128-randomstepconv-coverage-diagnostic-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "claim": {
            "log2_expected_bad_upper": aggregate / base.LOG2,
            "margin_bits": -aggregate / base.LOG2,
            "closes": aggregate < 0.0,
            "closes_40_bits": aggregate < -40.0 * base.LOG2,
            "dominant_occupation": dominant + 1,
            "dominant_pointwise_log2_upper": float(best[dominant]) / base.LOG2,
        },
        "parameters": {
            "outer_code": "extended BCH [128,64,22]",
            "outer_rows": base.L,
            "active_outer_rows": base.ACTIVE_ROWS,
            "message_bits": base.K * base.ACTIVE_ROWS,
            "output_bits": base.N,
            "distance_cutoff": base.D,
            "memory_bits": args.memory_bits,
        },
        "probability_space": (
            "one fixed repeated EBCH constituent; independent uniform outer-row "
            "coordinate permutations and region permutations; independent "
            "RandomStepConv maps sampled once and shared by every message"
        ),
        "method": (
            "exact Q=1 spectrum sum; for Q>=2, a spectrum-derived support-"
            "containment union bound plus monotone deletion to an optimized "
            "number of singly-active coordinate regions"
        ),
        "top_occupations": [
            {
                "active_outer_rows": int(q + 1),
                "pointwise_log2_upper": float(best[q]) / base.LOG2,
                "log_surprisal": float(witnesses[q]),
                "coverage_threshold": int(thresholds[q]),
            }
            for q in top
        ],
        "occupation_rows": [
            {
                "active_outer_rows": int(q + 1),
                "pointwise_log2_upper": float(best[q]) / base.LOG2,
                "log_surprisal": float(witnesses[q]),
                "coverage_threshold": int(thresholds[q]),
            }
            for q in range(base.ACTIVE_ROWS)
        ],
        "containment_log2": [
            None if not math.isfinite(value) else float(value / base.LOG2)
            for value in log_containment
        ],
        "limitations": [
            "Nearest binary64 arithmetic is not an outward certificate.",
            "The RandomStepConv monotone-deletion and coverage lemmas still require prose proof and independent review.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectrum", type=Path, default=base.SPECTRUM)
    parser.add_argument("--memory-bits", type=int, default=base.MEMORY)
    parser.add_argument("--grid-min", type=float, default=-12.0)
    parser.add_argument("--grid-max", type=float, default=1.0)
    parser.add_argument("--grid-step", type=float, default=0.5)
    parser.add_argument("--refine-min", type=float, default=0.65)
    parser.add_argument("--refine-max", type=float, default=0.80)
    parser.add_argument("--refine-step", type=float, default=0.025)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    payload = evaluate(args)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
