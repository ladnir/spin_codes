#!/usr/bin/env python3
"""Parity-pivot diagnostic for repeated EBCH128 and RandomStepConv.

The ordinary EBCH words are dominated by a uniform-even reference.  Each
reference row samples an independent parity pivot.  After deleting the pivot
input, the remaining 127 coordinates are independent fair candidates.  The
pivot loads across the 128 coordinate regions are multinomial.

For occupations Q>=65, two positive coefficient fugacities bound the exact
multinomial matrix transfer.  Occupation one uses the exact BCH spectrum.
Occupations 2 through 64 use the exact fixed-occupation region recurrence and
the valid aligned-pivot deletion bound.  Binary64 arithmetic makes the result
a diagnostic rather than an outward certificate.
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
OUTPUT = WORKSTREAM / "ebch128_randomstepconv_g1_s30_parity_pivot_d11_diagnostic.json"
LOW_MAX_Q = 64


def pivot_terms(q: int, log_r: float) -> np.ndarray:
    """Logs of r^m/(m! binom(L,q-m)), for 0<=m<=q."""
    terms = np.empty(q + 1, dtype=np.float64)
    terms[0] = -(
        math.lgamma(base.L + 1)
        - math.lgamma(q + 1)
        - math.lgamma(base.L - q + 1)
    )
    if q:
        m = np.arange(q, dtype=np.float64)
        increments = (
            log_r
            + np.log(base.L - q + m + 1.0)
            - np.log(m + 1.0)
            - np.log(q - m)
        )
        terms[1:] = terms[0] + np.cumsum(increments)
    return terms


def pivot_saddle(q: int) -> tuple[float, float]:
    """Choose a positive pivot fugacity and return its valid log bound."""
    mean_target = q / float(base.B)
    initial_r = (
        (mean_target + 1.0)
        * (q - mean_target)
        / (base.L - q + mean_target + 1.0)
    )
    log_r = math.log(max(initial_r, 1e-300))

    # The elementary value is sufficient unless the endpoint m=q overtakes
    # the central mode. The crossover depends on L, so detect it from the
    # coefficient terms instead of using an L-specific occupation cutoff.
    # Bisection is stable across this crossover.
    initial_terms = pivot_terms(q, log_r)
    if int(np.argmax(initial_terms)) > q // 2:
        low = log_r - 4.0
        high = log_r + 1.0
        indices = np.arange(q + 1, dtype=np.float64)
        for _ in range(40):
            middle = (low + high) / 2.0
            terms = pivot_terms(q, middle)
            normalized = np.exp(terms - float(logsumexp(terms)))
            mean = float(np.dot(normalized, indices))
            if mean < mean_target:
                low = middle
            else:
                high = middle
        log_r = (low + high) / 2.0

    terms = pivot_terms(q, log_r)
    log_g = float(logsumexp(terms))
    cost = base.B * log_g - q * log_r
    return log_r, cost


def precompute_pivot_costs() -> tuple[np.ndarray, np.ndarray]:
    log_r = np.full(base.ACTIVE_ROWS + 1, math.nan)
    cost = np.full(base.ACTIVE_ROWS + 1, math.nan)
    for q in range(LOW_MAX_Q + 1, base.ACTIVE_ROWS + 1):
        log_r[q], cost[q] = pivot_saddle(q)
        if q % 1024 == 0 or q == base.ACTIVE_ROWS:
            print(f"pivot,{q},{base.ACTIVE_ROWS}", flush=True)
    return log_r, cost


def q1_exact_moment(
    zero: np.ndarray,
    active: np.ndarray,
    spectrum: dict[int, int],
) -> float:
    bit_regions = uniform_coefficients(
        base.log_entries(zero), base.log_entries(active), base.L, 1
    )
    coordinate_products = uniform_coefficients(
        bit_regions[0], bit_regions[1], base.B, base.B
    )
    moments = np.logaddexp(
        coordinate_products[:, 0, 0], coordinate_products[:, 0, 1]
    )
    return float(
        logsumexp(
            [
                math.log(count) + moments[weight]
                for weight, count in spectrum.items()
                if weight > 0 and count
            ]
        )
    )


def dense_pivot_moments(
    zero: np.ndarray,
    candidate: np.ndarray,
    occupations: np.ndarray,
    pivot_cost: np.ndarray,
    v_offset: float,
) -> tuple[np.ndarray, np.ndarray]:
    q = occupations.astype(np.float64)
    mean_candidates = q * (base.B - 1.0) / base.B
    base_x = mean_candidates / (base.L - mean_candidates)
    scheduled_v = np.maximum(0.0, 0.60 * (1.0 - q / 16_000.0)) + v_offset
    x = base_x * np.exp(scheduled_v)
    polynomial = zero[None, :, :] + x[:, None, None] * candidate[None, :, :]
    region_power = base.log_power(base.log_entries(polynomial), base.L)
    all_regions = base.log_power(region_power, base.B)
    matrix_moment = np.logaddexp(all_regions[:, 0, 0], all_regions[:, 0, 1])
    conditional = (
        np.asarray([math.lgamma(int(value) + 1) for value in occupations])
        - q * math.log(base.B)
        + pivot_cost[occupations]
        + matrix_moment
        - q * (base.B - 1.0) * np.log(x)
    )
    return np.minimum(conditional, 0.0), scheduled_v


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    spectrum = base.load_spectrum(args.spectrum)
    envelope, maximizing_weight = base.body_envelope(spectrum)
    effective_mass = envelope + 1
    log_effective_mass = math.log(effective_mass.numerator) - math.log(
        effective_mass.denominator
    )
    occupations = np.arange(1, base.ACTIVE_ROWS + 1, dtype=np.int64)
    outer_logs = base.log_choose(base.ACTIVE_ROWS, occupations) + occupations * log_effective_mass
    dense_q = occupations[LOW_MAX_Q:]
    pivot_log_r, pivot_cost = precompute_pivot_costs()

    best = np.full(base.ACTIVE_ROWS, math.inf)
    witnesses = np.full(base.ACTIVE_ROWS, math.nan)
    witness_v = np.full(base.ACTIVE_ROWS, math.nan)
    u_values = np.arange(
        args.grid_min, args.grid_max + args.grid_step / 2.0, args.grid_step
    )
    for index, u in enumerate(u_values):
        surprisal = math.exp(float(u))
        z = math.exp(-surprisal)
        zero, active = base.step_matrices(z, args.memory_bits)
        candidate = 0.5 * (zero + active)

        # Q=1: exact BCH spectrum, including the unique all-one word.
        one_value = (
            math.log(base.ACTIVE_ROWS)
            + q1_exact_moment(zero, active, spectrum)
            + base.D * surprisal
        )
        if one_value < best[0]:
            best[0] = one_value
            witnesses[0] = float(u)
            witness_v[0] = math.nan

        # Q=2..64: exact region coefficients under the aligned-pivot
        # deletion bound.  Alignment is harmless in this sparse range.
        zero_region = base.final_zero_region(zero)
        sparse_regions = base.exact_region_coefficients(zero, candidate)[2 : LOW_MAX_Q + 1]
        sparse_moments = base.moments_from_regions(sparse_regions, zero_region)
        sparse_inner = np.minimum(0.0, sparse_moments + base.D * surprisal)
        sparse_values = outer_logs[1:LOW_MAX_Q] + sparse_inner
        sparse_improved = sparse_values < best[1:LOW_MAX_Q]
        best[1:LOW_MAX_Q][sparse_improved] = sparse_values[sparse_improved]
        witnesses[1:LOW_MAX_Q][sparse_improved] = float(u)

        # Q>=65: dispersed parity pivots.  Several fixed region-fugacity
        # schedules are evaluated; taking their minimum preserves validity.
        for v_offset in args.v_offsets:
            dense_conditional, scheduled_v = dense_pivot_moments(
                zero, candidate, dense_q, pivot_cost, v_offset
            )
            dense_inner = np.minimum(0.0, dense_conditional + base.D * surprisal)
            dense_values = outer_logs[LOW_MAX_Q:] + dense_inner
            improved = dense_values < best[LOW_MAX_Q:]
            best[LOW_MAX_Q:][improved] = dense_values[improved]
            witnesses[LOW_MAX_Q:][improved] = float(u)
            witness_v[LOW_MAX_Q:][improved] = scheduled_v[improved]

        print(f"tilt,{index + 1},{len(u_values)},u,{u:.6f}", flush=True)

    aggregate = float(logsumexp(best))
    dominant = int(np.argmax(best))
    top = np.argsort(best)[-20:][::-1]
    return {
        "schema": "ebch128-randomstepconv-parity-pivot-diagnostic-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "claim": {
            "log2_expected_bad_upper": aggregate / base.LOG2,
            "margin_bits": -aggregate / base.LOG2,
            "closes": aggregate < 0.0,
            "closes_20_bits": aggregate < -20.0 * base.LOG2,
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
            "exact_low_occupation_max": LOW_MAX_Q,
            "region_fugacity_offsets": args.v_offsets,
        },
        "outer_envelope": {
            "body_factor_log2": math.log2(envelope.numerator) - math.log2(envelope.denominator),
            "maximizing_body_weight": maximizing_weight,
            "effective_mass_log2": log_effective_mass / base.LOG2,
        },
        "probability_space": (
            "one fixed repeated EBCH constituent; independent uniform outer-row "
            "coordinate permutations and region permutations; independent "
            "RandomStepConv maps sampled once and shared by every message"
        ),
        "method": (
            "exact spectrum at Q=1; exact aligned-pivot region recurrence for "
            "2<=Q<=64; uniform-even pointwise majorant with independent parity "
            "pivots and two positive coefficient fugacities for Q>=65"
        ),
        "top_occupations": [
            {
                "active_outer_rows": int(q + 1),
                "pointwise_log2_upper": float(best[q]) / base.LOG2,
                "log_surprisal": float(witnesses[q]),
                "surprisal_hex": float(math.exp(witnesses[q])).hex(),
                "region_log_fugacity_offset": None if math.isnan(witness_v[q]) else float(witness_v[q]),
                "region_fugacity_hex": (
                    None
                    if math.isnan(witness_v[q])
                    else float(
                        ((q + 1) * (base.B - 1) / base.B)
                        / (base.L - ((q + 1) * (base.B - 1) / base.B))
                        * math.exp(witness_v[q])
                    ).hex()
                ),
                "pivot_fugacity_hex": (
                    None
                    if q + 1 <= LOW_MAX_Q
                    else float(math.exp(pivot_log_r[q + 1])).hex()
                ),
            }
            for q in top
        ],
        "occupation_rows": [
            {
                "active_outer_rows": int(q + 1),
                "pointwise_log2_upper": float(best[q]) / base.LOG2,
                "log_surprisal": float(witnesses[q]),
                "surprisal_hex": float(math.exp(witnesses[q])).hex(),
                "region_log_fugacity_offset": None if math.isnan(witness_v[q]) else float(witness_v[q]),
                "region_fugacity_hex": (
                    None
                    if math.isnan(witness_v[q])
                    else float(
                        ((q + 1) * (base.B - 1) / base.B)
                        / (base.L - ((q + 1) * (base.B - 1) / base.B))
                        * math.exp(witness_v[q])
                    ).hex()
                ),
                "pivot_fugacity_hex": (
                    None
                    if q + 1 <= LOW_MAX_Q
                    else float(math.exp(pivot_log_r[q + 1])).hex()
                ),
            }
            for q in range(base.ACTIVE_ROWS)
        ],
        "limitations": [
            "Nearest binary64 arithmetic is not an outward certificate.",
            "The parity-pivot and two-fugacity inequalities require independent proof review.",
            "The fixed fugacity schedules are witnesses, not claims of global numerical optimality.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectrum", type=Path, default=base.SPECTRUM)
    parser.add_argument("--outer-rows", type=int, default=base.L)
    parser.add_argument("--memory-bits", type=int, default=base.MEMORY)
    parser.add_argument("--distance-numerator", type=int, default=11)
    parser.add_argument("--distance-denominator", type=int, default=100)
    parser.add_argument("--grid-min", type=float, default=-10.0)
    parser.add_argument("--grid-max", type=float, default=1.0)
    parser.add_argument("--grid-step", type=float, default=0.25)
    parser.add_argument(
        "--v-offsets",
        type=float,
        nargs="+",
        default=[-0.10, 0.0, 0.10],
    )
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    base.configure_outer_rows(
        args.outer_rows,
        args.distance_numerator,
        args.distance_denominator,
    )
    payload = evaluate(args)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
