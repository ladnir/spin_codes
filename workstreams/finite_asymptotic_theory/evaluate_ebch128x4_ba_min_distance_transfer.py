#!/usr/bin/env python3
"""Minimum-distance-only transfer for EBCH128x4 + BA-t.

For each BA stage, Markov's inequality turns its expected spectrum into a
high-probability minimum-distance event.  Conditional on that event, this
program deliberately forgets the rest of the spectrum.  It deletes each
nonzero row to a uniform minimum-weight subset and applies a Bernoulli change
of measure before the structured transpose and RandomStepConv.

This is a nearest-binary64 diagnostic.  It tests whether minimum distance
alone is strong enough to avoid the unresolved fixed-spectrum concentration
problem.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

from evaluate_repeated_random512_randomstepconv_g1 import (
    log_choose,
    log_matmul_batch,
    log_matrix_entries,
    step_matrices,
)


WORKSTREAM = Path(__file__).resolve().parent
SPECTRA = WORKSTREAM / "ebch128x4_ba0_6_B512_expected_spectra.json"
OUTPUT = WORKSTREAM / "ebch128x4_ba0_6_B512_min_distance_transfer_m22.json"
B = 512
K = 256
L = 4096
N = 1 << 21
D = (109 * N + 999) // 1000
LOG2 = math.log(2.0)
MESSAGE_LOG = K * LOG2 + math.log1p(-math.ldexp(1.0, -K))


def region_coefficients(zero: np.ndarray, candidate: np.ndarray) -> np.ndarray:
    log_zero = log_matrix_entries(tuple(zero))
    log_candidate = log_matrix_entries(tuple(candidate))
    current = np.full((L + 1, 2, 2), -math.inf)
    current[0] = np.eye(2)
    current[0, 0, 0] = 0.0
    current[0, 1, 1] = 0.0
    current[0, 0, 1] = -math.inf
    current[0, 1, 0] = -math.inf
    for completed in range(L):
        maximum = completed + 1
        updated = np.full_like(current, -math.inf)
        zero_products = log_matmul_batch(current[:maximum], log_zero)
        degrees = np.arange(maximum, dtype=np.float64)
        updated[:maximum] = zero_products + np.log(
            (completed + 1.0 - degrees) / (completed + 1.0)
        )[:, None, None]
        candidate_products = log_matmul_batch(current[:maximum], log_candidate)
        selected = np.arange(1, maximum + 1, dtype=np.float64)
        candidate_terms = candidate_products + np.log(
            selected / (completed + 1.0)
        )[:, None, None]
        updated[1 : maximum + 1] = np.logaddexp(
            updated[1 : maximum + 1], candidate_terms
        )
        current = updated
    return current


def logdiffexp_array(high: np.ndarray, low: np.ndarray) -> np.ndarray:
    if np.any(high <= low):
        raise ArithmeticError("positive activation mass was lost")
    return high + np.log(-np.expm1(low - high))


def nonzero_moments(regions: np.ndarray, probability: float) -> np.ndarray:
    """Power all region matrices while excluding the global zero input."""
    occupations = np.arange(1, L + 1, dtype=np.float64)
    all_zero_region = occupations * math.log1p(-probability)
    activate_zero = logdiffexp_array(regions[1:, 0, 0], all_zero_region)
    activate_live = regions[1:, 0, 1]
    inactive = np.zeros(L)
    active_zero = np.full(L, -math.inf)
    active_live = np.full(L, -math.inf)
    for _ in range(B):
        next_zero = np.logaddexp.reduce(
            np.stack(
                (
                    inactive + activate_zero,
                    active_zero + regions[1:, 0, 0],
                    active_live + regions[1:, 1, 0],
                )
            ),
            axis=0,
        )
        next_live = np.logaddexp.reduce(
            np.stack(
                (
                    inactive + activate_live,
                    active_zero + regions[1:, 0, 1],
                    active_live + regions[1:, 1, 1],
                )
            ),
            axis=0,
        )
        inactive += all_zero_region
        active_zero, active_live = next_zero, next_live
    return np.logaddexp(active_zero, active_live)


def stage_spectra() -> dict[int, np.ndarray]:
    source = json.loads(SPECTRA.read_text(encoding="utf-8"))
    result = {}
    for stage in source["stages"]:
        values = np.full(B + 1, -math.inf)
        for row in stage["spectrum"]:
            if row["log2_expected_multiplicity"] is not None:
                values[int(row["weight"])] = float(row["log2_expected_multiplicity"]) * LOG2
        result[int(stage["accumulators"])] = values
    return result


def log_cumulative(values: np.ndarray, cutoff: int) -> float:
    selected = values[1 : cutoff + 1]
    finite = np.isfinite(selected)
    if not finite.any():
        return -math.inf
    return float(logsumexp(selected[finite]))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--memory-bits", type=int, default=22)
    parser.add_argument("--event-failure-bits", type=float, default=43.0)
    parser.add_argument("--stages", type=int, nargs="+", default=[0, 1, 2, 3, 4, 5, 6])
    parser.add_argument("--p-values", type=float, nargs="+", default=[0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.3, 0.4, 0.5])
    parser.add_argument("--grid-min", type=float, default=-12.0)
    parser.add_argument("--grid-max", type=float, default=1.0)
    parser.add_argument("--grid-step", type=float, default=0.5)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    spectra = stage_spectra()
    stage_rows = {}
    for stage in args.stages:
        spectrum = spectra[stage]
        cutoff = 0
        for candidate in range(1, B + 1):
            if log_cumulative(spectrum, candidate) / LOG2 <= -args.event_failure_bits:
                cutoff = candidate
        minimum_distance = cutoff + 1
        failure_log = log_cumulative(spectrum, cutoff)
        stage_rows[stage] = {
            "accumulators": stage,
            "minimum_distance_event": minimum_distance,
            "event_failure_log2_upper": failure_log / LOG2,
            "best": np.full(L + 1, math.inf),
            "witnesses": [None] * (L + 1),
        }

    u_values = np.arange(args.grid_min, args.grid_max + args.grid_step / 2, args.grid_step)
    total_pairs = len(args.p_values) * len(u_values)
    completed = 0
    for probability in args.p_values:
        if not 0 < probability < 1:
            parser.error("each Bernoulli probability must lie in (0,1)")
        for u_raw in u_values:
            completed += 1
            u = float(u_raw)
            surprisal = math.exp(u)
            z = math.exp(-surprisal)
            zero, active = step_matrices(z, args.memory_bits)
            candidate = (1.0 - probability) * zero + probability * active
            regions = region_coefficients(zero, candidate)
            moments = nonzero_moments(regions, probability)
            inner = np.minimum(0.0, moments + D * surprisal)
            occupations = np.arange(1, L + 1)
            row_choices = np.asarray([log_choose(L, int(q)) for q in occupations])
            for stage, row in stage_rows.items():
                minimum_distance = int(row["minimum_distance_event"])
                log_conditioning = (
                    log_choose(B, minimum_distance)
                    + minimum_distance * math.log(probability)
                    + (B - minimum_distance) * math.log1p(-probability)
                )
                per_row = MESSAGE_LOG - log_conditioning
                values = row_choices + occupations * per_row + inner
                best = row["best"]
                improved = values < best[1:]
                best[1:][improved] = values[improved]
                for q in occupations[improved]:
                    row["witnesses"][int(q)] = {
                        "candidate_probability": probability,
                        "log_surprisal": u,
                    }
            print(f"pair,{completed},{total_pairs},p,{probability:.6f},u,{u:.6f}", flush=True)

    rows = []
    for stage in args.stages:
        row = stage_rows[stage]
        best = row.pop("best")
        witnesses = row.pop("witnesses")
        aggregate = float(logsumexp(best[1:]))
        dominant_q = int(np.argmax(best[1:])) + 1
        event_log2 = float(row["event_failure_log2_upper"])
        joint_log2 = float(np.logaddexp(aggregate, event_log2 * LOG2) / LOG2)
        rows.append(
            {
                **row,
                "conditional_distance_log2_upper": aggregate / LOG2,
                "conditional_distance_margin_bits": -aggregate / LOG2,
                "combined_log2_upper": joint_log2,
                "combined_margin_bits": -joint_log2,
                "dominant_occupation": dominant_q,
                "dominant_pointwise_log2_upper": float(best[dominant_q] / LOG2),
                "dominant_witness": witnesses[dominant_q],
            }
        )

    payload = {
        "schema": "ebch128x4-ba-min-distance-randomstepconv-v1",
        "status": "BINARY64_MINIMUM_DISTANCE_WRAPPER_DIAGNOSTIC",
        "parameters": {
            "outer": "four EBCH [128,64,22] blocks followed by the stated number of uniform-interleaved accumulators",
            "outer_block_bits": B,
            "outer_dimension": K,
            "outer_rows": L,
            "output_bits": N,
            "distance_cutoff": D,
            "memory_bits": args.memory_bits,
            "event_failure_target_bits": args.event_failure_bits,
        },
        "rows": rows,
        "proof_model": {
            "event": "the sampled BA constituent has no nonzero word below the displayed minimum distance",
            "event_bound": "Markov applied to the expected cumulative BA spectrum",
            "conditional_transfer": "delete every nonzero row to a uniform minimum-weight subset, then use a Bernoulli change of measure and the exact structured region permutation",
        },
        "limitations": [
            "Nearest binary64 arithmetic makes all distance-transfer values diagnostic.",
            "The event probability uses the binary64 expected BA spectrum and must be outward-certified if this wrapper closes.",
            "A failure here rejects only the minimum-distance-only wrapper; it does not reject the BA constituent.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(rows, indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
