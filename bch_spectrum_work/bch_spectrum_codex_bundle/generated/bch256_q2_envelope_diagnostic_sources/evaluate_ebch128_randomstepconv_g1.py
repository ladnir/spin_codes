#!/usr/bin/env python3
"""Diagnostic for repeated EBCH[128,64,22] plus RandomStepConv.

The ordinary BCH words are pointwise dominated by a uniform-even reference.
The unique all-one word is treated as a separate row type. Monotonicity then
maps both row types to the same reference: 127 fair coordinates followed by
one zero coordinate. The last coordinate still emits a complete output
region; only its input bits are deleted.

Occupations through 64 use the exact uniform-subset region coefficient.
Larger occupations use a positive Bernoulli-conditioning upper bound.
Nearest binary64 arithmetic makes this a diagnostic, not a certificate.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp


WORKSTREAM = Path(__file__).resolve().parent
REPO_ROOT = WORKSTREAM.parents[1]
SPECTRUM = REPO_ROOT / "scripts" / "EBCH128_64.wd"
OUTPUT = WORKSTREAM / "ebch128_randomstepconv_g1_s30_allq_d11_diagnostic.json"

B = 128
K = 64
FAIR_REGIONS = 127
L = 16_560
ACTIVE_ROWS = 16_384
N = B * L
D = (11 * N + 99) // 100
MEMORY = 30
EXACT_MAX_Q = 64
LOG2 = math.log(2.0)


def configure_outer_rows(
    outer_rows: int,
    distance_numerator: int = 11,
    distance_denominator: int = 100,
) -> None:
    """Set the routed row count and exact relative-distance target."""

    if outer_rows < ACTIVE_ROWS:
        raise ValueError("outer row count cannot be below the message row count")
    if not 0 < distance_numerator < distance_denominator:
        raise ValueError("distance target must lie strictly between zero and one")
    global L, N, D
    L = outer_rows
    N = B * L
    D = (
        distance_numerator * N + distance_denominator - 1
    ) // distance_denominator


def load_spectrum(path: Path) -> dict[int, int]:
    spectrum: dict[int, int] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        weight_text, count_text = stripped.split()
        spectrum[int(weight_text)] = int(count_text)
    if sum(spectrum.values()) != 1 << K:
        raise ValueError("EBCH spectrum mass is not 2^64")
    if spectrum.get(0) != 1 or spectrum.get(B) != 1:
        raise ValueError("expected unique zero and all-one words")
    if any(spectrum.get(weight, 0) != spectrum.get(B - weight, 0) for weight in spectrum):
        raise ValueError("EBCH spectrum is not complement symmetric")
    if any(weight & 1 for weight, count in spectrum.items() if count):
        raise ValueError("EBCH spectrum is not even")
    return spectrum


def body_envelope(spectrum: dict[int, int]) -> tuple[Fraction, int]:
    candidates = [
        (Fraction(count * (1 << (B - 1)), math.comb(B, weight)), weight)
        for weight, count in spectrum.items()
        if 0 < weight < B and count
    ]
    return max(candidates)


def step_matrices(z: float, memory_bits: int) -> tuple[np.ndarray, np.ndarray]:
    q = math.ldexp(1.0, -memory_bits)
    b = (1.0 + z) / 2.0
    terminate = q * b
    survive = (1.0 - q) * b
    return (
        np.asarray(((1.0, 0.0), (terminate, survive))),
        np.asarray(((terminate, survive), (terminate, survive))),
    )


def log_entries(matrix: np.ndarray) -> np.ndarray:
    result = np.full_like(matrix, -math.inf, dtype=np.float64)
    positive = matrix > 0.0
    result[positive] = np.log(matrix[positive])
    return result


def log_matmul(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    result = np.empty(np.broadcast_shapes(left.shape, right.shape), dtype=np.float64)
    result[..., 0, 0] = np.logaddexp(
        left[..., 0, 0] + right[..., 0, 0],
        left[..., 0, 1] + right[..., 1, 0],
    )
    result[..., 0, 1] = np.logaddexp(
        left[..., 0, 0] + right[..., 0, 1],
        left[..., 0, 1] + right[..., 1, 1],
    )
    result[..., 1, 0] = np.logaddexp(
        left[..., 1, 0] + right[..., 0, 0],
        left[..., 1, 1] + right[..., 1, 0],
    )
    result[..., 1, 1] = np.logaddexp(
        left[..., 1, 0] + right[..., 0, 1],
        left[..., 1, 1] + right[..., 1, 1],
    )
    return result


def log_identity(count: int | None = None) -> np.ndarray:
    shape = (2, 2) if count is None else (count, 2, 2)
    result = np.full(shape, -math.inf)
    result[..., 0, 0] = 0.0
    result[..., 1, 1] = 0.0
    return result


def log_power(matrix: np.ndarray, exponent: int) -> np.ndarray:
    count = None if matrix.ndim == 2 else len(matrix)
    result = log_identity(count)
    power = matrix.copy()
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = log_matmul(result, power)
        remaining >>= 1
        if remaining:
            power = log_matmul(power, power)
    return result


def exact_region_coefficients(
    zero: np.ndarray, candidate: np.ndarray
) -> np.ndarray:
    log_zero = log_entries(zero)
    log_candidate = log_entries(candidate)
    current = np.full((EXACT_MAX_Q + 1, 2, 2), -math.inf)
    current[0] = log_identity()
    for completed in range(L):
        maximum = min(completed + 1, EXACT_MAX_Q)
        old_maximum = min(completed, EXACT_MAX_Q)
        updated = np.full_like(current, -math.inf)

        zero_products = log_matmul(current[: old_maximum + 1], log_zero)
        degrees = np.arange(old_maximum + 1, dtype=np.float64)
        weights = ((completed + 1.0) - degrees) / (completed + 1.0)
        updated[: old_maximum + 1] = (
            zero_products + np.log(weights)[:, None, None]
        )

        candidate_products = log_matmul(current[:maximum], log_candidate)
        selected = np.arange(1, maximum + 1, dtype=np.float64)
        terms = candidate_products + np.log(
            selected / (completed + 1.0)
        )[:, None, None]
        updated[1 : maximum + 1] = np.logaddexp(
            updated[1 : maximum + 1], terms
        )
        current = updated
    return current


def final_zero_region(zero: np.ndarray) -> np.ndarray:
    return log_power(log_entries(zero), L)


def moments_from_regions(
    regions: np.ndarray, zero_region: np.ndarray
) -> np.ndarray:
    powered = log_power(regions, FAIR_REGIONS)
    complete = log_matmul(powered, zero_region)
    return np.logaddexp(complete[..., 0, 0], complete[..., 0, 1])


def dense_conditioned_regions(
    zero: np.ndarray, candidate: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    occupations = np.arange(EXACT_MAX_Q + 1, ACTIVE_ROWS + 1, dtype=np.float64)
    probabilities = occupations / float(L)
    mixed = (
        (1.0 - probabilities)[:, None, None] * zero[None, :, :]
        + probabilities[:, None, None] * candidate[None, :, :]
    )
    regions = log_power(log_entries(mixed), L)
    log_point_probability = (
        np.asarray(
            [
                math.lgamma(L + 1)
                - math.lgamma(int(q) + 1)
                - math.lgamma(L - int(q) + 1)
                for q in occupations
            ]
        )
        + occupations * np.log(probabilities)
        + (L - occupations) * np.log1p(-probabilities)
    )
    regions -= log_point_probability[:, None, None]
    return occupations.astype(np.int64), regions


def log_choose(n: int, q: np.ndarray) -> np.ndarray:
    q_float = q.astype(np.float64)
    return (
        math.lgamma(n + 1)
        - np.asarray([math.lgamma(int(value) + 1) for value in q])
        - np.asarray([math.lgamma(n - int(value) + 1) for value in q])
    )


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    spectrum = load_spectrum(args.spectrum)
    envelope, maximizing_weight = body_envelope(spectrum)
    effective_mass = envelope + 1
    log_effective_mass = math.log(effective_mass.numerator) - math.log(
        effective_mass.denominator
    )
    occupations = np.arange(1, ACTIVE_ROWS + 1, dtype=np.int64)
    outer_logs = log_choose(ACTIVE_ROWS, occupations) + occupations * log_effective_mass
    best = np.full(ACTIVE_ROWS, math.inf)
    witnesses = np.full(ACTIVE_ROWS, math.nan)

    u_values = np.arange(
        args.grid_min,
        args.grid_max + args.grid_step / 2.0,
        args.grid_step,
    )
    u_values = np.asarray(
        sorted(
            set(float(value) for value in u_values)
            | set(
                float(value)
                for value in np.arange(
                    args.refine_min,
                    args.refine_max + args.refine_step / 2.0,
                    args.refine_step,
                )
            )
        )
    )
    for index, u in enumerate(u_values):
        surprisal = math.exp(float(u))
        z = math.exp(-surprisal)
        zero, active = step_matrices(z, args.memory_bits)
        candidate = 0.5 * (zero + active)
        zero_region = final_zero_region(zero)

        sparse_regions = exact_region_coefficients(zero, candidate)[1:]
        sparse_moments = moments_from_regions(sparse_regions, zero_region)
        dense_q, dense_regions = dense_conditioned_regions(zero, candidate)
        dense_moments = moments_from_regions(dense_regions, zero_region)
        moments = np.concatenate((sparse_moments, dense_moments))
        if not np.array_equal(dense_q, occupations[EXACT_MAX_Q:]):
            raise AssertionError("dense occupation indexing failed")
        inner = np.minimum(0.0, moments + D * surprisal)
        values = outer_logs + inner
        improved = values < best
        best[improved] = values[improved]
        witnesses[improved] = float(u)
        print(f"tilt,{index + 1},{len(u_values)},u,{u:.6f}", flush=True)

    aggregate = float(logsumexp(best))
    dominant = int(np.argmax(best))
    top = np.argsort(best)[-20:][::-1]
    source_bytes = args.spectrum.read_bytes()
    return {
        "schema": "ebch128-randomstepconv-g1-allq-diagnostic-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "claim": {
            "log2_expected_bad_upper": aggregate / LOG2,
            "margin_bits": -aggregate / LOG2,
            "closes_40_bits": aggregate < -40.0 * LOG2,
            "dominant_occupation": dominant + 1,
            "dominant_pointwise_log2_upper": float(best[dominant]) / LOG2,
        },
        "parameters": {
            "outer_code": "extended BCH [128,64,22]",
            "outer_bits": B,
            "outer_dimension": K,
            "outer_rows": L,
            "active_outer_rows": ACTIVE_ROWS,
            "shortened_outer_rows": L - ACTIVE_ROWS,
            "message_bits": K * ACTIVE_ROWS,
            "output_bits": N,
            "distance_cutoff": D,
            "memory_bits": args.memory_bits,
            "fair_reference_regions": FAIR_REGIONS,
            "deleted_input_regions": 1,
            "deleted_output_regions": 0,
            "exact_region_max_occupation": EXACT_MAX_Q,
        },
        "outer_envelope": {
            "body_factor_numerator": str(envelope.numerator),
            "body_factor_denominator": str(envelope.denominator),
            "body_factor_log2": math.log2(envelope.numerator) - math.log2(envelope.denominator),
            "maximizing_body_weight": maximizing_weight,
            "effective_mass_log2": log_effective_mass / LOG2,
            "description": (
                "each active row is either the unique all-one word with unit "
                "mass or an ordinary word dominated by body_factor times a "
                "uniform-even reference"
            ),
        },
        "probability_space": (
            "the EBCH constituent is fixed and repeated; row-coordinate and "
            "region permutations are independent; all RandomStepConv maps "
            "are independent, sampled once, and shared by all messages"
        ),
        "dense_relaxation": (
            "for Q>64, each fixed-Q region coefficient is bounded by the "
            "iid Bernoulli(Q/L) transfer divided by Pr[Bin(L,Q/L)=Q], then "
            "the factor is paid independently in all 127 fair regions"
        ),
        "top_occupations": [
            {
                "active_outer_rows": int(q + 1),
                "pointwise_log2_upper": float(best[q]) / LOG2,
                "log_surprisal": float(witnesses[q]),
            }
            for q in top
        ],
        "occupation_rows": [
            {
                "active_outer_rows": int(q + 1),
                "pointwise_log2_upper": float(best[q]) / LOG2,
                "log_surprisal": float(witnesses[q]),
            }
            for q in range(ACTIVE_ROWS)
        ],
        "spectrum": {
            "path": str(args.spectrum),
            "sha256": hashlib.sha256(source_bytes).hexdigest(),
            "mass": sum(spectrum.values()),
            "minimum_nonzero_weight": min(weight for weight in spectrum if weight),
        },
        "limitations": [
            "Nearest binary64 arithmetic is not an outward certificate.",
            "The Q>64 Bernoulli-conditioning step is rigorous algebraically but numerically diagnostic.",
            "The final parity input is deleted by RandomStepConv monotonicity; no output coordinate is deleted.",
        ],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectrum", type=Path, default=SPECTRUM)
    parser.add_argument("--outer-rows", type=int, default=L)
    parser.add_argument("--memory-bits", type=int, default=MEMORY)
    parser.add_argument("--distance-numerator", type=int, default=11)
    parser.add_argument("--distance-denominator", type=int, default=100)
    parser.add_argument("--grid-min", type=float, default=-12.0)
    parser.add_argument("--grid-max", type=float, default=1.0)
    parser.add_argument("--grid-step", type=float, default=0.5)
    parser.add_argument("--refine-min", type=float, default=0.65)
    parser.add_argument("--refine-max", type=float, default=0.80)
    parser.add_argument("--refine-step", type=float, default=0.025)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    configure_outer_rows(
        args.outer_rows,
        args.distance_numerator,
        args.distance_denominator,
    )
    payload = evaluate(args)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(json.dumps(payload["outer_envelope"], indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
