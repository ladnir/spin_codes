#!/usr/bin/env python3
"""One-active outer-shell diagnostic for Riffle LDPCSplitState."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

from analyze_riffle_ldpcsplitstate_nested_ensemble import (
    distinct_sampling_probability,
    expected_spectrum,
)


DEFAULT_OUTPUT = Path(
    "constructions/riffle_ldpcsplitstate_g4_t256_s64/"
    "receipts/one_active_outer_model.json"
)


def log_choose(n: int, k: int) -> float:
    if not 0 <= k <= n:
        return -math.inf
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def outer_log_spectrum(minimum_distance: int) -> dict[int, float]:
    interior = list(range(minimum_distance, 256 - minimum_distance + 1, 2))
    raw = np.asarray([log_choose(256, weight) for weight in interior])
    normalization = math.log((1 << 128) - 2) - float(logsumexp(raw))
    result = {
        weight: float(value + normalization)
        for weight, value in zip(interior, raw)
    }
    result[256] = 0.0
    return result


def conditioned_spectrum(
    cutoff: int, *, apply_conditioning_scale: bool
) -> tuple[np.ndarray, float]:
    spectrum = expected_spectrum(
        p_column_weight=6, parity_column_weight=3, accumulator_depth=2
    )
    expected_bad = float(np.sum(spectrum[1 : cutoff + 1]))
    distinct_probability = distinct_sampling_probability(math.comb(64, 3), 192)
    joint_probability = distinct_probability - expected_bad
    if joint_probability <= 0.0:
        raise RuntimeError("the requested expurgation has no positive mass bound")
    upper = spectrum.copy()
    upper[1 : cutoff + 1] = 0.0
    if apply_conditioning_scale:
        upper[cutoff + 1 :] /= joint_probability
    return upper, joint_probability


def epoch_matrices(
    *,
    z: float,
    upper_spectrum: np.ndarray,
    constituent_distance: int,
    deterministic_secant: bool,
) -> tuple[np.ndarray, np.ndarray, dict[str, float]]:
    denominator = (1 << 64) - 1
    weights = np.arange(257)
    spectrum_moment = float(upper_spectrum[1:] @ (z ** weights[1:])) / denominator
    if deterministic_secant:
        mean_weight = 256 * (1 << 63) / denominator
        low_fraction = (256 - mean_weight) / (256 - constituent_distance)
        nonzero_moment = (
            low_fraction * z**constituent_distance
            + (1.0 - low_fraction) * z**256
        )
    else:
        mean_weight = math.nan
        low_fraction = math.nan
        nonzero_moment = spectrum_moment
    empty = np.asarray([[1.0, 0.0], [0.0, nonzero_moment]], dtype=np.longdouble)

    sparse_floor = z ** (constituent_distance - 1)
    full_coset_upper = nonzero_moment + (1.0 - z) / denominator
    live_total = min(sparse_floor, full_coset_upper)
    termination = sparse_floor / denominator
    impulse = np.asarray(
        [[0.0, z], [termination, live_total]], dtype=np.longdouble
    )
    return empty, impulse, {
        "empty_live_moment_log2": math.log2(nonzero_moment),
        "impulse_live_total_log2": math.log2(live_total),
        "impulse_termination_log2": math.log2(termination),
        "impulse_sparse_floor_log2": math.log2(sparse_floor),
        "impulse_full_coset_upper_log2": math.log2(full_coset_upper),
        "spectrum_live_moment_log2": math.log2(spectrum_moment),
        "deterministic_secant": deterministic_secant,
        "nonzero_codeword_mean_weight": mean_weight,
        "secant_low_endpoint_fraction": low_fraction,
    }


def region_log_matrices(
    *, z: float, epoch_details: dict[str, float]
) -> tuple[np.ndarray, np.ndarray]:
    log_m = epoch_details["empty_live_moment_log2"] * math.log(2.0)
    log_live = epoch_details["impulse_live_total_log2"] * math.log(2.0)
    log_termination = epoch_details["impulse_termination_log2"] * math.log(2.0)
    m = math.exp(log_m)
    geometric = sum(m**power for power in range(32))
    log_geometric_average = math.log(geometric) - math.log(32.0)

    empty_region = np.full((2, 2), -math.inf)
    empty_region[0, 0] = 0.0
    empty_region[1, 1] = 32 * log_m

    active_region = np.full((2, 2), -math.inf)
    active_region[0, 1] = math.log(z) + log_geometric_average
    active_region[1, 0] = log_termination + log_geometric_average
    active_region[1, 1] = log_live + 31 * log_m
    return empty_region, active_region


def all_region_weight_log_moments(
    empty_region: np.ndarray, active_region: np.ndarray
) -> np.ndarray:
    vectors = np.full((257, 2), -math.inf)
    vectors[0, 0] = 0.0
    for region in range(256):
        following = np.full_like(vectors, -math.inf)
        for active_count in range(region + 1):
            for start_state in range(2):
                start = vectors[active_count, start_state]
                if not math.isfinite(start):
                    continue
                for end_state in range(2):
                    following[active_count, end_state] = np.logaddexp(
                        following[active_count, end_state],
                        start + empty_region[start_state, end_state],
                    )
                    following[active_count + 1, end_state] = np.logaddexp(
                        following[active_count + 1, end_state],
                        start + active_region[start_state, end_state],
                    )
        vectors = following
    result = np.full(257, -math.inf)
    for weight in range(257):
        result[weight] = (
            float(logsumexp(vectors[weight])) - log_choose(256, weight)
        )
    return result


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    upper_spectrum, joint_probability = conditioned_spectrum(
        args.cutoff,
        apply_conditioning_scale=not args.omit_conditioning_scale,
    )
    outer_spectrum = outer_log_spectrum(args.outer_minimum_distance)
    target_distance = math.floor(args.relative_distance * (1 << 21))
    rows = []
    best: dict[str, object] | None = None
    for log_surprisal in np.linspace(-10.0, 2.0, 145):
        surprisal = math.exp(float(log_surprisal))
        z = math.exp(-surprisal)
        empty, impulse, epoch_details = epoch_matrices(
            z=z,
            upper_spectrum=upper_spectrum,
            constituent_distance=args.cutoff + 1,
            deterministic_secant=args.deterministic_secant,
        )
        empty_region, active_region = region_log_matrices(
            z=z, epoch_details=epoch_details
        )
        log_moments = all_region_weight_log_moments(
            empty_region, active_region
        )
        terms = [
            math.log(args.outer_blocks)
            + log_multiplicity
            + log_moments[weight]
            + target_distance * surprisal
            for weight, log_multiplicity in outer_spectrum.items()
        ]
        log_bound = float(logsumexp(terms))
        dominant_weight = max(
            outer_spectrum,
            key=lambda weight: (
                outer_spectrum[weight] + log_moments[weight]
            ),
        )
        row = {
            "log_surprisal": float(log_surprisal),
            "z": z,
            "aggregate_log2_bound": log_bound / math.log(2.0),
            "aggregate_margin_bits": -log_bound / math.log(2.0),
            "dominant_outer_weight": dominant_weight,
            "epoch": epoch_details,
        }
        rows.append(row)
        if best is None or row["aggregate_log2_bound"] < best["aggregate_log2_bound"]:
            best = row
    assert best is not None
    return {
        "schema": "riffle-ldpcsplitstate-one-active-v1",
        "candidate": "Riffle LDPCSplitState g=4 t=256 s=64",
        "parameters": {
            "outer_blocks": args.outer_blocks,
            "outer_model": (
                "real-valued complement-symmetric even [256,128,38]-shaped"
            ),
            "outer_minimum_distance": args.outer_minimum_distance,
            "constituent_expurgation_cutoff": args.cutoff,
            "constituent_joint_event_probability_lower_bound": joint_probability,
            "conditioning_scale_applied": not args.omit_conditioning_scale,
            "deterministic_secant_moment": args.deterministic_secant,
            "target_distance": target_distance,
            "relative_distance": args.relative_distance,
        },
        "best": best,
        "grid": rows,
        "scope": (
            "Floating-point one-active diagnostic. The zero/live recurrence "
            "and outer coordinate placement are exact for one active block. "
            "The live moments use the conditional expected-spectrum existence "
            "bound at each tested tilt; one common fixed pair across the full "
            "tilt grid is not certified. If the conditioning scale is omitted, "
            "the result is an ideal expurgated-shape diagnostic rather than an "
            "existence bound. The outer spectrum is modeled."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--outer-blocks", type=int, default=8192)
    parser.add_argument("--outer-minimum-distance", type=int, default=38)
    parser.add_argument("--cutoff", type=int, default=39)
    parser.add_argument("--relative-distance", type=float, default=0.09)
    parser.add_argument("--omit-conditioning-scale", action="store_true")
    parser.add_argument("--deterministic-secant", action="store_true")
    args = parser.parse_args()
    payload = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    best = payload["best"]
    print(
        f"margin_bits,{best['aggregate_margin_bits']:.6f},"
        f"dominant_outer_weight,{best['dominant_outer_weight']},"
        f"log_surprisal,{best['log_surprisal']:.6f},z,{best['z']:.9e}"
    )
    print(
        "empty_live_moment_log2,"
        f"{best['epoch']['empty_live_moment_log2']:.6f},"
        "impulse_live_total_log2,"
        f"{best['epoch']['impulse_live_total_log2']:.6f},"
        "termination_log2,"
        f"{best['epoch']['impulse_termination_log2']:.6f}"
    )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
