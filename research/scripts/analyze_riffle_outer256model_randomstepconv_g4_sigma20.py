#!/usr/bin/env python3
"""Evaluate a fixed-rate [256,128,38]-shaped outer model.

This script is a model diagnostic, not a construction certificate. It pairs
the 16384 GF(2^64) data symbols, applies one modeled rate-one-half binary
constituent per pair, and retains the one-lap RandomStepConv inner. The model
uses a complement-symmetric random-like even weight spectrum with minimum
nonzero weight 38.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import gammaln, logsumexp
from scipy.stats import binom

from analyze_riffle_randomstepconv_g4_sigma20_goal02 import (
    TARGET_G,
    TARGET_PACKETS,
    TARGET_SIGMA,
    log_binomial,
)
from analyze_riffle_randomstepconv_g4_sigma20_goal03 import (
    DEFAULT_INNER_RECEIPT,
    load_inner_caps,
)
from analyze_riffle_randomstepconv_g4_sigma20_goal04 import (
    inner_parameters,
    log_binomial_tail_from_two,
    refined_segments,
)


MODEL_BLOCK_BITS = 256
MODEL_INPUT_BITS = 128
DEFAULT_MODEL_DISTANCE = 38
PACKET_BITS = 4
PACKETS_PER_BLOCK = MODEL_BLOCK_BITS // PACKET_BITS
DATA_PAIRS = 1 << 13
LOG2 = math.log(2.0)
DEFAULT_INFLATIONS = (0.0, 10.0, 20.0, 40.0)
DEFAULT_OUTPUT = Path(
    "constructions/riffle_outer256model_randomstepconv_g4_sigma20/receipts/"
    "fixed_rate_outer256_model.json"
)


def log_choose(n: int, k: int) -> float:
    if not 0 <= k <= n:
        return -math.inf
    return float(gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1))


def packet_support_count(weight: int, support: int) -> int:
    coefficient = 0
    for selected in range(support + 1):
        available = PACKET_BITS * selected
        if available < weight:
            continue
        coefficient += (
            (-1) ** (support - selected)
            * math.comb(support, selected)
            * math.comb(available, weight)
        )
    return math.comb(PACKETS_PER_BLOCK, support) * coefficient


class ModeledLocalMoments:
    """Packet-support moments for the real-valued spectrum model."""

    def __init__(self, minimum_distance: int, low_weight_inflation_bits: float) -> None:
        if minimum_distance < 2 or minimum_distance >= MODEL_BLOCK_BITS // 2:
            raise ValueError("modeled distance must lie in [2,127]")
        if minimum_distance % 2:
            raise ValueError("the even-weight model requires an even distance")
        interior_weights = np.asarray(
            list(range(minimum_distance, MODEL_BLOCK_BITS - minimum_distance + 1, 2)),
            dtype=np.int64,
        )
        raw_logs = np.asarray(
            [log_choose(MODEL_BLOCK_BITS, int(weight)) for weight in interior_weights],
            dtype=np.float64,
        )
        sensitive = (interior_weights <= minimum_distance + 32) | (
            interior_weights >= MODEL_BLOCK_BITS - minimum_distance - 32
        )
        raw_logs = raw_logs + sensitive * (low_weight_inflation_bits * LOG2)
        target_interior_mass = float((1 << MODEL_INPUT_BITS) - 2)
        normalized_logs = raw_logs + math.log(target_interior_mass) - float(logsumexp(raw_logs))

        self.weights = np.concatenate(
            (np.asarray([0]), interior_weights, np.asarray([MODEL_BLOCK_BITS]))
        )
        self.log_multiplicities = np.concatenate(
            (np.asarray([0.0]), normalized_logs, np.asarray([0.0]))
        )
        self.nonzero = self.weights != 0
        self.supports = np.arange(PACKETS_PER_BLOCK + 1, dtype=np.float64)
        self.log_support_probabilities = np.full(
            (len(self.weights), PACKETS_PER_BLOCK + 1),
            -math.inf,
            dtype=np.float64,
        )
        for row, weight_value in enumerate(self.weights):
            weight = int(weight_value)
            if weight == 0:
                self.log_support_probabilities[row, 0] = 0.0
                continue
            denominator = math.comb(MODEL_BLOCK_BITS, weight)
            total = 0
            for support in range(1, PACKETS_PER_BLOCK + 1):
                count = packet_support_count(weight, support)
                if count == 0:
                    continue
                total += count
                self.log_support_probabilities[row, support] = (
                    math.log(count) - math.log(denominator)
                )
            if total != denominator:
                raise AssertionError(f"packet-support law failed at weight {weight}")

        observed_log_mass = float(logsumexp(self.log_multiplicities))
        if not math.isclose(observed_log_mass, MODEL_INPUT_BITS * LOG2, abs_tol=2e-13):
            raise AssertionError("modeled spectrum mass is not 2^128")
        self.low_weight_inflation_bits = low_weight_inflation_bits
        self.minimum_distance = minimum_distance

    def log_q(self, log_t: float) -> np.ndarray:
        return logsumexp(
            self.log_support_probabilities + self.supports * log_t,
            axis=1,
        )

    def log_moments(self, log_t: float) -> tuple[float, float]:
        log_q = self.log_q(log_t)
        log_m1 = float(
            logsumexp(self.log_multiplicities[self.nonzero] + log_q[self.nonzero])
        )
        log_m2 = float(
            logsumexp(
                self.log_multiplicities[self.nonzero] + 2.0 * log_q[self.nonzero]
            )
        )
        return log_m1, log_m2

    def active_block_profile(self, log_t: float) -> dict[str, float | int]:
        """Return the local saddle profile conditioned on a nonzero block."""
        joint = (
            self.log_multiplicities[:, None]
            + self.log_support_probabilities
            + self.supports[None, :] * log_t
        )
        joint = joint[self.nonzero]
        normalizer = float(logsumexp(joint))
        probabilities = np.exp(joint - normalizer)
        active_weights = self.weights[self.nonzero]
        weight_probabilities = probabilities.sum(axis=1)
        support_probabilities = probabilities.sum(axis=0)
        joint_mode = np.unravel_index(int(np.argmax(joint)), joint.shape)
        return {
            "mean_binary_weight": float(weight_probabilities @ active_weights),
            "modal_binary_weight": int(active_weights[int(np.argmax(weight_probabilities))]),
            "mean_packet_support": float(support_probabilities @ self.supports),
            "modal_packet_support": int(np.argmax(support_probabilities)),
            "joint_modal_binary_weight": int(active_weights[joint_mode[0]]),
            "joint_modal_packet_support": int(joint_mode[1]),
            "combined_packet_support_tilt": math.exp(log_t),
        }


def canonical_late_suffix_event(packet_support: int, distance: int) -> dict[str, object]:
    """Optimize an explicit low-weight event with no post-activation reset.

    All active packets land in a suffix of length R. After the first active
    packet, the sigma-bit state remains nonzero through the suffix. This event
    needs no inner termination.
    """
    if not 1 <= packet_support <= TARGET_PACKETS:
        raise ValueError("packet support must lie in [1,N]")

    log_no_reset_step = math.log1p(-(2.0 ** -TARGET_SIGMA))

    def terms(suffix_packets: int) -> tuple[float, float, float, float]:
        placement = (
            log_choose(suffix_packets, packet_support)
            - log_choose(TARGET_PACKETS, packet_support)
        )
        output_tail = float(binom.logcdf(distance, TARGET_G * suffix_packets, 0.5))
        no_reset = suffix_packets * log_no_reset_step
        return placement + output_tail + no_reset, placement, output_tail, no_reset

    coarse_step = 128
    candidates = range(packet_support, TARGET_PACKETS + 1, coarse_step)
    best_suffix = max(candidates, key=lambda suffix: terms(suffix)[0])
    fine_left = max(packet_support, best_suffix - 2 * coarse_step)
    fine_right = min(TARGET_PACKETS, best_suffix + 2 * coarse_step)
    best_suffix = max(range(fine_left, fine_right + 1), key=lambda suffix: terms(suffix)[0])
    total, placement, output_tail, no_reset = terms(best_suffix)
    return {
        "event": (
            "all active packets lie in the final R positions; after the first "
            "active packet, the inner state stays nonzero; the emitted binary "
            "weight is at most D"
        ),
        "probability_space": "global packet permutation and independent inner maps",
        "packet_support_h": packet_support,
        "suffix_packets_R": best_suffix,
        "forced_zero_prefix_packets": TARGET_PACKETS - best_suffix,
        "distance_threshold_D": distance,
        "post_activation_termination_count": 0,
        "log2_probability_lower": total / LOG2,
        "log2_placement_probability": placement / LOG2,
        "log2_output_tail_probability_lower": output_tail / LOG2,
        "log2_no_termination_probability_lower": no_reset / LOG2,
        "mean_relative_binary_weight_if_all_R_steps_live": (
            best_suffix / (2.0 * TARGET_PACKETS)
        ),
        "scope": (
            "This is a per-word probability lower bound for any fixed outer "
            "word with packet support h. The modeled outer saddle identifies "
            "a representative h but does not certify that an explicit code "
            "contains the representative local profile."
        ),
    }


def all_occupation_log_moment(
    model: ModeledLocalMoments, log_t: float
) -> tuple[float, dict[str, float]]:
    log_m1, log_m2 = model.log_moments(log_t)

    # With one active data pair, the two field parity equations give an
    # invertible map from that pair to the parity pair. Cauchy bounds the
    # product of their local moments by M2.
    log_one_pair = math.log(DATA_PAIRS) + log_m2

    # For two or more active data pairs, discard the parity block when t<=1.
    # For t>1, its packet-support moment is at most t^64.
    log_many_pairs = log_binomial_tail_from_two(log_m1, DATA_PAIRS)
    log_many_pairs += max(0.0, PACKETS_PER_BLOCK * log_t)
    total = float(np.logaddexp(log_one_pair, log_many_pairs))

    if log_m1 > 40.0:
        mode = float(DATA_PAIRS)
    else:
        value = math.exp(log_m1)
        mode = max(2.0, math.floor((DATA_PAIRS + 1) * value / (1.0 + value)))
    return total, {
        "log2_one_pair_component": log_one_pair / LOG2,
        "log2_many_pair_component": log_many_pairs / LOG2,
        "many_pair_occupation_mode": mode,
        "log2_M1": log_m1 / LOG2,
        "log2_M2": log_m2 / LOG2,
    }


def optimize_weighted_interval(
    model: ModeledLocalMoments,
    left: int,
    right: int,
    log_support_base: float,
) -> dict[str, object]:
    unconditioned, unconditioned_meta = all_occupation_log_moment(
        model, log_support_base
    )

    def optimize_side(side: str) -> tuple[float, float, bool, str]:
        def objective(log_surprisal: float) -> float:
            magnitude = math.exp(log_surprisal)
            if side == "lower":
                log_s = -magnitude
                cutoff = right
            else:
                log_s = magnitude
                cutoff = left
            moment, _meta = all_occupation_log_moment(
                model, log_support_base + log_s
            )
            return moment - cutoff * log_s

        grid = np.linspace(-16.0, 3.0, 77)
        values = np.asarray([objective(float(point)) for point in grid])
        index = int(np.argmin(values))
        low = float(grid[max(0, index - 1)])
        high = float(grid[min(len(grid) - 1, index + 1)])
        result = minimize_scalar(
            objective,
            bounds=(low, high),
            method="bounded",
            options={"xatol": 1e-11, "maxiter": 500},
        )
        return float(result.fun), math.exp(float(result.x)), bool(result.success), str(result.message)

    lower, lower_magnitude, lower_success, lower_message = optimize_side("lower")
    upper, upper_magnitude, upper_success, upper_message = optimize_side("upper")
    candidates = (
        ("unconditioned", unconditioned, 0.0, True, "no support tilt"),
        ("lower_tail", lower, -lower_magnitude, lower_success, lower_message),
        ("upper_tail", upper, upper_magnitude, upper_success, upper_message),
    )
    selected = min(candidates, key=lambda item: item[1])
    _selected_value, selected_meta = all_occupation_log_moment(
        model, log_support_base + selected[2]
    )
    return {
        "selected_tail": selected[0],
        "weighted_log_moment_upper": selected[1],
        "support_log_tilt": selected[2],
        "support_tilt": math.exp(selected[2]),
        "optimizer_success": selected[3],
        "optimizer_message": selected[4],
        "one_pair_component_log2": selected_meta["log2_one_pair_component"],
        "many_pair_component_log2": selected_meta["log2_many_pair_component"],
        "many_pair_occupation_mode": selected_meta["many_pair_occupation_mode"],
        "unconditioned_many_pair_occupation_mode": unconditioned_meta[
            "many_pair_occupation_mode"
        ],
    }


def evaluate_distance(
    model: ModeledLocalMoments, rows: list[dict[str, object]]
) -> dict[str, object]:
    results: list[dict[str, object]] = []
    for left, right, source_row in refined_segments(rows):
        parameters = inner_parameters(source_row)
        negative_log_binomial_left = -log_binomial(TARGET_PACKETS, left)
        negative_log_binomial_right = -log_binomial(TARGET_PACKETS, right)
        if left == right:
            secant_slope = 0.0
            secant_intercept = negative_log_binomial_left
        else:
            secant_slope = (
                negative_log_binomial_right - negative_log_binomial_left
            ) / (right - left)
            secant_intercept = negative_log_binomial_left - secant_slope * left
        weighted = optimize_weighted_interval(
            model,
            left,
            right,
            parameters["log_support_base"] + secant_slope,
        )
        contribution = (
            parameters["log_constant"]
            + secant_intercept
            + float(weighted["weighted_log_moment_upper"])
        )
        result: dict[str, object] = {
            "support_left": left,
            "support_right": right,
            "source_certificate_left": int(source_row["support_left"]),
            "source_certificate_right": int(source_row["support_right"]),
            "reciprocal_binomial_secant_intercept_log2": secant_intercept / LOG2,
            "reciprocal_binomial_secant_slope_log2": secant_slope / LOG2,
            "log2_contribution_upper": contribution / LOG2,
        }
        result.update(parameters)
        result.update(weighted)
        results.append(result)
        print(
            f"inflation,{model.low_weight_inflation_bits:.0f},"
            f"distance,{model.minimum_distance},"
            f"segment,{left}:{right},log2,{contribution / LOG2:.6f},"
            f"tail,{weighted['selected_tail']},"
            f"mode,{weighted['many_pair_occupation_mode']:.0f}",
            flush=True,
        )
    total = float(
        logsumexp([float(row["log2_contribution_upper"]) * LOG2 for row in results])
    ) / LOG2
    dominant = max(results, key=lambda row: float(row["log2_contribution_upper"]))
    combined_log_t = (
        float(dominant["log_support_base"])
        + float(dominant["reciprocal_binomial_secant_slope_log2"]) * LOG2
        + float(dominant["support_log_tilt"])
    )
    local_profile = model.active_block_profile(combined_log_t)
    active_blocks = int(dominant["many_pair_occupation_mode"])
    representative_support = round(
        active_blocks * float(local_profile["mean_packet_support"])
    )
    representative_support = min(
        int(dominant["support_right"]),
        max(int(dominant["support_left"]), representative_support),
    )
    canonical_bad_event = {
        "diagnostic_role": (
            "a concrete low-weight mechanism for comparison with the dominant "
            "upper bound; it is not an explicit outer-codeword witness or a "
            "construction-level failure witness"
        ),
        "outer_saddle_profile": {
            "active_data_pair_blocks": active_blocks,
            **local_profile,
            "representative_total_packet_support": representative_support,
        },
        "inner_event": canonical_late_suffix_event(
            representative_support, int(rows[0]["distance"])
        ),
        "interpretation": (
            "The global permutation delays activation until a late suffix. "
            "The inner remains live after activation, so accidental inner "
            "terminations are not part of this event. The per-word event must "
            "still be combined with an outer coefficient before drawing a "
            "construction-level conclusion."
        ),
    }
    return {
        "relative_binary_weight": float(rows[0]["relative_binary_weight"]),
        "distance": int(rows[0]["distance"]),
        "low_weight_inflation_bits": model.low_weight_inflation_bits,
        "modeled_minimum_distance": model.minimum_distance,
        "total_log2_first_moment_upper": total,
        "dominant_support_interval": [
            int(dominant["support_left"]),
            int(dominant["support_right"]),
        ],
        "dominant_occupation_mode": dominant["many_pair_occupation_mode"],
        "canonical_bad_event": canonical_bad_event,
        "segments": results,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inner-receipt", type=Path, default=DEFAULT_INNER_RECEIPT)
    parser.add_argument("--relative-distance", type=float, default=0.09)
    parser.add_argument(
        "--modeled-distances", nargs="+", type=int, default=[DEFAULT_MODEL_DISTANCE]
    )
    parser.add_argument(
        "--inflation-bits", nargs="+", type=float, default=list(DEFAULT_INFLATIONS)
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--no-write", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    grouped = load_inner_caps(args.inner_receipt)
    if args.relative_distance not in grouped:
        raise ValueError("inner receipt does not contain the requested distance")
    results = [
        evaluate_distance(
            ModeledLocalMoments(distance, inflation),
            grouped[args.relative_distance],
        )
        for distance in args.modeled_distances
        for inflation in args.inflation_bits
    ]
    for result in results:
        print(
            f"result,distance,{result['modeled_minimum_distance']},"
            f"inflation,{result['low_weight_inflation_bits']:.0f},"
            f"total_log2,{result['total_log2_first_moment_upper']:.6f},"
            f"dominant,{result['dominant_support_interval']},"
            f"mode,{result['dominant_occupation_mode']:.0f}",
            flush=True,
        )
    payload = {
        "schema": "riffle-outer256model-randomstepconv-v2",
        "model": "Riffle Outer256Model-RandomStepConv g=4 sigma=20",
        "evidence": {
            "local_weight_spectrum": "MODELED_REAL_VALUED_NOT_A_CODE_SPECTRUM",
            "packet_support_law_given_weight": "EXACT_INTEGER_VALIDATED",
            "one_pair_parity_cauchy": "EXACT",
            "many_pair_parity_treatment": "RIGOROUS_WORST_CASE_RELATIVE_TO_MODEL",
            "inner_caps": "OUTWARD_ROUNDED_FROM_ONE_LAP_GOAL02",
            "optimization_and_final_sum": "FLOATING_DIAGNOSTIC",
        },
        "parameters": {
            "modeled_local_code_lengths": [MODEL_BLOCK_BITS, MODEL_INPUT_BITS],
            "modeled_minimum_distances": list(args.modeled_distances),
            "packet_bits": PACKET_BITS,
            "packets_per_local_block": PACKETS_PER_BLOCK,
            "data_pairs": DATA_PAIRS,
            "parity_pairs": 1,
            "packet_positions": TARGET_PACKETS,
            "relative_binary_weight": args.relative_distance,
            "low_weight_inflation_bits": list(args.inflation_bits),
        },
        "spectrum_model": {
            "support": "zero, all-one, and even weights from d through 256-d",
            "interior_shape": "binomial coefficients normalized to total mass 2^128-2",
            "sensitivity_window": "weights 38..70 and 186..218",
            "complement_symmetric": True,
        },
        "results": results,
        "scope": (
            "This receipt evaluates a guessed local spectrum. It does not assert "
            "that a [256,128,38] code with this spectrum exists. The reported "
            "first-moment values are floating diagnostics, not certificates."
        ),
    }
    if not args.no_write:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(f"receipt,{args.output}")


if __name__ == "__main__":
    main()
