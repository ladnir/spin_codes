#!/usr/bin/env python3
"""Evaluate exact RM(4,9) [512,256,32] blocks with RandomStepConv.

The 4096 data groups use the exact RM(4,9) weight spectrum. The final
128 parity-input bits retain one 256-bit rate-half block, which preserves the
existing total of 524352 four-bit packets. That parity block is treated
worst-case in the many-group moment.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import expit, gammaln, logsumexp
from scipy.stats import binom

from analyze_riffle_outer256model_randomstepconv_g4_sigma20 import (
    LOG2,
    canonical_late_suffix_event,
)
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
from analyze_riffle_randomstepconv_g4_sigma20_goal04 import inner_parameters


LOCAL_BITS = 512
LOCAL_INPUT_BITS = 256
PACKET_BITS = 4
PACKETS_PER_LOCAL = LOCAL_BITS // PACKET_BITS
DATA_GROUPS = 1 << 12
PARITY_BLOCK_PACKETS = 64
DEFAULT_SPECTRUM = Path(__file__).with_name("rm512_256_spectrum.csv")
DEFAULT_OUTPUT = Path(
    "constructions/riffle_rm512_randomstepconv_g4_sigma20/receipts/"
    "rm512_randomstepconv_distance09.json"
)


def log_choose(n: int, k: int) -> float:
    if not 0 <= k <= n:
        return -math.inf
    return float(gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1))


def weighted_quantiles(
    values: np.ndarray, probabilities: np.ndarray
) -> dict[str, int]:
    levels = (0.01, 0.10, 0.25, 0.50, 0.75, 0.90, 0.99)
    cumulative = np.cumsum(probabilities)
    return {
        f"q{round(100 * level):02d}": int(
            values[min(int(np.searchsorted(cumulative, level)), len(values) - 1)]
        )
        for level in levels
    }


def read_spectrum(path: Path) -> tuple[np.ndarray, list[int], dict[str, object]]:
    rows: list[tuple[int, int]] = []
    with path.open(newline="", encoding="utf-8") as source:
        for row in csv.DictReader(source):
            weight = int(row["weight"])
            count = int(row["count"])
            if not 0 <= weight <= LOCAL_BITS or count <= 0:
                raise ValueError("invalid RM spectrum row")
            rows.append((weight, count))
    rows.sort()
    if len({weight for weight, _count in rows}) != len(rows):
        raise ValueError("duplicate RM spectrum weight")
    counts = dict(rows)
    if sum(counts.values()) != 1 << LOCAL_INPUT_BITS:
        raise ValueError("RM spectrum mass is not 2^256")
    if counts.get(0) != 1 or counts.get(LOCAL_BITS) != 1:
        raise ValueError("RM spectrum endpoints changed")
    if any(counts.get(weight, 0) != counts.get(LOCAL_BITS - weight, 0) for weight in counts):
        raise ValueError("RM spectrum lost complement symmetry")
    minimum_distance = min(weight for weight, count in rows if weight > 0 and count)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    return (
        np.asarray([weight for weight, _count in rows], dtype=np.int64),
        [count for _weight, count in rows],
        {
            "path": str(path),
            "sha256": digest,
            "total_mass": str(1 << LOCAL_INPUT_BITS),
            "minimum_distance": minimum_distance,
            "complement_symmetric": True,
        },
    )


def packet_support_coefficients() -> list[list[int]]:
    """Return coefficients of ((1+x)^4-1)^s for all s."""
    polynomials: list[list[int]] = [[1]]
    for _support in range(1, PACKETS_PER_LOCAL + 1):
        previous = polynomials[-1]
        current = [0] * (len(previous) + PACKET_BITS)
        for degree, coefficient in enumerate(previous):
            current[degree + 1] += 4 * coefficient
            current[degree + 2] += 6 * coefficient
            current[degree + 3] += 4 * coefficient
            current[degree + 4] += coefficient
        polynomials.append(current)
    return polynomials


class ExactRMMoments:
    """Exact local spectrum and packet-support moments for RM(4,9)."""

    def __init__(self, spectrum_path: Path) -> None:
        self.weights, counts, self.spectrum_receipt = read_spectrum(spectrum_path)
        self.log_multiplicities = np.asarray(
            [math.log(count) for count in counts], dtype=np.float64
        )
        self.nonzero = self.weights != 0
        self.supports = np.arange(PACKETS_PER_LOCAL + 1, dtype=np.float64)
        self.log_support_probabilities = np.full(
            (len(self.weights), PACKETS_PER_LOCAL + 1),
            -math.inf,
            dtype=np.float64,
        )
        support_polynomials = packet_support_coefficients()
        for row, weight_value in enumerate(self.weights):
            weight = int(weight_value)
            if weight == 0:
                self.log_support_probabilities[row, 0] = 0.0
                continue
            denominator = math.comb(LOCAL_BITS, weight)
            observed = 0
            for support in range(1, PACKETS_PER_LOCAL + 1):
                polynomial = support_polynomials[support]
                if weight >= len(polynomial):
                    continue
                coefficient = polynomial[weight]
                if coefficient == 0:
                    continue
                count = math.comb(PACKETS_PER_LOCAL, support) * coefficient
                observed += count
                self.log_support_probabilities[row, support] = (
                    math.log(count) - math.log(denominator)
                )
            if observed != denominator:
                raise AssertionError(f"packet-support law failed at weight {weight}")

    def log_q(self, log_t: float) -> np.ndarray:
        return logsumexp(
            self.log_support_probabilities + self.supports * log_t,
            axis=1,
        )

    def log_moment(self, log_t: float) -> float:
        log_q = self.log_q(log_t)
        return float(
            logsumexp(self.log_multiplicities[self.nonzero] + log_q[self.nonzero])
        )

    def active_profile(self, log_t: float) -> dict[str, object]:
        joint_logs = (
            self.log_multiplicities[:, None]
            + self.log_support_probabilities
            + self.supports[None, :] * log_t
        )[self.nonzero]
        joint = np.exp(joint_logs - float(logsumexp(joint_logs)))
        weights = self.weights[self.nonzero]
        weight_probabilities = joint.sum(axis=1)
        support_probabilities = joint.sum(axis=0)
        joint_mode = np.unravel_index(int(np.argmax(joint_logs)), joint_logs.shape)
        return {
            "mean_binary_weight": float(weight_probabilities @ weights),
            "binary_weight_quantiles": weighted_quantiles(
                weights, weight_probabilities
            ),
            "mean_packet_support": float(
                support_probabilities @ self.supports
            ),
            "packet_support_quantiles": weighted_quantiles(
                self.supports.astype(np.int64), support_probabilities
            ),
            "joint_modal_binary_weight": int(weights[joint_mode[0]]),
            "joint_modal_packet_support": int(joint_mode[1]),
            "combined_packet_support_tilt": math.exp(log_t),
        }

    def tilted_block_distribution(
        self, log_t: float
    ) -> tuple[float, np.ndarray]:
        active_by_support = logsumexp(
            self.log_multiplicities[self.nonzero, None]
            + self.log_support_probabilities[self.nonzero],
            axis=0,
        )
        logs = active_by_support + self.supports * log_t
        logs[0] = np.logaddexp(logs[0], 0.0)
        normalizer = float(logsumexp(logs))
        return normalizer, np.exp(logs - normalizer)


def nonempty_power_log(log_x: float, positions: int) -> float:
    total = positions * float(np.logaddexp(0.0, log_x))
    if total > 40.0:
        return total + math.log1p(-math.exp(-total))
    return math.log(math.expm1(total))


def all_occupation_log_moment(
    model: ExactRMMoments, log_t: float
) -> tuple[float, dict[str, float]]:
    log_m1 = model.log_moment(log_t)
    data_groups = nonempty_power_log(log_m1, DATA_GROUPS)
    parity_charge = max(0.0, PARITY_BLOCK_PACKETS * log_t)
    total = data_groups + parity_charge
    active_probability = float(expit(log_m1))
    mode = min(
        float(DATA_GROUPS),
        max(1.0, math.floor((DATA_GROUPS + 1) * active_probability)),
    )
    return total, {
        "log2_data_group_component": data_groups / LOG2,
        "parity_charge_log2": parity_charge / LOG2,
        "occupation_mode": mode,
        "active_probability": active_probability,
        "log2_M1": log_m1 / LOG2,
    }


def refined_segments_with_low_support(
    rows: list[dict[str, object]],
) -> list[tuple[int, int, dict[str, object]]]:
    boundaries = {1, 18, TARGET_PACKETS // 2 + 1, TARGET_PACKETS + 1}
    for row in rows:
        boundaries.add(int(row["support_left"]))
        boundaries.add(int(row["support_right"]) + 1)
    for boundary in tuple(boundaries):
        reflected = TARGET_PACKETS - boundary + 1
        if 1 <= reflected <= TARGET_PACKETS:
            boundaries.add(reflected)
    ordered = sorted(boundary for boundary in boundaries if 1 <= boundary <= TARGET_PACKETS + 1)
    segments: list[tuple[int, int, dict[str, object]]] = []
    row_index = 0
    for left, stop in zip(ordered, ordered[1:]):
        right = stop - 1
        if right < 18:
            source = rows[0]
        else:
            while left > int(rows[row_index]["support_right"]):
                row_index += 1
            source = rows[row_index]
            if right > int(source["support_right"]):
                raise AssertionError("segment crosses an inner certificate boundary")
        segments.append((left, right, source))
    if segments[0][0] != 1 or segments[-1][1] != TARGET_PACKETS:
        raise AssertionError("segments do not cover support 1 through N")
    return segments


def optimize_weighted_interval(
    model: ExactRMMoments,
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
        return (
            float(result.fun),
            math.exp(float(result.x)),
            bool(result.success),
            str(result.message),
        )

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
        "occupation_mode": selected_meta["occupation_mode"],
        "active_probability": selected_meta["active_probability"],
        "parity_charge_log2": selected_meta["parity_charge_log2"],
        "unconditioned_occupation_mode": unconditioned_meta["occupation_mode"],
    }


def convolved_support_distribution(block: np.ndarray) -> np.ndarray:
    maximum_degree = DATA_GROUPS * PACKETS_PER_LOCAL
    transform_size = 1 << maximum_degree.bit_length()
    padded = np.pad(block, (0, transform_size - len(block)))
    distribution = np.fft.irfft(
        np.fft.rfft(padded) ** DATA_GROUPS, transform_size
    )[: maximum_degree + 1]
    minimum = float(distribution.min())
    if minimum < -1e-12:
        raise ArithmeticError("FFT support distribution has material negative mass")
    distribution = np.maximum(distribution, 0.0)
    distribution /= distribution.sum()
    return distribution


def evaluate(
    model: ExactRMMoments, rows: list[dict[str, object]]
) -> dict[str, object]:
    segments: list[dict[str, object]] = []
    for left, right, source_row in refined_segments_with_low_support(rows):
        parameters = inner_parameters(source_row)
        negative_left = -log_binomial(TARGET_PACKETS, left)
        negative_right = -log_binomial(TARGET_PACKETS, right)
        if left == right:
            secant_slope = 0.0
            secant_intercept = negative_left
        else:
            secant_slope = (negative_right - negative_left) / (right - left)
            secant_intercept = negative_left - secant_slope * left
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
        segments.append(result)
        print(
            f"segment,{left}:{right},log2,{contribution / LOG2:.6f},"
            f"tail,{weighted['selected_tail']},"
            f"mode,{weighted['occupation_mode']:.0f}",
            flush=True,
        )

    total = float(
        logsumexp(
            [float(row["log2_contribution_upper"]) * LOG2 for row in segments]
        )
    ) / LOG2
    dominant = max(segments, key=lambda row: float(row["log2_contribution_upper"]))
    combined_log_t = (
        float(dominant["log_support_base"])
        + float(dominant["reciprocal_binomial_secant_slope_log2"]) * LOG2
        + float(dominant["support_log_tilt"])
    )
    conditioned_profile = model.active_profile(combined_log_t)
    active_probability = float(dominant["active_probability"])
    mean_active_groups = DATA_GROUPS * active_probability
    representative_support = round(
        mean_active_groups * float(conditioned_profile["mean_packet_support"])
    )
    representative_support = min(
        int(dominant["support_right"]),
        max(int(dominant["support_left"]), representative_support),
    )

    log_block_normalizer, block_distribution = model.tilted_block_distribution(
        combined_log_t
    )
    total_distribution = convolved_support_distribution(block_distribution)
    exact_support_probability = float(total_distribution[representative_support])
    relaxed_outer_coefficient_log2 = (
        DATA_GROUPS * log_block_normalizer
        + math.log(exact_support_probability)
        - representative_support * combined_log_t
    ) / LOG2
    inner_upper_log2 = (
        float(dominant["log_constant"])
        + representative_support * float(dominant["log_support_base"])
        - log_binomial(TARGET_PACKETS, representative_support)
    ) / LOG2
    pointwise_upper_log2 = relaxed_outer_coefficient_log2 + inner_upper_log2

    secant_intercept = (
        float(dominant["reciprocal_binomial_secant_intercept_log2"]) * LOG2
    )
    secant_slope = (
        float(dominant["reciprocal_binomial_secant_slope_log2"]) * LOG2
    )
    exact_reciprocal = -log_binomial(TARGET_PACKETS, representative_support)
    secant_gap_bits = (
        secant_intercept
        + secant_slope * representative_support
        - exact_reciprocal
    ) / LOG2

    left = int(dominant["support_left"])
    right = int(dominant["support_right"])
    support_indices = np.arange(left, right + 1)
    support_tilt = float(dominant["support_tilt"])
    selected_tail = str(dominant["selected_tail"])
    if selected_tail == "lower_tail":
        exponent = right - support_indices
    elif selected_tail == "upper_tail":
        exponent = left - support_indices
    else:
        exponent = np.zeros_like(support_indices)
    chernoff_ratio = float(
        np.sum(
            total_distribution[left : right + 1]
            * np.exp(exponent * math.log(support_tilt))
        )
    )

    full_support_upper = min(
        TARGET_PACKETS, representative_support + PARITY_BLOCK_PACKETS
    )
    explicit_event = canonical_late_suffix_event(
        full_support_upper, int(rows[0]["distance"])
    )
    explicit_event["scope"] = (
        "This per-word probability lower bound uses the maximum possible full "
        "support for every exact RM data configuration counted at the "
        "representative data support."
    )
    explicit_expected_count_log2 = (
        relaxed_outer_coefficient_log2
        + float(explicit_event["log2_probability_lower"])
    )

    return {
        "relative_binary_weight": float(rows[0]["relative_binary_weight"]),
        "distance": int(rows[0]["distance"]),
        "total_log2_first_moment_upper": total,
        "dominant_support_interval": [left, right],
        "dominant_occupation_mode": dominant["occupation_mode"],
        "unconditioned_active_block": model.active_profile(0.0),
        "dominant_conditioned_active_block": conditioned_profile,
        "dominant_occupation": {
            "mean_active_groups": mean_active_groups,
            "standard_deviation_active_groups": math.sqrt(
                DATA_GROUPS * active_probability * (1.0 - active_probability)
            ),
            "active_group_quantiles": {
                "q01": int(binom.ppf(0.01, DATA_GROUPS, active_probability)),
                "q10": int(binom.ppf(0.10, DATA_GROUPS, active_probability)),
                "q50": int(binom.ppf(0.50, DATA_GROUPS, active_probability)),
                "q90": int(binom.ppf(0.90, DATA_GROUPS, active_probability)),
                "q99": int(binom.ppf(0.99, DATA_GROUPS, active_probability)),
            },
        },
        "representative_support_audit": {
            "packet_support": representative_support,
            "relaxed_data_outer_coefficient_log2": relaxed_outer_coefficient_log2,
            "inner_probability_upper_log2": inner_upper_log2,
            "pointwise_first_moment_upper_log2": pointwise_upper_log2,
            "dominant_segment_upper_log2": float(
                dominant["log2_contribution_upper"]
            ),
            "segment_minus_pointwise_upper_bits": (
                float(dominant["log2_contribution_upper"])
                - pointwise_upper_log2
            ),
            "reciprocal_binomial_secant_gap_bits": secant_gap_bits,
            "support_interval_chernoff_loss_bits": -math.log2(chernoff_ratio),
            "exact_support_probability_under_saddle": exact_support_probability,
        },
        "canonical_late_placement_test": {
            "data_packet_support": representative_support,
            "full_packet_support_charged": full_support_upper,
            "parity_packets_charged": full_support_upper - representative_support,
            "inner_event": explicit_event,
            "modeled_expected_count_log2": explicit_expected_count_log2,
            "is_refutation_witness": explicit_expected_count_log2 >= 0.0,
            "validity": (
                "The event charges the maximum 64 packets from the determined "
                "parity block. It is therefore uniform over the parity value "
                "for every data configuration counted at the representative "
                "data support."
            ),
        },
        "segments": segments,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--spectrum", type=Path, default=DEFAULT_SPECTRUM)
    parser.add_argument("--inner-receipt", type=Path, default=DEFAULT_INNER_RECEIPT)
    parser.add_argument("--relative-distance", type=float, default=0.09)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    grouped = load_inner_caps(args.inner_receipt)
    if args.relative_distance not in grouped:
        raise ValueError("inner receipt does not contain the requested distance")
    model = ExactRMMoments(args.spectrum)
    result = evaluate(model, grouped[args.relative_distance])
    payload = {
        "schema": "riffle-rm512-randomstepconv-v1",
        "model": "Riffle RM512-RandomStepConv g=4 sigma=20",
        "parameters": {
            "data_local_code": [LOCAL_BITS, LOCAL_INPUT_BITS, 32],
            "data_groups": DATA_GROUPS,
            "data_packets_per_group": PACKETS_PER_LOCAL,
            "parity_input_bits": 128,
            "parity_output_bits": 256,
            "parity_block_packets": PARITY_BLOCK_PACKETS,
            "packet_positions": TARGET_PACKETS,
            "packet_bits": PACKET_BITS,
            "relative_binary_weight": args.relative_distance,
        },
        "spectrum": model.spectrum_receipt,
        "evidence": {
            "data_local_weight_spectrum": "EXACT_INTEGER_RM_4_9",
            "packet_support_law_given_weight": "EXACT_INTEGER_VALIDATED",
            "data_configuration_sum": "EXACT_RELATIVE_TO_SPECTRUM",
            "parity_block": "RIGOROUS_WORST_CASE_MOMENT",
            "inner_caps": "OUTWARD_ROUNDED_FROM_ONE_LAP_GOAL02",
            "optimization_fft_and_final_sum": "FLOATING_DIAGNOSTIC",
        },
        "result": result,
        "scope": (
            "The RM data-block spectrum and packetization law are exact. The "
            "parity block is pessimistic. Inner caps are certified, while the "
            "outer optimizations, FFT audit, and final sum are floating "
            "diagnostics rather than a proof certificate."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"total_log2,{result['total_log2_first_moment_upper']:.6f}")
    print(f"dominant_support,{result['dominant_support_interval']}")
    print(f"dominant_mode,{result['dominant_occupation_mode']:.0f}")
    print(
        "pointwise_upper_log2,"
        f"{result['representative_support_audit']['pointwise_first_moment_upper_log2']:.6f}"
    )
    print(
        "explicit_event_expected_log2,"
        f"{result['canonical_late_placement_test']['modeled_expected_count_log2']:.6f}"
    )
    print(f"receipt,{args.output}")


if __name__ == "__main__":
    main()
