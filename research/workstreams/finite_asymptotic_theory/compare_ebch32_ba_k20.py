#!/usr/bin/env python3
"""Compare Golay--BA-3/B=240 and EBCH32--BA-3/B=256 at k=2^20.

This program is a nearest-binary64 diagnostic.  It computes each exact
direct-sum spectrum, propagates it through two uniform-interleaver
accumulators, sweeps symmetric tail-free conditioning windows, and evaluates
the finite RM2Sub occupation-one transfer.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.special import gammaln, logsumexp


WORKSTREAM = Path(__file__).resolve().parent
sys.path.insert(0, str(WORKSTREAM))

from evaluate_golay_ba_rm2sub_finite import one_active  # noqa: E402
from evaluate_random_outer_rm2sub_one_active import (  # noqa: E402
    DEFAULT_ACTIVATION,
    DEFAULT_LIVE_SPECTRUM,
)


MESSAGE_BITS = 1 << 20
STEP_BITS = 128
STATE_BITS = 19
CONSTITUENT_DISTANCE = 48
OUTPUT = WORKSTREAM / "ebch32_ba3_vs_golay_ba3_k20_diagnostic.json"


def logadd(values) -> float:
    finite = [float(value) for value in values if math.isfinite(float(value))]
    if not finite:
        return -math.inf
    return float(logsumexp(np.asarray(finite, dtype=np.float64)))


def log_choose(total: int, selected: np.ndarray | int) -> np.ndarray:
    values = np.asarray(selected)
    answer = np.full(values.shape, -math.inf, dtype=np.float64)
    valid = (0 <= values) & (values <= total)
    answer[valid] = (
        gammaln(total + 1)
        - gammaln(values[valid] + 1)
        - gammaln(total - values[valid] + 1)
    )
    return answer


def direct_sum_log_counts(
    length: int, local_length: int, local_dimension: int, local_spectrum: dict[int, int]
) -> np.ndarray:
    if length % local_length:
        raise ValueError("outer length must be divisible by constituent length")
    distribution = np.asarray([1.0], dtype=np.float64)
    polynomial = np.zeros(local_length + 1, dtype=np.float64)
    for weight, count in local_spectrum.items():
        polynomial[weight] = count / (1 << local_dimension)
    for _ in range(length // local_length):
        distribution = np.convolve(distribution, polynomial)
    result = np.full(length + 1, -math.inf, dtype=np.float64)
    positive = distribution > 0.0
    result[positive] = np.log(distribution[positive]) + (length // 2) * math.log(2.0)
    return result


def accumulator_transition(length: int) -> np.ndarray:
    transition = np.full((length + 1, length + 1), -math.inf, dtype=np.float64)
    transition[0, 0] = 0.0
    outputs = np.arange(length + 1)
    for input_weight in range(1, length + 1):
        runs = (input_weight + 1) // 2
        left = log_choose(length, input_weight)
        first = np.full(length + 1, -math.inf, dtype=np.float64)
        valid_first = outputs >= runs
        first[valid_first] = (
            gammaln(outputs[valid_first])
            - gammaln(runs)
            - gammaln(outputs[valid_first] - runs + 1)
        )
        second_total = length - outputs
        second_selected = input_weight - runs
        second = np.full(length + 1, -math.inf, dtype=np.float64)
        valid_second = second_total >= second_selected
        second[valid_second] = (
            gammaln(second_total[valid_second] + 1)
            - gammaln(second_selected + 1)
            - gammaln(second_total[valid_second] - second_selected + 1)
        )
        transition[input_weight] = first + second - float(left)
    return transition


def expected_ba_log_spectrum(
    length: int,
    local_length: int,
    local_dimension: int,
    local_spectrum: dict[int, int],
    parity_fanout: tuple[int, int] | None = None,
) -> np.ndarray:
    initial = direct_sum_log_counts(
        length, local_length, local_dimension, local_spectrum
    )
    if parity_fanout is not None:
        initial = parity_fanout_log_spectrum(
            initial, length, parity_fanout[0], parity_fanout[1]
        )
    transition = accumulator_transition(length)
    after_one = logsumexp(initial[:, None] + transition, axis=0)
    return logsumexp(after_one[:, None] + transition, axis=0)


def parity_fanout_log_spectrum(
    spectrum: np.ndarray, length: int, source_size: int, target_size: int
) -> np.ndarray:
    """Average a disjoint random parity fanout exactly in combinatorial form."""
    if source_size <= 0 or target_size <= 0 or source_size + target_size > length:
        raise ValueError("invalid parity-fanout sizes")
    buckets: list[list[float]] = [[] for _ in range(length + 1)]
    source_denominator = float(log_choose(length, source_size))
    target_denominator = float(log_choose(length - source_size, target_size))
    for weight in range(length + 1):
        log_count = float(spectrum[weight])
        if not math.isfinite(log_count):
            continue
        minimum_source = max(0, source_size - (length - weight))
        maximum_source = min(source_size, weight)
        for source_overlap in range(minimum_source, maximum_source + 1):
            log_source_probability = (
                float(log_choose(weight, source_overlap))
                + float(log_choose(length - weight, source_size - source_overlap))
                - source_denominator
            )
            if source_overlap % 2 == 0:
                buckets[weight].append(log_count + log_source_probability)
                continue
            remaining_weight = weight - source_overlap
            remaining_length = length - source_size
            minimum_target = max(0, target_size - (remaining_length - remaining_weight))
            maximum_target = min(target_size, remaining_weight)
            for target_overlap in range(minimum_target, maximum_target + 1):
                output_weight = weight + target_size - 2 * target_overlap
                log_target_probability = (
                    float(log_choose(remaining_weight, target_overlap))
                    + float(
                        log_choose(
                            remaining_length - remaining_weight,
                            target_size - target_overlap,
                        )
                    )
                    - target_denominator
                )
                buckets[output_weight].append(
                    log_count + log_source_probability + log_target_probability
                )
    result = np.full(length + 1, -math.inf, dtype=np.float64)
    for weight, terms in enumerate(buckets):
        if terms:
            result[weight] = logadd(terms)
    before = logadd(spectrum)
    after = logadd(result)
    if abs(before - after) > 1e-9:
        raise ArithmeticError("parity fanout did not preserve spectrum mass")
    return result


def sweep_window(spectrum: np.ndarray, active: dict[str, object], length: int) -> list[dict[str, object]]:
    pointwise = {
        int(row["outer_weight"]): float(row["pointwise_log2_upper"]) * math.log(2.0)
        for row in active["weight_rows"]
    }
    rows = []
    for lower in range(1, length // 2 + 1):
        upper = length - lower
        tail = logadd(
            list(spectrum[1:lower]) + list(spectrum[upper + 1 :])
        )
        tail_mass = 0.0 if tail == -math.inf else math.exp(tail)
        if tail_mass >= 1.0:
            continue
        good = 1.0 - tail_mass
        inside = logadd(
            pointwise.get(weight, -math.inf)
            for weight in range(lower, upper + 1)
        )
        conditional = inside - math.log(good)
        rows.append(
            {
                "window": [lower, upper],
                "tail_expected_words": tail_mass,
                "good_event_probability_lower": good,
                "expected_trials_upper": 1.0 / good,
                "conditioned_q1_margin_bits": -conditional / math.log(2.0),
            }
        )
    return rows


def random_spectrum_excess(spectrum: np.ndarray, length: int) -> dict[str, object]:
    dimension = length // 2
    normalization = math.log((2**dimension - 1) / (2**length - 1))
    rows = []
    for weight in range(1, length + 1):
        if not math.isfinite(float(spectrum[weight])):
            continue
        random_log = float(log_choose(length, weight)) + normalization
        rows.append((weight, (float(spectrum[weight]) - random_log) / math.log(2.0)))
    worst = max(rows, key=lambda row: row[1])
    return {
        "maximum_excess_over_uniform_random_bits": worst[1],
        "maximum_excess_weight": worst[0],
        "central_excess_at_half_weight_bits": next(
            excess for weight, excess in rows if weight == length // 2
        ),
    }


def evaluate(configuration: dict[str, object]) -> dict[str, object]:
    outer_bits = int(configuration["outer_bits"])
    outer_rows = int(configuration["outer_rows"])
    spectrum = expected_ba_log_spectrum(
        outer_bits,
        int(configuration["local_length"]),
        int(configuration["local_dimension"]),
        configuration["local_spectrum"],
        configuration.get("parity_fanout"),
    )
    results = []
    for relative_distance in (0.11, 0.109):
        distance = math.floor(relative_distance * outer_bits * outer_rows)
        active = one_active(
            spectrum=spectrum,
            outer_bits=outer_bits,
            outer_blocks=outer_rows,
            distance=distance,
            step_bits=STEP_BITS,
            state_bits=STATE_BITS,
            constituent_distance=CONSTITUENT_DISTANCE,
            activation=DEFAULT_ACTIVATION,
            live_spectrum_path=DEFAULT_LIVE_SPECTRUM,
            tilt_minimum=-10.0,
            tilt_maximum=-6.0,
            tilt_spacing=0.01,
        )
        windows = sweep_window(spectrum, active, outer_bits)
        best = max(windows, key=lambda row: float(row["conditioned_q1_margin_bits"]))
        results.append(
            {
                "relative_distance": relative_distance,
                "distance": distance,
                "unconditioned_q1_margin_bits": active["aggregate_margin_bits"],
                "best_symmetric_conditioning_window": best,
                "selected_windows": [
                    row
                    for row in windows
                    if int(row["window"][0]) in range(20, 33)
                ],
            }
        )
    return {
        "name": configuration["name"],
        "outer_bits": outer_bits,
        "outer_rows": outer_rows,
        "output_bits": outer_bits * outer_rows,
        "parent_dimension": outer_bits * outer_rows // 2,
        "shortened_coordinates": outer_bits * outer_rows // 2 - MESSAGE_BITS,
        "inner_epochs": outer_bits * outer_rows // STEP_BITS,
        "constituent": {
            "length": configuration["local_length"],
            "dimension": configuration["local_dimension"],
            "spectrum": configuration["local_spectrum"],
        },
        "parity_fanout": configuration.get("parity_fanout"),
        "spectrum_comparison": random_spectrum_excess(spectrum, outer_bits),
        "distance_results": results,
    }


def main() -> None:
    configurations = [
        {
            "name": "Golay--BA-3/B=240",
            "outer_bits": 240,
            "outer_rows": 8832,
            "local_length": 24,
            "local_dimension": 12,
            "local_spectrum": {0: 1, 8: 759, 12: 2576, 16: 759, 24: 1},
        },
        {
            "name": "EBCH32--BA-3/B=256",
            "outer_bits": 256,
            "outer_rows": 8192,
            "local_length": 32,
            "local_dimension": 16,
            "local_spectrum": {
                0: 1,
                8: 620,
                12: 13888,
                16: 36518,
                20: 13888,
                24: 620,
                32: 1,
            },
        },
        {
            "name": "EBCH32--ParityFanout31x33--BA-3/B=256",
            "outer_bits": 256,
            "outer_rows": 8192,
            "local_length": 32,
            "local_dimension": 16,
            "local_spectrum": {
                0: 1,
                8: 620,
                12: 13888,
                16: 36518,
                20: 13888,
                24: 620,
                32: 1,
            },
            "parity_fanout": (31, 33),
        },
    ]
    payload = {
        "schema": "ebch32-ba3-vs-golay-ba3-k20-diagnostic-v1",
        "status": "NEAREST_BINARY64_DIAGNOSTIC",
        "message_bits": MESSAGE_BITS,
        "configurations": [evaluate(configuration) for configuration in configurations],
        "limitations": [
            "The BA expectations and RM2Sub transfer are evaluated in nearest binary64 arithmetic.",
            "The comparison covers occupation one and simple spectrum diagnostics, not the dense occupation cover.",
            "A conditioning window requires an efficient membership test before it defines a practical setup algorithm.",
        ],
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                row["name"]: {
                    "geometry": [row["outer_bits"], row["outer_rows"], row["output_bits"]],
                    "spectrum": row["spectrum_comparison"],
                    "distance_results": row["distance_results"],
                }
                for row in payload["configurations"]
            },
            indent=2,
        )
    )
    print(f"output={OUTPUT}")


if __name__ == "__main__":
    main()
