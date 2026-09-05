#!/usr/bin/env python3
"""Pointwise global audit for FrozenRandom512 with field parities.

The audit retains the exact random-linear packet spectrum and the exact
parity distribution. It uses a Chernoff coefficient upper for the global
outer enumerator and the termination-separated pointwise inner upper.
"""

from __future__ import annotations

import argparse
import json
import math
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import gammaln, logsumexp

from analyze_riffle_rm512_randomstepconv_g4_sigma20 import LOG2
from audit_riffle_frozenrandom512_p2_onegroup import (
    DEFAULT_SPECTRUM,
    ebch_packet_support_distribution,
)
from analyze_riffle_randomstepconv_g4_sigma20_goal01 import (
    TARGET_G,
    optimize_tilts,
)


TOTAL_INPUT_BITS = 1 << 20
DATA_PACKETS = 1 << 19
RELATIVE_DISTANCE = 0.09
DEFAULT_OUTPUT = Path(
    "constructions/riffle_frozenrandom512_p2_randomstepconv_g4_sigma20/"
    "receipts/global_pointwise_p2.json"
)


def termination_upper_profile(
    support: int, sigma: int, packet_positions: int, distance: int
) -> dict[str, object]:
    optimized = optimize_tilts(
        packet_positions, TARGET_G, sigma, support, distance
    )
    return {
        "termination_separated_log2_upper": optimized["diagnostic_log2_upper"],
        "dominant_termination_count": optimized["dominant_termination_count"],
    }


def log_two_power_minus_one(exponent: int) -> float:
    return exponent * LOG2 + math.log1p(-math.ldexp(1.0, -exponent))


class GlobalOuterMoment:
    """Expected global support enumerator under independently frozen maps."""

    def __init__(
        self,
        ebch_spectrum: Path,
        parity_count: int = 2,
        outer_bits: int = 512,
    ) -> None:
        if not 1 <= parity_count <= 8:
            raise ValueError("current parity audit requires 1 <= P <= 8")
        if outer_bits < 128 or outer_bits > 2048 or outer_bits % 128:
            raise ValueError("outer bits must be a multiple of 128 in [128, 2048]")
        self.local_output_bits = outer_bits
        self.local_input_bits = outer_bits // 2
        if TOTAL_INPUT_BITS % self.local_input_bits:
            raise ValueError("local input size must divide the total input size")
        self.data_groups = TOTAL_INPUT_BITS // self.local_input_bits
        self.local_output_packets = outer_bits // TARGET_G
        self.group_input_count = (1 << self.local_input_bits) - 1
        self.field_input_symbols = self.local_input_bits // 64
        self.parity_count = parity_count
        self.single_parity = ebch_packet_support_distribution(ebch_spectrum)
        self.uniform_pair = np.asarray([1.0])
        for _index in range(parity_count):
            self.uniform_pair = np.convolve(
                self.uniform_pair, self.single_parity
            )
        if parity_count <= self.field_input_symbols:
            self.one_group_pair = self.uniform_pair.copy()
            if self.local_input_bits <= 512:
                total_messages = float(1 << self.local_input_bits)
                self.one_group_pair[0] = (
                    total_messages * self.one_group_pair[0] - 1.0
                ) / (total_messages - 1.0)
                self.one_group_pair[1:] *= (
                    total_messages / (total_messages - 1.0)
                )
            # Above 512 input bits, conditioning away the all-zero input
            # changes double-precision probabilities by less than 2^-512.
            self.one_group_pair /= self.one_group_pair.sum()
        else:
            self.one_group_pair = None
        self.uniform_pair_supports = np.arange(len(self.uniform_pair))
        self.one_group_pair_supports = (
            np.arange(len(self.one_group_pair))
            if self.one_group_pair is not None
            else np.asarray([], dtype=np.int64)
        )
        self.log_uniform_pair = np.full(len(self.uniform_pair), -math.inf)
        self.log_one_group_pair = np.full(
            len(self.one_group_pair) if self.one_group_pair is not None else 0,
            -math.inf,
        )
        uniform_positive = self.uniform_pair > 0.0
        one_positive = (
            self.one_group_pair > 0.0
            if self.one_group_pair is not None
            else np.asarray([], dtype=bool)
        )
        self.log_uniform_pair[uniform_positive] = np.log(
            self.uniform_pair[uniform_positive]
        )
        if self.one_group_pair is not None:
            self.log_one_group_pair[one_positive] = np.log(
                self.one_group_pair[one_positive]
            )
        self.log_group_input_count = log_two_power_minus_one(
            self.local_input_bits
        )
        self.one_group_mds_weights: list[tuple[int, int]] = []
        if parity_count > self.field_input_symbols:
            distance = parity_count - self.field_input_symbols + 1
            q = 1 << 64
            for weight in range(distance, parity_count + 1):
                alternating = sum(
                    (-1) ** index
                    * math.comb(weight - 1, index)
                    * q ** (weight - distance - index)
                    for index in range(weight - distance + 1)
                )
                count = (
                    math.comb(parity_count, weight)
                    * (q - 1)
                    * alternating
                )
                if count <= 0:
                    raise AssertionError("invalid MDS weight count")
                self.one_group_mds_weights.append((weight, count))
            if (
                sum(count for _weight, count in self.one_group_mds_weights)
                != self.group_input_count
            ):
                raise AssertionError("MDS weight enumerator mass changed")
        self.log_output_nonzero_ratio = (
            self.log_group_input_count
            - log_two_power_minus_one(self.local_output_bits)
        )
        occupations = np.arange(self.data_groups + 1, dtype=np.float64)
        self.occupations = occupations
        self.log_choose_groups = (
            gammaln(self.data_groups + 1)
            - gammaln(occupations + 1)
            - gammaln(self.data_groups - occupations + 1)
        )

    def log_active_data_moment(self, log_t: float) -> float:
        t = math.exp(log_t)
        exponent = self.local_output_packets * math.log1p(15.0 * t)
        log_expm1 = (
            exponent + math.log1p(-math.exp(-exponent))
            if exponent > 50.0
            else math.log(math.expm1(exponent))
        )
        return self.log_output_nonzero_ratio + log_expm1

    @staticmethod
    def log_probability_moment(
        log_probabilities: np.ndarray,
        supports: np.ndarray,
        log_t: float,
    ) -> float:
        return float(logsumexp(log_probabilities + supports * log_t))

    def log_moment(self, log_t: float, minimum_occupation: int = 1) -> float:
        log_a = self.log_active_data_moment(log_t)
        if self.one_group_pair is not None:
            log_one_parity = self.log_probability_moment(
                self.log_one_group_pair, self.one_group_pair_supports, log_t
            )
        elif log_t <= 0.0:
            log_one_parity = float(
                logsumexp(
                    [
                        math.log(count)
                        - self.log_group_input_count
                        + 6.0 * weight * log_t
                        for weight, count in self.one_group_mds_weights
                    ]
                )
            )
        else:
            log_one_parity = 32.0 * self.parity_count * log_t
        log_uniform_parity = self.log_probability_moment(
            self.log_uniform_pair, self.uniform_pair_supports, log_t
        )

        one_group = (
            self.log_choose_groups[1] + log_a + log_one_parity
            if minimum_occupation <= 1
            else -math.inf
        )
        start = max(2, minimum_occupation)
        occupation_logs = (
            self.log_choose_groups[start:]
            + self.occupations[start:] * log_a
        )
        two_or_more = float(logsumexp(occupation_logs)) + log_uniform_parity

        # For r>=2, the exact parity moment differs from the uniform moment by
        # at most 2^-512 relative to the nonuniform component. This adjustment
        # is below double precision and is recorded in the receipt.
        return float(np.logaddexp(one_group, two_or_more))

    def tilted_profile(self, log_t: float, support: int) -> dict[str, object]:
        log_a = self.log_active_data_moment(log_t)
        if self.one_group_pair is not None:
            log_one_parity = self.log_probability_moment(
                self.log_one_group_pair, self.one_group_pair_supports, log_t
            )
        elif log_t <= 0.0:
            log_one_parity = float(
                logsumexp(
                    [
                        math.log(count)
                        - self.log_group_input_count
                        + 6.0 * weight * log_t
                        for weight, count in self.one_group_mds_weights
                    ]
                )
            )
        else:
            log_one_parity = 32.0 * self.parity_count * log_t
        log_uniform_parity = self.log_probability_moment(
            self.log_uniform_pair, self.uniform_pair_supports, log_t
        )
        occupation_logs = np.full(self.data_groups + 1, -math.inf)
        minimum_occupation = max(
            1,
            math.ceil(
                (support - 32 * self.parity_count)
                / self.local_output_packets
            ),
        )
        if minimum_occupation <= 1:
            occupation_logs[1] = self.log_choose_groups[1] + log_a + log_one_parity
        start = max(2, minimum_occupation)
        occupation_logs[start:] = (
            self.log_choose_groups[start:]
            + self.occupations[start:] * log_a
            + log_uniform_parity
        )
        occupation_probabilities = np.exp(
            occupation_logs - float(logsumexp(occupation_logs))
        )
        occupation_cumulative = np.cumsum(occupation_probabilities)

        t = math.exp(log_t)
        exponent = self.local_output_packets * math.log1p(15.0 * t)
        mean_data_support = (
            self.local_output_packets
            * 15.0
            * t
            / (1.0 + 15.0 * t)
            / (1.0 - math.exp(-exponent))
        )
        uniform_parity_probabilities = np.exp(
            self.log_uniform_pair
            + self.uniform_pair_supports * log_t
            - log_uniform_parity
        )
        if self.one_group_pair is not None:
            one_group_parity_mean = float(
                np.exp(
                    self.log_one_group_pair
                    + self.one_group_pair_supports * log_t
                    - log_one_parity
                )
                @ self.one_group_pair_supports
            )
        else:
            one_group_parity_mean = float("nan")
        return {
            "occupation_mode": int(np.argmax(occupation_logs)),
            "mean_active_groups": float(
                occupation_probabilities @ self.occupations
            ),
            "active_group_quantiles": {
                "q10": int(np.searchsorted(occupation_cumulative, 0.10)),
                "q50": int(np.searchsorted(occupation_cumulative, 0.50)),
                "q90": int(np.searchsorted(occupation_cumulative, 0.90)),
            },
            "mean_data_support_per_active_group": mean_data_support,
            "mean_parity_support_for_one_group": one_group_parity_mean,
            "mean_parity_support_for_multiple_groups": float(
                uniform_parity_probabilities @ self.uniform_pair_supports
            ),
        }

    @lru_cache(maxsize=None)
    def coefficient_chernoff_log2(self, support: int) -> tuple[float, float]:
        minimum_occupation = max(
            1,
            math.ceil(
                (support - 32 * self.parity_count)
                / self.local_output_packets
            ),
        )

        def objective(log_t: float) -> float:
            return (
                self.log_moment(log_t, minimum_occupation)
                - support * log_t
            )

        grid = np.linspace(-12.0, 2.5, 59)
        values = np.asarray([objective(float(point)) for point in grid])
        index = int(np.argmin(values))
        low = float(grid[max(0, index - 1)])
        high = float(grid[min(len(grid) - 1, index + 1)])
        result = minimize_scalar(
            objective,
            bounds=(low, high),
            method="bounded",
            options={"xatol": 2e-9, "maxiter": 200},
        )
        if not result.success:
            raise ArithmeticError(str(result.message))
        return float(result.fun) / LOG2, float(result.x)


def support_schedule(packet_positions: int) -> list[int]:
    supports = set(range(1, 513))
    supports.update(range(520, 4097, 8))
    supports.update(range(4352, 8193, 256))
    supports.update((9216, 10240, 12288, 14336, 16384, 20480, 24576, 28672, 32767))
    return sorted(support for support in supports if support <= packet_positions)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ebch-spectrum", type=Path, default=DEFAULT_SPECTRUM)
    parser.add_argument("--parity-count", type=int, default=2)
    parser.add_argument(
        "--outer-bits",
        type=int,
        default=512,
        help="rate-half random outer constituent output size",
    )
    parser.add_argument("--sigma", type=int, default=20)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    packet_positions = DATA_PACKETS + 32 * args.parity_count
    distance = math.floor(RELATIVE_DISTANCE * TARGET_G * packet_positions)
    outer = GlobalOuterMoment(
        args.ebch_spectrum, args.parity_count, args.outer_bits
    )
    rows: list[dict[str, object]] = []
    schedule = support_schedule(packet_positions)
    for index, support in enumerate(schedule):
        outer_log2, log_t = outer.coefficient_chernoff_log2(support)
        inner = termination_upper_profile(
            support, args.sigma, packet_positions, distance
        )
        inner_log2 = float(inner["termination_separated_log2_upper"])
        rows.append(
            {
                "packet_support": support,
                "outer_coefficient_chernoff_log2": outer_log2,
                "outer_log_tilt": log_t,
                "inner_probability_upper_log2": inner_log2,
                "pointwise_expected_count_upper_log2": outer_log2 + inner_log2,
                "dominant_termination_count": int(
                    inner["dominant_termination_count"]
                ),
            }
        )
        if index % 100 == 0:
            print(
                f"progress,{index + 1}/{len(schedule)},support,{support},"
                f"combined,{outer_log2 + inner_log2:.6f}",
                flush=True,
            )

    dominant = sorted(
        rows,
        key=lambda row: float(row["pointwise_expected_count_upper_log2"]),
        reverse=True,
    )
    payload = {
        "schema": "riffle-randomouter-global-pointwise-v2",
        "model": (
            f"Riffle FrozenRandom{args.outer_bits}-P{args.parity_count}-"
            f"RandomStepConv g=4 sigma={args.sigma}"
        ),
        "outer_constituent": {
            "length_bits": outer.local_output_bits,
            "dimension_bits": outer.local_input_bits,
            "data_groups": outer.data_groups,
            "output_packets_per_group": outer.local_output_packets,
        },
        "state_bits": args.sigma,
        "field_parity_symbols": args.parity_count,
        "packet_positions": packet_positions,
        "distance": distance,
        "support_schedule_count": len(schedule),
        "support_schedule_scope": (
            "Every support through 512, every eighth support through 4096, "
            "then checkpoints through 32767. This is a diagnostic scan, not "
            "a complete support certificate."
        ),
        "dominant_sampled_profiles": dominant[:30],
        "dominant_sampled_outer_profile": outer.tilted_profile(
            float(dominant[0]["outer_log_tilt"]),
            int(dominant[0]["packet_support"]),
        ),
        "rows": rows,
        "evidence": {
            "data_packet_spectrum": "exact random-linear ensemble expectation",
            "one_group_parity_distribution": (
                "exact surjective distribution when P does not exceed the "
                "number of input field symbols; otherwise an MDS field-weight "
                "enumerator with the EBCH packet-support lower bound"
            ),
            "multi_group_parity_distribution": "uniform pair moment with relative correction below 2^-512",
            "outer_coefficient": "rigorous Chernoff upper apart from recorded 2^-512 parity correction",
            "inner": "termination-separated pointwise floating upper",
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    top = dominant[0]
    print(
        "dominant_sampled,"
        f"support={top['packet_support']},"
        f"outer={top['outer_coefficient_chernoff_log2']:.6f},"
        f"inner={top['inner_probability_upper_log2']:.6f},"
        f"combined={top['pointwise_expected_count_upper_log2']:.6f}"
    )
    print(f"receipt,{args.output}")


if __name__ == "__main__":
    main()
