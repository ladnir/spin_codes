#!/usr/bin/env python3
"""Diagnose the random-outer RM2Sub dense-occupation exponent.

The derivation is recorded in RM2SUB_DENSE_OCCUPATION.md.  The 7-by-7
association matrix is computed with exact integer Walsh transforms.  The
continuous saddle and Perron roots use binary64 arithmetic and are diagnostic.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize


DEFAULT_SELECTION = Path(
    "constructions/"
    "riffle_bchperm_transpose_bitshuffle_splitstate_preaddmul_rm2sub_t128_s20/"
    "receipts/min_state/s19_rm2sub_selection.json"
)


def fwht_in_place(values: np.ndarray) -> None:
    """Apply the unnormalized Walsh-Hadamard transform in place."""
    width = 1
    while width < values.size:
        view = values.reshape(-1, 2 * width)
        left = view[:, :width].copy()
        right = view[:, width:].copy()
        view[:, :width] = left + right
        view[:, width:] = left - right
        width *= 2


def state_code_weights(
    generator_words: list[int], output_bits: int
) -> np.ndarray:
    state_count = 1 << len(generator_words)
    weights = np.empty(state_count, dtype=np.uint8)
    weights[0] = 0
    codeword = 0
    previous_gray = 0
    for counter in range(1, state_count):
        gray = counter ^ (counter >> 1)
        changed = gray ^ previous_gray
        codeword ^= generator_words[changed.bit_length() - 1]
        weights[gray] = codeword.bit_count()
        previous_gray = gray
    if np.any(weights > output_bits):
        raise AssertionError("state-code weight exceeds the output length")
    return weights


def verify_transpose_columns(
    *, generator_words: list[int], columns: list[int], output_bits: int
) -> None:
    if len(columns) != output_bits:
        raise ValueError("the transpose map has the wrong number of columns")
    for position in range(output_bits):
        derived = sum(
            ((word >> position) & 1) << state_bit
            for state_bit, word in enumerate(generator_words)
        )
        if derived != columns[position]:
            raise ValueError(
                f"B column {position} differs from the transpose of A"
            )


def association_matrix(
    weights: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    classes = np.unique(weights)
    masks = [weights == weight for weight in classes]
    counts = np.asarray([int(np.count_nonzero(mask)) for mask in masks])
    matrix = np.empty((classes.size, classes.size), dtype=np.int64)
    for column, mask in enumerate(masks):
        transformed = mask.astype(np.int64)
        fwht_in_place(transformed)
        for row, source_mask in enumerate(masks):
            matrix[row, column] = int(np.sum(transformed[source_mask]))

    state_count = weights.size
    expected_axis_sum = np.zeros(classes.size, dtype=np.int64)
    expected_axis_sum[0] = state_count
    if not np.array_equal(np.sum(matrix, axis=0), expected_axis_sum):
        raise AssertionError("association column sums failed")
    if not np.array_equal(np.sum(matrix, axis=1), expected_axis_sum):
        raise AssertionError("association row sums failed")
    if not np.array_equal(matrix[0], counts):
        raise AssertionError("association zero row differs from the spectrum")
    if not np.array_equal(matrix[:, 0], counts):
        raise AssertionError("association zero column differs from the spectrum")
    return classes, counts, matrix


def kl_bits(alpha: float, probability: float) -> float:
    if alpha == 0.0:
        return -math.log2(1.0 - probability)
    if alpha == 1.0:
        return -math.log2(probability)
    return (
        alpha * math.log2(alpha / probability)
        + (1.0 - alpha)
        * math.log2((1.0 - alpha) / (1.0 - probability))
    )


class DenseEnvelope:
    def __init__(
        self,
        *,
        classes: np.ndarray,
        counts: np.ndarray,
        association: np.ndarray,
        state_bits: int,
        step_bits: int,
        delta: float,
    ) -> None:
        self.classes = classes.astype(np.float64)
        self.counts = counts.astype(np.float64)
        self.association = association.astype(np.float64)
        self.state_bits = state_bits
        self.step_bits = step_bits
        self.delta = delta
        self.state_space = 1 << state_bits
        self.live_states = self.state_space - 1
        self.punctured_factor = self.live_states / (self.live_states - 1)

    def transfer(
        self, *, candidate_probability: float, surprisal: float
    ) -> tuple[np.ndarray, dict[str, float]]:
        z = math.exp(-surprisal)
        active_bit_probability = candidate_probability / 2.0
        zero_coordinate_moment = (
            1.0 - active_bit_probability + active_bit_probability * z
        )
        one_coordinate_moment = (
            active_bit_probability
            + (1.0 - active_bit_probability) * z
        )
        signed_coordinate_moment = (
            1.0 - active_bit_probability - active_bit_probability * z
        )

        live_output_by_weight = (
            zero_coordinate_moment ** (self.step_bits - self.classes)
            * one_coordinate_moment**self.classes
        )
        syndrome_fourier_by_weight = (
            zero_coordinate_moment ** (self.step_bits - self.classes)
            * signed_coordinate_moment**self.classes
        )

        zero_to_zero = float(
            np.dot(self.counts, syndrome_fourier_by_weight)
            / self.state_space
        )
        zero_total = zero_coordinate_moment**self.step_bits
        zero_to_deterministic = max(0.0, zero_total - zero_to_zero)

        paired_total = float(
            syndrome_fourier_by_weight
            @ self.association
            @ live_output_by_weight
            / self.state_space
        )
        paired_nonzero = max(
            0.0,
            paired_total
            - zero_to_zero * float(live_output_by_weight[0]),
        )
        if zero_to_deterministic > 1e-300:
            deterministic_moment = paired_nonzero / zero_to_deterministic
        else:
            deterministic_moment = float(np.max(live_output_by_weight[1:]))

        uniform_live_moment = float(
            (
                np.dot(self.counts, live_output_by_weight)
                - live_output_by_weight[0]
            )
            / self.live_states
        )
        punctured_live_moment = (
            self.punctured_factor * uniform_live_moment
        )

        transfer = np.asarray(
            [
                [zero_to_zero, zero_to_deterministic, 0.0],
                [
                    deterministic_moment / self.live_states,
                    0.0,
                    deterministic_moment,
                ],
                [
                    punctured_live_moment / self.live_states,
                    0.0,
                    punctured_live_moment,
                ],
            ],
            dtype=np.float64,
        )
        details = {
            "z": z,
            "zero_to_zero": zero_to_zero,
            "zero_to_deterministic": zero_to_deterministic,
            "deterministic_moment": deterministic_moment,
            "uniform_live_moment": uniform_live_moment,
            "punctured_live_moment": punctured_live_moment,
        }
        return transfer, details

    def transfer_small_density(
        self, *, candidate_probability: float, surprisal: float
    ) -> tuple[np.ndarray, dict[str, float]]:
        """Return a four-state envelope whose puncturing loss is transient.

        The states are zero, deterministic, uniform-live, and
        punctured-live.  The last class remembers that one nonzero field
        element can be absent, but it does not charge that likelihood ratio
        on every later uniform-live epoch.
        """
        three_state, details = self.transfer(
            candidate_probability=candidate_probability,
            surprisal=surprisal,
        )
        deterministic_moment = details["deterministic_moment"]
        uniform_live_moment = details["uniform_live_moment"]
        z = details["z"]

        # A candidate position carries a fair bit.  Hence its actual input
        # bit is Bernoulli(candidate_probability / 2).  Fourier inversion
        # over im(A)=row(C) gives the exact probability that C(X)=0.
        kernel_probability = float(
            np.dot(
                self.counts,
                (1.0 - candidate_probability) ** self.classes,
            )
            / self.state_space
        )
        activation_probability = max(0.0, 1.0 - kernel_probability)

        active_lower_bound = z**self.step_bits * activation_probability
        deterministic_inactive = max(
            0.0, deterministic_moment - active_lower_bound
        )
        uniform_inactive = max(
            0.0, uniform_live_moment - active_lower_bound
        )
        survival_probability = (self.live_states - 1) / self.live_states
        deterministic_active = min(
            deterministic_moment, activation_probability
        ) * survival_probability
        uniform_active = min(
            uniform_live_moment, activation_probability
        ) * survival_probability
        punctured_inactive = self.punctured_factor * uniform_inactive
        punctured_active = self.punctured_factor * uniform_active
        deterministic_termination = (
            min(deterministic_moment, activation_probability)
            / self.live_states
        )
        uniform_termination = (
            min(uniform_live_moment, activation_probability)
            / self.live_states
        )
        punctured_termination = (
            self.punctured_factor * uniform_termination
        )

        zero_to_zero = details["zero_to_zero"]
        zero_to_deterministic = details["zero_to_deterministic"]
        transfer = np.asarray(
            [
                [zero_to_zero, zero_to_deterministic, 0.0, 0.0],
                [
                    deterministic_termination,
                    0.0,
                    deterministic_inactive,
                    deterministic_active,
                ],
                [
                    uniform_termination,
                    0.0,
                    uniform_inactive,
                    uniform_active,
                ],
                [
                    punctured_termination,
                    0.0,
                    punctured_inactive,
                    punctured_active,
                ],
            ],
            dtype=np.float64,
        )
        details.update(
            {
                "kernel_probability": kernel_probability,
                "activation_probability": activation_probability,
                "active_lower_bound": active_lower_bound,
                "deterministic_inactive_bound": deterministic_inactive,
                "uniform_inactive_bound": uniform_inactive,
                "punctured_inactive_bound": punctured_inactive,
                "deterministic_active_bound": deterministic_active,
                "uniform_active_bound": uniform_active,
                "punctured_active_bound": punctured_active,
            }
        )
        return transfer, details

    def objective(
        self,
        *,
        alpha: float,
        candidate_probability: float,
        surprisal: float,
        envelope: str = "three_state",
    ) -> tuple[float, dict[str, float]]:
        if envelope == "three_state":
            transfer, details = self.transfer(
                candidate_probability=candidate_probability,
                surprisal=surprisal,
            )
        elif envelope == "four_state":
            transfer, details = self.transfer_small_density(
                candidate_probability=candidate_probability,
                surprisal=surprisal,
            )
        else:
            raise ValueError(f"unknown envelope {envelope!r}")
        eigenvalues = np.linalg.eigvals(transfer)
        radius = float(np.max(np.abs(eigenvalues)))
        value = (
            alpha / 2.0
            + kl_bits(alpha, candidate_probability)
            + math.log2(radius) / self.step_bits
            + self.delta * surprisal / math.log(2.0)
        )
        details.update(
            {
                "perron_radius": radius,
                "objective_bits_per_output_bit": value,
            }
        )
        return value, details

    def optimize_variant(
        self, alpha: float, envelope: str
    ) -> dict[str, float | bool | str]:
        epsilon = 1e-10
        alpha_start = min(1.0 - epsilon, max(epsilon, alpha))
        scaled_probability = min(
            1.0 - epsilon, max(epsilon, 1.64 * alpha)
        )
        scaled_surprisal = max(epsilon, 1.65 * alpha)
        starts = [
            (scaled_probability, scaled_surprisal),
            (
                min(1.0 - epsilon, max(epsilon, 1.2 * alpha)),
                max(epsilon, alpha),
            ),
            (
                min(1.0 - epsilon, max(epsilon, 2.0 * alpha)),
                max(epsilon, 2.0 * alpha),
            ),
            (alpha_start, 0.02),
            (alpha_start, 0.2),
            (alpha_start, 0.8),
            (alpha_start, 1.5),
            (0.1, 0.2),
            (0.25, 0.8),
            (0.5, 1.2),
            (0.75, 1.8),
            (0.95, 2.2),
        ]
        best = None
        for start in starts:
            result = minimize(
                lambda point: self.objective(
                    alpha=alpha,
                    candidate_probability=float(point[0]),
                    surprisal=float(point[1]),
                    envelope=envelope,
                )[0],
                x0=np.asarray(start, dtype=np.float64),
                bounds=((epsilon, 1.0 - epsilon), (epsilon, 8.0)),
                method="Nelder-Mead",
                options={"maxiter": 500, "xatol": 1e-10, "fatol": 1e-12},
            )
            if best is None or float(result.fun) < float(best.fun):
                best = result
        assert best is not None
        probability = float(best.x[0])
        surprisal = float(best.x[1])
        value, details = self.objective(
            alpha=alpha,
            candidate_probability=probability,
            surprisal=surprisal,
            envelope=envelope,
        )
        return {
            "envelope": envelope,
            "alpha": alpha,
            "candidate_probability": probability,
            "surprisal": surprisal,
            "z": details["z"],
            "objective_bits_per_output_bit": value,
            "margin_bits_per_output_bit": -value,
            "optimizer_success": bool(best.success),
            "perron_radius": details["perron_radius"],
            "zero_to_zero": details["zero_to_zero"],
            "zero_to_deterministic": details["zero_to_deterministic"],
            "deterministic_moment": details["deterministic_moment"],
            "punctured_live_moment": details["punctured_live_moment"],
        }

    def optimize(self, alpha: float) -> dict[str, object]:
        variants = [
            self.optimize_variant(alpha, "three_state"),
            self.optimize_variant(alpha, "four_state"),
        ]
        best = min(
            variants,
            key=lambda row: float(row["objective_bits_per_output_bit"]),
        )
        return {
            **best,
            "variant_rows": variants,
        }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--selection", type=Path, default=DEFAULT_SELECTION)
    parser.add_argument("--delta", type=float, default=0.11)
    parser.add_argument(
        "--alphas",
        type=float,
        nargs="+",
        default=[
            0.0001,
            0.001,
            0.01,
            0.05,
            0.1,
            0.2,
            0.3,
            0.4,
            0.5,
            0.6,
            0.7,
            0.8,
            0.9,
            0.99,
            1.0,
        ],
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    payload = json.loads(args.selection.read_text(encoding="utf-8"))
    selected = payload["selected"]
    generator_words = [
        int(value, 16) for value in selected["A_generator_words_hex"]
    ]
    columns = [int(value, 16) for value in selected["B_columns_hex"]]
    step_bits = len(columns)
    state_bits = len(generator_words)
    verify_transpose_columns(
        generator_words=generator_words,
        columns=columns,
        output_bits=step_bits,
    )
    weights = state_code_weights(generator_words, step_bits)
    classes, counts, association = association_matrix(weights)

    envelope = DenseEnvelope(
        classes=classes,
        counts=counts,
        association=association,
        state_bits=state_bits,
        step_bits=step_bits,
        delta=args.delta,
    )
    rows = [envelope.optimize(alpha) for alpha in args.alphas]
    maximum = max(rows, key=lambda row: float(row["objective_bits_per_output_bit"]))
    result = {
        "schema": "rm2sub-dense-occupation-v2",
        "status": "EXACT_INTEGER_ASSOCIATION_BINARY64_SADDLE_DIAGNOSTIC",
        "selection": str(args.selection),
        "delta": args.delta,
        "state_bits": state_bits,
        "step_bits": step_bits,
        "checks": {
            "B_is_transpose_of_A": True,
            "association_axis_sums": True,
            "spectrum_mass": int(np.sum(counts)),
        },
        "weight_classes": [int(value) for value in classes],
        "weight_counts": [int(value) for value in counts],
        "association_matrix": association.tolist(),
        "rows": rows,
        "maximum_sampled_objective": maximum,
        "all_sampled_closed": all(
            float(row["objective_bits_per_output_bit"]) < 0.0 for row in rows
        ),
        "limitations": [
            "The association matrix and transpose checks use exact integers.",
            "The continuous saddle and Perron roots use binary64 arithmetic.",
            "A sampled alpha list is not a uniform interval certificate.",
            "The three-state and four-state transfers are entrywise envelopes, not equalities.",
            "The four-state active-event split uses a rigorous but loose probability bound.",
        ],
    }
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()
