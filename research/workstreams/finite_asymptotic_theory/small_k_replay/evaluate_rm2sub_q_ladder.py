#!/usr/bin/env python3
"""Screen fixed occupations with an exact RM spectrum envelope and RM2Sub.

The evaluator covers every message with a fixed number Q of active outer
rows.  An exact Holder change of measure handles every non-all-one RM
codeword.  The unique all-one codeword is handled by an exact mixture over
its row count.  Each mixture term may choose its own Chernoff tilt.

Nearest binary64 arithmetic makes the output diagnostic rather than outward
certified.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.special import logsumexp


HERE = Path(__file__).resolve().parent
WORKSTREAM = HERE.parent
REPOSITORY = WORKSTREAM.parent.parent
SCRIPTS = REPOSITORY / "scripts"
sys.path.insert(0, str(WORKSTREAM))
sys.path.insert(0, str(SCRIPTS))

from analyze_riffle_bitshuffle_splitstate_regular_bulk import (  # noqa: E402
    LOG2,
    binomial_transform,
    candidate_epoch_logs,
    load_nonzero_spectrum,
    load_uniform_nonactivation,
    log_choose,
    log_matmul_batch,
    log_matrix_power_moments_batch,
    regular_region_log_matrices,
    spectrum_bernoulli_envelope_log,
    splitstate_impulse_matrices,
)
from small_k_replay.evaluate_exact_spectra_q1_phase import (  # noqa: E402
    CONSTITUENTS,
    load_spectrum,
)
from small_k_replay.evaluate_rm2sub_primary_tranche import configuration  # noqa: E402


DEFAULT_OUTPUT = HERE / "rm2sub_q_ladder_d100.json"


def holder_grid() -> np.ndarray:
    """Orders dense near one and extending to the essential supremum."""
    return 1.0 + np.exp2(np.linspace(-20.0, 20.0, 1281))


def counting_density_statistics(
    spectrum_logs: np.ndarray, block_bits: int, probability: float
) -> tuple[np.ndarray, np.ndarray, float, int]:
    """Return exact moments of the counting density under Bernoulli(p)."""
    reference_logs = []
    density_logs = []
    weights = []
    for weight in range(1, block_bits):
        multiplicity = float(spectrum_logs[weight])
        if not math.isfinite(multiplicity):
            continue
        shell = log_choose(block_bits, weight)
        point_log = (
            weight * math.log(probability)
            + (block_bits - weight) * math.log1p(-probability)
        )
        reference_logs.append(shell + point_log)
        density_logs.append(multiplicity - shell - point_log)
        weights.append(weight)
    reference = np.asarray(reference_logs, dtype=np.float64)
    density = np.asarray(density_logs, dtype=np.float64)
    maximum_index = int(np.argmax(density))
    return reference, density, float(density[maximum_index]), weights[maximum_index]


def log_multinomial_rows(total_rows: int, regular: int, all_one: int) -> float:
    return log_choose(total_rows, all_one) + log_choose(
        total_rows - all_one, regular
    )


def log_position_multinomial(total: int, regular: int, forced: int) -> float:
    """Log number of disjoint regular and forced positions."""
    if regular < 0 or forced < 0 or regular + forced > total:
        return -math.inf
    return (
        math.lgamma(total + 1)
        - math.lgamma(regular + 1)
        - math.lgamma(forced + 1)
        - math.lgamma(total - regular - forced + 1)
    )


def mixed_epoch_log_matrices(
    impulses: np.ndarray, probability: float, maximum_occupation: int
) -> np.ndarray:
    """Epoch transfers with regular Bernoulli and forced-one positions."""
    result = np.full(
        (maximum_occupation + 1, maximum_occupation + 1, 2, 2),
        -math.inf,
    )
    for regular in range(maximum_occupation + 1):
        for forced in range(maximum_occupation - regular + 1):
            matrix = np.zeros((2, 2), dtype=np.float64)
            for live_regular in range(regular + 1):
                mass = (
                    math.comb(regular, live_regular)
                    * probability**live_regular
                    * (1.0 - probability) ** (regular - live_regular)
                )
                matrix += mass * impulses[forced + live_regular]
            with np.errstate(divide="ignore"):
                result[regular, forced] = np.log(matrix)
    return result


def mixed_region_log_matrices(
    epoch_logs: np.ndarray,
    epoch_bits: int,
    epochs_per_region: int,
    maximum_occupation: int,
) -> np.ndarray:
    """Average products with two disjoint position types in one region."""
    current = epoch_logs.copy()
    for completed_epochs in range(1, epochs_per_region):
        old_slots = completed_epochs * epoch_bits
        new_slots = old_slots + epoch_bits
        updated = np.full_like(current, -math.inf)
        for total_regular in range(maximum_occupation + 1):
            for total_forced in range(maximum_occupation - total_regular + 1):
                denominator = log_position_multinomial(
                    new_slots, total_regular, total_forced
                )
                for epoch_regular in range(total_regular + 1):
                    previous_regular = total_regular - epoch_regular
                    for epoch_forced in range(total_forced + 1):
                        previous_forced = total_forced - epoch_forced
                        if epoch_regular + epoch_forced > epoch_bits:
                            continue
                        if previous_regular + previous_forced > old_slots:
                            continue
                        weight = (
                            log_position_multinomial(
                                old_slots, previous_regular, previous_forced
                            )
                            + log_position_multinomial(
                                epoch_bits, epoch_regular, epoch_forced
                            )
                            - denominator
                        )
                        product = log_matmul_batch(
                            current[previous_regular, previous_forced][None, ...],
                            epoch_logs[epoch_regular, epoch_forced],
                        )[0]
                        updated[total_regular, total_forced] = np.logaddexp(
                            updated[total_regular, total_forced], product + weight
                        )
        current = updated
    return current


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step-bits", type=int, required=True)
    parser.add_argument("--state-bits", type=int, required=True)
    parser.add_argument("--message-exponent", type=int, default=16)
    parser.add_argument("--occupations", type=int, nargs="+", default=[3, 4, 8])
    parser.add_argument("--distance-numerator", type=int, default=100)
    parser.add_argument("--distance-denominator", type=int, default=1000)
    parser.add_argument("--outer-weight-lower", type=int, default=1)
    parser.add_argument("--outer-weight-upper", type=int)
    parser.add_argument(
        "--log-surprisals",
        type=float,
        nargs="+",
        default=[-8.0, -7.5, -7.0, -6.5, -6.0, -5.5, -5.0, -4.5,
                 -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0],
    )
    parser.add_argument(
        "--candidate-probabilities", type=float, nargs="+", default=[0.5]
    )
    parser.add_argument(
        "--regular-only",
        action="store_true",
        help="screen only mixtures with no all-one outer row",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    if args.step_bits not in (64, 128, 256):
        parser.error("this tranche supports epoch lengths 64, 128, and 256")
    expected_persistence = args.state_bits + int(math.log2(args.step_bits))
    config = configuration(args.step_bits, expected_persistence)
    if int(config["state_bits"]) != args.state_bits:
        raise AssertionError("configuration lookup returned the wrong state size")

    constituent = CONSTITUENTS["rm49"]
    message_bits = 1 << args.message_exponent
    if message_bits % constituent.dimension:
        parser.error("outer dimension does not divide the message length")
    outer_rows = message_bits // constituent.dimension
    if outer_rows % args.step_bits:
        parser.error("epoch length does not divide the transposed region")
    if any(not 3 <= occupation <= outer_rows for occupation in args.occupations):
        parser.error("occupations must lie in [3,L]")
    if any(not 0.0 < probability < 1.0 for probability in args.candidate_probabilities):
        parser.error("candidate probabilities must lie in (0,1)")

    output_bits = 2 * message_bits
    bad_weight = math.floor(
        args.distance_numerator * output_bits / args.distance_denominator
    )
    maximum_occupation = max(args.occupations)
    outer_weight_upper = (
        constituent.block_bits - 1
        if args.outer_weight_upper is None
        else args.outer_weight_upper
    )
    if not 1 <= args.outer_weight_lower <= outer_weight_upper < constituent.block_bits:
        parser.error("outer weight band must lie in [1,B-1]")
    if (
        (args.outer_weight_lower != 1 or outer_weight_upper != constituent.block_bits - 1)
        and not args.regular_only
    ):
        parser.error("a restricted outer band requires --regular-only")
    outer_spectrum = np.full(constituent.block_bits + 1, -math.inf)
    for weight, count in load_spectrum(constituent).items():
        if count and (
            weight == 0
            or args.outer_weight_lower <= weight <= outer_weight_upper
        ):
            outer_spectrum[weight] = math.log(count)
    envelopes = {
        probability: spectrum_bernoulli_envelope_log(
            constituent.block_bits, outer_spectrum, probability, 0
        )
        for probability in args.candidate_probabilities
    }
    holder_orders = holder_grid()
    holder_statistics = {}
    for probability in args.candidate_probabilities:
        (
            density_reference_logs,
            density_logs,
            maximum_density_log,
            maximum_density_weight,
        ) = counting_density_statistics(
            outer_spectrum, constituent.block_bits, probability
        )
        holder_statistics[probability] = {
            "log_moments": np.asarray(
                [
                    float(
                        logsumexp(
                            density_reference_logs + order * density_logs
                        )
                    )
                    for order in holder_orders
                ]
            ),
            "maximum_density_log": maximum_density_log,
            "maximum_density_weight": maximum_density_weight,
        }
    nonactivation = load_uniform_nonactivation(Path(config["b"]), args.step_bits)
    live_spectrum = load_nonzero_spectrum(
        Path(config["a"]), args.step_bits, args.state_bits
    )
    epochs_per_region = outer_rows // args.step_bits

    keys = [
        (occupation, all_one)
        for occupation in args.occupations
        for all_one in (
            [0] if args.regular_only else range(occupation + 1)
        )
    ]
    best = {key: math.inf for key in keys}
    witnesses = {
        key: {
            "log_surprisal": math.nan,
            "candidate_probability": math.nan,
            "holder_order": math.nan,
        }
        for key in keys
    }

    for log_surprisal in args.log_surprisals:
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        impulses = splitstate_impulse_matrices(
            z=z,
            step_bits=args.step_bits,
            state_bits=args.state_bits,
            constituent_distance=int(config["a_distance"]),
            live_moment_order=3,
            live_model="factorial-moment",
            nonactivation=nonactivation,
            live_spectrum=live_spectrum,
        )
        for probability in args.candidate_probabilities:
            holder = holder_statistics[probability]
            density_log_moments = holder["log_moments"]
            maximum_density_log = float(holder["maximum_density_log"])
            candidate_logs = candidate_epoch_logs(
                impulses, binomial_transform(args.step_bits, probability)
            )
            regular_reference = regular_region_log_matrices(
                candidate_logs,
                args.step_bits,
                epochs_per_region,
                maximum_occupation,
            )[: maximum_occupation + 1]
            if args.regular_only:
                regular_moments = log_matrix_power_moments_batch(
                    regular_reference, constituent.block_bits
                )
                moments = regular_moments[:, None]
            else:
                epoch_logs = mixed_epoch_log_matrices(
                    impulses, probability, maximum_occupation
                )
                region_logs = mixed_region_log_matrices(
                    epoch_logs,
                    args.step_bits,
                    epochs_per_region,
                    maximum_occupation,
                )
                flat_region_logs = region_logs.reshape(
                    (-1, *region_logs.shape[-2:])
                )
                moments = log_matrix_power_moments_batch(
                    flat_region_logs, constituent.block_bits
                ).reshape(region_logs.shape[:2])

                # The zero-forced slice is the established one-colour recurrence.
                finite = np.isfinite(regular_reference) & np.isfinite(
                    region_logs[:, 0]
                )
                if np.any(finite):
                    error = float(
                        np.max(
                            np.abs(
                                regular_reference[finite]
                                - region_logs[:, 0][finite]
                            )
                        )
                    )
                    if error > 2e-10:
                        raise ArithmeticError(
                            "two-colour recurrence disagrees with regular "
                            f"recurrence: {error}"
                        )
            for occupation, all_one in keys:
                regular = occupation - all_one
                inner_moment = (
                    moments[regular, 0]
                    if args.regular_only
                    else moments[regular, all_one]
                )
                inner = min(
                    0.0,
                    float(inner_moment) + bad_weight * surprisal,
                )
                if regular == 0:
                    candidate = inner
                    holder_order: float | str = 1.0
                else:
                    corrections = (
                        regular * density_log_moments / holder_orders
                        + (1.0 - 1.0 / holder_orders) * inner
                    )
                    holder_index = int(np.argmin(corrections))
                    candidate = float(corrections[holder_index])
                    holder_order = float(holder_orders[holder_index])
                    infinity_candidate = regular * maximum_density_log + inner
                    if infinity_candidate < candidate:
                        candidate = infinity_candidate
                        holder_order = "infinity"
                key = (occupation, all_one)
                if candidate < best[key]:
                    best[key] = candidate
                    witnesses[key] = {
                        "log_surprisal": log_surprisal,
                        "candidate_probability": probability,
                        "holder_order": holder_order,
                    }

    occupation_rows = []
    for occupation in args.occupations:
        mixture_rows = []
        mixture_logs = []
        for all_one in ([0] if args.regular_only else range(occupation + 1)):
            regular = occupation - all_one
            outer_log = log_multinomial_rows(outer_rows, regular, all_one)
            contribution = outer_log + best[(occupation, all_one)]
            mixture_logs.append(contribution)
            mixture_rows.append(
                {
                    "regular_rows": regular,
                    "all_one_rows": all_one,
                    "outer_log2_row_choices": outer_log / LOG2,
                    "inner_plus_change_of_measure_log2_upper": best[(occupation, all_one)]
                    / LOG2,
                    "pointwise_log2_upper": contribution / LOG2,
                    "pointwise_margin_bits": -contribution / LOG2,
                    **witnesses[(occupation, all_one)],
                }
            )
        aggregate = float(logsumexp(np.asarray(mixture_logs)))
        dominant = max(mixture_rows, key=lambda row: row["pointwise_log2_upper"])
        occupation_rows.append(
            {
                "occupation": occupation,
                "log2_upper_diagnostic": aggregate / LOG2,
                "margin_bits_diagnostic": -aggregate / LOG2,
                "dominant_all_one_rows": dominant["all_one_rows"],
                "dominant_regular_rows": dominant["regular_rows"],
                "mixture_rows": mixture_rows,
            }
        )

    aggregate = float(
        logsumexp(
            np.asarray(
                [-row["margin_bits_diagnostic"] * LOG2 for row in occupation_rows]
            )
        )
    )
    payload = {
        "schema": "rm2sub-fixed-occupation-rm-holder-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "construction": (
            "one fixed RM(4,9) constituent repeated in every row, uniform "
            "routing, and one fixed audited RM2Sub A/B pair"
        ),
        "probability_space": {
            "outer": "one fixed authenticated RM(4,9) constituent",
            "routing": "independent uniform row-coordinate and region permutations",
            "inner": "independent nonzero field scalar in every RM2Sub epoch",
        },
        "parameters": {
            "message_exponent": args.message_exponent,
            "message_bits": message_bits,
            "output_bits": output_bits,
            "outer_rows": outer_rows,
            "outer_block_bits": constituent.block_bits,
            "outer_dimension": constituent.dimension,
            "outer_weight_lower": args.outer_weight_lower,
            "outer_weight_upper": outer_weight_upper,
            "bad_weight": bad_weight,
            "distance_numerator": args.distance_numerator,
            "distance_denominator": args.distance_denominator,
            "step_bits": args.step_bits,
            "state_bits": args.state_bits,
            "persistence_exponent": expected_persistence,
            "epochs_per_region": epochs_per_region,
            "occupations": args.occupations,
            "log_surprisals": args.log_surprisals,
            "candidate_probabilities": args.candidate_probabilities,
            "a_minimum_distance": config["a_distance"],
            "kernel_minimum_distance": config["kernel_distance"],
            "kernel_weight_four": config["kernel_weight_four"],
        },
        "comparison_outer_envelopes_not_used_by_holder_bound": [
            {
                "candidate_probability": probability,
                "log2_mass": envelope[0] / LOG2,
                "maximizing_weight": envelope[1],
            }
            for probability, envelope in envelopes.items()
        ],
        "holder_change_of_measure": {
            "reference": "the listed Bernoulli product measure on F_2^512 for each regular row",
            "density": "exact nonzero, non-all-one RM(4,9) counting measure after a uniform coordinate permutation",
            "order_grid_size": int(len(holder_orders)),
            "order_grid_minimum": float(holder_orders[0]),
            "order_grid_maximum": float(holder_orders[-1]),
            "by_candidate_probability": [
                {
                    "candidate_probability": probability,
                    "maximum_log2_density": float(
                        holder_statistics[probability]["maximum_density_log"]
                    )
                    / LOG2,
                    "maximum_density_weight": holder_statistics[probability][
                        "maximum_density_weight"
                    ],
                }
                for probability in args.candidate_probabilities
            ],
        },
        "proof_reduction": [
            "Separate the unique all-one codeword from the other nonzero RM codewords.",
            "Express each routed non-all-one RM counting measure by its exact density relative to the uniform product measure.",
            "For fixed Q, sum exactly over the number of rows carrying the all-one codeword.",
            "A positive two-colour hypergeometric recurrence places ordinary candidate rows and forced all-one rows without overlap in every transposed region.",
            "Apply Holder across the independent regular-row density factors and the uniform-reference bad event.",
            "Optimize the Chernoff and Holder witnesses separately for every mixture term before summing positive terms.",
        ],
        "occupation_rows": occupation_rows,
        "partial_log2_upper": aggregate / LOG2,
        "partial_margin_bits": -aggregate / LOG2,
        "limitations": [
            "Nearest binary64 arithmetic is not outward rounded.",
            "The finite tilt and Bernoulli grids are not asserted optimal.",
            "The finite Holder-order grid is not asserted optimal.",
            "Only the displayed occupations are covered.",
            *(
                ["Only the no-all-one-row face is covered."]
                if args.regular_only
                else []
            ),
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    for row in occupation_rows:
        print(
            f"Q,{row['occupation']},margin,{row['margin_bits_diagnostic']:.6f},"
            f"dominant_all_one,{row['dominant_all_one_rows']}"
        )
    print(f"partial_margin,{payload['partial_margin_bits']:.6f}")
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
