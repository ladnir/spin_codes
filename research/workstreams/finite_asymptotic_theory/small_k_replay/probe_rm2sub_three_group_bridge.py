#!/usr/bin/env python3
"""Cover all low/central/high RM compositions with the sparse RM2Sub bound.

The exact nonzero RM(4,9) spectrum is partitioned into three groups.  Each
group is dominated by a Bernoulli product measure with probability 1/4, 1/2,
or 3/4.  For a fixed group composition, independent thinning followed by the
uniform region permutation leaves a uniform support conditional on its total
size.  A three-binomial convolution therefore feeds the established positive
forced-support RM2Sub recurrence.

Nearest-binary64 arithmetic and a finite Chernoff grid make the output a
diagnostic rather than an outward certificate.
"""

from __future__ import annotations

import argparse
import hashlib
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
    load_nonzero_spectrum,
    load_uniform_nonactivation,
    regular_region_log_matrices,
    splitstate_impulse_matrices,
)
from small_k_replay.evaluate_exact_spectra_q1_phase import (  # noqa: E402
    CONSTITUENTS,
    load_spectrum,
)
from small_k_replay.evaluate_rm2sub_primary_tranche import configuration  # noqa: E402
from small_k_replay.probe_rm2sub_two_band_bridge import (  # noqa: E402
    band_envelope_at_probability,
    band_density_statistics,
    holder_grid,
    minimize_band_envelope,
    parse_band_range,
    verify_exact_nonactivation,
)


GROUPS = (
    ("low", 1, 47, 0.25),
    ("central", 48, 416, 0.50),
    ("high", 417, 512, 0.75),
)
DEFAULT_OUTPUT = HERE / "rm2sub_three_group_bridge_probe_d100.json"


def compositions_three(total: int):
    for low in range(total + 1):
        for central in range(total - low + 1):
            yield low, central, total - low - central


def log_multinomial(total: int, counts: tuple[int, ...]) -> float:
    if sum(counts) > total or min(counts) < 0:
        return -math.inf
    return (
        math.lgamma(total + 1)
        - math.lgamma(total - sum(counts) + 1)
        - sum(math.lgamma(count + 1) for count in counts)
    )


def three_binomial_distributions(
    compositions: list[tuple[int, int, int]],
    maximum_occupation: int,
    probabilities: tuple[float, float, float] | None = None,
) -> np.ndarray:
    selected_probabilities = (
        tuple(group[3] for group in GROUPS)
        if probabilities is None
        else probabilities
    )
    transforms = [
        binomial_transform(maximum_occupation, probability)
        for probability in selected_probabilities
    ]
    result = np.zeros(
        (len(compositions), maximum_occupation + 1), dtype=np.float64
    )
    for index, counts in enumerate(compositions):
        distribution = np.asarray([1.0])
        for count, transform in zip(counts, transforms):
            distribution = np.convolve(
                distribution, transform[count, : count + 1]
            )
        result[index, : distribution.size] = distribution
    error = float(np.max(np.abs(np.sum(result, axis=1) - 1.0)))
    if error > 2e-13:
        raise ArithmeticError(f"three-binomial mass error: {error}")
    return result


def log_matmul_batch(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    result = np.empty_like(left)
    result[:, 0, 0] = np.logaddexp(
        left[:, 0, 0] + right[:, 0, 0],
        left[:, 0, 1] + right[:, 1, 0],
    )
    result[:, 0, 1] = np.logaddexp(
        left[:, 0, 0] + right[:, 0, 1],
        left[:, 0, 1] + right[:, 1, 1],
    )
    result[:, 1, 0] = np.logaddexp(
        left[:, 1, 0] + right[:, 0, 0],
        left[:, 1, 1] + right[:, 1, 0],
    )
    result[:, 1, 1] = np.logaddexp(
        left[:, 1, 0] + right[:, 0, 1],
        left[:, 1, 1] + right[:, 1, 1],
    )
    return result


def log_row_matmul_batch(rows: np.ndarray, matrices: np.ndarray) -> np.ndarray:
    result = np.empty_like(rows)
    result[:, 0] = np.logaddexp(
        rows[:, 0] + matrices[:, 0, 0],
        rows[:, 1] + matrices[:, 1, 0],
    )
    result[:, 1] = np.logaddexp(
        rows[:, 0] + matrices[:, 0, 1],
        rows[:, 1] + matrices[:, 1, 1],
    )
    return result


def log_matrix_power_moments_binary(
    matrices: np.ndarray, power: int
) -> np.ndarray:
    """Return log(e_0^T M^power 1) for a batch of positive 2-by-2 matrices."""
    rows = np.full((matrices.shape[0], 2), -math.inf, dtype=np.float64)
    rows[:, 0] = 0.0
    powers = matrices.copy()
    exponent = power
    while exponent:
        if exponent & 1:
            rows = log_row_matmul_batch(rows, powers)
        exponent >>= 1
        if exponent:
            powers = log_matmul_batch(powers, powers)
    return np.logaddexp(rows[:, 0], rows[:, 1])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step-bits", type=int, default=64)
    parser.add_argument("--state-bits", type=int, default=14)
    parser.add_argument("--message-exponent", type=int, default=16)
    parser.add_argument("--minimum-occupation", type=int, default=30)
    parser.add_argument("--maximum-occupation", type=int, default=52)
    parser.add_argument("--distance-numerator", type=int, default=100)
    parser.add_argument("--distance-denominator", type=int, default=1000)
    parser.add_argument(
        "--log-surprisals",
        type=float,
        nargs="+",
        default=[value / 10.0 for value in range(-40, 1)],
    )
    parser.add_argument(
        "--change-of-measure",
        choices=["envelope", "holder"],
        default="envelope",
    )
    parser.add_argument("--holder-chunk-size", type=int, default=2048)
    parser.add_argument(
        "--reference-probabilities",
        type=float,
        nargs=3,
        help=(
            "override the three Bernoulli references; omit to retain the "
            "standard 1/4, 1/2, 3/4 schedule"
        ),
    )
    parser.add_argument(
        "--group-ranges",
        type=parse_band_range,
        nargs=3,
        help="override the three contiguous spectrum groups",
    )
    parser.add_argument(
        "--optimize-reference-probabilities",
        action="store_true",
        help="minimize each group's pointwise density envelope",
    )
    parser.add_argument(
        "--common-witness-only",
        action="store_true",
        help=(
            "store only one Chernoff witness and aggregate bound per "
            "occupation; this avoids materializing the per-composition table"
        ),
    )
    parser.add_argument(
        "--packed-witness-output",
        type=Path,
        help=(
            "write one signed-byte, exact-tenth Chernoff witness per "
            "composition in canonical occupation/low/central order"
        ),
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    if args.step_bits not in (64, 128, 256):
        parser.error("supported epoch lengths are 64, 128, and 256")
    if not 1 <= args.minimum_occupation <= args.maximum_occupation:
        parser.error("the occupation interval must be nonempty")
    if args.reference_probabilities is not None and args.optimize_reference_probabilities:
        parser.error("choose explicit or optimized reference probabilities, not both")
    if args.reference_probabilities is not None and any(
        not 0.0 < probability < 1.0
        for probability in args.reference_probabilities
    ):
        parser.error("reference probabilities must lie in (0,1)")

    persistence = args.state_bits + int(math.log2(args.step_bits))
    config = configuration(args.step_bits, persistence)
    if int(config["state_bits"]) != args.state_bits:
        raise AssertionError("configuration lookup returned the wrong state size")
    constituent = CONSTITUENTS["rm49"]
    message_bits = 1 << args.message_exponent
    if message_bits % constituent.dimension:
        parser.error("outer dimension does not divide the message length")
    outer_rows = message_bits // constituent.dimension
    if outer_rows % args.step_bits:
        parser.error("epoch length does not divide one transposed region")
    if args.maximum_occupation > outer_rows:
        parser.error("maximum occupation exceeds the number of outer rows")

    output_bits = 2 * message_bits
    bad_weight = math.floor(
        args.distance_numerator * output_bits / args.distance_denominator
    )
    spectrum = load_spectrum(constituent)
    group_rows = []
    envelope_logs = []
    holder_orders = holder_grid()
    holder_log_moments = []
    covered_mass = 0
    base_groups = GROUPS
    if args.group_ranges is not None:
        if (
            args.group_ranges[0][0] != 1
            or args.group_ranges[-1][1] != constituent.block_bits
            or any(
                args.group_ranges[index][1] + 1
                != args.group_ranges[index + 1][0]
                for index in range(2)
            )
        ):
            parser.error("group ranges must partition weights 1 through 512")
        base_groups = tuple(
            (GROUPS[index][0], lower, upper, GROUPS[index][3])
            for index, (lower, upper) in enumerate(args.group_ranges)
        )
    selected_groups = []
    for index, (name, lower, upper, default_probability) in enumerate(base_groups):
        if args.optimize_reference_probabilities:
            probability, _, _ = minimize_band_envelope(
                spectrum, constituent.block_bits, lower, upper
            )
        elif args.reference_probabilities is not None:
            probability = args.reference_probabilities[index]
        else:
            probability = default_probability
        selected_groups.append((name, lower, upper, probability))
    for name, lower, upper, probability in selected_groups:
        envelope, maximizing_weight = band_envelope_at_probability(
            spectrum, constituent.block_bits, lower, upper, probability
        )
        mass = sum(
            count
            for weight, count in spectrum.items()
            if lower <= weight <= upper
        )
        covered_mass += mass
        envelope_logs.append(envelope)
        moments, maximum_density, maximum_weight = band_density_statistics(
            spectrum,
            constituent.block_bits,
            lower,
            upper,
            probability,
            holder_orders,
        )
        if abs(maximum_density - envelope) > 2e-12:
            raise ArithmeticError("Holder endpoint disagrees with envelope")
        if maximum_weight != maximizing_weight:
            raise ArithmeticError("Holder endpoint has the wrong maximizing weight")
        holder_log_moments.append(moments)
        group_rows.append(
            {
                "name": name,
                "lower": lower,
                "upper": upper,
                "reference_probability": probability,
                "maximum_density_log2": envelope / LOG2,
                "maximum_density_weight": maximizing_weight,
                "exact_mass": str(mass),
            }
        )
    if covered_mass != (1 << constituent.dimension) - 1:
        raise AssertionError("three groups do not partition the nonzero spectrum")

    compositions = [
        composition
        for occupation in range(
            args.minimum_occupation, args.maximum_occupation + 1
        )
        for composition in compositions_three(occupation)
    ]
    composition_array = np.asarray(compositions, dtype=np.int64)
    distributions = three_binomial_distributions(
        compositions,
        args.maximum_occupation,
        tuple(group[3] for group in selected_groups),
    )
    location_logs = np.asarray(
        [
            log_multinomial(outer_rows, composition)
            for composition in compositions
        ],
        dtype=np.float64,
    )
    envelope_corrections = composition_array @ np.asarray(
        envelope_logs, dtype=np.float64
    )
    holder_per_row = (
        np.asarray(holder_log_moments, dtype=np.float64)
        / holder_orders[np.newaxis, :]
    )

    nonactivation_path = Path(config["b"])
    nonactivation = load_uniform_nonactivation(nonactivation_path, args.step_bits)
    verify_exact_nonactivation(nonactivation_path, nonactivation)
    live_spectrum = load_nonzero_spectrum(
        Path(config["a"]), args.step_bits, args.state_bits
    )
    epochs_per_region = outer_rows // args.step_bits
    best = np.full(len(compositions), math.inf, dtype=np.float64)
    witnesses = np.full(len(compositions), math.nan, dtype=np.float64)
    holder_witnesses = np.full(len(compositions), math.nan, dtype=np.float64)
    occupation_slices: dict[int, slice] = {}
    cursor = 0
    for occupation in range(args.minimum_occupation, args.maximum_occupation + 1):
        count = (occupation + 1) * (occupation + 2) // 2
        occupation_slices[occupation] = slice(cursor, cursor + count)
        cursor += count
    if cursor != len(compositions):
        raise AssertionError("occupation slices do not cover the compositions")
    common_best = {
        occupation: math.inf for occupation in occupation_slices
    }
    common_witness = {
        occupation: math.nan for occupation in occupation_slices
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
            live_model="support-averaged-preaddmul",
            nonactivation=nonactivation,
            live_spectrum=live_spectrum,
        )
        for input_weight in range(1, args.step_bits + 1):
            impulses[input_weight, 0, 1] = max(
                0.0, 1.0 - float(nonactivation[input_weight])
            ) * z**input_weight
        with np.errstate(divide="ignore"):
            forced_epoch_logs = np.log(impulses)
        forced_region_logs = regular_region_log_matrices(
            forced_epoch_logs,
            args.step_bits,
            epochs_per_region,
            args.maximum_occupation,
        )[: args.maximum_occupation + 1]
        forced_region = np.exp(forced_region_logs).reshape(
            (args.maximum_occupation + 1, 4)
        )
        reference_regions = (distributions @ forced_region).reshape((-1, 2, 2))
        with np.errstate(divide="ignore"):
            reference_logs = np.log(reference_regions)
        moments = log_matrix_power_moments_binary(
            reference_logs, constituent.block_bits
        )
        inner = np.minimum(0.0, moments + bad_weight * surprisal)
        if args.change_of_measure == "envelope":
            corrections = envelope_corrections + inner
            selected_holder_orders = np.full(len(compositions), math.inf)
        else:
            corrections = np.empty(len(compositions), dtype=np.float64)
            selected_holder_orders = np.empty(len(compositions), dtype=np.float64)
            holder_coefficients = 1.0 - 1.0 / holder_orders
            for start in range(0, len(compositions), args.holder_chunk_size):
                stop = min(start + args.holder_chunk_size, len(compositions))
                finite = (
                    composition_array[start:stop] @ holder_per_row
                    + inner[start:stop, np.newaxis] * holder_coefficients
                )
                indexes = np.argmin(finite, axis=1)
                chunk_rows = np.arange(stop - start)
                finite_best = finite[chunk_rows, indexes]
                infinity = envelope_corrections[start:stop] + inner[start:stop]
                use_infinity = infinity < finite_best
                corrections[start:stop] = np.where(
                    use_infinity, infinity, finite_best
                )
                selected_holder_orders[start:stop] = np.where(
                    use_infinity, math.inf, holder_orders[indexes]
                )
        candidates = location_logs + corrections
        for occupation, selected in occupation_slices.items():
            common_candidate = float(logsumexp(candidates[selected]))
            if common_candidate < common_best[occupation]:
                common_best[occupation] = common_candidate
                common_witness[occupation] = log_surprisal
        improved = candidates < best
        best[improved] = candidates[improved]
        witnesses[improved] = log_surprisal
        holder_witnesses[improved] = selected_holder_orders[improved]
        print(
            f"tilt,{log_surprisal:.2f},"
            f"current_min_margin,{(-float(np.max(best)) / LOG2):.6f}",
            flush=True,
        )

    if args.common_witness_only:
        packed_witnesses = None
        if args.packed_witness_output is not None:
            witness_tenths = np.rint(10.0 * witnesses).astype(np.int8)
            if not np.allclose(
                witnesses, witness_tenths.astype(np.float64) / 10.0,
                rtol=0.0, atol=1e-12,
            ):
                raise AssertionError("a selected witness is not on the tenth grid")
            raw_witnesses = witness_tenths.tobytes(order="C")
            args.packed_witness_output.write_bytes(raw_witnesses)
            packed_witnesses = {
                "path": str(args.packed_witness_output),
                "sha256": hashlib.sha256(raw_witnesses).hexdigest(),
                "encoding": "signed int8 storing 10*log_surprisal",
                "composition_order": "occupation ascending, then low ascending, then central ascending",
                "count": len(witness_tenths),
            }
        common_rows = [
            {
                "occupation": occupation,
                "composition_count": (occupation + 1) * (occupation + 2) // 2,
                "log2_upper_diagnostic": common_best[occupation] / LOG2,
                "margin_bits_diagnostic": -common_best[occupation] / LOG2,
                "log_surprisal": common_witness[occupation],
            }
            for occupation in occupation_slices
        ]
        common_aggregate = float(
            logsumexp(np.asarray([common_best[q] for q in occupation_slices]))
        )
        optimized_rows = [
            {
                "occupation": occupation,
                "composition_count": (occupation + 1) * (occupation + 2) // 2,
                "log2_upper_diagnostic": float(logsumexp(best[selected])) / LOG2,
                "margin_bits_diagnostic": -float(logsumexp(best[selected])) / LOG2,
            }
            for occupation, selected in occupation_slices.items()
        ]
        optimized_aggregate = float(
            logsumexp(
                np.asarray(
                    [
                        float(logsumexp(best[selected]))
                        for selected in occupation_slices.values()
                    ]
                )
            )
        )
        payload = {
            "schema": "rm2sub-fixed-rm-three-group-common-witness-v1",
            "status": "BINARY64_DIAGNOSTIC",
            "construction": (
                "one fixed RM(4,9) constituent repeated in every row, uniform "
                "routing, and one fixed audited RM2Sub A/B pair"
            ),
            "parameters": {
                "message_exponent": args.message_exponent,
                "message_bits": message_bits,
                "output_bits": output_bits,
                "outer_rows": outer_rows,
                "outer_block_bits": constituent.block_bits,
                "step_bits": args.step_bits,
                "state_bits": args.state_bits,
                "bad_weight": bad_weight,
                "distance_numerator": args.distance_numerator,
                "distance_denominator": args.distance_denominator,
                "occupation_interval": [
                    args.minimum_occupation,
                    args.maximum_occupation,
                ],
                "log_surprisals": args.log_surprisals,
                "live_model": "support-averaged-preaddmul",
                "change_of_measure": (
                    "pointwise-density-envelope"
                    if args.change_of_measure == "envelope"
                    else "joint-holder"
                ),
                "exact_zero_state_destination_split": True,
                "exact_nonactivation_verified_from_kernel_shell_counts": True,
            },
            "groups": group_rows,
            "occupation_rows": common_rows,
            "interval_log2_upper_diagnostic": common_aggregate / LOG2,
            "interval_margin_bits_diagnostic": -common_aggregate / LOG2,
            "per_composition_witness_rows": optimized_rows,
            "per_composition_witness_interval_log2_upper_diagnostic": optimized_aggregate / LOG2,
            "per_composition_witness_interval_margin_bits_diagnostic": -optimized_aggregate / LOG2,
            "packed_witnesses": packed_witnesses,
            "limitations": [
                "Nearest binary64 arithmetic is not outward rounded.",
                "Each occupation uses one displayed finite-grid Chernoff witness for every composition.",
                "Only the displayed occupation interval is covered.",
            ],
        }
        args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        print(f"common_interval_margin,{payload['interval_margin_bits_diagnostic']:.6f}")
        print(
            "per_composition_interval_margin,"
            f"{payload['per_composition_witness_interval_margin_bits_diagnostic']:.6f}"
        )
        print(f"output={args.output}")
        return

    occupation_rows = []
    occupation_values = []
    composition_rows = []
    for index, composition in enumerate(compositions):
        value = float(best[index])
        composition_rows.append(
            {
                "low_rows": composition[0],
                "central_rows": composition[1],
                "high_rows": composition[2],
                "occupation": sum(composition),
                "log2_upper_diagnostic": value / LOG2,
                "margin_bits_diagnostic": -value / LOG2,
                "gap_to_40_bits": -value / LOG2 - 40.0,
                "log_surprisal": float(witnesses[index]),
                "holder_order": (
                    "infinity"
                    if math.isinf(float(holder_witnesses[index]))
                    else float(holder_witnesses[index])
                ),
            }
        )
    for occupation in range(
        args.minimum_occupation, args.maximum_occupation + 1
    ):
        indexes = np.flatnonzero(np.sum(composition_array, axis=1) == occupation)
        aggregate = float(logsumexp(best[indexes]))
        occupation_values.append(aggregate)
        dominant_index = int(indexes[int(np.argmax(best[indexes]))])
        dominant = composition_rows[dominant_index]
        occupation_rows.append(
            {
                "occupation": occupation,
                "composition_count": int(indexes.size),
                "log2_upper_diagnostic": aggregate / LOG2,
                "margin_bits_diagnostic": -aggregate / LOG2,
                "gap_to_40_bits": -aggregate / LOG2 - 40.0,
                "dominant_composition": {
                    "low_rows": dominant["low_rows"],
                    "central_rows": dominant["central_rows"],
                    "high_rows": dominant["high_rows"],
                    "pointwise_margin_bits_diagnostic": dominant[
                        "margin_bits_diagnostic"
                    ],
                },
            }
        )
    aggregate = float(logsumexp(np.asarray(occupation_values)))
    payload = {
        "schema": "rm2sub-fixed-rm-three-group-sparse-bridge-v1",
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
            "step_bits": args.step_bits,
            "state_bits": args.state_bits,
            "epochs_per_region": epochs_per_region,
            "bad_weight": bad_weight,
            "distance_numerator": args.distance_numerator,
            "distance_denominator": args.distance_denominator,
            "occupation_interval": [
                args.minimum_occupation,
                args.maximum_occupation,
            ],
            "log_surprisals": args.log_surprisals,
            "live_model": "support-averaged-preaddmul",
            "change_of_measure": (
                "pointwise-density-envelope"
                if args.change_of_measure == "envelope"
                else "joint-holder"
            ),
            "holder_order_grid": (
                holder_orders.tolist()
                if args.change_of_measure == "holder"
                else None
            ),
            "exact_zero_state_destination_split": True,
            "exact_nonactivation_verified_from_kernel_shell_counts": True,
        },
        "groups": group_rows,
        "proof_reduction": [
            "Partition the exact nonzero outer spectrum into the displayed low, central, and high groups.",
            "Express each group as a density relative to its displayed Bernoulli product measure.",
            "For a fixed three-group composition, uniformly permute the three disjoint candidate-position sets in every region.",
            "Thin each position independently with its group's reference probability.",
            "Conditional on total live weight, permutation symmetry makes the live support uniform.",
            "Convolve the three binomial laws and apply the positive forced-support RM2Sub recurrence.",
            (
                "Apply one joint Holder inequality to all group-density factors and the reference bad event."
                if args.change_of_measure == "holder"
                else "Apply the pointwise density envelope to every group-density factor."
            ),
            "Optimize the finite Chernoff witness separately for every composition before summing positive terms.",
        ],
        "occupation_rows": occupation_rows,
        "composition_rows": composition_rows,
        "interval_log2_upper_diagnostic": aggregate / LOG2,
        "interval_margin_bits_diagnostic": -aggregate / LOG2,
        "interval_gap_to_40_bits": -aggregate / LOG2 - 40.0,
        "limitations": [
            "Nearest binary64 arithmetic is not outward rounded.",
            "The finite Chernoff grid is not asserted optimal.",
            "Only the displayed occupation interval is covered.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"interval_margin,{payload['interval_margin_bits_diagnostic']:.6f}")
    print(f"gap40,{payload['interval_gap_to_40_bits']:.6f}")
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
