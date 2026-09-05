#!/usr/bin/env python3
"""Probe mixed RM spectrum bands with the sparse RM2Sub recurrence.

For each selected pair of outer-weight bands, the probe fixes the number of
rows from each band.  Each band is dominated by its own Bernoulli product
measure.  The region permutation makes the candidate positions exchangeable.
After the two Bernoulli thinnings, the live support is therefore uniform
conditional on its total size.  A convolution of two binomial laws can then
be applied to the established forced-support region recurrence.

All reductions are positive.  Nearest-binary64 arithmetic and a finite tilt
grid make the resulting margins diagnostic rather than outward certificates.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
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
    log_choose,
    log_matrix_power_moments_batch,
    regular_region_log_matrices,
    splitstate_impulse_matrices,
)
from small_k_replay.evaluate_exact_spectra_q1_phase import (  # noqa: E402
    CONSTITUENTS,
    load_spectrum,
)
from small_k_replay.evaluate_rm2sub_primary_tranche import configuration  # noqa: E402
from small_k_replay.probe_rm2sub_fixed_rm_dense_bands import BANDS  # noqa: E402


DEFAULT_OUTPUT = HERE / "rm2sub_two_band_bridge_probe_d100.json"


def holder_grid() -> np.ndarray:
    return 1.0 + np.exp2(np.linspace(-20.0, 20.0, 1281))


def band_density_statistics(
    spectrum: dict[int, int],
    block_bits: int,
    lower: int,
    upper: int,
    probability: float,
    orders: np.ndarray,
) -> tuple[np.ndarray, float, int]:
    """Return Renyi moments of one band's counting density."""
    reference_logs = []
    density_logs = []
    weights = []
    for weight, count in spectrum.items():
        if not count or not lower <= weight <= upper:
            continue
        shell = log_choose(block_bits, weight)
        if probability == 0.0:
            point_log = 0.0 if weight == 0 else -math.inf
        elif probability == 1.0:
            point_log = 0.0 if weight == block_bits else -math.inf
        else:
            point_log = (
                weight * math.log(probability)
                + (block_bits - weight) * math.log1p(-probability)
            )
        if not math.isfinite(point_log):
            raise ValueError(
                "endpoint reference probability does not support its band"
            )
        reference_logs.append(shell + point_log)
        density_logs.append(math.log(count) - shell - point_log)
        weights.append(weight)
    reference = np.asarray(reference_logs, dtype=np.float64)
    density = np.asarray(density_logs, dtype=np.float64)
    moments = np.asarray(
        [float(logsumexp(reference + order * density)) for order in orders]
    )
    maximum_index = int(np.argmax(density))
    return moments, float(density[maximum_index]), weights[maximum_index]


def parse_band_range(value: str) -> tuple[int, int]:
    try:
        lower_raw, upper_raw = value.split(":", maxsplit=1)
        lower = int(lower_raw)
        upper = int(upper_raw)
    except (ValueError, TypeError) as error:
        raise argparse.ArgumentTypeError(
            "band ranges must have the form lower:upper"
        ) from error
    if not 1 <= lower <= upper <= 512:
        raise argparse.ArgumentTypeError("band range must lie in [1,512]")
    return lower, upper


def band_envelope_at_probability(
    spectrum: dict[int, int],
    block_bits: int,
    lower: int,
    upper: int,
    probability: float,
) -> tuple[float, int]:
    """Return the logarithmic Bernoulli majorant for one spectrum band."""
    def point_log(weight: int) -> float:
        if probability == 0.0:
            return 0.0 if weight == 0 else -math.inf
        if probability == 1.0:
            return 0.0 if weight == block_bits else -math.inf
        return (
            weight * math.log(probability)
            + (block_bits - weight) * math.log1p(-probability)
        )

    rows = [
        (
            math.log(count)
            - log_choose(block_bits, weight)
            - point_log(weight),
            weight,
        )
        for weight, count in spectrum.items()
        if count and lower <= weight <= upper
    ]
    if not rows:
        raise ValueError(f"empty spectrum band {lower}:{upper}")
    return max(rows)


def minimize_band_envelope(
    spectrum: dict[int, int], block_bits: int, lower: int, upper: int
) -> tuple[float, float, int]:
    """Choose a Bernoulli probability that minimizes the band majorant."""
    epsilon = 1e-10
    result = minimize_scalar(
        lambda probability: band_envelope_at_probability(
            spectrum, block_bits, lower, upper, float(probability)
        )[0],
        bounds=(epsilon, 1.0 - epsilon),
        method="bounded",
        options={"xatol": 1e-13},
    )
    probability = float(result.x)
    envelope, weight = band_envelope_at_probability(
        spectrum, block_bits, lower, upper, probability
    )
    return probability, envelope, weight


def two_binomial_region_log_matrices(
    forced_region_logs: np.ndarray,
    first_probability: float,
    second_probability: float,
    maximum_occupation: int,
) -> np.ndarray:
    """Average forced-support transfers after two Bernoulli thinnings.

    Entry ``[a,b]`` is the region transfer for ``a`` candidate positions of
    the first type and ``b`` candidate positions of the second type.  The two
    position sets are disjoint and uniformly permuted inside the region.
    """
    def transform(probability: float) -> np.ndarray:
        if probability == 0.0:
            result = np.zeros(
                (maximum_occupation + 1, maximum_occupation + 1),
                dtype=np.float64,
            )
            result[:, 0] = 1.0
            return result
        if probability == 1.0:
            return np.eye(maximum_occupation + 1, dtype=np.float64)
        return binomial_transform(maximum_occupation, probability)

    first = transform(first_probability)
    second = transform(second_probability)
    result = np.full(
        (maximum_occupation + 1, maximum_occupation + 1, 2, 2),
        -math.inf,
        dtype=np.float64,
    )
    for first_count in range(maximum_occupation + 1):
        for second_count in range(maximum_occupation - first_count + 1):
            distribution = np.convolve(
                first[first_count, : first_count + 1],
                second[second_count, : second_count + 1],
            )
            with np.errstate(divide="ignore"):
                log_distribution = np.log(distribution)
            result[first_count, second_count] = logsumexp(
                forced_region_logs[: distribution.size]
                + log_distribution[:, None, None],
                axis=0,
            )
    return result


def log_multinomial_rows(total: int, first: int, second: int) -> float:
    if min(first, second, total - first - second) < 0:
        return -math.inf
    return (
        math.lgamma(total + 1)
        - math.lgamma(first + 1)
        - math.lgamma(second + 1)
        - math.lgamma(total - first - second + 1)
    )


def verify_exact_nonactivation(path: Path, values: np.ndarray) -> None:
    """Check that the loaded table equals the exact kernel shell ratios."""
    payload = json.loads(path.read_text(encoding="utf-8"))
    seen = set()
    for row in payload["by_total_weight"]:
        weight = int(row["total_weight"])
        if weight >= values.size:
            continue
        exact = int(row["kernel_words"]) / int(row["shell_size"])
        if not math.isclose(
            float(values[weight]), exact, rel_tol=2e-15, abs_tol=2e-15
        ):
            raise ValueError(
                f"nonactivation entry {weight} is not its exact kernel ratio"
            )
        seen.add(weight)
    if seen != set(range(values.size)):
        raise ValueError("exact nonactivation table is incomplete")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step-bits", type=int, default=64)
    parser.add_argument("--state-bits", type=int, default=14)
    parser.add_argument("--message-exponent", type=int, default=16)
    parser.add_argument("--minimum-occupation", type=int, default=30)
    parser.add_argument("--maximum-occupation", type=int, default=52)
    parser.add_argument("--anchor-band-index", type=int, default=0)
    parser.add_argument(
        "--band-ranges",
        type=parse_band_range,
        nargs="+",
        help="replace the default bands with lower:upper ranges",
    )
    parser.add_argument(
        "--partner-band-indexes",
        type=int,
        nargs="+",
    )
    parser.add_argument("--distance-numerator", type=int, default=100)
    parser.add_argument("--distance-denominator", type=int, default=1000)
    parser.add_argument(
        "--live-model",
        choices=["factorial-moment", "support-averaged-preaddmul"],
        default="factorial-moment",
    )
    parser.add_argument(
        "--change-of-measure",
        choices=["envelope", "holder"],
        default="envelope",
    )
    parser.add_argument(
        "--common-reference-probability",
        type=float,
        help="override every band's envelope-minimizing Bernoulli probability",
    )
    parser.add_argument(
        "--reference-probabilities",
        type=float,
        nargs="+",
        help="one Bernoulli reference probability for every displayed band",
    )
    parser.add_argument(
        "--log-surprisals",
        type=float,
        nargs="+",
        default=[value / 4.0 for value in range(-32, 1)],
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    bands = tuple(args.band_ranges) if args.band_ranges is not None else BANDS
    if args.partner_band_indexes is None:
        args.partner_band_indexes = list(range(1, len(bands)))

    if args.step_bits not in (64, 128, 256):
        parser.error("supported epoch lengths are 64, 128, and 256")
    if (
        args.common_reference_probability is not None
        and not 0.0 < args.common_reference_probability < 1.0
    ):
        parser.error("common reference probability must lie in (0,1)")
    if (
        args.common_reference_probability is not None
        and args.reference_probabilities is not None
    ):
        parser.error("use either a common or per-band reference probability")
    if args.reference_probabilities is not None and (
        len(args.reference_probabilities) != len(bands)
        or any(not 0.0 <= value <= 1.0 for value in args.reference_probabilities)
    ):
        parser.error("reference-probabilities must give one value per band in [0,1]")
    if not 1 <= args.minimum_occupation <= args.maximum_occupation:
        parser.error("the occupation interval must be nonempty")
    if not 0 <= args.anchor_band_index < len(bands):
        parser.error("anchor band index is out of range")
    if any(
        not 0 <= index < len(bands) or index == args.anchor_band_index
        for index in args.partner_band_indexes
    ):
        parser.error("partner band indexes must be distinct from the anchor")

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
    selected_indexes = sorted(
        set([args.anchor_band_index, *args.partner_band_indexes])
    )
    band_data: dict[int, dict[str, float | int]] = {}
    holder_orders = holder_grid()
    holder_data: dict[int, dict[str, np.ndarray | float | int]] = {}
    for index in selected_indexes:
        lower, upper = bands[index]
        if args.reference_probabilities is not None:
            probability = args.reference_probabilities[index]
            envelope, maximizing_weight = band_envelope_at_probability(
                spectrum,
                constituent.block_bits,
                lower,
                upper,
                probability,
            )
        elif args.common_reference_probability is None:
            probability, envelope, maximizing_weight = minimize_band_envelope(
                spectrum, constituent.block_bits, lower, upper
            )
        else:
            probability = args.common_reference_probability
            envelope, maximizing_weight = band_envelope_at_probability(
                spectrum,
                constituent.block_bits,
                lower,
                upper,
                probability,
            )
        band_data[index] = {
            "index": index,
            "lower": lower,
            "upper": upper,
            "probability": probability,
            "envelope_log2": envelope / LOG2,
            "maximizing_weight": maximizing_weight,
            "exact_mass": str(
                sum(
                    count
                    for weight, count in spectrum.items()
                    if lower <= weight <= upper
                )
            ),
        }
        moments, maximum_density, maximum_weight = band_density_statistics(
            spectrum,
            constituent.block_bits,
            lower,
            upper,
            probability,
            holder_orders,
        )
        holder_data[index] = {
            "log_moments": moments,
            "maximum_density_log": maximum_density,
            "maximum_density_weight": maximum_weight,
        }

    nonactivation_path = Path(config["b"])
    nonactivation = load_uniform_nonactivation(nonactivation_path, args.step_bits)
    verify_exact_nonactivation(nonactivation_path, nonactivation)
    live_spectrum = load_nonzero_spectrum(
        Path(config["a"]), args.step_bits, args.state_bits
    )
    epochs_per_region = outer_rows // args.step_bits
    maximum_occupation = args.maximum_occupation

    pairs: dict[int, dict[str, object]] = {}
    for partner in args.partner_band_indexes:
        pairs[partner] = {
            "best": np.full(
                (maximum_occupation + 1, maximum_occupation + 1),
                math.inf,
                dtype=np.float64,
            ),
            "witness": np.full(
                (maximum_occupation + 1, maximum_occupation + 1),
                math.nan,
                dtype=np.float64,
            ),
            "holder_witness": np.full(
                (maximum_occupation + 1, maximum_occupation + 1),
                None,
                dtype=object,
            ),
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
            live_model=args.live_model,
            nonactivation=nonactivation,
            live_spectrum=live_spectrum,
        )
        # When the incoming state is zero, the emitted word is the epoch
        # input itself.  Its moment is z^w, and B(X)=0 has the exact averaged
        # probability nonactivation[w].  Split the zero and live destinations
        # instead of charging the complete moment to both entries.
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
            maximum_occupation,
        )[: maximum_occupation + 1]

        anchor = band_data[args.anchor_band_index]
        for partner_index in args.partner_band_indexes:
            partner = band_data[partner_index]
            reference_regions = two_binomial_region_log_matrices(
                forced_region_logs,
                float(anchor["probability"]),
                float(partner["probability"]),
                maximum_occupation,
            )
            flat = reference_regions.reshape((-1, 2, 2))
            moments = log_matrix_power_moments_batch(
                flat, constituent.block_bits
            ).reshape(reference_regions.shape[:2])
            best = pairs[partner_index]["best"]
            witness = pairs[partner_index]["witness"]
            holder_witness = pairs[partner_index]["holder_witness"]
            assert isinstance(best, np.ndarray)
            assert isinstance(witness, np.ndarray)
            assert isinstance(holder_witness, np.ndarray)
            for occupation in range(
                args.minimum_occupation, maximum_occupation + 1
            ):
                # Include pure faces as internal overlap checks.  The reported
                # mixed union below uses only 1 <= first_count < occupation.
                for first_count in range(occupation + 1):
                    second_count = occupation - first_count
                    inner = min(
                        0.0,
                        float(moments[first_count, second_count])
                        + bad_weight * surprisal,
                    )
                    location = log_multinomial_rows(
                        outer_rows, first_count, second_count
                    )
                    if args.change_of_measure == "envelope":
                        correction = (
                            first_count
                            * float(anchor["envelope_log2"])
                            * LOG2
                            + second_count
                            * float(partner["envelope_log2"])
                            * LOG2
                            + inner
                        )
                        holder_order: float | str = "infinity"
                    else:
                        first_holder = holder_data[args.anchor_band_index]
                        second_holder = holder_data[partner_index]
                        finite = (
                            first_count
                            * np.asarray(first_holder["log_moments"])
                            / holder_orders
                            + second_count
                            * np.asarray(second_holder["log_moments"])
                            / holder_orders
                            + (1.0 - 1.0 / holder_orders) * inner
                        )
                        holder_index = int(np.argmin(finite))
                        correction = float(finite[holder_index])
                        holder_order = float(holder_orders[holder_index])
                        infinity = (
                            first_count
                            * float(first_holder["maximum_density_log"])
                            + second_count
                            * float(second_holder["maximum_density_log"])
                            + inner
                        )
                        if infinity < correction:
                            correction = infinity
                            holder_order = "infinity"
                    candidate = location + correction
                    if candidate < best[first_count, second_count]:
                        best[first_count, second_count] = candidate
                        witness[first_count, second_count] = log_surprisal
                        holder_witness[first_count, second_count] = holder_order

    pair_rows = []
    global_mixed_values = []
    for partner_index in args.partner_band_indexes:
        best = pairs[partner_index]["best"]
        witness = pairs[partner_index]["witness"]
        holder_witness = pairs[partner_index]["holder_witness"]
        assert isinstance(best, np.ndarray)
        assert isinstance(witness, np.ndarray)
        assert isinstance(holder_witness, np.ndarray)
        occupation_rows = []
        for occupation in range(
            args.minimum_occupation, maximum_occupation + 1
        ):
            split_rows = []
            mixed_values = []
            for first_count in range(occupation + 1):
                second_count = occupation - first_count
                value = float(best[first_count, second_count])
                split_rows.append(
                    {
                        "anchor_rows": first_count,
                        "partner_rows": second_count,
                        "is_mixed": 0 < first_count < occupation,
                        "log2_upper_diagnostic": value / LOG2,
                        "margin_bits_diagnostic": -value / LOG2,
                        "gap_to_40_bits": -value / LOG2 - 40.0,
                        "log_surprisal": float(
                            witness[first_count, second_count]
                        ),
                        "holder_order": holder_witness[
                            first_count, second_count
                        ],
                    }
                )
                if 0 < first_count < occupation:
                    mixed_values.append(value)
            aggregate = float(logsumexp(np.asarray(mixed_values)))
            global_mixed_values.append(aggregate)
            dominant = max(
                (row for row in split_rows if row["is_mixed"]),
                key=lambda row: float(row["log2_upper_diagnostic"]),
            )
            occupation_rows.append(
                {
                    "occupation": occupation,
                    "mixed_log2_upper_diagnostic": aggregate / LOG2,
                    "mixed_margin_bits_diagnostic": -aggregate / LOG2,
                    "mixed_gap_to_40_bits": -aggregate / LOG2 - 40.0,
                    "dominant_mixed_split": {
                        "anchor_rows": dominant["anchor_rows"],
                        "partner_rows": dominant["partner_rows"],
                        "margin_bits_diagnostic": dominant[
                            "margin_bits_diagnostic"
                        ],
                    },
                    "split_rows": split_rows,
                }
            )
        pair_aggregate = float(
            logsumexp(
                np.asarray(
                    [
                        row["mixed_log2_upper_diagnostic"] * LOG2
                        for row in occupation_rows
                    ]
                )
            )
        )
        pair_rows.append(
            {
                "partner_band_index": partner_index,
                "partner_band": band_data[partner_index],
                "mixed_interval_log2_upper_diagnostic": pair_aggregate / LOG2,
                "mixed_interval_margin_bits_diagnostic": -pair_aggregate / LOG2,
                "mixed_interval_gap_to_40_bits": -pair_aggregate / LOG2 - 40.0,
                "occupation_rows": occupation_rows,
            }
        )

    global_aggregate = float(logsumexp(np.asarray(global_mixed_values)))
    payload = {
        "schema": "rm2sub-fixed-rm-two-band-sparse-bridge-v1",
        "status": "BINARY64_PARTIAL_DIAGNOSTIC",
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
            "anchor_band_index": args.anchor_band_index,
            "partner_band_indexes": args.partner_band_indexes,
            "band_ranges": [list(band) for band in bands],
            "live_model": args.live_model,
            "exact_zero_state_destination_split": True,
            "exact_nonactivation_verified_from_kernel_shell_counts": True,
            "change_of_measure": args.change_of_measure,
            "common_reference_probability": args.common_reference_probability,
            "reference_probabilities": args.reference_probabilities,
            "log_surprisals": args.log_surprisals,
        },
        "anchor_band": band_data[args.anchor_band_index],
        "selected_bands": [band_data[index] for index in selected_indexes],
        "proof_reduction": [
            "Partition the exact outer spectrum into the displayed bands.",
            "Dominate each band by its displayed Bernoulli product measure and multiplicative envelope.",
            "For a fixed two-band composition, uniformly permute the two disjoint sets of candidate positions in every region.",
            "Thin each position with its band's Bernoulli probability.",
            "Conditional on the total live weight, permutation symmetry makes the live support uniform.",
            "Convolve the two binomial laws and apply the positive forced-support RM2Sub recurrence.",
            *(
                ["Apply one joint Holder inequality to the independent band-density factors and the reference bad event."]
                if args.change_of_measure == "holder"
                else ["Apply the pointwise Bernoulli density envelopes."]
            ),
            "Optimize the finite Chernoff witness separately for every composition before summing positive terms.",
        ],
        "pair_rows": pair_rows,
        "union_of_displayed_mixed_pairs": {
            "log2_upper_diagnostic": global_aggregate / LOG2,
            "margin_bits_diagnostic": -global_aggregate / LOG2,
            "gap_to_40_bits": -global_aggregate / LOG2 - 40.0,
        },
        "limitations": [
            "Only compositions supported on the anchor band and one displayed partner band are included.",
            "Compositions using three or more bands remain open.",
            "Pure-band faces are recorded as overlap checks but excluded from the displayed mixed unions.",
            "Each band uses its envelope-minimizing Bernoulli probability; joint probability optimization may improve the margins.",
            "Nearest binary64 arithmetic is not outward rounded.",
            "The finite Chernoff grid is not asserted optimal.",
            *(
                ["The finite Holder-order grid is not asserted optimal."]
                if args.change_of_measure == "holder"
                else []
            ),
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    for row in pair_rows:
        partner = row["partner_band"]
        print(
            f"pair,{args.anchor_band_index},{row['partner_band_index']},"
            f"weights,{band_data[args.anchor_band_index]['lower']}-"
            f"{band_data[args.anchor_band_index]['upper']},"
            f"{partner['lower']}-{partner['upper']},"
            f"margin,{row['mixed_interval_margin_bits_diagnostic']:.6f},"
            f"gap40,{row['mixed_interval_gap_to_40_bits']:.6f}",
            flush=True,
        )
    print(
        "displayed_pair_union_margin,"
        f"{payload['union_of_displayed_mixed_pairs']['margin_bits_diagnostic']:.6f}",
        flush=True,
    )
    print(f"output={args.output}", flush=True)


if __name__ == "__main__":
    main()
