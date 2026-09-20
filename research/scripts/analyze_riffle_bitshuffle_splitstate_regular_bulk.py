#!/usr/bin/env python3
"""Middle-density diagnostic for bit-shuffle plus SplitState epochs.

This deliberately mirrors the regular-occupation calculation used for
FieldCheckpointAccumulate.  The outer spectrum-density envelope and the
hypergeometric distribution of active outer blocks across one transposed
region are identical.  Only the epoch transfer is changed.

The epoch and state widths are parameters. Zero-state nonactivation uses the
uniform-support column of the selected activation receipt. An ensemble-average
receipt is a target profile; it is not automatically valid for one fixed map.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from import_wd_spectrum import parse_wd


LOG2 = math.log(2.0)


DEFAULT_ACTIVATION = Path(
    "constructions/riffle_ldpcsplitstate_g4_t256_s64/"
    "receipts/zero_state_activation_table.json"
)
DEFAULT_OUTPUT = Path(
    "constructions/riffle_bchperm_transpose_bitshuffle_structuredstepconv/"
    "receipts/splitstate_regular_bulk_comparison.json"
)


def log_choose(total: int, selected: int) -> float:
    if not 0 <= selected <= total:
        return -math.inf
    return (
        math.lgamma(total + 1.0)
        - math.lgamma(selected + 1.0)
        - math.lgamma(total - selected + 1.0)
    )


def log_two_power_minus_one(bits: int) -> float:
    return bits * LOG2 + math.log1p(-math.ldexp(1.0, -bits))


def logsumexp(values: np.ndarray) -> float:
    maximum = float(np.max(values))
    if maximum == -math.inf:
        return -math.inf
    return maximum + math.log(float(np.sum(np.exp(values - maximum))))


def modeled_even_floor_spectrum_logs(
    outer_bits: int, minimum_distance: int
) -> np.ndarray:
    dimension = outer_bits // 2
    weights = list(range(minimum_distance, outer_bits - minimum_distance + 1, 2))
    raw = np.asarray([log_choose(outer_bits, weight) for weight in weights])
    interior_log_mass = log_two_power_minus_one(dimension) + math.log1p(
        -math.exp(-log_two_power_minus_one(dimension))
    )
    normalizer = logsumexp(raw)
    values = np.full(outer_bits + 1, -math.inf)
    for weight, raw_log in zip(weights, raw):
        values[weight] = float(raw_log) + interior_log_mass - normalizer
    values[outer_bits] = 0.0
    return values


def load_outer_spectrum_logs(path: Path, outer_bits: int) -> np.ndarray:
    """Load an exact or ensemble-expected homogeneous outer spectrum.

    The B=1024 BA receipts store base-2 expected multiplicities.  An expected
    spectrum is valid here when the local BA schedules are sampled
    independently for each Riffle outer-block instance.
    """

    values = np.full(outer_bits + 1, -math.inf)
    if path.suffix.lower() == ".wd":
        rows = parse_wd(path.read_text(encoding="utf-8"))
        for weight_raw, count_raw in rows:
            count = int(count_raw)
            if count:
                values[int(weight_raw)] = math.log(count)
    else:
        payload = json.loads(path.read_text(encoding="utf-8"))
        if "spectrum" in payload:
            for row in payload["spectrum"]:
                weight = int(row["weight"])
                log2_multiplicity = row.get("log2_expected_multiplicity")
                if log2_multiplicity is not None:
                    values[weight] = float(log2_multiplicity) * LOG2
        elif "weight_counts" in payload:
            raw = payload["weight_counts"]
            rows = enumerate(raw) if isinstance(raw, list) else raw.items()
            for weight_raw, count_raw in rows:
                count = int(count_raw)
                if count:
                    values[int(weight_raw)] = math.log(count)
        else:
            raise ValueError("outer spectrum receipt has no supported spectrum field")
    if not math.isclose(values[0], 0.0, abs_tol=1e-12):
        raise ValueError("outer spectrum must contain the unique zero word")
    return values


def random_linear_extension_spectrum_logs(
    path: Path, outer_bits: int, extension_bits: int = 1
) -> np.ndarray:
    """Expected spectrum after adding random message-linear coordinates.

    The source has length outer_bits-extension_bits.  Independently for each
    outer-block instance, setup samples uniform linear functionals on its
    message space and appends their values.  For every fixed nonzero source
    word, the added weight is Binomial(extension_bits, 1/2).  The zero word
    remains the unique zero word.
    """
    if not 1 <= extension_bits < outer_bits:
        raise ValueError("extension_bits must lie in [1, outer_bits)")
    source = load_outer_spectrum_logs(path, outer_bits - extension_bits)
    result = np.full(outer_bits + 1, -math.inf)
    result[0] = 0.0
    for weight in range(1, outer_bits - extension_bits + 1):
        multiplicity = float(source[weight])
        if not math.isfinite(multiplicity):
            continue
        for added_weight in range(extension_bits + 1):
            contribution = (
                multiplicity
                + math.log(math.comb(extension_bits, added_weight))
                - extension_bits * LOG2
            )
            destination = weight + added_weight
            result[destination] = float(
                np.logaddexp(result[destination], contribution)
            )
    return result


def spectrum_density_envelope_log(
    outer_bits: int,
    dimension: int,
    spectrum: np.ndarray,
    tail_endpoint_width: int = 0,
) -> tuple[float, int]:
    random_log_density = log_two_power_minus_one(dimension) - outer_bits * LOG2
    best = -math.inf
    best_weight = -1
    for weight in range(1, outer_bits):
        if min(weight, outer_bits - weight) <= tail_endpoint_width:
            continue
        if not math.isfinite(float(spectrum[weight])):
            continue
        ratio = (
            float(spectrum[weight])
            - log_choose(outer_bits, weight)
            - random_log_density
        )
        if ratio > best:
            best = ratio
            best_weight = weight
    return best, best_weight


def binomial_transform(maximum: int, probability: float = 0.5) -> np.ndarray:
    if not 0.0 < probability < 1.0:
        raise ValueError("candidate probability must lie strictly between zero and one")
    transform = np.zeros((maximum + 1, maximum + 1), dtype=np.float64)
    for candidates in range(maximum + 1):
        for active in range(candidates + 1):
            transform[candidates, active] = (
                math.comb(candidates, active)
                * probability**active
                * (1.0 - probability) ** (candidates - active)
            )
    return transform


def spectrum_bernoulli_envelope_log(
    outer_bits: int,
    spectrum: np.ndarray,
    probability: float,
    tail_endpoint_width: int = 0,
) -> tuple[float, int]:
    """Return log mass of a Bernoulli envelope for the nonzero spectrum.

    For every support S of weight w, the permuted-codeword mass is
    A_w/C(B,w).  This routine returns the least C on the supplied Bernoulli
    product measure such that A_w/C(B,w) <= C p^w(1-p)^(B-w) for every
    nonzero, non-all-one weight.  The all-one word remains a separate class.
    """
    best = -math.inf
    best_weight = -1
    for weight in range(1, outer_bits):
        if min(weight, outer_bits - weight) <= tail_endpoint_width:
            continue
        if not math.isfinite(float(spectrum[weight])):
            continue
        candidate = (
            float(spectrum[weight])
            - log_choose(outer_bits, weight)
            - weight * math.log(probability)
            - (outer_bits - weight) * math.log1p(-probability)
        )
        if candidate > best:
            best = candidate
            best_weight = weight
    return best, best_weight


def deterministic_live_moment(
    z: float,
    distance: int,
    order: int,
    step_bits: int,
    state_bits: int,
) -> float:
    denominator = (1 << state_bits) - 1
    if order == 1:
        mean_weight = step_bits * (1 << (state_bits - 1)) / denominator
        low_fraction = (step_bits - mean_weight) / (step_bits - distance)
        return low_fraction * z**distance + (1.0 - low_fraction) * z**step_bits

    # SciPy is loaded only for the stronger factorial-moment experiment.
    # The mean-only parameter sweeps therefore avoid its large startup cost.
    from scipy.optimize import linprog

    weights = np.arange(distance, step_bits + 1, dtype=np.float64)
    equality_rows = [np.ones_like(weights)]
    equality_values = [1.0]
    falling = np.ones_like(weights)
    normalization = 1.0
    for degree in range(1, order + 1):
        falling *= weights - (degree - 1)
        normalization *= step_bits - (degree - 1)
        equality_rows.append(falling / normalization)
        equality_values.append((1 << (state_bits - degree)) / denominator)
    objective_scale = z**distance
    result = linprog(
        c=-(z**weights) / objective_scale,
        A_eq=np.asarray(equality_rows),
        b_eq=np.asarray(equality_values),
        bounds=(0.0, None),
        method="highs",
    )
    if not result.success:
        raise RuntimeError(f"moment LP failed: {result.message}")
    return -float(result.fun) * objective_scale


def support_averaged_kernel(
    *, z: float, step_bits: int, input_weight: int, codeword_weight: int
) -> float:
    """Return E[z^wt(X+c)] for X uniform on one weight shell."""
    denominator = math.comb(step_bits, input_weight)
    result = 0.0
    for intersection in range(
        max(0, input_weight - (step_bits - codeword_weight)),
        min(input_weight, codeword_weight) + 1,
    ):
        probability = (
            math.comb(codeword_weight, intersection)
            * math.comb(step_bits - codeword_weight, input_weight - intersection)
            / denominator
        )
        result += probability * z ** (
            input_weight + codeword_weight - 2 * intersection
        )
    return result


def mean_constrained_upper(values: np.ndarray, mean: float, minimum: int) -> float:
    """Maximize a discrete expectation with fixed mass and mean."""
    lower_maximum = min(int(math.floor(mean)), len(values) - 1)
    upper_minimum = max(int(math.ceil(mean)), minimum)
    best = -math.inf
    if mean.is_integer() and minimum <= int(mean) < len(values):
        best = float(values[int(mean)])
    for lower in range(minimum, lower_maximum + 1):
        for upper in range(upper_minimum, len(values)):
            if lower == upper:
                continue
            lower_mass = (upper - mean) / (upper - lower)
            candidate = (
                lower_mass * float(values[lower])
                + (1.0 - lower_mass) * float(values[upper])
            )
            best = max(best, candidate)
    if not math.isfinite(best):
        raise ArithmeticError("empty mean-constrained spectrum polytope")
    return best


def support_averaged_fixed_code_moments(
    *,
    z: float,
    step_bits: int,
    state_bits: int,
    distance: int,
    moment_order: int,
) -> np.ndarray:
    """Bound every uniform-support affine moment for one fixed state code.

    The bound uses only minimum distance and full output support.  Full output
    support fixes the mean weight of a uniform nonzero state codeword.
    """
    denominator = (1 << state_bits) - 1
    mean = step_bits * (1 << (state_bits - 1)) / denominator
    result = np.empty(step_bits + 1, dtype=np.float64)
    if moment_order > 1:
        from scipy.optimize import linprog

        weights = np.arange(distance, step_bits + 1, dtype=np.float64)
        equality_rows = [np.ones_like(weights)]
        equality_values = [1.0]
        falling = np.ones_like(weights)
        normalization = 1.0
        for degree in range(1, moment_order + 1):
            falling *= weights - (degree - 1)
            normalization *= step_bits - (degree - 1)
            equality_rows.append(falling / normalization)
            equality_values.append(
                (1 << (state_bits - degree)) / denominator
            )
        equality_matrix = np.asarray(equality_rows)
        equality_vector = np.asarray(equality_values)
    for input_weight in range(step_bits + 1):
        by_codeword_weight = np.zeros(step_bits + 1, dtype=np.float64)
        for codeword_weight in range(distance, step_bits + 1):
            by_codeword_weight[codeword_weight] = support_averaged_kernel(
                z=z,
                step_bits=step_bits,
                input_weight=input_weight,
                codeword_weight=codeword_weight,
            )
        if moment_order == 1:
            result[input_weight] = mean_constrained_upper(
                by_codeword_weight, mean, distance
            )
        else:
            objective = by_codeword_weight[distance:]
            scale = float(np.max(objective))
            if scale == 0.0:
                result[input_weight] = 0.0
                continue
            solution = linprog(
                c=-objective / scale,
                A_eq=equality_matrix,
                b_eq=equality_vector,
                bounds=(0.0, None),
                method="highs",
            )
            if not solution.success:
                raise RuntimeError(
                    "support-averaged moment LP failed: "
                    f"{solution.message}"
                )
            result[input_weight] = -float(solution.fun) * scale
    return result


def log_matmul_batch(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    result = np.empty_like(left)
    result[:, 0, 0] = np.logaddexp(
        left[:, 0, 0] + right[0, 0], left[:, 0, 1] + right[1, 0]
    )
    result[:, 0, 1] = np.logaddexp(
        left[:, 0, 0] + right[0, 1], left[:, 0, 1] + right[1, 1]
    )
    result[:, 1, 0] = np.logaddexp(
        left[:, 1, 0] + right[0, 0], left[:, 1, 1] + right[1, 0]
    )
    result[:, 1, 1] = np.logaddexp(
        left[:, 1, 0] + right[0, 1], left[:, 1, 1] + right[1, 1]
    )
    return result


def log_choose_vector(total: int) -> np.ndarray:
    return np.asarray([log_choose(total, count) for count in range(total + 1)])


def regular_region_log_matrices(
    candidate_logs: np.ndarray,
    epoch_bits: int,
    epochs_per_region: int,
    maximum_occupation: int,
    progress: bool = False,
) -> np.ndarray:
    current = candidate_logs.copy()
    current_maximum = epoch_bits
    epoch_choose = log_choose_vector(epoch_bits)
    for completed_epochs in range(1, epochs_per_region):
        next_maximum = min((completed_epochs + 1) * epoch_bits, maximum_occupation)
        updated = np.full((next_maximum + 1, 2, 2), -math.inf)
        current_choose = log_choose_vector(completed_epochs * epoch_bits)
        next_choose = log_choose_vector((completed_epochs + 1) * epoch_bits)
        for next_epoch_count in range(epoch_bits + 1):
            source_count = min(current_maximum, next_maximum - next_epoch_count) + 1
            if source_count <= 0:
                break
            products = log_matmul_batch(
                current[:source_count], candidate_logs[next_epoch_count]
            )
            destinations = slice(next_epoch_count, next_epoch_count + source_count)
            weights = (
                epoch_choose[next_epoch_count]
                + current_choose[:source_count]
                - next_choose[next_epoch_count : next_epoch_count + source_count]
            )
            updated[destinations] = np.logaddexp(
                updated[destinations], products + weights[:, None, None]
            )
        current = updated
        current_maximum = next_maximum
        if progress:
            print(
                f"region_epoch,{completed_epochs + 1},{epochs_per_region},"
                f"occupations,{current_maximum + 1}",
                flush=True,
            )
    return current


def log_matrix_power_moments_batch(
    matrices: np.ndarray, power: int
) -> np.ndarray:
    rows = np.full((matrices.shape[0], 2), -math.inf)
    rows[:, 0] = 0.0
    for _ in range(power):
        next_zero = np.logaddexp(
            rows[:, 0] + matrices[:, 0, 0],
            rows[:, 1] + matrices[:, 1, 0],
        )
        next_live = np.logaddexp(
            rows[:, 0] + matrices[:, 0, 1],
            rows[:, 1] + matrices[:, 1, 1],
        )
        rows[:, 0] = next_zero
        rows[:, 1] = next_live
    return np.logaddexp(rows[:, 0], rows[:, 1])


def one_forced_from_regular_logs(regular_logs: np.ndarray) -> np.ndarray:
    """Convert a fair candidate into one forced one in every region."""
    minuend = LOG2 + regular_logs[1:]
    subtrahend = regular_logs[:-1]
    both_zero = np.isneginf(subtrahend) & np.isneginf(minuend)
    with np.errstate(invalid="ignore"):
        gap = subtrahend - minuend
    gap[both_zero] = -math.inf
    if np.any(gap > 2e-11):
        raise ArithmeticError("one-forced finite difference became negative")
    gap = np.minimum(gap, 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        correction = np.log1p(-np.exp(gap))
    return minuend + correction


def load_uniform_nonactivation(path: Path, step_bits: int) -> np.ndarray:
    payload = json.loads(path.read_text(encoding="utf-8"))
    values = np.full(step_bits + 1, math.nan, dtype=np.float64)
    values[0] = 1.0
    for row in payload["by_total_weight"]:
        weight = int(row["total_weight"])
        if weight <= step_bits:
            key = "uniform_support_average_distinct_upper_bound"
            if key not in row:
                key = "uniform_256_support_average_distinct_upper_bound"
            values[weight] = float(row[key])
    missing = np.flatnonzero(np.isnan(values))
    if missing.size:
        raise ValueError(f"activation table is missing weight {int(missing[0])}")
    return values


def load_nonzero_spectrum(
    path: Path, step_bits: int, state_bits: int
) -> np.ndarray:
    payload = json.loads(path.read_text(encoding="utf-8"))
    counts = np.zeros(step_bits + 1, dtype=np.float64)
    for row in payload["spectrum"]:
        counts[int(row["weight"])] = int(row["count"])
    if int(np.sum(counts)) != 1 << state_bits or counts[0] != 1.0:
        raise ValueError("live spectrum has the wrong mass")
    counts[0] = 0.0
    return counts / ((1 << state_bits) - 1)


def support_averaged_exact_spectrum_moments(
    *, z: float, step_bits: int, spectrum: np.ndarray
) -> np.ndarray:
    result = np.zeros(step_bits + 1, dtype=np.float64)
    live_weights = np.flatnonzero(spectrum)
    for input_weight in range(step_bits + 1):
        result[input_weight] = sum(
            float(spectrum[codeword_weight])
            * support_averaged_kernel(
                z=z,
                step_bits=step_bits,
                input_weight=input_weight,
                codeword_weight=int(codeword_weight),
            )
            for codeword_weight in live_weights
        )
    return result


def splitstate_impulse_matrices(
    *,
    z: float,
    step_bits: int,
    state_bits: int,
    constituent_distance: int,
    live_moment_order: int,
    live_model: str,
    nonactivation: np.ndarray,
    live_spectrum: np.ndarray | None = None,
) -> np.ndarray:
    """Return entrywise upper epoch transfers conditioned on input weight."""
    denominator = math.ldexp(1.0, state_bits) - 1.0
    support_averaged_moments = None
    if live_model in (
        "support-averaged-fixed-code",
        "support-averaged-preaddmul",
    ):
        if live_spectrum is None:
            support_averaged_moments = support_averaged_fixed_code_moments(
                z=z,
                step_bits=step_bits,
                state_bits=state_bits,
                distance=constituent_distance,
                moment_order=live_moment_order,
            )
        else:
            support_averaged_moments = (
                support_averaged_exact_spectrum_moments(
                    z=z,
                    step_bits=step_bits,
                    spectrum=live_spectrum,
                )
            )
        live_moment = float(support_averaged_moments[0])
    elif live_model in (
        "factorial-moment",
        "factorial-live-averaged-coset",
        "combined-coset-envelope",
    ):
        live_moment = deterministic_live_moment(
            z,
            constituent_distance,
            live_moment_order,
            step_bits,
            state_bits,
        )
    elif live_model in ("ideal-random-spectrum", "ideal-random-coset"):
        # Expected spectrum of a uniform injective s-to-t binary linear map,
        # conditioned on a fixed nonzero input.  This is a design-landscape
        # model, not a certificate for the recorded degree-7 constituent.
        zero_mass = math.ldexp(1.0, -step_bits)
        live_moment = (
            math.ldexp((1.0 + z) ** step_bits, -step_bits) - zero_mass
        ) / (1.0 - zero_mass)
    else:
        raise ValueError(f"unknown live model: {live_model}")
    matrices = np.zeros((step_bits + 1, 2, 2), dtype=np.float64)
    matrices[0, 0, 0] = 1.0
    if live_model == "support-averaged-preaddmul":
        # A live state before an epoch is either uniform on F* or uniform on
        # F* with one additional value omitted.  The latter distribution is
        # bounded by this factor times the uniform-nonzero average.
        punctured_factor = denominator / (denominator - 1.0)
        matrices[0, 1, 1] = punctured_factor * live_moment
    else:
        punctured_factor = 1.0
        matrices[0, 1, 1] = live_moment
    ambient_zero_mass = math.ldexp(1.0, -step_bits)
    ambient_nonzero_mass = 1.0 - ambient_zero_mass
    ambient_termination = ambient_zero_mass / ambient_nonzero_mass
    for weight in range(1, step_bits + 1):
        output_factor = z**weight
        # These are separate entrywise upper bounds.  In particular the
        # zero-row entries need not sum to the exact row moment.
        matrices[weight, 0, 0] = nonactivation[weight] * output_factor
        matrices[weight, 0, 1] = output_factor
        if live_model == "support-averaged-preaddmul":
            # Variant update: Q' = alpha*Q + B(X).  Conditioned on Q != 0,
            # fresh uniform alpha makes alpha*Q uniform on F*.  If B(X) is
            # nonzero, termination has probability 1/(2^s-1) independently
            # of the emitted word X+A(Q).  Conditional on entering live,
            # the next state is uniform on F* with at most one value omitted.
            averaged = punctured_factor * float(
                support_averaged_moments[weight]
            )
            matrices[weight, 1, 0] = averaged / denominator
            matrices[weight, 1, 1] = averaged
        elif live_model == "support-averaged-fixed-code":
            # Conditioned on its weight, the bit-shuffled epoch input has
            # uniform support.  Averaging X+A(Q) over this support and over
            # uniform nonzero Q depends only on the ordinary spectrum of
            # im(A), not on an affine-coset maximum.
            matrices[weight, 1, 0] = z / denominator
            matrices[weight, 1, 1] = float(
                support_averaged_moments[weight]
            )
        elif live_model in (
            "ideal-random-coset",
            "factorial-live-averaged-coset",
            "combined-coset-envelope",
        ):
            # Ensemble average for C uniform over the nonzero ambient vectors:
            # X+C is uniform over every vector except X.  The one terminating
            # choice C=X is separated from the remaining live outputs.
            total_shifted_moment = (
                math.ldexp((1.0 + z) ** step_bits, -step_bits)
                - ambient_zero_mass * output_factor
            ) / ambient_nonzero_mass
            averaged_live = max(
                0.0, total_shifted_moment - ambient_termination
            )
            matrices[weight, 1, 0] = ambient_termination
            if live_model == "combined-coset-envelope":
                sparse_floor = z ** max(0, constituent_distance - weight)
                fixed_coset_bound = (
                    live_moment + (1.0 - output_factor) / denominator
                )
                matrices[weight, 1, 1] = min(
                    sparse_floor, fixed_coset_bound, averaged_live
                )
            else:
                matrices[weight, 1, 1] = averaged_live
        else:
            sparse_floor = z ** max(0, constituent_distance - weight)
            coset_bound = live_moment + (1.0 - output_factor) / denominator
            matrices[weight, 1, 0] = sparse_floor / denominator
            matrices[weight, 1, 1] = min(sparse_floor, coset_bound)
    return matrices


def candidate_epoch_logs(
    impulses: np.ndarray, transform: np.ndarray
) -> np.ndarray:
    """Average over fair values of the candidate bits in one epoch."""
    candidates = (transform @ impulses.reshape(impulses.shape[0], 4)).reshape(
        impulses.shape
    )
    # Matrix multiplication can create tiny negative roundoff from structural
    # zeros even though every mathematical coefficient is nonnegative.
    minimum = float(np.min(candidates))
    if minimum < -2e-15:
        raise ArithmeticError(f"candidate transfer became negative: {minimum}")
    candidates[candidates < 0.0] = 0.0
    with np.errstate(divide="ignore"):
        return np.log(candidates)


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    dimension = (
        args.outer_dimension
        if args.outer_dimension is not None
        else args.outer_bits // 2
    )
    if args.message_bits % dimension:
        raise ValueError("outer dimension must divide message bits")
    outer_blocks = args.message_bits // dimension
    output_bits = outer_blocks * args.outer_bits
    distance = math.floor(args.relative_distance * output_bits)
    region_bits = output_bits // args.outer_bits
    if region_bits % args.step_bits:
        raise ValueError("step bits must divide one transposed region")
    epochs_per_region = region_bits // args.step_bits
    maximum_occupation = max(args.active_blocks)
    if args.include_one_all_one:
        maximum_occupation += 1
    if maximum_occupation > region_bits:
        raise ValueError("active-block occupation exceeds one region")

    spectrum = (
        (
            random_linear_extension_spectrum_logs(
                args.outer_spectrum,
                args.outer_bits,
                args.random_linear_extension_bits,
            )
            if args.random_linear_extension
            else load_outer_spectrum_logs(args.outer_spectrum, args.outer_bits)
        )
        if args.outer_spectrum is not None
        else modeled_even_floor_spectrum_logs(
            args.outer_bits, args.modeled_minimum_distance
        )
    )
    log_eta, envelope_weight = spectrum_density_envelope_log(
        args.outer_bits,
        dimension,
        spectrum,
        args.tail_endpoint_width,
    )
    probability_envelopes = {
        probability: spectrum_bernoulli_envelope_log(
            args.outer_bits,
            spectrum,
            probability,
            args.tail_endpoint_width,
        )
        for probability in args.candidate_probabilities
    }
    if args.ideal_random_activation:
        nonactivation = np.ones(args.step_bits + 1, dtype=np.float64)
        nonactivation[1:] = math.ldexp(1.0, -args.state_bits)
    else:
        nonactivation = load_uniform_nonactivation(
            args.activation, args.step_bits
        )
    live_spectrum = (
        load_nonzero_spectrum(
            args.live_spectrum, args.step_bits, args.state_bits
        )
        if args.live_spectrum is not None
        else None
    )

    best = {occupation: math.inf for occupation in args.active_blocks}
    best_tilt = {occupation: math.nan for occupation in args.active_blocks}
    best_probability = {
        occupation: math.nan for occupation in args.active_blocks
    }
    one_all_one_best = {
        occupation: math.inf
        for occupation in args.active_blocks
        if occupation < outer_blocks
    }
    one_all_one_best_tilt = {
        occupation: math.nan for occupation in one_all_one_best
    }
    tilt_rows = []
    for tilt_index, log_surprisal in enumerate(args.log_surprisals):
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        impulses = splitstate_impulse_matrices(
            z=z,
            step_bits=args.step_bits,
            state_bits=args.state_bits,
            constituent_distance=args.constituent_distance,
            live_moment_order=args.live_moment_order,
            live_model=args.live_model,
            nonactivation=nonactivation,
            live_spectrum=live_spectrum,
        )
        if args.occupation_family != "regular":
            # Every selected all-one outer block contributes one forced one
            # to each transposed region.  The region permutation makes the
            # resulting support uniform without randomizing its values.
            with np.errstate(divide="ignore"):
                candidate_logs = np.log(impulses)
        probabilities = (
            args.candidate_probabilities
            if args.occupation_family == "regular"
            else (0.5,)
        )
        for probability in probabilities:
            if args.occupation_family == "regular":
                transform = binomial_transform(args.step_bits, probability)
                candidate_logs = candidate_epoch_logs(impulses, transform)
            regions = regular_region_log_matrices(
                candidate_logs,
                args.step_bits,
                epochs_per_region,
                maximum_occupation,
                progress=args.progress,
            )
            moments = log_matrix_power_moments_batch(regions, args.outer_bits)
            if args.include_one_all_one:
                if args.occupation_family != "regular":
                    raise ValueError(
                        "--include-one-all-one requires --occupation-family regular"
                    )
                forced_regions = one_forced_from_regular_logs(regions)
                forced_moments = log_matrix_power_moments_batch(
                    forced_regions, args.outer_bits
                )
            envelope_log = (
                probability_envelopes[probability][0]
                if args.occupation_family == "regular"
                else 0.0
            )
            for occupation in args.active_blocks:
                inner_candidate = (
                    float(moments[occupation]) + distance * surprisal
                )
                candidate = (
                    occupation * envelope_log + min(0.0, inner_candidate)
                )
                if not args.omit_tilt_rows:
                    tilt_rows.append(
                        {
                            "active_regular_outer_blocks": occupation,
                            "log_surprisal": log_surprisal,
                            "candidate_probability": probability,
                            "inner_plus_spectrum_log2_upper": candidate / LOG2,
                        }
                    )
                if candidate < best[occupation]:
                    best[occupation] = candidate
                    best_tilt[occupation] = log_surprisal
                    best_probability[occupation] = probability
                if args.include_one_all_one and occupation < outer_blocks:
                    forced_candidate = (
                        occupation * envelope_log
                        + min(
                            0.0,
                            float(forced_moments[occupation])
                            + distance * surprisal,
                        )
                    )
                    if forced_candidate < one_all_one_best[occupation]:
                        one_all_one_best[occupation] = forced_candidate
                        one_all_one_best_tilt[occupation] = log_surprisal
        print(
            f"tilt,{tilt_index + 1},{len(args.log_surprisals)},"
            f"log_surprisal,{log_surprisal:.6f}",
            flush=True,
        )

    rows = []
    for occupation in args.active_blocks:
        outer_log = log_choose(outer_blocks, occupation)
        inner_log = best[occupation]
        rows.append(
            {
                "active_regular_outer_blocks": occupation,
                "best_log_surprisal": best_tilt[occupation],
                "best_candidate_probability": best_probability[occupation],
                "outer_log2_envelope": outer_log / LOG2,
                "inner_log2_upper": inner_log / LOG2,
                "pointwise_log2_upper": (outer_log + inner_log) / LOG2,
                "pointwise_margin_bits": -(outer_log + inner_log) / LOG2,
            }
        )

    pointwise_logs = np.asarray(
        [-row["pointwise_margin_bits"] * LOG2 for row in rows],
        dtype=np.float64,
    )
    partial_log = float(logsumexp(pointwise_logs))
    dominant_rows = sorted(
        rows,
        key=lambda row: float(row["pointwise_log2_upper"]),
        reverse=True,
    )[:20]

    one_all_one_rows = []
    one_all_one_partial_log = -math.inf
    if args.include_one_all_one:
        forced_logs = []
        for occupation in args.active_blocks:
            if occupation >= outer_blocks:
                continue
            outer_log = (
                log_choose(outer_blocks, occupation)
                + math.log(outer_blocks - occupation)
            )
            inner_log = one_all_one_best[occupation]
            contribution = outer_log + inner_log
            forced_logs.append(contribution)
            one_all_one_rows.append(
                {
                    "regular_active_blocks": occupation,
                    "all_one_active_blocks": 1,
                    "best_log_surprisal": one_all_one_best_tilt[occupation],
                    "outer_log2_envelope": outer_log / LOG2,
                    "inner_log2_upper": inner_log / LOG2,
                    "pointwise_log2_upper": contribution / LOG2,
                    "pointwise_margin_bits": -contribution / LOG2,
                }
            )
        one_all_one_partial_log = float(logsumexp(np.asarray(forced_logs)))

    return {
        "schema": "riffle-bitshuffle-splitstate-regular-bulk-v1",
        "candidate": (
            "Riffle BCHPerm-TransposeBitShuffle-SplitState "
            f"t={args.step_bits} s={args.state_bits}"
        ),
        "method": {
            "outer": "regular spectrum-density envelope from the supplied exact/expected spectrum or the modeled fallback",
            "region_recurrence": "exact hypergeometric conditioning of candidate positions",
            "epoch_transfer": "entrywise SplitState upper matrix",
            "arithmetic": "nearest binary64 log semiring",
            "outward_rounded": False,
        },
        "parameters": {
            "message_bits": args.message_bits,
            "output_bits": output_bits,
            "relative_distance": args.relative_distance,
            "distance": distance,
            "outer_bits": args.outer_bits,
            "outer_dimension": dimension,
            "outer_blocks": outer_blocks,
            "modeled_minimum_distance": args.modeled_minimum_distance,
            "outer_spectrum": (
                str(args.outer_spectrum) if args.outer_spectrum is not None else None
            ),
            "random_linear_extension": args.random_linear_extension,
            "random_linear_extension_bits": (
                args.random_linear_extension_bits
                if args.random_linear_extension
                else 0
            ),
            "tail_endpoint_width": args.tail_endpoint_width,
            "step_bits": args.step_bits,
            "state_bits": args.state_bits,
            "constituent_distance": args.constituent_distance,
            "live_moment_order": args.live_moment_order,
            "live_model": args.live_model,
            "region_bits": region_bits,
            "epochs_per_region": epochs_per_region,
            "active_blocks": list(args.active_blocks),
            "occupation_family": args.occupation_family,
            "activation_model": (
                "ideal-random-linear-map"
                if args.ideal_random_activation
                else str(args.activation)
            ),
            "live_spectrum": (
                str(args.live_spectrum)
                if args.live_spectrum is not None
                else None
            ),
            "log_surprisals": list(args.log_surprisals),
            "candidate_probabilities": list(args.candidate_probabilities),
        },
        "envelope": {
            "log2_eta": log_eta / LOG2,
            "maximizing_regular_weight": envelope_weight,
            "bernoulli_envelopes": [
                {
                    "probability": probability,
                    "log2_mass": envelope[0] / LOG2,
                    "maximizing_weight": envelope[1],
                }
                for probability, envelope in probability_envelopes.items()
            ],
        },
        "occupation_rows": rows,
        "dominant_rows": dominant_rows,
        "partial_log2_upper": partial_log / LOG2,
        "partial_margin_bits": -partial_log / LOG2,
        "one_all_one_rows": one_all_one_rows,
        "one_all_one_partial_log2_upper": one_all_one_partial_log / LOG2,
        "one_all_one_partial_margin_bits": -one_all_one_partial_log / LOG2,
        "tilt_rows": tilt_rows,
        "scope": (
            "Floating-point diagnostic. The activation table is consumed as "
            "a fixed-map profile; an ensemble-average table requires a separate "
            "fixed-map existence or concentration argument. The factorial-moment live model is an upper "
            "bound for the audited constituent. The ideal-random spectrum "
            "and coset models are only constituent-code landscape assumptions. "
            "Numerical arithmetic is not outward rounded."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=1 << 20)
    parser.add_argument("--outer-bits", type=int, default=256)
    parser.add_argument(
        "--outer-dimension",
        type=int,
        help="outer message dimension (default: outer-bits // 2)",
    )
    parser.add_argument("--modeled-minimum-distance", type=int, default=38)
    parser.add_argument("--outer-spectrum", type=Path)
    parser.add_argument(
        "--random-linear-extension",
        action="store_true",
        help=(
            "treat outer-spectrum as a length-(outer-bits-1) source and "
            "append an independently sampled message-linear bit per block"
        ),
    )
    parser.add_argument(
        "--random-linear-extension-bits",
        type=int,
        default=1,
        help="number of independently sampled appended linear bits",
    )
    parser.add_argument(
        "--tail-endpoint-width",
        type=int,
        default=0,
        help="exclude these endpoint shells from the regular envelope",
    )
    parser.add_argument("--relative-distance", type=float, default=0.09)
    parser.add_argument("--step-bits", type=int, default=256)
    parser.add_argument("--state-bits", type=int, default=64)
    parser.add_argument("--constituent-distance", type=int, default=40)
    parser.add_argument("--live-moment-order", type=int, default=3)
    parser.add_argument(
        "--occupation-family",
        choices=("regular", "all-one-only"),
        default="regular",
    )
    parser.add_argument(
        "--live-model",
        choices=(
            "factorial-moment",
            "ideal-random-spectrum",
            "ideal-random-coset",
            "factorial-live-averaged-coset",
            "combined-coset-envelope",
            "support-averaged-fixed-code",
            "support-averaged-preaddmul",
        ),
        default="factorial-moment",
    )
    parser.add_argument("--active-blocks", type=int, nargs="+", default=(512, 768, 1024))
    parser.add_argument(
        "--active-range",
        type=int,
        nargs=2,
        metavar=("MINIMUM", "MAXIMUM"),
        help="replace --active-blocks with an inclusive occupation range",
    )
    parser.add_argument(
        "--log-surprisals",
        type=float,
        nargs="+",
        default=(-2.5, -2.25, -2.0, -1.75, -1.5, -1.25),
    )
    parser.add_argument(
        "--candidate-probabilities",
        type=float,
        nargs="+",
        default=(0.5,),
        help=(
            "Bernoulli reference probabilities for the spectrum-density "
            "envelope; every value is optimized independently"
        ),
    )
    parser.add_argument("--activation", type=Path, default=DEFAULT_ACTIVATION)
    parser.add_argument("--live-spectrum", type=Path)
    parser.add_argument("--ideal-random-activation", action="store_true")
    parser.add_argument("--progress", action="store_true")
    parser.add_argument("--omit-tilt-rows", action="store_true")
    parser.add_argument("--include-one-all-one", action="store_true")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.random_linear_extension and args.outer_spectrum is None:
        parser.error("--random-linear-extension requires --outer-spectrum")
    if args.random_linear_extension_bits < 1:
        parser.error("--random-linear-extension-bits must be positive")
    if not 0 <= args.tail_endpoint_width < args.outer_bits // 2:
        parser.error("--tail-endpoint-width must lie in [0, outer-bits/2)")
    if any(not 0.0 < probability < 1.0 for probability in args.candidate_probabilities):
        parser.error("--candidate-probabilities must lie strictly between zero and one")
    if args.active_range is not None:
        minimum, maximum = args.active_range
        if not 1 <= minimum <= maximum:
            parser.error("--active-range must satisfy 1 <= MINIMUM <= MAXIMUM")
        args.active_blocks = tuple(range(minimum, maximum + 1))
    payload = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    for row in payload["occupation_rows"]:
        if len(payload["occupation_rows"]) > 32 and row not in payload["dominant_rows"]:
            continue
        print(
            f"occupation,{row['active_regular_outer_blocks']},"
            f"tilt,{row['best_log_surprisal']:.6f},"
            f"pointwise_log2,{row['pointwise_log2_upper']:.6f},"
            f"margin,{row['pointwise_margin_bits']:.6f}",
            flush=True,
        )
    print(f"partial_margin,{payload['partial_margin_bits']:.6f}", flush=True)
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
