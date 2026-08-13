#!/usr/bin/env python3
"""65-state matrix probes for the identity Riffle inner.

This removes the Singer-style episode decomposition.  After exposing the
first nonzero input block, a nonlinear positive operator tracks the exact
state weight through every remaining block.  It includes zero-state turnoff,
input/state cancellation, restart, and terminal-state discard automatically.

For occupancy pole x, output pole z, first-block weight r, remaining input
weight H, and suffix length L, Cauchy's bound is applied as

    C(64,r) / C(N,h) * x^-H * a_r^T F_x,z^L(1).

The BCH input slices near zero and 64 are exact.  Every other slice uses the
ordinary BCH spectrum caps adversarially at each operator application.  This
is deliberately pessimistic.  Long-double/log-gamma arithmetic makes this a
design diagnostic, not yet a theorem-facing certificate.

At occupancy pole x=1 the suffix operator uses the exact global identity that
state XOR uniform input is uniform, avoiding inconsistent per-slice maxima.
Complement mode writes each input block as all-ones XOR a defect block and
uses the same operator with reversed XOR rows; it is intended for h>N/2.
"""

from __future__ import annotations

import argparse
import bisect
import math
from pathlib import Path

import numpy as np

from analyze_nosinger_bch_kernel import (
    B,
    EXACT_SLICES,
    SPECTRUM,
    accumulator_weight_entries,
    exact_slice,
    generator_rows,
    load_exact_slices,
)
from certificate_spectra import load_csv_spectrum, load_ebch128_spectrum
from punctured_ebch_outer import optimize_heterogeneous_log2_many


ROOT = Path(__file__).resolve().parent
OUTER_SPECTRUM = ROOT / "ebch128_64_spectrum.csv"
N = 2_097_408
INNER_BLOCKS = 32_772
OUTER_BLOCKS = 16_386
DISTANCE = 9 * N // 100


def log2_binomial(n: int, k: int) -> float:
    if k < 0 or k > n:
        return -math.inf
    return (
        math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
    ) / math.log(2)


def build_input_xor_matrix(pole: float) -> np.ndarray:
    matrix = np.zeros((B + 1, B + 1), dtype=np.longdouble)
    powers = [np.longdouble(pole) ** weight for weight in range(B + 1)]
    for state_weight in range(B + 1):
        for input_weight in range(B + 1):
            for intersection in range(
                max(0, state_weight + input_weight - B),
                min(state_weight, input_weight) + 1,
            ):
                xor_weight = state_weight + input_weight - 2 * intersection
                matrix[state_weight, xor_weight] += (
                    math.comb(state_weight, intersection)
                    * math.comb(B - state_weight, input_weight - intersection)
                    * powers[input_weight]
                )
    return matrix


def build_split_matrix(weights: tuple[int, ...], pole: float) -> np.ndarray:
    matrix = np.zeros((len(weights), B + 1), dtype=np.longdouble)
    for index, codeword_weight in enumerate(weights):
        denominator = math.comb(2 * B, codeword_weight)
        for state_weight in range(
            max(0, codeword_weight - B), min(B, codeword_weight) + 1
        ):
            output_weight = codeword_weight - state_weight
            matrix[index, state_weight] = (
                np.longdouble(math.comb(B, output_weight))
                * math.comb(B, state_weight)
                / denominator
                * np.longdouble(pole) ** output_weight
            )
    return matrix


def optimize_outer_log2(
    weight: int, spectrum: list[tuple[int, int]]
) -> float:
    low, high = 1e-10, 1.0 - 1e-10
    for _ in range(90):
        pole = (low + high) / 2
        terms = [count * pole**local_weight for local_weight, count in spectrum]
        mean = sum(
            local_weight * term
            for (local_weight, _count), term in zip(spectrum, terms)
        ) / sum(terms)
        if OUTER_BLOCKS * mean < weight:
            low = pole
        else:
            high = pole
    pole = (low + high) / 2
    enumerator = sum(
        count * pole**local_weight for local_weight, count in spectrum
    )
    return OUTER_BLOCKS * math.log2(enumerator) - weight * math.log2(pole)


def optimize_outer_log2_many(
    message_weights: np.ndarray, spectrum: list[tuple[int, int]]
) -> np.ndarray:
    """Vectorized Chernoff optimization for the symmetric h<=N/2 side."""
    local_weights = np.array([weight for weight, _count in spectrum], dtype=float)
    log_counts = np.log(
        np.array([count for _weight, count in spectrum], dtype=float)
    )
    result = np.empty(len(message_weights), dtype=float)
    for start in range(0, len(message_weights), 65536):
        stop = min(start + 65536, len(message_weights))
        weights = np.asarray(message_weights[start:stop], dtype=float)
        targets = weights / OUTER_BLOCKS
        low = np.full(len(weights), math.log(1e-10))
        high = np.full(len(weights), math.log(1.0 - 1e-10))
        for _ in range(70):
            log_poles = (low + high) / 2.0
            exponents = log_counts[:, None] + local_weights[:, None] * log_poles
            maxima = np.max(exponents, axis=0)
            terms = np.exp(exponents - maxima)
            normalizers = np.sum(terms, axis=0)
            means = np.sum(local_weights[:, None] * terms, axis=0) / normalizers
            go_right = means < targets
            low = np.where(go_right, log_poles, low)
            high = np.where(go_right, high, log_poles)
        log_poles = (low + high) / 2.0
        exponents = log_counts[:, None] + local_weights[:, None] * log_poles
        maxima = np.max(exponents, axis=0)
        log_enumerators = maxima + np.log(
            np.sum(np.exp(exponents - maxima), axis=0)
        )
        result[start:stop] = (
            OUTER_BLOCKS * log_enumerators - weights * log_poles
        ) / math.log(2.0)
    return result


def main() -> None:
    global N, INNER_BLOCKS, OUTER_BLOCKS, DISTANCE
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mode",
        choices=("first-active", "complement"),
        default="first-active",
        help="sparse first-active ledger or all-one/complement-defect ledger",
    )
    parser.add_argument(
        "--suffix-envelope",
        choices=("slice", "transport"),
        default="slice",
        help=(
            "unknown BCH slices maximized independently, or a per-transition "
            "transportation envelope enforcing all slice/weight marginals"
        ),
    )
    parser.add_argument("--h-min", type=int)
    parser.add_argument("--h-max", type=int)
    parser.add_argument(
        "--h-values",
        type=str,
        help="comma-separated diagnostic weights instead of a contiguous interval",
    )
    parser.add_argument(
        "--term-output",
        type=Path,
        help=(
            "write per-suffix log2 coefficient terms as a NumPy array; "
            "restricted to explicit h-values for paired-pole diagnostics"
        ),
    )
    parser.add_argument("--h-step", type=int, default=2)
    parser.add_argument("--output-pole", type=float, required=True)
    parser.add_argument("--occupancy-pole", type=float, required=True)
    parser.add_argument("--exact-radius", type=int, default=5)
    parser.add_argument(
        "--max-suffix-blocks",
        type=int,
        default=None,
        help=(
            "diagnostic truncation of the suffix recurrence; values below "
            f"{INNER_BLOCKS} do not bound the full placement sum"
        ),
    )
    parser.add_argument(
        "--power2-punctured",
        action="store_true",
        help="use 256 punctured EBCH blocks and N=2^21",
    )
    parser.add_argument("--spectrum", type=Path, default=SPECTRUM)
    parser.add_argument("--exact-slices", type=Path, default=EXACT_SLICES)
    parser.add_argument("--outer-spectrum", type=Path, default=OUTER_SPECTRUM)
    args = parser.parse_args()
    if args.power2_punctured:
        N = 1 << 21
        INNER_BLOCKS = N // B
        OUTER_BLOCKS = 16_386
        DISTANCE = 9 * N // 100
    if args.max_suffix_blocks is None:
        args.max_suffix_blocks = INNER_BLOCKS
    if args.h_values:
        selected_weights = sorted({int(value) for value in args.h_values.split(",")})
        if not selected_weights or selected_weights[0] < 1 or selected_weights[-1] > N:
            raise SystemExit("identity matrix probe: invalid h-values")
    else:
        if args.h_min is None or args.h_max is None:
            raise SystemExit("identity matrix probe: provide h-values or h-min/h-max")
        if not 1 <= args.h_min <= args.h_max <= N:
            raise SystemExit("identity matrix probe: invalid weight interval")
    if args.term_output and (not args.h_values or len(selected_weights) > 64):
        raise SystemExit("identity matrix probe: term-output needs at most 64 h-values")
    if args.term_output and args.mode != "first-active":
        raise SystemExit("identity matrix probe: term-output is first-active only")
    if args.h_step < 1:
        raise SystemExit("identity matrix probe: h-step must be positive")
    if not 1 <= args.max_suffix_blocks <= INNER_BLOCKS:
        raise SystemExit("identity matrix probe: invalid max-suffix-blocks")
    if not 0 < args.output_pole < 1 or not 0 < args.occupancy_pole:
        raise SystemExit("identity matrix probe: invalid pole")

    generator = generator_rows()
    exact_weights = tuple(range(1, args.exact_radius + 1)) + tuple(
        range(B - args.exact_radius, B + 1)
    )
    exact = {weight: exact_slice(generator, weight) for weight in exact_weights}
    for weight, histogram in load_exact_slices(args.exact_slices).items():
        exact[weight] = histogram

    loaded = load_ebch128_spectrum(args.spectrum)
    spectrum = {weight: count for weight, count in loaded if weight > 0 and count > 0}
    weights = tuple(sorted(spectrum))
    weight_index = {weight: index for index, weight in enumerate(weights)}
    populations = [math.comb(B, weight) for weight in range(B + 1)]
    known = frozenset(exact)

    exact_matrix = np.zeros((B + 1, len(weights)), dtype=np.longdouble)
    for input_weight, histogram in exact.items():
        for output_weight, count in histogram.items():
            exact_matrix[input_weight, weight_index[output_weight]] = (
                np.longdouble(count) / populations[input_weight]
            )

    xor_matrix = build_input_xor_matrix(args.occupancy_pole)
    split_matrix = build_split_matrix(weights, args.output_pole)

    def robust_rows(value: np.ndarray) -> np.ndarray:
        scores = split_matrix @ value
        order = np.argsort(-scores)
        ordered_counts = [spectrum[weights[index]] for index in order]
        cumulative_counts: list[int] = []
        running = 0
        for count in ordered_counts:
            running += count
            cumulative_counts.append(running)
        cumulative_values = np.cumsum(
            np.array(
                [
                    np.longdouble(count) * scores[index]
                    for count, index in zip(ordered_counts, order)
                ]
            )
        )

        result = exact_matrix @ scores
        result[0] = value[0]
        for input_weight in range(1, B + 1):
            if input_weight in known:
                continue
            population = populations[input_weight]
            index = bisect.bisect_left(cumulative_counts, population)
            previous_count = 0 if index == 0 else cumulative_counts[index - 1]
            previous_value = (
                np.longdouble(0) if index == 0 else cumulative_values[index - 1]
            )
            result[input_weight] = (
                previous_value
                + (population - previous_count) * scores[order[index]]
            ) / population
        return result

    def uniform_input_suffix(value: np.ndarray) -> np.ndarray:
        """Apply H_1 G exactly, using that state XOR uniform input is uniform."""
        scores = split_matrix @ value
        total = np.longdouble(value[0])
        for index, codeword_weight in enumerate(weights):
            total += np.longdouble(spectrum[codeword_weight]) * scores[index]
        return np.full(B + 1, total, dtype=np.longdouble)

    residual_spectrum = dict(spectrum)
    for histogram in exact.values():
        for codeword_weight, count in histogram.items():
            residual_spectrum[codeword_weight] -= count
            if residual_spectrum[codeword_weight] < 0:
                raise SystemExit("identity matrix probe: exact slices exceed spectrum")
    unknown_weights = tuple(weight for weight in range(1, B + 1) if weight not in known)
    if sum(populations[weight] for weight in unknown_weights) != sum(
        residual_spectrum.values()
    ):
        raise SystemExit("identity matrix probe: transportation marginals disagree")

    def transport_suffix(value: np.ndarray, *, complement: bool = False) -> np.ndarray:
        """Apply H_x G with one globally consistent unknown slice table."""
        scores = split_matrix @ value
        column_order = sorted(
            range(len(weights)), key=lambda index: scores[index], reverse=True
        )
        exact_averages = exact_matrix @ scores
        exact_averages[0] = value[0]
        result = np.zeros(B + 1, dtype=np.longdouble)
        for state_weight in range(B + 1):
            xor_row = xor_matrix[B - state_weight if complement else state_weight]
            exact_total = np.longdouble(0)
            for input_weight in (0, *sorted(known)):
                exact_total += xor_row[input_weight] * exact_averages[input_weight]

            rows = sorted(
                (
                    (
                        xor_row[input_weight] / populations[input_weight],
                        populations[input_weight],
                    )
                    for input_weight in unknown_weights
                ),
                key=lambda item: item[0],
                reverse=True,
            )
            row_index = 0
            row_remaining = rows[0][1]
            unknown_total = np.longdouble(0)
            for column_index in column_order:
                column_remaining = residual_spectrum[weights[column_index]]
                while column_remaining:
                    amount = min(row_remaining, column_remaining)
                    unknown_total += (
                        np.longdouble(amount)
                        * rows[row_index][0]
                        * scores[column_index]
                    )
                    row_remaining -= amount
                    column_remaining -= amount
                    if row_remaining == 0:
                        row_index += 1
                        if row_index < len(rows):
                            row_remaining = rows[row_index][1]
            result[state_weight] = exact_total + unknown_total
        return result

    if args.h_values:
        message_weights = np.array(selected_weights)
    else:
        message_weights = np.arange(args.h_min, args.h_max + 1, args.h_step)
    bounds = np.full(len(message_weights), -np.inf)
    saved_terms: list[np.ndarray] = []
    largest_term = -math.inf
    largest_term_suffix = -1
    value = np.ones(B + 1, dtype=np.longdouble)
    log_scale = 0.0
    log_occupancy = math.log2(args.occupancy_pole)
    if args.mode == "first-active":
        first_factors = np.array(
            [
                math.comb(B, weight) * args.occupancy_pole**weight
                for weight in range(B + 1)
            ]
        )
        for suffix_blocks in range(args.max_suffix_blocks):
            activation = robust_rows(value)
            prefix = np.cumsum(first_factors * np.asarray(activation, dtype=float))
            low = np.maximum(1, message_weights - B * suffix_blocks)
            high = np.minimum(B, message_weights)
            valid = low <= high
            segment = np.zeros(len(message_weights))
            indices = np.where(valid)[0]
            segment[indices] = prefix[high[indices]] - np.where(
                low[indices] > 0, prefix[low[indices] - 1], 0
            )
            term = np.full(len(message_weights), -np.inf)
            positive = segment > 0
            term[positive] = (
                log_scale
                + np.log2(segment[positive])
                - message_weights[positive] * log_occupancy
            )
            term_maximum = float(np.max(term))
            if term_maximum > largest_term:
                largest_term = term_maximum
                largest_term_suffix = suffix_blocks
            bounds = np.logaddexp2(bounds, term)
            if args.term_output:
                saved_terms.append(term.copy())

            if args.occupancy_pole == 1.0:
                value = uniform_input_suffix(value)
            elif args.suffix_envelope == "transport":
                value = transport_suffix(value)
            else:
                value = xor_matrix @ activation
            maximum = float(np.max(value))
            value /= maximum
            log_scale += math.log2(maximum)
    else:
        defect_weights = N - message_weights
        for block_index in range(args.max_suffix_blocks):
            activation = robust_rows(value)
            if args.occupancy_pole == 1.0:
                value = uniform_input_suffix(value)
            elif args.suffix_envelope == "transport":
                value = transport_suffix(value, complement=True)
            else:
                value = xor_matrix[::-1] @ activation
            maximum = float(np.max(value))
            value /= maximum
            log_scale += math.log2(maximum)
        bounds = (
            log_scale
            + np.log2(float(value[0]))
            - defect_weights * log_occupancy
        )
        largest_term = float(np.max(bounds))
        largest_term_suffix = args.max_suffix_blocks

    coefficient_weights = np.minimum(message_weights, N - message_weights)
    if args.power2_punctured:
        outer_bounds = optimize_heterogeneous_log2_many(
            coefficient_weights,
            full_blocks=16_130,
            punctured_blocks=256,
            spectrum_path=args.outer_spectrum,
        )
    else:
        outer = list(
            load_csv_spectrum(
                args.outer_spectrum,
                name="extended BCH [128,64,22] outer spectrum",
                length=128,
                dimension=64,
                minimum_distance=22,
            )
        )
        outer_by_weight = dict(outer)
        if any(
            outer_by_weight.get(weight, 0) != outer_by_weight.get(128 - weight, 0)
            for weight in range(129)
        ):
            raise SystemExit("identity matrix probe: outer spectrum is not symmetric")
        outer_bounds = optimize_outer_log2_many(coefficient_weights, outer)
    volume_logs = np.array(
        [log2_binomial(N, int(weight)) for weight in message_weights]
    )
    row_bounds = (
        outer_bounds
        - volume_logs
        - DISTANCE * math.log2(args.output_pole)
        + bounds
    )
    peak_index = max(range(len(row_bounds)), key=row_bounds.__getitem__)
    maximum = row_bounds[peak_index]
    total = maximum + math.log2(
        sum(2 ** (value - maximum) for value in row_bounds)
    )
    print(f"identity {args.mode} matrix ledger probe")
    print(
        f"construction={'power2-punctured' if args.power2_punctured else 'unpunctured'} "
        f"N={N} distance={DISTANCE} "
        f"mode={args.mode} "
        f"suffix_envelope={args.suffix_envelope} "
        f"h_values={','.join(str(int(weight)) for weight in message_weights) if args.h_values else '-'} "
        f"h_min={args.h_min} h_max={args.h_max} h_step={args.h_step} "
        f"output_pole={args.output_pole} occupancy_pole={args.occupancy_pole} "
        f"max_suffix_blocks={args.max_suffix_blocks}"
    )
    print(
        f"peak_h={int(message_weights[peak_index])} "
        f"peak_log2_upper_approx={maximum:.12f}"
    )
    if args.h_values:
        for weight, bound in zip(message_weights, row_bounds):
            print(f"sample_h={int(weight)} log2_upper_approx={bound:.12f}")
    print(f"interval_log2_upper_approx={total:.12f}")
    print(f"after_graph_r24_log2_upper_approx={total - 24:.12f}")
    print(
        f"recurrence_log2_scale={log_scale:.12f} "
        f"largest_coefficient_term_log2={largest_term:.12f} "
        f"at_suffix_blocks={largest_term_suffix}"
    )
    if args.max_suffix_blocks < INNER_BLOCKS:
        print("status=TRUNCATED_DIAGNOSTIC_NOT_A_FULL_PLACEMENT_BOUND")
    else:
        print("status=DIAGNOSTIC_ONLY_NEEDS_EXACT_OR_OUTWARD_ARITHMETIC")
    if args.term_output:
        np.save(args.term_output, np.vstack(saved_terms))
        print(f"term_output={args.term_output}")


if __name__ == "__main__":
    main()
