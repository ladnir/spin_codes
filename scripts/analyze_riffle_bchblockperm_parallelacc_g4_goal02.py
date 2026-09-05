#!/usr/bin/env python3
"""Transfer a three-block BCH slice to one uniform 384-bit superblock.

The real minimum outer profile contains three independently permuted BCH
words, each of weight 22.  Equivalently, it is a uniform weight-66 vector in
384 bits conditioned on every consecutive 128-bit part having weight 22.
For any event E after the independent packet and gap permutations,

    Pr[E | balanced] <= Pr[E] / Pr[balanced].

This script computes the conditioning probability and two bounds on Pr[E] for
the unconditioned superblock.  The first bound sums zero-state gaps exactly.
The second also uses the Hamming weights of positive accumulator states.  A
reproducible Monte Carlo search tests the real balanced law for active packet
orders with unusually many zero prefix states.  The script prints one JSON
receipt and does not modify the workspace.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math

import numpy as np
from scipy.optimize import minimize
from scipy.special import gammaln, logsumexp


PART_BITS = 128
PART_WEIGHT = 22
PART_COUNT = 3
PACKET_BITS = 4
SUPERBLOCK_BITS = PART_BITS * PART_COUNT
SUPERBLOCK_WEIGHT = PART_WEIGHT * PART_COUNT
PACKETS_PER_SUPERBLOCK = SUPERBLOCK_BITS // PACKET_BITS
STATE_COUNT = 1 << PACKET_BITS


def log_binom(n: int, k: int) -> float:
    if not 0 <= k <= n:
        return -math.inf
    return float(gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1))


def balance_log_probability() -> float:
    return (
        PART_COUNT * log_binom(PART_BITS, PART_WEIGHT)
        - log_binom(SUPERBLOCK_BITS, SUPERBLOCK_WEIGHT)
    )


def packet_support_count(input_weight: int, support: int) -> int:
    """Count weight-h vectors in 96 packets with exactly H nonzero packets."""
    coefficient = 0
    for selected in range(support + 1):
        available = PACKET_BITS * selected
        if available < input_weight:
            continue
        coefficient += (
            (-1) ** (support - selected)
            * math.comb(support, selected)
            * math.comb(available, input_weight)
        )
    return math.comb(PACKETS_PER_SUPERBLOCK, support) * coefficient


def weighted_path_log_counts(state_factor: np.ndarray) -> np.ndarray:
    """Return log ordered-state sums, rescaling after every packet slot."""
    shape = (SUPERBLOCK_WEIGHT + 1, PACKETS_PER_SUPERBLOCK + 1, STATE_COUNT)
    current = np.zeros(shape, dtype=np.float64)
    current[0, 0, 0] = 1.0
    states = np.arange(STATE_COUNT)
    log_scale = 0.0

    for slot in range(PACKETS_PER_SUPERBLOCK):
        following = np.zeros_like(current)
        max_weight = min(SUPERBLOCK_WEIGHT, PACKET_BITS * slot)
        max_support = min(slot, SUPERBLOCK_WEIGHT)
        for value in range(STATE_COUNT):
            weight = value.bit_count()
            support_increment = int(value != 0)
            if max_weight + weight > SUPERBLOCK_WEIGHT:
                source_weight_stop = SUPERBLOCK_WEIGHT + 1 - weight
            else:
                source_weight_stop = max_weight + 1
            source_support_stop = min(
                max_support + 1,
                PACKETS_PER_SUPERBLOCK + 1 - support_increment,
            )
            if source_weight_stop <= 0 or source_support_stop <= 0:
                continue
            source = current[
                :source_weight_stop,
                :source_support_stop,
                states ^ value,
            ]
            if value:
                source = source * state_factor[states]
            following[
                weight : weight + source_weight_stop,
                support_increment : support_increment + source_support_stop,
                states,
            ] += source
        scale = float(np.max(following))
        if not math.isfinite(scale) or scale <= 0.0:
            return np.full(PACKETS_PER_SUPERBLOCK + 1, math.inf)
        following /= scale
        log_scale += math.log(scale)
        current = following

    totals = np.sum(current[SUPERBLOCK_WEIGHT, :, :], axis=1)
    result = np.full(PACKETS_PER_SUPERBLOCK + 1, -math.inf)
    positive = totals > 0.0
    result[positive] = np.log(totals[positive]) + log_scale
    return result


def expected_log_path_products(
    u: float, scaled_v: float, packet_positions: int
) -> np.ndarray:
    """Log-average each ordered path product, separated by active support H."""
    v_gap = scaled_v / packet_positions
    state_factor = np.empty(STATE_COUNT, dtype=np.float64)
    for state in range(STATE_COUNT):
        state_weight = state.bit_count()
        exponent = v_gap + state_weight * u
        denominator = -math.expm1(-exponent)
        state_factor[state] = math.exp(-state_weight * u) / denominator

    normalization_log = log_binom(SUPERBLOCK_BITS, SUPERBLOCK_WEIGHT)
    return weighted_path_log_counts(state_factor) - normalization_log


def moment_log_bound(
    point: np.ndarray, distance: int, packet_positions: int
) -> float:
    u = float(point[0])
    scaled_v = float(point[1])
    log_averages = expected_log_path_products(u, scaled_v, packet_positions)
    v_gap = scaled_v / packet_positions
    terms = []
    for support, log_average in enumerate(log_averages):
        if not math.isfinite(log_average):
            continue
        terms.append(
            (packet_positions - support) * v_gap
            - log_binom(packet_positions, support)
            + log_average
        )
    if not terms:
        return math.inf
    return distance * u + float(logsumexp(terms))


def optimize_superblock_moment(
    distance: int, packet_positions: int, maxiter: int
) -> dict[str, object]:
    start_u = max(1e-3, math.log(packet_positions / max(1, distance)) / 4.0)
    starts = (np.asarray([start_u, float(SUPERBLOCK_WEIGHT)]),)
    results = []
    for start in starts:
        results.append(
            minimize(
                moment_log_bound,
                start,
                args=(distance, packet_positions),
                method="Nelder-Mead",
                bounds=((1e-9, 30.0), (1e-6, float(packet_positions))),
                options={"xatol": 2e-7, "fatol": 2e-9, "maxiter": maxiter},
            )
        )
    result = min(results, key=lambda item: float(item.fun))
    u = float(result.x[0])
    scaled_v = float(result.x[1])
    raw_log_bound = moment_log_bound(
        np.asarray([u, scaled_v]), distance, packet_positions
    )
    unconditional_log_bound = min(0.0, raw_log_bound)
    conditioning_log = balance_log_probability()
    conditional_log_bound = min(0.0, unconditional_log_bound - conditioning_log)
    return {
        "unconditional_log2_bound": unconditional_log_bound / math.log(2.0),
        "conditioning_penalty_bits": -conditioning_log / math.log(2.0),
        "balanced_three_block_log2_bound": conditional_log_bound / math.log(2.0),
        "balanced_effective_root_per_input_bit": math.exp(
            conditional_log_bound / SUPERBLOCK_WEIGHT
        ),
        "z": math.exp(-u),
        "tau": math.exp(-scaled_v / packet_positions),
        "scaled_gap_parameter": scaled_v,
        "optimizer_success": bool(result.success),
        "optimizer_message": str(result.message),
        "optimizer_evaluations_selected_start": int(result.nfev),
        "optimizer_total_evaluations": int(sum(item.nfev for item in results)),
    }


def zero_return_joint_counts(
    positive_state_factor: np.ndarray | None = None,
) -> tuple[np.ndarray, float, float]:
    """Count paths by (H,m0), optionally weighting positive prefix states."""
    max_returns = SUPERBLOCK_WEIGHT // 2
    shape = (
        SUPERBLOCK_WEIGHT + 1,
        PACKETS_PER_SUPERBLOCK + 1,
        STATE_COUNT,
        max_returns + 1,
    )
    current = np.zeros(shape, dtype=np.float64)
    current[0, 0, 0, 0] = 1.0
    log_scale = 0.0

    for slot in range(PACKETS_PER_SUPERBLOCK):
        following = np.zeros_like(current)
        max_weight = min(SUPERBLOCK_WEIGHT, PACKET_BITS * slot)
        max_support = min(slot, SUPERBLOCK_WEIGHT)
        max_observed_returns = min(max_returns, max_support // 2)
        for value in range(STATE_COUNT):
            value_weight = value.bit_count()
            support_increment = int(value != 0)
            source_weight_stop = min(
                max_weight + 1,
                SUPERBLOCK_WEIGHT + 1 - value_weight,
            )
            source_support_stop = min(
                max_support + 1,
                PACKETS_PER_SUPERBLOCK + 1 - support_increment,
            )
            if source_weight_stop <= 0 or source_support_stop <= 0:
                continue
            for next_state in range(STATE_COUNT):
                return_increment = int(value != 0 and next_state == 0)
                source_return_stop = min(
                    max_observed_returns + 1,
                    max_returns + 1 - return_increment,
                )
                following[
                    value_weight : value_weight + source_weight_stop,
                    support_increment : support_increment + source_support_stop,
                    next_state,
                    return_increment : return_increment + source_return_stop,
                ] += current[
                    :source_weight_stop,
                    :source_support_stop,
                    next_state ^ value,
                    :source_return_stop,
                ] * (
                    positive_state_factor[next_state]
                    if value != 0
                    and next_state != 0
                    and positive_state_factor is not None
                    else 1.0
                )
        scale = float(np.max(following))
        if not math.isfinite(scale) or scale <= 0.0:
            raise FloatingPointError("zero-return DP produced an invalid scale")
        following /= scale
        log_scale += math.log(scale)
        current = following

    normalization_log = log_binom(SUPERBLOCK_BITS, SUPERBLOCK_WEIGHT)
    joint_counts = np.sum(current[SUPERBLOCK_WEIGHT, :, :, :], axis=1)
    total_log = math.log(float(np.sum(joint_counts))) + log_scale
    if positive_state_factor is None and not math.isclose(
        total_log, normalization_log, rel_tol=0.0, abs_tol=3e-10
    ):
        raise RuntimeError("zero-return DP failed the uniform-slice normalization gate")
    if positive_state_factor is None:
        support_counts = np.sum(joint_counts, axis=1)
        for support in range(1, SUPERBLOCK_WEIGHT + 1):
            exact_count = packet_support_count(SUPERBLOCK_WEIGHT, support)
            observed_count = float(support_counts[support])
            if exact_count == 0 and observed_count == 0.0:
                continue
            observed_log = math.log(observed_count) + log_scale
            if not math.isclose(
                observed_log, math.log(exact_count), rel_tol=0.0, abs_tol=4e-10
            ):
                raise RuntimeError("zero-return DP failed the packet-support gate")
    return joint_counts, log_scale, normalization_log


def exact_zero_gap_statistics(
    distance: int, packet_positions: int
) -> dict[str, object]:
    """Count (H,m0) exactly and sum every zero-state gap without cost."""
    max_returns = SUPERBLOCK_WEIGHT // 2
    joint_counts, log_scale, normalization_log = zero_return_joint_counts()
    unconditioned_bound = 0.0
    contributions = []
    for support in range(1, SUPERBLOCK_WEIGHT + 1):
        for zero_returns in range(max_returns + 1):
            count = float(joint_counts[support, zero_returns])
            if count <= 0.0:
                continue
            slice_log_probability = math.log(count) + log_scale - normalization_log
            positive_gaps = support - zero_returns
            if distance < positive_gaps:
                conditional_log_bound = -math.inf
            else:
                conditional_log_bound = min(
                    0.0,
                    log_binom(distance, positive_gaps)
                    + log_binom(
                        packet_positions - support + zero_returns,
                        zero_returns,
                    )
                    - log_binom(packet_positions, support),
                )
            contribution_log = slice_log_probability + conditional_log_bound
            contribution = math.exp(contribution_log)
            unconditioned_bound += contribution
            contributions.append(
                {
                    "packet_support": support,
                    "zero_prefix_states": zero_returns,
                    "slice_log2_probability": slice_log_probability / math.log(2.0),
                    "conditional_gap_log2_bound": conditional_log_bound
                    / math.log(2.0),
                    "weighted_log2_contribution": contribution_log / math.log(2.0),
                }
            )

    unconditional_log_bound = min(0.0, math.log(unconditioned_bound))
    conditioning_log = balance_log_probability()
    balanced_log_bound = min(0.0, unconditional_log_bound - conditioning_log)
    return_counts = np.sum(joint_counts, axis=0)
    return_probabilities = return_counts * math.exp(log_scale - normalization_log)
    return {
        "unconditional_log2_bound": unconditional_log_bound / math.log(2.0),
        "balanced_three_block_log2_bound": balanced_log_bound / math.log(2.0),
        "balanced_effective_root_per_input_bit": math.exp(
            balanced_log_bound / SUPERBLOCK_WEIGHT
        ),
        "probability_any_zero_prefix_state_unconditioned": float(
            1.0 - return_probabilities[0]
        ),
        "mean_zero_prefix_states_unconditioned": float(
            sum(
                index * probability
                for index, probability in enumerate(return_probabilities)
            )
        ),
        "top_weighted_contributions": sorted(
            contributions,
            key=lambda item: item["weighted_log2_contribution"],
            reverse=True,
        )[:15],
    }


def weighted_positive_gap_bound(
    distance: int, packet_positions: int, scaled_cost: float
) -> dict[str, object]:
    """Use actual positive state weights and sum zero-state gaps exactly."""
    if distance <= 0:
        raise ValueError("weighted positive-gap bound requires positive distance")
    u = scaled_cost / distance
    positive_state_factor = np.ones(STATE_COUNT, dtype=np.float64)
    for state in range(1, STATE_COUNT):
        state_weight = state.bit_count()
        exponent = u * state_weight
        positive_state_factor[state] = math.exp(-exponent) / -math.expm1(-exponent)

    joint_counts, log_scale, normalization_log = zero_return_joint_counts(
        positive_state_factor
    )
    terms = []
    contributions = []
    max_returns = SUPERBLOCK_WEIGHT // 2
    for support in range(1, SUPERBLOCK_WEIGHT + 1):
        for zero_returns in range(max_returns + 1):
            count = float(joint_counts[support, zero_returns])
            if count <= 0.0:
                continue
            term = (
                math.log(count)
                + log_scale
                - normalization_log
                + log_binom(
                    packet_positions - support + zero_returns,
                    zero_returns,
                )
                - log_binom(packet_positions, support)
            )
            terms.append(term)
            contributions.append(
                {
                    "packet_support": support,
                    "zero_prefix_states": zero_returns,
                    "pre_chernoff_log2_contribution": term / math.log(2.0),
                }
            )
    raw_log_bound = distance * u + float(logsumexp(terms))
    unconditional_log_bound = min(0.0, raw_log_bound)
    conditioning_log = balance_log_probability()
    balanced_log_bound = min(0.0, unconditional_log_bound - conditioning_log)
    return {
        "scaled_cost_parameter": scaled_cost,
        "z": math.exp(-u),
        "unconditional_log2_bound": unconditional_log_bound / math.log(2.0),
        "balanced_three_block_log2_bound": balanced_log_bound / math.log(2.0),
        "balanced_effective_root_per_input_bit": math.exp(
            balanced_log_bound / SUPERBLOCK_WEIGHT
        ),
        "top_pre_chernoff_contributions": sorted(
            contributions,
            key=lambda item: item["pre_chernoff_log2_contribution"],
            reverse=True,
        )[:15],
        "interpretation": (
            "Every zero-state gap is free. Each positive-state gap uses its actual "
            "Hamming weight through z^w/(1-z^w)."
        ),
    }


def validate_gap_inequalities() -> str:
    """Brute-force the two per-path gap inequalities on small instances."""
    u = 0.37
    cases = 0
    for active_count in range(1, 5):
        for state_weights in itertools.product(range(PACKET_BITS + 1), repeat=active_count):
            zero_count = state_weights.count(0)
            positive_count = active_count - zero_count
            positive_factor = 1.0
            for weight in state_weights:
                if weight:
                    positive_factor *= math.exp(-u * weight) / -math.expm1(
                        -u * weight
                    )
            for packet_positions in range(active_count, 10):
                for distance in range(0, 4 * packet_positions + 1):
                    exact = 0
                    for positions in itertools.combinations(
                        range(1, packet_positions + 1), active_count
                    ):
                        durations = [
                            positions[index + 1] - positions[index]
                            for index in range(active_count - 1)
                        ]
                        durations.append(packet_positions + 1 - positions[-1])
                        output_weight = sum(
                            weight * duration
                            for weight, duration in zip(state_weights, durations)
                        )
                        exact += int(output_weight <= distance)

                    zero_allocations = math.comb(
                        packet_positions - active_count + zero_count,
                        zero_count,
                    )
                    coarse = (
                        math.comb(distance, positive_count) * zero_allocations
                        if distance >= positive_count
                        else 0
                    )
                    weighted = (
                        math.exp(distance * u)
                        * positive_factor
                        * zero_allocations
                    )
                    if exact > coarse or exact > weighted + 2e-10:
                        raise RuntimeError("small-instance gap inequality gate failed")
                    cases += 1
    return f"PASS ({cases} exhaustive small cases)"


def sample_weight_slice(rng: np.random.Generator, bits: int, weight: int) -> np.ndarray:
    vector = np.zeros(bits, dtype=np.uint8)
    vector[rng.choice(bits, size=weight, replace=False)] = 1
    return vector


def packets_from_bits(bits: np.ndarray) -> np.ndarray:
    reshaped = bits.reshape(-1, PACKET_BITS)
    powers = np.asarray([1, 2, 4, 8], dtype=np.uint8)
    return reshaped @ powers


def refutation_search(samples: int, seed: int) -> dict[str, object]:
    """Sample the real balanced law and random active packet interleavings."""
    rng = np.random.default_rng(seed)
    zero_returns = np.empty(samples, dtype=np.int16)
    supports = np.empty(samples, dtype=np.int16)
    positive_state_weight_sums = np.empty(samples, dtype=np.int16)
    best_zero_returns = -1
    best_record: dict[str, object] | None = None

    for sample_index in range(samples):
        packets = []
        for _ in range(PART_COUNT):
            packets.append(
                packets_from_bits(sample_weight_slice(rng, PART_BITS, PART_WEIGHT))
            )
        active = np.concatenate(packets)
        active = active[active != 0]
        rng.shuffle(active)

        state = 0
        returns = 0
        state_weight_sum = 0
        for value in active:
            state ^= int(value)
            if state == 0:
                returns += 1
            else:
                state_weight_sum += state.bit_count()

        support = int(active.size)
        zero_returns[sample_index] = returns
        supports[sample_index] = support
        positive_state_weight_sums[sample_index] = state_weight_sum
        if returns > best_zero_returns:
            best_zero_returns = returns
            best_record = {
                "sample_index": sample_index,
                "packet_support": support,
                "zero_prefix_states": returns,
                "positive_prefix_state_weight_sum": state_weight_sum,
                "active_packet_values": [int(value) for value in active],
            }

    assert best_record is not None
    quantiles = [0.5, 0.9, 0.99, 0.999]
    return {
        "law": "three independent weight-22 block slices; uniform active-packet order",
        "samples": samples,
        "seed": seed,
        "mean_packet_support": float(np.mean(supports)),
        "mean_zero_prefix_states": float(np.mean(zero_returns)),
        "probability_any_zero_prefix_state": float(np.mean(zero_returns > 0)),
        "zero_return_quantiles": {
            str(quantile): float(np.quantile(zero_returns, quantile, method="higher"))
            for quantile in quantiles
        },
        "maximum_zero_prefix_states_observed": int(np.max(zero_returns)),
        "minimum_positive_prefix_state_weight_sum_observed": int(
            np.min(positive_state_weight_sums)
        ),
        "highest_zero_return_sample": best_record,
        "interpretation": (
            "Numerical obstruction search only. It neither upper-bounds the bad "
            "event nor proves that the recorded active order has appreciable law."
        ),
    }


def main() -> None:
    global PART_WEIGHT, SUPERBLOCK_WEIGHT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--part-weight", type=int, default=PART_WEIGHT)
    parser.add_argument("--packet-positions", type=int, default=524_352)
    parser.add_argument("--distance", type=int, default=188_766)
    parser.add_argument("--search-samples", type=int, default=20_000)
    parser.add_argument("--seed", type=int, default=20_260_823)
    parser.add_argument("--optimizer-maxiter", type=int, default=8)
    parser.add_argument("--include-joint-moment", action="store_true")
    parser.add_argument("--skip-exact-zero-gap", action="store_true")
    parser.add_argument("--positive-scaled-cost", type=float, default=33.7)
    parser.add_argument("--skip-positive-weight-bound", action="store_true")
    parser.add_argument("--skip-refutation-search", action="store_true")
    args = parser.parse_args()

    if not 1 <= args.part_weight <= PART_BITS:
        raise ValueError("part weight lies outside a BCH block")
    PART_WEIGHT = args.part_weight
    SUPERBLOCK_WEIGHT = PART_WEIGHT * PART_COUNT
    if args.packet_positions < PACKETS_PER_SUPERBLOCK:
        raise ValueError("global packet count must contain the three BCH blocks")
    if not 0 <= args.distance <= PACKET_BITS * args.packet_positions:
        raise ValueError("distance lies outside the output block")
    if args.search_samples < 1:
        raise ValueError("search sample count must be positive")
    if args.optimizer_maxiter < 1:
        raise ValueError("optimizer iteration count must be positive")
    if args.positive_scaled_cost <= 0.0:
        raise ValueError("positive scaled cost must be positive")

    conditioning_log = balance_log_probability()
    moment = (
        optimize_superblock_moment(
            args.distance, args.packet_positions, args.optimizer_maxiter
        )
        if args.include_joint_moment
        else {
            "status": "SKIPPED_BY_DEFAULT",
            "reason": "Goal 01 established the zero-state singularity; exact zero-gap bounds replace it.",
        }
    )
    zero_gap = (
        None
        if args.skip_exact_zero_gap
        else exact_zero_gap_statistics(args.distance, args.packet_positions)
    )
    positive_weight_bound = (
        None
        if args.skip_positive_weight_bound
        else weighted_positive_gap_bound(
            args.distance, args.packet_positions, args.positive_scaled_cost
        )
    )
    search = (
        None
        if args.skip_refutation_search
        else refutation_search(args.search_samples, args.seed)
    )
    payload = {
        "schema": "riffle-bchblockperm-parallelacc-g4-goal02-v1",
        "evidence_label": "RIGOROUS_SUPERBLOCK_TRANSFER_WITH_NUMERICAL_OBSTRUCTION_SEARCH",
        "parameters": {
            "bch_block_bits": PART_BITS,
            "bch_block_weight_profile": [PART_WEIGHT] * PART_COUNT,
            "packet_bits": PACKET_BITS,
            "packets_in_three_blocks": PACKETS_PER_SUPERBLOCK,
            "global_packet_positions": args.packet_positions,
            "global_binary_length": PACKET_BITS * args.packet_positions,
            "failure_weight_inclusive": args.distance,
        },
        "superblock_transfer": {
            "unconditioned_law": (
                f"uniform binary weight-{SUPERBLOCK_WEIGHT} slice in 384 coordinates"
            ),
            "balanced_event": (
                f"each labeled 128-bit part has binary weight {PART_WEIGHT}"
            ),
            "balance_probability": math.exp(conditioning_log),
            "balance_log2_probability": conditioning_log / math.log(2.0),
            "radon_nikodym_upper_bound": math.exp(-conditioning_log),
            "inequality": "Pr_real[E] <= Pr_superblock[E] / Pr_superblock[balanced]",
        },
        "joint_lane_moment": moment,
        "exact_zero_gap": zero_gap,
        "positive_state_weight_bound": positive_weight_bound,
        "refutation_search": search,
        "validation": {
            "uniform_slice_normalization_gate": (
                "PASS" if zero_gap is not None else "NOT_RUN"
            ),
            "packet_support_marginal_gate": (
                "PASS" if zero_gap is not None else "NOT_RUN"
            ),
            "small_gap_inequality_gate": validate_gap_inequalities(),
        },
        "scope": (
            "The rigorous bound covers only the declared equal-weight profile "
            f"({PART_WEIGHT},{PART_WEIGHT},{PART_WEIGHT}). It does not cover other BCH weights or sum the field-outer "
            "spectrum. The refutation search is numerical evidence only."
        ),
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
