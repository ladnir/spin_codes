#!/usr/bin/env python3
"""Collapse shared packet profiles through an active-packet envelope.

For a fixed Chernoff point, let U_k be the one-region transfer conditioned on
exactly k nonzero packets after PacketMul.  For 0 < theta <= 1, define the
entrywise matrix envelope

    M(theta)[a,b] = max_k U_k[a,b] / theta**k.

If a shared group contains r active outer blocks, its packet survives in one
region with probability q_r = 1 - (1-p)**r.  Consequently every profile n
satisfies

    R_n <= M(theta) * product_r (1-q_r+q_r*theta)**n_r.

The scalar factor commutes with the 256-region matrix product.  Summing all
profiles of occupation h then becomes one coefficient of a degree-four
polynomial raised to the 2048th power.  The calculation below tests the loss
of this profile-free bound against the recorded exact profile sum.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from analyze_riffle_ldpcsplitstate_occupation_ladder import (
    epoch_transfers,
    matrix_power_moment,
    regular_outer_density_log2,
)
from analyze_riffle_ldpcsplitstate_packetmul_shared_groups import (
    universal_nonzero_region_transfers,
)


CONSTRUCTION = Path(
    "constructions/riffle_ldpcsplitstate_packetmul_g4_t256_s64"
)
DEFAULT_ACTIVATION = Path(
    "constructions/riffle_ldpcsplitstate_g4_t256_s64/"
    "receipts/zero_state_activation_table.json"
)
DEFAULT_REFERENCE = CONSTRUCTION / "receipts/occupation128_full_profiles_fixed_point.json"
DEFAULT_OUTPUT = CONSTRUCTION / "receipts/active_count_envelope_occupation128.json"


def logadd(left: float, right: float) -> float:
    """Return log(exp(left)+exp(right)) without underflow."""
    if left == -math.inf:
        return right
    if right == -math.inf:
        return left
    if left < right:
        left, right = right, left
    return left + math.log1p(math.exp(right - left))


def weighted_occupation_coefficient_log(
    occupation: int, log_group_weights: list[float], groups: int = 2048
) -> float:
    """Return log [x^occupation] (sum_r C(4,r) a_r x^r)^groups."""
    if len(log_group_weights) != 5:
        raise ValueError("expected group weights for widths zero through four")
    log_terms = [
        math.log(math.comb(4, width)) + log_group_weights[width]
        for width in range(5)
    ]
    coefficients = [-math.inf] * (occupation + 1)
    coefficients[0] = 0.0
    for _ in range(groups):
        following = [-math.inf] * (occupation + 1)
        for degree, value in enumerate(coefficients):
            if value == -math.inf:
                continue
            for width, term in enumerate(log_terms):
                target = degree + width
                if target > occupation:
                    break
                following[target] = logadd(following[target], value + term)
        coefficients = following
    return coefficients[occupation]


def entrywise_exponential_envelope(
    universal: list[np.ndarray], theta: float
) -> tuple[np.ndarray, list[int]]:
    """Return M(theta) and the maximizing k for each matrix entry."""
    if not 0.0 < theta <= 1.0:
        raise ValueError("theta must lie in (0,1]")
    matrix = np.zeros((2, 2), dtype=np.float64)
    witnesses = [0, 0, 0, 0]
    log_theta = math.log(theta)
    for row in range(2):
        for column in range(2):
            best_log = -math.inf
            best_k = 0
            for packet_count, transfer in enumerate(universal):
                value = float(transfer[row, column])
                if value <= 0.0:
                    continue
                candidate = math.log(value) - packet_count * log_theta
                if candidate > best_log:
                    best_log = candidate
                    best_k = packet_count
            matrix[row, column] = 0.0 if best_log == -math.inf else math.exp(best_log)
            witnesses[2 * row + column] = best_k
    return matrix, witnesses


def log_entrywise_exponential_envelope(
    universal: list[np.ndarray], theta: float, row: int, column: int
) -> tuple[float, int]:
    """Return log M_ij(theta) and a maximizing packet count."""
    log_theta = math.log(theta)
    best = -math.inf
    witness = 0
    for packet_count, transfer in enumerate(universal):
        value = float(transfer[row, column])
        if value <= 0.0:
            continue
        candidate = math.log(value) - packet_count * log_theta
        if candidate > best:
            best = candidate
            witness = packet_count
    return best, witness


def optimize_entry_thetas(
    universal: list[np.ndarray], packed_groups: int
) -> tuple[list[float], list[float], list[int]]:
    """Tune each transition against the persistent all-width-four profile."""
    thetas: list[float] = []
    log_envelopes: list[float] = []
    witnesses: list[int] = []
    coarse_logs = np.linspace(math.log(0.001), math.log(0.999), 1001)
    for row in range(2):
        for column in range(2):
            candidates = []
            for log_theta in coarse_logs:
                theta = math.exp(float(log_theta))
                log_envelope, witness = log_entrywise_exponential_envelope(
                    universal, theta, row, column
                )
                q4 = 15.0 / 16.0
                log_scalar = math.log1p(q4 * (theta - 1.0))
                candidates.append(
                    (log_envelope + packed_groups * log_scalar, theta, witness)
                )
            _, coarse_theta, _ = min(candidates)
            step = float(coarse_logs[1] - coarse_logs[0])
            refined_logs = np.linspace(
                math.log(coarse_theta) - step,
                math.log(coarse_theta) + step,
                401,
            )
            refined = []
            for log_theta in refined_logs:
                theta = min(0.999999, max(0.000001, math.exp(float(log_theta))))
                log_envelope, witness = log_entrywise_exponential_envelope(
                    universal, theta, row, column
                )
                q4 = 15.0 / 16.0
                log_scalar = math.log1p(q4 * (theta - 1.0))
                refined.append(
                    (log_envelope + packed_groups * log_scalar, theta, witness)
                )
            _, theta, witness = min(refined)
            log_envelope, witness = log_entrywise_exponential_envelope(
                universal, theta, row, column
            )
            thetas.append(theta)
            log_envelopes.append(log_envelope)
            witnesses.append(witness)
    return thetas, log_envelopes, witnesses


def polynomial_power_coefficient_log(
    occupation: int, log_terms: list[float], power: int = 2048
) -> float:
    """Return log [x^occupation] P(x)^power for P(0)=1, deg(P)=4."""
    values = [-math.inf] * (occupation + 1)
    values[0] = 0.0
    for degree in range(1, occupation + 1):
        terms = []
        for width in range(1, min(4, degree) + 1):
            multiplier = (power + 1) * width - degree
            terms.append(
                math.log(multiplier)
                + log_terms[width]
                + values[degree - width]
            )
        total = -math.inf
        for term in terms:
            total = logadd(total, term)
        values[degree] = total - math.log(degree)
    return values[occupation]


def state_path_types(regions: int):
    """Yield transition counts and multiplicities for paths starting in Z."""
    yield (regions, 0, 0, 0), 1
    for crossings in range(1, regions // 2 + 1):
        remaining = regions - 2 * crossings
        for count00 in range(remaining + 1):
            count11 = remaining - count00
            multiplicity = math.comb(count00 + crossings, crossings)
            multiplicity *= math.comb(count11 + crossings - 1, crossings - 1)
            yield (count00, crossings, crossings, count11), multiplicity
    for returns in range((regions - 1) // 2 + 1):
        count01 = returns + 1
        count10 = returns
        remaining = regions - count01 - count10
        for count00 in range(remaining + 1):
            count11 = remaining - count00
            multiplicity = math.comb(count00 + returns, returns)
            multiplicity *= math.comb(count11 + returns, returns)
            yield (count00, count01, count10, count11), multiplicity


def per_transition_path_envelope(
    *,
    universal: list[np.ndarray],
    occupation: int,
    probability: float,
    regions: int,
    conditioning_bits: float,
    target_bits: float,
) -> dict[str, object]:
    """Sum a four-rate active-count envelope over all binary state paths."""
    if occupation % 4:
        raise ValueError("the current theta tuner expects an all-packed endpoint")
    thetas, log_envelopes, witnesses = optimize_entry_thetas(
        universal, occupation // 4
    )
    log_factors = [[0.0] * 4 for _ in range(5)]
    for width in range(1, 5):
        survival = 1.0 - (1.0 - probability) ** width
        for transition, theta in enumerate(thetas):
            log_factors[width][transition] = math.log1p(
                survival * (theta - 1.0)
            )

    path_total = -math.inf
    path_type_count = 0
    multiplicity_total = 0
    largest = None
    for counts, multiplicity in state_path_types(regions):
        path_type_count += 1
        multiplicity_total += multiplicity
        log_matrix_weight = math.log(multiplicity) + sum(
            counts[index] * log_envelopes[index] for index in range(4)
        )
        log_terms = [0.0]
        for width in range(1, 5):
            profile_scalar = sum(
                counts[index] * log_factors[width][index]
                for index in range(4)
            )
            log_terms.append(math.log(math.comb(4, width)) + profile_scalar)
        coefficient = polynomial_power_coefficient_log(occupation, log_terms)
        contribution = log_matrix_weight + coefficient
        path_total = logadd(path_total, contribution)
        if largest is None or contribution > largest[0]:
            largest = (contribution, counts, multiplicity, coefficient)

    if multiplicity_total != 1 << regions:
        raise RuntimeError("state-path multiplicities do not sum to 2^regions")
    aggregate = path_total / math.log(2.0) + conditioning_bits + target_bits
    assert largest is not None
    return {
        "aggregate_log2_bound": aggregate,
        "aggregate_margin_bits": -aggregate,
        "theta_by_transition_00_01_10_11": thetas,
        "log2_envelope_by_transition_00_01_10_11": [
            value / math.log(2.0) for value in log_envelopes
        ],
        "maximizing_packet_count_by_transition_00_01_10_11": witnesses,
        "state_path_type_count": path_type_count,
        "state_path_multiplicity_audit": "PASS",
        "largest_path_type": {
            "transition_counts_00_01_10_11": list(largest[1]),
            "state_sequence_multiplicity_log2": math.log2(largest[2]),
            "weighted_profile_coefficient_log2": largest[3] / math.log(2.0),
            "path_contribution_log2": largest[0] / math.log(2.0),
        },
    }


def load_activation(path: Path, maximum: int) -> list[float]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    activation = [0.0] * (maximum + 1)
    activation[0] = 1.0
    for row in payload["by_total_weight"]:
        weight = int(row["total_weight"])
        if weight > maximum:
            break
        activation[weight] = float(row["maximum_distinct_conditioned_upper_bound"])
    return activation


def evaluate_theta(
    *,
    theta: float,
    universal: list[np.ndarray],
    occupation: int,
    probability: float,
    regions: int,
    conditioning_bits: float,
    target_bits: float,
) -> dict[str, object]:
    envelope, witnesses = entrywise_exponential_envelope(universal, theta)
    matrix_log2 = matrix_power_moment(envelope, regions) / math.log(2.0)

    log_group_weights = [0.0]
    for width in range(1, 5):
        survival = 1.0 - (1.0 - probability) ** width
        scalar = 1.0 - survival + survival * theta
        log_group_weights.append(regions * math.log(scalar))
    coefficient_log2 = (
        weighted_occupation_coefficient_log(occupation, log_group_weights)
        / math.log(2.0)
    )
    aggregate = matrix_log2 + coefficient_log2 + conditioning_bits + target_bits
    return {
        "theta": theta,
        "aggregate_log2_bound": aggregate,
        "aggregate_margin_bits": -aggregate,
        "matrix_power_log2": matrix_log2,
        "weighted_profile_coefficient_log2": coefficient_log2,
        "envelope_matrix": envelope.tolist(),
        "entrywise_maximizing_packet_counts": witnesses,
    }


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    reference = json.loads(args.reference.read_text(encoding="utf-8"))
    rows = {
        int(row["occupation"]): row for row in reference["occupation_rows"]
    }
    if args.occupation not in rows:
        raise ValueError(f"reference has no occupation {args.occupation}")
    reference_row = rows[args.occupation]
    log_surprisal = float(reference_row["log_surprisal"])
    probability = float(reference_row["placement_probability"])
    surprisal = math.exp(log_surprisal)
    z = math.exp(-surprisal)
    activation = load_activation(args.activation, 256)
    epoch = epoch_transfers(
        z=z,
        distance=args.constituent_distance,
        moment_order=args.live_moment_order,
        activation_upper=activation,
        maximum=256,
    )
    universal = universal_nonzero_region_transfers(epoch, args.occupation)

    conditioning_bits = args.occupation * regular_outer_density_log2(probability)
    target_distance = math.floor(args.relative_distance * (1 << 21))
    target_bits = target_distance * surprisal / math.log(2.0)

    coarse_thetas = np.exp(
        np.linspace(math.log(args.minimum_theta), math.log(args.maximum_theta), args.grid)
    )
    coarse = [
        evaluate_theta(
            theta=float(theta),
            universal=universal,
            occupation=args.occupation,
            probability=probability,
            regions=args.regions,
            conditioning_bits=conditioning_bits,
            target_bits=target_bits,
        )
        for theta in coarse_thetas
    ]
    best_coarse = min(coarse, key=lambda row: row["aggregate_log2_bound"])
    center = math.log(float(best_coarse["theta"]))
    coarse_step = (
        math.log(args.maximum_theta) - math.log(args.minimum_theta)
    ) / max(1, args.grid - 1)
    refined_thetas = np.exp(np.linspace(center - coarse_step, center + coarse_step, args.refine_grid))
    refined = [
        evaluate_theta(
            theta=float(min(1.0, max(args.minimum_theta, theta))),
            universal=universal,
            occupation=args.occupation,
            probability=probability,
            regions=args.regions,
            conditioning_bits=conditioning_bits,
            target_bits=target_bits,
        )
        for theta in refined_thetas
    ]
    best = min([*coarse, *refined], key=lambda row: row["aggregate_log2_bound"])
    exact_margin = float(reference_row["profile_sum_margin_bits"])
    best["loss_from_exact_profile_sum_bits"] = exact_margin - float(
        best["aggregate_margin_bits"]
    )
    path_envelope = per_transition_path_envelope(
        universal=universal,
        occupation=args.occupation,
        probability=probability,
        regions=args.regions,
        conditioning_bits=conditioning_bits,
        target_bits=target_bits,
    )
    path_envelope["loss_from_exact_profile_sum_bits"] = exact_margin - float(
        path_envelope["aggregate_margin_bits"]
    )
    print(
        f"occupation,{args.occupation},theta,{best['theta']:.12g},"
        f"envelope_margin,{best['aggregate_margin_bits']:.6f},"
        f"path_envelope_margin,{path_envelope['aggregate_margin_bits']:.6f},"
        f"exact_profile_sum_margin,{exact_margin:.6f},"
        f"path_loss,{path_envelope['loss_from_exact_profile_sum_bits']:.6f}",
        flush=True,
    )
    return {
        "schema": "riffle-ldpcsplitstate-packetmul-active-count-envelope-v1",
        "candidate": "Riffle LDPCSplitState PacketMul g=4 t=256 s=64",
        "parameters": {
            "occupation": args.occupation,
            "regions": args.regions,
            "placement_probability": probability,
            "log_surprisal": log_surprisal,
            "constituent_distance": args.constituent_distance,
            "live_moment_order": args.live_moment_order,
            "relative_distance": args.relative_distance,
        },
        "reference_exact_profile_sum_margin_bits": exact_margin,
        "best_active_count_envelope": best,
        "per_transition_path_envelope": path_envelope,
        "coarse_grid": coarse,
        "scope": (
            "Float64 diagnostic of a rigorous algebraic envelope. The bound "
            "uses only the number of surviving nonzero packets inside the "
            "inner transfer and sums all shared-group profiles with one "
            "univariate occupation coefficient. It is not outward-rounded."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--activation", type=Path, default=DEFAULT_ACTIVATION)
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--occupation", type=int, default=128)
    parser.add_argument("--regions", type=int, default=256)
    parser.add_argument("--constituent-distance", type=int, default=40)
    parser.add_argument("--live-moment-order", type=int, choices=(1, 2, 3), default=3)
    parser.add_argument("--relative-distance", type=float, default=0.09)
    parser.add_argument("--minimum-theta", type=float, default=0.005)
    parser.add_argument("--maximum-theta", type=float, default=0.999)
    parser.add_argument("--grid", type=int, default=121)
    parser.add_argument("--refine-grid", type=int, default=101)
    args = parser.parse_args()
    payload = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
