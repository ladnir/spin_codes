#!/usr/bin/env python3
"""Sample the full packet-weight simplex against fixed proof witnesses.

Each outer/inner witness consists of:

* valid asymmetric linear Brascamp--Lieb band coefficients and positive outer
  profile fugacities; and
* a fixed robust 65-state inner pole/fugacity witness.

For any exact profile, keeping those parameters fixed gives a rigorous
inequality apart from the scripts' current binary64 arithmetic.  Its profile
dependence is explicit and convex: a linear Cauchy charge minus the exact
multinomial log normalization.  This probe samples the complete simplex,
including low-dimensional faces, to discover anchor profiles needed by a
later adaptive outward-certified cover.  Sampling is diagnostic and cannot
certify coverage by itself.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.special import gammaln, logsumexp

from probe_packet8_repeated_value_state_transfer import (
    EXACT_SLICES,
    SPECTRUM,
    build_split_caps,
    load_ebch128_spectrum,
    load_exact,
    witness,
)
from probe_packet8_outer_profile_enumerator import (
    optimize_outer as optimize_spectrum_outer,
)
from probe_packet8_three_band_linear_bl2 import optimize_outer
from probe_packet8_three_band_profile_enumerator import (
    BINOMIAL_LOG_PROBABILITIES,
    D,
    K,
    N,
    group_log_moments,
)
from probe_packet8_weight_profile_scalar import CLASSES
from probe_packet8_weight_profile_state_transfer import (
    INNER_BLOCKS,
    RobustProfileKernel,
    block_histograms,
)


M = N // 8
PROFILE_COUNT_LOG2 = math.log2(math.comb(M + 8, 8))
HAMMING_BALL_LOG2_UPPER = (
    gammaln(N + 1)
    - gammaln(D + 1)
    - gammaln(N - D + 1)
    + math.log(N - D + 1)
    - math.log(N - 2 * D + 1)
) / math.log(2.0)


@dataclass(frozen=True)
class Anchor:
    name: str
    profile: tuple[int, ...]
    pole: float
    fugacities: tuple[float, ...]
    outer_log_bound: float = 40.0


INTERPOLATION_ANCHORS = (
    Anchor(
        "a10",
        (197325, 819, 2867, 5735, 7168, 5734, 2867, 819, 38810),
        0.7,
        (1, .00867998254, .00963178625, .005209731038, .003042882698,
         .002742188002, .002742188002, .003013344862, .6268151541),
    ),
    Anchor(
        "a20",
        (175514, 1638, 5735, 11469, 14336, 11469, 5734, 1638, 34611),
        0.5,
        (1, .1209987819, .02658537044, .01569923696, .01206414251,
         .01569923696, .01396503608, .01569923696, .8571492948),
    ),
    Anchor(
        "a30",
        (153702, 2458, 8602, 17203, 21504, 17203, 8602, 2457, 30413),
        0.4,
        (1, .06193202725, .05358743245, .0458101931, .03531996808,
         .03775393496, .03746169431, .03463853702, .8831595108),
    ),
    Anchor(
        "a40",
        (131891, 3277, 11469, 22938, 28672, 22937, 11469, 3277, 26214),
        0.35,
        (1, .04568371862, .08516703941, .07675091576, .05955769122,
         .05029298109, .04071170448, .05297870902, .8571492948),
    ),
    Anchor(
        "a50",
        (110080, 4096, 14336, 28672, 35840, 28672, 14336, 4096, 22016),
        0.15,
        (1, .2934843373, .3594729332, .3075276847, .2844416768,
         .323950169, .2771380993, .1645230549, 1.021671305),
    ),
    Anchor(
        "a60",
        (88269, 4915, 17203, 34407, 43008, 34406, 17203, 4915, 17818),
        0.15,
        (1, .52232165, .14745139, .48401486, .65593856, .76255769,
         .26213042, .23924421, .86534489),
    ),
)


def pure_anchor(
    weight: int,
    pole: float,
    floor: float = 0.1,
    outer_log_bound: float = math.log(10.0),
) -> Anchor:
    """Return a finite-neighborhood witness for one nonzero simplex vertex."""
    profile = [0] * 9
    profile[weight] = M
    fugacities = [floor] * 9
    fugacities[weight] = 1.0
    return Anchor(
        f"pure{weight}",
        tuple(profile),
        pole,
        tuple(fugacities),
        outer_log_bound,
    )


# Poles and the 0.1 inactive floor were checked by
# probe_packet8_pure_class_frontier.py.  The positive floor is deliberate:
# zero would cover the vertex but give an infinite charge on every adjacent
# face, while the much smaller initial 1e-12 floor produced unusably narrow
# neighborhoods despite enormous vertex margin.
PURE_ANCHORS = tuple(
    pure_anchor(weight, pole, outer_log_bound=4.0 if weight == 8 else math.log(10.0))
    for weight, pole in enumerate(
        (0.0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
    )
    if weight
)


# First adaptive edge anchor.  The active-class ratio is the normalized output
# of the robust coordinate tuner at the 55/45 weight-7/weight-8 profile; the
# inactive 0.1 floor and bounded outer fugacities keep the witness finite on a
# neighborhood of that edge.
FACE_ANCHORS = (
    Anchor(
        "edge78_55_45",
        (0, 0, 0, 0, 0, 0, 0, 144179, 117965),
        0.1,
        (0.37568, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 1.0, 2.5016),
        math.log(10.0),
    ),
    Anchor(
        "adaptive_0_1_5_7",
        (144190, 47344, 0, 0, 0, 12538, 0, 58072, 0),
        0.2,
        (1.0, 1.10459274, 0.1, 0.1, 0.1, 0.05955769122,
         0.1, 1.10459274, 0.1),
        6.0,
    ),
    Anchor(
        "adaptive_0_1_2_6",
        (129738, 62442, 333, 0, 0, 0, 69631, 0, 0),
        0.2,
        (1.0, 1.10459274, 0.005209731038, 0.1, 0.1, 0.1,
         0.625288432, 0.1, 0.1),
        6.0,
    ),
    Anchor(
        "adaptive_0_6_7_8",
        (80130, 0, 0, 0, 0, 0, 85046, 84178, 12790),
        0.2,
        (1.0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.473910989,
         0.9267177055, 0.5258777108),
        6.0,
    ),
    Anchor(
        "adaptive_0_5_6_7_8",
        (21895, 0, 0, 0, 0, 12888, 92752, 82937, 51672),
        0.2,
        (1.0, 0.1, 0.1, 0.1, 0.1, 0.05955769122,
         0.625288432, 1.028336959, 1.021671305),
        6.0,
    ),
    Anchor(
        "adaptive_1_2_4_6_7_8",
        (0, 65738, 66406, 0, 2172, 0, 35982, 58398, 33448),
        0.2,
        (1.0, 2.154734535, 1.375690028, 0.1, 0.02519664178,
         0.1, 0.8136965599, 2.045501395, 2.154734535),
        6.0,
    ),
)


def midpoint_anchor(left: int, right: int) -> Anchor:
    profile = [0] * 9
    profile[left] = M // 2
    profile[right] = M - profile[left]
    fugacities = [0.1] * 9
    fugacities[left] = 1.0
    fugacities[right] = 1.0
    return Anchor(
        f"mid{left}{right}",
        tuple(profile),
        0.1,
        tuple(fugacities),
        math.log(10.0),
    )


EDGE_MIDPOINT_ANCHORS = tuple(
    midpoint_anchor(left, right)
    for left in range(9)
    for right in range(left + 1, 9)
)


ANCHORS = (
    INTERPOLATION_ANCHORS
    + PURE_ANCHORS
    + FACE_ANCHORS
    + EDGE_MIDPOINT_ANCHORS
)


@dataclass
class FixedWitness:
    name: str
    linear_charge: np.ndarray
    constant_log2: float


def normalization_log2(profiles: np.ndarray) -> np.ndarray:
    return (
        gammaln(M + 1)
        - np.sum(gammaln(profiles + 1), axis=1)
        + profiles @ np.log(np.asarray(CLASSES, dtype=np.float64))
    ) / math.log(2.0)


def stars_bars_log2(profiles: np.ndarray) -> np.ndarray:
    classes = np.asarray(CLASSES, dtype=np.float64)
    return np.sum(
        gammaln(profiles + classes)
        - gammaln(profiles + 1.0)
        - gammaln(classes),
        axis=1,
    ) / math.log(2.0)


def load_anchor_file(path: Path) -> tuple[Anchor, ...]:
    with path.open(encoding="utf-8") as handle:
        rows = json.load(handle)
    return tuple(
        Anchor(
            str(row["name"]),
            tuple(int(value) for value in row["profile"]),
            float(row["pole"]),
            tuple(float(value) for value in row["fugacities"]),
            float(row["outer_log_bound"]),
        )
        for row in rows
    )


def build_fixed_witness(anchor: Anchor, split_caps) -> FixedWitness:
    profile = np.asarray(anchor.profile, dtype=np.float64)
    coefficient_log, outer_result = optimize_outer(
        list(anchor.profile),
        optimize_band_coefficients=True,
        log_fugacity_bound=anchor.outer_log_bound,
    )
    log_t = np.concatenate(([0.0], outer_result.x[:8]))
    outer_mgf_log2 = (
        coefficient_log + float(np.dot(profile, log_t))
    ) / math.log(2.0)

    fugacities = np.asarray(anchor.fugacities, dtype=np.float64)
    kernel = RobustProfileKernel(
        block_histograms(fugacities),
        split_caps,
        anchor.pole,
        float(np.max(fugacities)) ** 8,
    )
    eigenvalue, domination, values, _worst_state = witness(kernel, 64)
    inner_mgf_log2 = (
        math.log2(domination)
        + INNER_BLOCKS * math.log2(eigenvalue)
        + math.log2(float(values[0]))
    )
    constant = outer_mgf_log2 + inner_mgf_log2 - D * math.log2(anchor.pole)
    charge = (log_t + np.log(fugacities)) / math.log(2.0)
    return FixedWitness(anchor.name, charge, constant)


def build_fixed_spectrum_witness(anchor: Anchor, split_caps) -> FixedWitness:
    """Combine the same inner witness with the exact total-spectrum outer bound."""
    profile = np.asarray(anchor.profile, dtype=np.float64)
    coefficient_log, outer_result, _components = optimize_spectrum_outer(
        list(anchor.profile)
    )
    log_t = np.concatenate(([0.0], outer_result.x[:8]))
    outer_mgf_log2 = (
        coefficient_log + float(np.dot(profile, log_t))
    ) / math.log(2.0)

    fugacities = np.asarray(anchor.fugacities, dtype=np.float64)
    kernel = RobustProfileKernel(
        block_histograms(fugacities),
        split_caps,
        anchor.pole,
        float(np.max(fugacities)) ** 8,
    )
    eigenvalue, domination, values, _worst_state = witness(kernel, 64)
    inner_mgf_log2 = (
        math.log2(domination)
        + INNER_BLOCKS * math.log2(eigenvalue)
        + math.log2(float(values[0]))
    )
    constant = outer_mgf_log2 + inner_mgf_log2 - D * math.log2(anchor.pole)
    charge = (log_t + np.log(fugacities)) / math.log(2.0)
    return FixedWitness(anchor.name + "_spectrum", charge, constant)


def build_witnesses(anchors=ANCHORS) -> list[FixedWitness]:
    spectrum = {
        weight: count
        for weight, count in load_ebch128_spectrum(SPECTRUM)
        if weight and count
    }
    split_caps = build_split_caps(spectrum, load_exact(EXACT_SLICES))
    result = []
    for anchor in anchors:
        result.append(build_fixed_witness(anchor, split_caps))
        result.append(build_fixed_spectrum_witness(anchor, split_caps))
    return result


def sample_profiles(rng: np.random.Generator, per_family: int) -> np.ndarray:
    rows = []
    for concentration in (0.03, 0.1, 0.3, 1.0, 3.0):
        rows.append(rng.dirichlet(np.full(9, concentration), per_family))
    for support in range(1, 9):
        for concentration in (0.1, 1.0):
            values = np.zeros((per_family, 9), dtype=np.float64)
            for row in range(per_family):
                active = rng.choice(9, support, replace=False)
                values[row, active] = rng.dirichlet(
                    np.full(support, concentration)
                )
            rows.append(values)
    # Deterministic edges catch narrow boundary failures that random Dirichlet
    # draws can miss.
    edge_rows = []
    for left in range(9):
        for right in range(left + 1, 9):
            for step in range(21):
                value = np.zeros(9)
                value[left] = step / 20
                value[right] = 1.0 - value[left]
                edge_rows.append(value)
    rows.append(np.asarray(edge_rows))
    return M * np.vstack(rows)


def integral_profile(profile: np.ndarray) -> np.ndarray:
    counts = np.floor(profile).astype(np.int64)
    remainder = M - int(np.sum(counts))
    if remainder:
        fractions = profile - counts
        counts[np.argsort(fractions)[-remainder:]] += 1
    return counts


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples-per-family", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260810)
    parser.add_argument("--show-worst", type=int, default=20)
    parser.add_argument("--extra-anchors", type=Path)
    args = parser.parse_args()
    if args.samples_per_family <= 0 or args.show_worst <= 0:
        raise SystemExit("profile landscape: invalid sample count")

    anchors = ANCHORS
    if args.extra_anchors is not None:
        anchors += load_anchor_file(args.extra_anchors)
    witnesses = build_witnesses(anchors)
    profiles = sample_profiles(
        np.random.default_rng(args.seed), args.samples_per_family
    )
    normalizations = normalization_log2(profiles)
    values = np.full(len(profiles), math.inf)
    leaders = np.full(len(profiles), "", dtype=object)

    # Exact full-message bijection moment with the trivial 2^K outer family.
    high_constant = N * math.log2(11) - (N - D) * math.log2(10) + K
    high_values = high_constant - normalizations
    improve = high_values < values
    values[improve] = high_values[improve]
    leaders[improve] = "full_bijection"

    # Continuous diagnostic extension of the exact integer inverse-orbit
    # branch in certify_packet8_profile_high_branches.py.  At integral sample
    # profiles this formula is exactly its logged combinatorial expression.
    orbit_values = (
        K
        + HAMMING_BALL_LOG2_UPPER
        - normalizations
        + stars_bars_log2(profiles)
    )
    improve = orbit_values < values
    values[improve] = orbit_values[improve]
    leaders[improve] = "inverse_orbit"

    for fixed in witnesses:
        candidate = (
            fixed.constant_log2
            - profiles @ fixed.linear_charge
            - normalizations
        )
        improve = candidate < values
        values[improve] = candidate[improve]
        leaders[improve] = fixed.name

    order = np.argsort(values)[::-1]
    print("packet-8 full profile-simplex witness landscape")
    print(f"samples={len(profiles)} seed={args.seed}")
    print(f"fixed_transfer_witnesses={len(witnesses)}")
    print(f"profile_count_log2={PROFILE_COUNT_LOG2:.12f}")
    print(f"per_profile_union_target_log2={-40-PROFILE_COUNT_LOG2:.12f}")
    print(f"worst_sample_log2={values[order[0]]:.12f}")
    print(f"uncovered_above_union_target={int(np.sum(values > -40-PROFILE_COUNT_LOG2))}")
    print(f"positive_samples={int(np.sum(values > 0))}")
    for rank, index in enumerate(order[: args.show_worst], 1):
        proportions = profiles[index] / M
        rounded = integral_profile(profiles[index])
        print(
            f"rank={rank} exponent={values[index]:.9f} "
            f"leader={leaders[index]} profile="
            + ",".join(f"{value:.9f}" for value in proportions)
            + " rounded_counts="
            + ",".join(
                map(
                    str,
                    rounded,
                )
            )
        )
    print("status=DIAGNOSTIC_SAMPLED_SIMPLEX_NOT_A_COVER_CERTIFICATE")


if __name__ == "__main__":
    main()
