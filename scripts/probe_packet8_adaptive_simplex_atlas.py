#!/usr/bin/env python3
"""Greedy cutting-plane atlas for the packet-weight profile simplex.

This diagnostic starts from the committed fixed witnesses, samples boundary
faces and the interior, then repeatedly:

1. selects the currently worst sampled profile;
2. tunes a positive-floor robust inner witness at that profile;
3. selects the broadest outer fugacity bound that leaves anchor margin; and
4. inserts the resulting fixed witness and rescans every sample.

The generated JSON can be passed to
``probe_packet8_profile_simplex_landscape.py --extra-anchors``.  It is a
numerical atlas-discovery artifact, not a finite cover certificate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from probe_packet8_profile_simplex_landscape import (
    ANCHORS,
    HAMMING_BALL_LOG2_UPPER,
    M,
    PROFILE_COUNT_LOG2,
    Anchor,
    FixedWitness,
    build_fixed_witness,
    build_fixed_spectrum_witness,
    build_witnesses,
    integral_profile,
    load_anchor_file,
    normalization_log2,
    sample_profiles,
    stars_bars_log2,
)
from probe_packet8_repeated_value_state_transfer import (
    EXACT_SLICES,
    SPECTRUM,
    build_split_caps,
    load_ebch128_spectrum,
    load_exact,
    witness,
)
from probe_packet8_small_alphabet import log2_multinomial
from probe_packet8_three_band_linear_bl2 import optimize_outer
from probe_packet8_outer_profile_enumerator import (
    optimize_outer as optimize_spectrum_outer,
)
from probe_packet8_three_band_profile_enumerator import D, K, N
from probe_packet8_weight_profile_scalar import CLASSES, optimize_profile
from probe_packet8_weight_profile_scalar import parse_profile
from probe_packet8_weight_profile_state_transfer import (
    INNER_BLOCKS,
    RobustProfileKernel,
    block_histograms,
)


def split_cap_table():
    spectrum = {
        weight: count
        for weight, count in load_ebch128_spectrum(SPECTRUM)
        if weight and count
    }
    return build_split_caps(spectrum, load_exact(EXACT_SLICES))


def profile_normalization(profile: list[int]) -> float:
    return log2_multinomial(profile) + sum(
        count * math.log2(classes)
        for count, classes in zip(profile, CLASSES)
        if count
    )


def robust_probability(
    profile: list[int],
    pole: float,
    fugacities: np.ndarray,
    split_caps,
    iterations: int,
):
    kernel = RobustProfileKernel(
        block_histograms(fugacities),
        split_caps,
        pole,
        float(np.max(fugacities)) ** 8,
    )
    eigenvalue, domination, values, worst_state = witness(kernel, iterations)
    mgf = (
        math.log2(domination)
        + INNER_BLOCKS * math.log2(eigenvalue)
        + math.log2(float(values[0]))
    )
    charge = sum(
        count * math.log2(float(fugacity))
        for count, fugacity in zip(profile, fugacities)
        if count
    )
    probability = (
        mgf
        - charge
        - profile_normalization(profile)
        - D * math.log2(pole)
    )
    return probability, eigenvalue, domination, worst_state


def tune_inner(
    profile: list[int],
    split_caps,
    poles: list[float],
    floor: float,
    passes: int,
    coordinate_iterations: int,
    witness_iterations: int,
):
    active = [weight for weight, count in enumerate(profile) if count]
    reference = active[0]
    variables = [weight for weight in active if weight != reference]

    starts = []
    for pole in poles:
        _scalar, _sequence, _turnoff, _live, initial = optimize_profile(
            profile, pole, passes=4, iterations=24
        )
        reference_value = float(initial[reference])
        if reference_value <= 0.0:
            reference_value = 1.0
        fugacities = np.full(9, floor, dtype=np.float64)
        fugacities[reference] = 1.0
        for weight in variables:
            value = float(initial[weight]) / reference_value
            fugacities[weight] = max(math.exp(-8.0), min(math.exp(3.0), value))
        row = robust_probability(
            profile, pole, fugacities, split_caps, witness_iterations
        )
        starts.append((row[0], pole, fugacities, row))
    _initial_value, pole, fugacities, initial_row = min(
        starts, key=lambda row: row[0]
    )

    point = np.log(fugacities[variables])
    evaluations = len(starts)

    def objective(candidate: np.ndarray):
        nonlocal evaluations
        evaluations += 1
        trial = np.full(9, floor, dtype=np.float64)
        trial[reference] = 1.0
        trial[variables] = np.exp(candidate)
        row = robust_probability(
            profile, pole, trial, split_caps, witness_iterations
        )
        return row[0], trial, row

    best = (initial_row[0], fugacities, initial_row)
    for _pass in range(passes):
        old = point.copy()
        for coordinate in range(len(point)):
            low = -8.0
            high = 3.0
            for _ in range(coordinate_iterations):
                left = (2.0 * low + high) / 3.0
                right = (low + 2.0 * high) / 3.0
                left_point = point.copy()
                right_point = point.copy()
                left_point[coordinate] = left
                right_point[coordinate] = right
                if objective(left_point)[0] <= objective(right_point)[0]:
                    high = right
                else:
                    low = left
            point[coordinate] = (low + high) / 2.0
            candidate = objective(point)
            if candidate[0] < best[0]:
                best = candidate
        if not len(point) or float(np.max(np.abs(point - old))) < 1e-5:
            break
    final = objective(point)
    if best[0] < final[0]:
        final = best
    return pole, final[1], final[0], evaluations


def make_anchor(
    profile: list[int],
    name: str,
    split_caps,
    args,
    inactive_floor: float,
):
    pole, fugacities, inner_probability, evaluations = tune_inner(
        profile,
        split_caps,
        args.poles,
        inactive_floor,
        args.tuning_passes,
        args.coordinate_iterations,
        args.witness_iterations,
    )
    choices = []
    spectrum_log, _spectrum_result, _spectrum_components = optimize_spectrum_outer(
        profile
    )
    spectrum_outer = spectrum_log / math.log(2.0)
    for bound in args.outer_bounds:
        coefficient_log, _result = optimize_outer(
            profile,
            optimize_band_coefficients=True,
            log_fugacity_bound=bound,
        )
        outer_probability = coefficient_log / math.log(2.0)
        combined = min(
            outer_probability + inner_probability,
            spectrum_outer + inner_probability,
        )
        choices.append((bound, combined, min(outer_probability, spectrum_outer)))
        if combined <= args.anchor_target:
            break
    bound, combined, outer_probability = choices[-1]
    anchor = Anchor(
        name,
        tuple(profile),
        pole,
        tuple(float(value) for value in fugacities),
        bound,
    )
    return anchor, combined, outer_probability, inner_probability, evaluations


def write_anchors(path: Path, anchors: list[Anchor]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = [
        {
            "name": anchor.name,
            "profile": list(anchor.profile),
            "pole": anchor.pole,
            "fugacities": list(anchor.fugacities),
            "outer_log_bound": anchor.outer_log_bound,
        }
        for anchor in anchors
    ]
    with path.open("w", encoding="utf-8") as handle:
        json.dump(rows, handle, indent=2)
        handle.write("\n")


def witness_cache_signature(anchors) -> str:
    digest = hashlib.sha256(repr(tuple(anchors)).encode("utf-8"))
    source_names = (
        "probe_packet8_profile_simplex_landscape.py",
        "probe_packet8_three_band_linear_bl2.py",
        "probe_packet8_outer_profile_enumerator.py",
        "probe_packet8_weight_profile_state_transfer.py",
        "probe_packet8_repeated_value_state_transfer.py",
    )
    root = Path(__file__).resolve().parent
    for name in source_names:
        digest.update((root / name).read_bytes())
    return digest.hexdigest()


def load_witness_cache(path: Path, anchors) -> list[FixedWitness] | None:
    if not path.exists():
        return None
    expected = witness_cache_signature(anchors)
    with np.load(path, allow_pickle=False) as cache:
        if str(cache["signature"].item()) != expected:
            return None
        charges = cache["charges"]
        constants = cache["constants"]
        names = cache["names"]
    if charges.shape != (len(constants), 9) or len(names) != len(constants):
        return None
    return [
        FixedWitness(str(name), charge, float(constant))
        for name, charge, constant in zip(names, charges, constants)
    ]


def save_witness_cache(path: Path, anchors, witnesses: list[FixedWitness]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path,
        signature=np.asarray(witness_cache_signature(anchors)),
        names=np.asarray([fixed.name for fixed in witnesses]),
        charges=np.vstack([fixed.linear_charge for fixed in witnesses]),
        constants=np.asarray([fixed.constant_log2 for fixed in witnesses]),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--iterations", type=int, default=10)
    parser.add_argument("--samples-per-family", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260810)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--poles", default="0.1,0.15,0.2,0.3")
    parser.add_argument("--inactive-floor", type=float, default=0.1)
    parser.add_argument(
        "--sharp-companion",
        action="store_true",
        help="also tune a support-local witness with a negligible inactive floor",
    )
    parser.add_argument("--sharp-floor", type=float, default=1e-30)
    parser.add_argument("--outer-bounds", default="6,10,40")
    parser.add_argument("--anchor-target", type=float, default=-1000.0)
    parser.add_argument("--tuning-passes", type=int, default=1)
    parser.add_argument("--coordinate-iterations", type=int, default=8)
    parser.add_argument("--witness-iterations", type=int, default=48)
    parser.add_argument(
        "--forced-profile",
        help="tune this exact nine-count profile instead of selecting the sampled worst row",
    )
    args = parser.parse_args()
    args.poles = [float(value) for value in args.poles.split(",")]
    args.outer_bounds = [float(value) for value in args.outer_bounds.split(",")]
    if args.iterations <= 0 or args.samples_per_family <= 0:
        raise SystemExit("adaptive atlas: invalid iteration or sample count")
    forced_profile = None
    if args.forced_profile is not None:
        if args.iterations != 1:
            raise SystemExit("adaptive atlas: a forced profile requires --iterations 1")
        try:
            forced_profile = parse_profile(args.forced_profile)
        except ValueError as error:
            raise SystemExit(f"adaptive atlas: {error}") from error

    extras = (
        list(load_anchor_file(args.output))
        if args.resume and args.output.exists()
        else []
    )
    split_caps = split_cap_table()
    all_anchors = ANCHORS + tuple(extras)
    cache_path = args.output.with_suffix(".witnesses.npz")
    base_witnesses = load_witness_cache(cache_path, all_anchors)
    cache_status = "HIT"
    if base_witnesses is None:
        cache_status = "MISS"
        base_witnesses = build_witnesses(all_anchors)
        save_witness_cache(cache_path, all_anchors, base_witnesses)
    profiles = sample_profiles(
        np.random.default_rng(args.seed), args.samples_per_family
    )
    normalizations = normalization_log2(profiles)
    values = (
        N * math.log2(11)
        - (N - D) * math.log2(10)
        + K
        - normalizations
    )
    orbit = (
        K
        + HAMMING_BALL_LOG2_UPPER
        - normalizations
        + stars_bars_log2(profiles)
    )
    values = np.minimum(values, orbit)
    for fixed in base_witnesses:
        values = np.minimum(
            values,
            fixed.constant_log2
            - profiles @ fixed.linear_charge
            - normalizations,
        )

    target = -40.0 - PROFILE_COUNT_LOG2
    print("packet-8 greedy adaptive simplex atlas", flush=True)
    print(
        f"samples={len(profiles)} base_witnesses={len(base_witnesses)} "
        f"resumed_anchors={len(extras)} target={target:.9f} "
        f"witness_cache={cache_status}",
        flush=True,
    )
    for iteration in range(len(extras), len(extras) + args.iterations):
        if forced_profile is None:
            worst_index = int(np.argmax(values))
            profile = integral_profile(profiles[worst_index]).tolist()
            before = float(values[worst_index])
        else:
            profile = list(forced_profile)
            point = np.asarray(profile, dtype=np.float64)
            point_normalization = float(normalization_log2(point[None, :])[0])
            before = min(
                fixed.constant_log2
                - float(np.dot(point, fixed.linear_charge))
                - point_normalization
                for fixed in base_witnesses
            )
        anchor, self_value, outer, inner, evaluations = make_anchor(
            profile,
            f"auto_{iteration:03d}_broad" if args.sharp_companion else f"auto_{iteration:03d}",
            split_caps,
            args,
            args.inactive_floor,
        )
        new_anchors = [anchor]
        rows = [
            (
                anchor,
                self_value,
                outer,
                inner,
                evaluations,
                "broad",
            )
        ]
        if args.sharp_companion:
            sharp = make_anchor(
                profile,
                f"auto_{iteration:03d}_sharp",
                split_caps,
                args,
                args.sharp_floor,
            )
            new_anchors.append(sharp[0])
            rows.append((*sharp, "sharp"))

        fixed_rows = []
        for new_anchor in new_anchors:
            fixed_rows.extend(
                (
                    build_fixed_witness(new_anchor, split_caps),
                    build_fixed_spectrum_witness(new_anchor, split_caps),
                )
            )
        base_witnesses.extend(fixed_rows)
        candidate = np.full(len(profiles), math.inf)
        for fixed in fixed_rows:
            candidate = np.minimum(
                candidate,
                fixed.constant_log2
                - profiles @ fixed.linear_charge
                - normalizations,
            )
        improved = int(np.sum(candidate < values))
        values = np.minimum(values, candidate)
        extras.extend(new_anchors)
        write_anchors(args.output, extras)
        save_witness_cache(
            cache_path,
            ANCHORS + tuple(extras),
            base_witnesses,
        )
        uncovered = int(np.sum(values > target))
        positive = int(np.sum(values > 0.0))
        row_summary = ";".join(
            f"{kind}:self={row_self:.6f},outer={row_outer:.6f},"
            f"inner={row_inner:.6f},pole={row_anchor.pole:.4f},"
            f"bound={row_anchor.outer_log_bound:.1f},evals={row_evaluations}"
            for row_anchor, row_self, row_outer, row_inner, row_evaluations, kind in rows
        )
        print(
            f"iteration={iteration} selected_before={before:.9f} "
            f"anchors={row_summary} improved_samples={improved} "
            f"worst_after={float(np.max(values)):.9f} "
            f"uncovered={uncovered} positive={positive} profile="
            + ",".join(map(str, profile)),
            flush=True,
        )
        if uncovered == 0:
            print("sampled_landscape_covered=PASS", flush=True)
            break
    print(f"atlas_path={args.output}", flush=True)
    print("status=DIAGNOSTIC_ADAPTIVE_SAMPLED_ATLAS", flush=True)


if __name__ == "__main__":
    main()
