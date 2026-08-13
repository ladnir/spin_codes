#!/usr/bin/env python3
"""Cover every uniform-subset Hasse edge by unions of witness intervals.

The ordered-chamber probe only credits an edge when one fixed witness is safe
at both endpoints.  That sufficient test throws away valid handoffs between
witnesses.  Along a segment ``p(t) = (1-t) p0 + t p1``, every fixed-witness
exponent is convex: its only nonlinear term is minus the concave multinomial
log normalization.  Consequently the safe set of one witness is an interval.

This diagnostic computes those intervals by derivative and value bisection,
unions them, and checks all 2295 Hasse edges.  It is sample-independent but is
still binary64, not an outward-rounded certificate.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import gammaln, psi

from probe_packet8_adaptive_simplex_atlas import load_witness_cache
from probe_packet8_profile_ordered_chambers import (
    add_full_bijection,
    uniform_subset,
)
from probe_packet8_profile_simplex_landscape import (
    ANCHORS,
    HAMMING_BALL_LOG2_UPPER,
    K,
    M,
    PROFILE_COUNT_LOG2,
    integral_profile,
    load_anchor_file,
)
from probe_packet8_shared_witness_io import load_shared_witness_arrays
from probe_packet8_weight_profile_scalar import CLASSES


LOG2 = math.log(2.0)
LOG_CLASSES = np.log(np.asarray(CLASSES, dtype=np.float64))
LOG_FACTORIAL_M = float(gammaln(M + 1))
ORBIT_CONSTANT = (
    K
    + HAMMING_BALL_LOG2_UPPER
    - LOG_FACTORIAL_M / LOG2
    - float(np.sum(gammaln(np.asarray(CLASSES, dtype=np.float64)))) / LOG2
)


def witness_value(
    t: float,
    left: np.ndarray,
    delta: np.ndarray,
    constant: float,
    charge: np.ndarray,
    target: float,
) -> float:
    point = left + t * delta
    normalization = (
        LOG_FACTORIAL_M
        - float(np.sum(gammaln(point + 1.0)))
        + float(point @ LOG_CLASSES)
    ) / LOG2
    return constant - float(point @ charge) - normalization - target


def normalization_value(point: np.ndarray) -> float:
    return (
        LOG_FACTORIAL_M
        - float(np.sum(gammaln(point + 1.0)))
        + float(point @ LOG_CLASSES)
    ) / LOG2


def witness_derivative(
    t: float,
    left: np.ndarray,
    delta: np.ndarray,
    charge: np.ndarray,
) -> float:
    point = left + t * delta
    normalization_derivative = (
        -float(psi(point + 1.0) @ delta) + float(LOG_CLASSES @ delta)
    ) / LOG2
    return -float(delta @ charge) - normalization_derivative


def safe_interval(
    left: np.ndarray,
    delta: np.ndarray,
    constant: float,
    charge: np.ndarray,
    target: float,
) -> tuple[float, float] | None:
    derivative0 = witness_derivative(0.0, left, delta, charge)
    derivative1 = witness_derivative(1.0, left, delta, charge)
    if derivative0 >= 0.0:
        minimum = 0.0
    elif derivative1 <= 0.0:
        minimum = 1.0
    else:
        minimum = brentq(
            witness_derivative,
            0.0,
            1.0,
            args=(left, delta, charge),
            xtol=2e-14,
            rtol=4 * np.finfo(float).eps,
        )

    value_minimum = witness_value(
        minimum, left, delta, constant, charge, target
    )
    if value_minimum > 0.0:
        return None

    value0 = witness_value(0.0, left, delta, constant, charge, target)
    if value0 <= 0.0:
        lower = 0.0
    else:
        lower = brentq(
            witness_value,
            0.0,
            minimum,
            args=(left, delta, constant, charge, target),
            xtol=2e-14,
            rtol=4 * np.finfo(float).eps,
        )

    value1 = witness_value(1.0, left, delta, constant, charge, target)
    if value1 <= 0.0:
        upper = 1.0
    else:
        upper = brentq(
            witness_value,
            minimum,
            1.0,
            args=(left, delta, constant, charge, target),
            xtol=2e-14,
            rtol=4 * np.finfo(float).eps,
        )
    return lower, upper


def orbit_value(
    t: float,
    left: np.ndarray,
    delta: np.ndarray,
    target: float,
) -> float:
    point = left + t * delta
    return (
        ORBIT_CONSTANT
        + float(np.sum(gammaln(point + np.asarray(CLASSES)))) / LOG2
        - float(point @ LOG_CLASSES) / LOG2
        - target
    )


def orbit_derivative(
    t: float,
    left: np.ndarray,
    delta: np.ndarray,
) -> float:
    point = left + t * delta
    return float((psi(point + np.asarray(CLASSES)) - LOG_CLASSES) @ delta) / LOG2


def orbit_safe_interval(
    left: np.ndarray,
    delta: np.ndarray,
    target: float,
) -> tuple[float, float] | None:
    derivative0 = orbit_derivative(0.0, left, delta)
    derivative1 = orbit_derivative(1.0, left, delta)
    if derivative0 >= 0.0:
        minimum = 0.0
    elif derivative1 <= 0.0:
        minimum = 1.0
    else:
        minimum = brentq(
            orbit_derivative,
            0.0,
            1.0,
            args=(left, delta),
            xtol=2e-14,
            rtol=4 * np.finfo(float).eps,
        )
    if orbit_value(minimum, left, delta, target) > 0.0:
        return None
    if orbit_value(0.0, left, delta, target) <= 0.0:
        lower = 0.0
    else:
        lower = brentq(
            orbit_value,
            0.0,
            minimum,
            args=(left, delta, target),
            xtol=2e-14,
            rtol=4 * np.finfo(float).eps,
        )
    if orbit_value(1.0, left, delta, target) <= 0.0:
        upper = 1.0
    else:
        upper = brentq(
            orbit_value,
            minimum,
            1.0,
            args=(left, delta, target),
            xtol=2e-14,
            rtol=4 * np.finfo(float).eps,
        )
    return lower, upper


def edge_cover(
    left: np.ndarray,
    right: np.ndarray,
    constants: np.ndarray,
    charges: np.ndarray,
    target: float,
) -> tuple[bool, float, float, int]:
    delta = right - left
    endpoint0 = constants - charges @ left - normalization_value(left) - target
    endpoint1 = constants - charges @ right - normalization_value(right) - target
    shared = np.flatnonzero((endpoint0 <= 0.0) & (endpoint1 <= 0.0))
    if len(shared):
        return True, 1.0, 1.0, int(shared[0])
    orbit0 = orbit_value(0.0, left, delta, target)
    orbit1 = orbit_value(1.0, left, delta, target)
    if orbit0 <= 0.0 and orbit1 <= 0.0:
        return True, 1.0, 1.0, -2

    intervals = []
    for witness_index, (constant, charge) in enumerate(zip(constants, charges)):
        interval = safe_interval(left, delta, constant, charge, target)
        if interval is not None:
            intervals.append((interval[0], interval[1], witness_index))
    orbit_interval = orbit_safe_interval(left, delta, target)
    if orbit_interval is not None:
        intervals.append((orbit_interval[0], orbit_interval[1], -2))
    intervals.sort()

    frontier = 0.0
    frontier_witness = -1
    for lower, upper, witness_index in intervals:
        if lower > frontier + 2e-12:
            return False, frontier, lower, frontier_witness
        if upper > frontier:
            frontier = upper
            frontier_witness = witness_index
        if frontier >= 1.0 - 2e-12:
            return True, 1.0, 1.0, frontier_witness
    return False, frontier, 1.0, frontier_witness


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", type=Path, required=True)
    parser.add_argument(
        "--edge-family", choices=("hasse", "comparable"), default="hasse"
    )
    parser.add_argument("--show-uncovered", type=int, default=20)
    parser.add_argument("--upgraded-cache", type=Path)
    parser.add_argument(
        "--extra-shared-report", type=Path, action="append", default=[]
    )
    args = parser.parse_args()
    if args.show_uncovered <= 0:
        raise SystemExit("edge interval cover: invalid uncovered-row limit")

    if args.upgraded_cache:
        names, constants, charges = load_shared_witness_arrays(
            args.upgraded_cache, args.extra_shared_report
        )
        full = add_full_bijection([])[0]
        names.append(full.name)
        constants = np.append(constants, full.constant_log2)
        charges = np.vstack((charges, full.linear_charge))
    else:
        if args.extra_shared_report:
            raise SystemExit(
                "edge interval cover: --extra-shared-report requires --upgraded-cache"
            )
        anchors = ANCHORS + load_anchor_file(args.atlas)
        witnesses = load_witness_cache(
            args.atlas.with_suffix(".witnesses.npz"), anchors
        )
        if witnesses is None:
            raise SystemExit("edge interval cover: matching witness cache missing")
        witnesses = add_full_bijection(witnesses)
        names = [row.name for row in witnesses]
        constants = np.asarray([row.constant_log2 for row in witnesses])
        charges = np.vstack([row.linear_charge for row in witnesses])
    target = -40.0 - PROFILE_COUNT_LOG2

    total = 0
    covered = 0
    uncovered = []
    edge_rows = []
    for subset in range(2, 1 << 9):
        complement = ((1 << 9) - 1) ^ subset
        if args.edge_family == "hasse":
            additions = [
                1 << packet_class
                for packet_class in range(9)
                if complement >> packet_class & 1
            ]
        else:
            additions = []
            addition = complement
            while addition:
                additions.append(addition)
                addition = (addition - 1) & complement
        edge_rows.extend((subset, subset | addition) for addition in additions)

    for subset, larger in edge_rows:
        left = uniform_subset(subset)
        total += 1
        is_covered, gap_left, gap_right, witness_index = edge_cover(
            left,
            uniform_subset(larger),
            constants,
            charges,
            target,
        )
        if is_covered:
            covered += 1
        else:
            midpoint = integral_profile(
                left
                + 0.5
                * (gap_left + gap_right)
                * (uniform_subset(larger) - left)
            )
            uncovered.append(
                (
                    subset,
                    larger,
                    gap_left,
                    gap_right,
                    witness_index,
                    midpoint,
                )
            )

    uncovered.sort(key=lambda row: row[3] - row[2], reverse=True)

    print(f"packet-8 uniform {args.edge_family}-edge interval-union cover probe")
    print(
        f"atlas={args.atlas} witnesses={len(constants)} "
        f"target={target:.12f}"
    )
    print(f"total_edges={total}")
    print("excluded_zero_incident_edges=8")
    print(f"covered_edges={covered}")
    print(f"uncovered_edges={total-covered}")
    for rank, (
        left,
        right,
        gap_left,
        gap_right,
        witness_index,
        midpoint,
    ) in enumerate(uncovered[: args.show_uncovered], 1):
        leader = (
            names[witness_index]
            if witness_index >= 0
            else "inverse_orbit" if witness_index == -2 else "none"
        )
        print(
            f"rank={rank} edge={left:03x}->{right:03x} "
            f"first_gap={gap_left:.12f},{gap_right:.12f} "
            f"gap_width={gap_right-gap_left:.12f} "
            f"frontier_leader={leader} midpoint_profile="
            + ",".join(map(str, midpoint))
        )
    print(f"complete_edge_interval_cover={'PASS' if covered == total else 'NO'}")
    print("status=DIAGNOSTIC_BINARY64_SAMPLE_INDEPENDENT_EDGE_INTERVALS")


if __name__ == "__main__":
    main()
