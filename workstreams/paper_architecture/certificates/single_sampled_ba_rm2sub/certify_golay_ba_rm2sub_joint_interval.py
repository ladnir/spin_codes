#!/usr/bin/env python3
"""Outward interval certificate for the BA-majorant/RM2Sub joint exponent."""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

import numpy as np
from mpmath import iv
from scipy.optimize import minimize

from analyze_golay_ba_rm2sub_joint import (
    joint_objective,
    load_envelope,
    logistic,
    logit,
)
from analyze_rm2sub_dense_occupation import DEFAULT_SELECTION


iv.dps = 100
ZERO = Fraction(0)
ONE = Fraction(1)
HALF = Fraction(1, 2)
DELTA = Fraction(11, 100)


def as_interval(value: Fraction):
    return iv.mpf(value.numerator) / value.denominator


def upper_float(value) -> float:
    return math.nextafter(float(value.b), math.inf)


def point_fraction(value: float, digits: int = 17) -> Fraction:
    return Fraction(format(value, f".{digits}g"))


def fast_witness(envelope, alpha: float, row_density: float) -> dict[str, float]:
    epsilon = 1.0e-9
    starts = [
        (alpha, row_density, max(0.02, 1.6 * alpha * row_density)),
        (alpha, 0.5, max(0.02, 0.8 * alpha)),
    ]
    best = None
    for probability, value_probability, surprisal in starts:
        probability = min(1 - epsilon, max(epsilon, probability))
        value_probability = min(1 - epsilon, max(epsilon, value_probability))
        result = minimize(
            lambda point: joint_objective(
                envelope,
                alpha=alpha,
                row_density=row_density,
                outer_exponent_natural=0.0,
                candidate_probability=logistic(float(point[0])),
                value_probability=logistic(float(point[1])),
                surprisal=math.exp(float(point[2])),
                transfer_name="three-state",
            )[0],
            np.asarray(
                [
                    logit(probability),
                    logit(value_probability),
                    math.log(surprisal),
                ]
            ),
            method="Nelder-Mead",
            options={"maxiter": 300, "xatol": 1.0e-8, "fatol": 1.0e-11},
        )
        if best is None or float(result.fun) < float(best.fun):
            best = result
    assert best is not None
    p = logistic(float(best.x[0]))
    y = logistic(float(best.x[1]))
    surprisal = math.exp(float(best.x[2]))
    return {"candidate_probability": p, "value_probability": y, "z": math.exp(-surprisal)}


@dataclass(frozen=True)
class SpectrumSegment:
    index: int
    lower: Fraction
    upper: Fraction
    slope: Fraction
    intercept: Fraction


def load_segments(path: Path) -> list[SpectrumSegment]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("status") != "proved":
        raise ValueError("the BA majorant receipt is not proved")
    left = []
    for index, record in enumerate(payload["left_segments"]):
        left.append(
            SpectrumSegment(
                index,
                Fraction(record["weight"][0]),
                Fraction(record["weight"][1]),
                Fraction(record["slope"]),
                Fraction(record["intercept"]),
            )
        )
    central_lower = left[-1].upper
    central_upper = ONE - central_lower
    central = SpectrumSegment(
        len(left),
        central_lower,
        central_upper,
        ZERO,
        Fraction(payload["central_constant_upper"]),
    )
    right = []
    for segment in reversed(left):
        right.append(
            SpectrumSegment(
                2 * len(left) - segment.index,
                ONE - segment.upper,
                ONE - segment.lower,
                -segment.slope,
                segment.slope + segment.intercept,
            )
        )
    segments = left + [central] + right
    for previous, following in zip(segments, segments[1:], strict=False):
        if previous.upper != following.lower:
            raise AssertionError("spectrum segments do not form a partition")
    return segments


def interval_three_state_transfer(
    *,
    actual_bit_probability: Fraction,
    z: Fraction,
    classes: list[int],
    counts: list[int],
    association: list[list[int]],
    state_space: int,
    step_bits: int,
) -> list[list[object]]:
    beta = as_interval(actual_bit_probability)
    z_iv = as_interval(z)
    u = 1 - beta + beta * z_iv
    v = beta + (1 - beta) * z_iv
    signed = 1 - beta - beta * z_iv
    live = [u ** (step_bits - weight) * v**weight for weight in classes]
    syndrome = [u ** (step_bits - weight) * signed**weight for weight in classes]
    zero_to_zero = sum(
        count * moment for count, moment in zip(counts, syndrome, strict=True)
    ) / state_space
    zero_total = u**step_bits
    zero_to_deterministic = zero_total - zero_to_zero
    paired_total = iv.mpf(0)
    for row, syndrome_moment in enumerate(syndrome):
        paired_total += syndrome_moment * sum(
            association[row][column] * live[column]
            for column in range(len(classes))
        )
    paired_total /= state_space
    paired_nonzero = paired_total - zero_to_zero * live[0]
    deterministic = paired_nonzero / zero_to_deterministic
    live_states = state_space - 1
    uniform_live = (
        sum(count * moment for count, moment in zip(counts, live, strict=True))
        - live[0]
    ) / live_states
    punctured = as_interval(Fraction(live_states, live_states - 1)) * uniform_live
    return [
        [zero_to_zero, zero_to_deterministic, iv.mpf(0)],
        [deterministic / live_states, iv.mpf(0), deterministic],
        [punctured / live_states, iv.mpf(0), punctured],
    ]


def collatz_upper(transfer: list[list[object]], vector: list[Fraction]) -> Fraction:
    ratios = []
    for row, denominator in zip(transfer, vector, strict=True):
        numerator = sum(
            entry * as_interval(coordinate)
            for entry, coordinate in zip(row, vector, strict=True)
        )
        ratios.append(numerator / as_interval(denominator))
    return point_fraction(
        math.nextafter(max(upper_float(value) for value in ratios), math.inf)
    )


def binary_kl(value: Fraction, probability: Fraction):
    x = as_interval(value)
    p = as_interval(probability)
    result = iv.mpf(0)
    if value > 0:
        result += x * iv.log(x / p)
    if value < 1:
        result += (1 - x) * iv.log((1 - x) / (1 - p))
    return result


def exponent_upper(
    *,
    alpha: Fraction,
    row_density: Fraction,
    segment: SpectrumSegment,
    candidate_probability: Fraction,
    value_probability: Fraction,
    z: Fraction,
    radius_upper: Fraction,
) -> float:
    aa = as_interval(alpha)
    xx = as_interval(row_density)
    outer = aa * (
        as_interval(segment.slope) * xx + as_interval(segment.intercept)
    )
    result = outer
    result += binary_kl(alpha, candidate_probability)
    result += aa * binary_kl(row_density, value_probability)
    result += iv.log(as_interval(radius_upper)) / 128
    result -= as_interval(DELTA) * iv.log(as_interval(z))
    return upper_float(result)


@dataclass(frozen=True)
class JointBox:
    segment: SpectrumSegment
    alpha_lo: Fraction
    alpha_hi: Fraction
    x_lo: Fraction
    x_hi: Fraction

    def vertices(self) -> list[tuple[Fraction, Fraction]]:
        return [
            (alpha, density)
            for alpha in (self.alpha_lo, self.alpha_hi)
            for density in (self.x_lo, self.x_hi)
        ]

    def split(self) -> tuple["JointBox", "JointBox"]:
        alpha_width = self.alpha_hi - self.alpha_lo
        x_width = self.x_hi - self.x_lo
        # The optimizer changes most rapidly near the occupation endpoints.
        # Refine occupation first.  Refine row density only after occupation
        # is local enough that a failed vertex test diagnoses the x-width.
        if alpha_width > Fraction(1, 65536) or x_width == 0:
            midpoint = (self.alpha_lo + self.alpha_hi) * HALF
            return (
                JointBox(self.segment, self.alpha_lo, midpoint, self.x_lo, self.x_hi),
                JointBox(self.segment, midpoint, self.alpha_hi, self.x_lo, self.x_hi),
            )
        midpoint = (self.x_lo + self.x_hi) * HALF
        return (
            JointBox(self.segment, self.alpha_lo, self.alpha_hi, self.x_lo, midpoint),
            JointBox(self.segment, self.alpha_lo, self.alpha_hi, midpoint, self.x_hi),
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--selection", type=Path, default=DEFAULT_SELECTION)
    parser.add_argument(
        "--majorant",
        type=Path,
        default=Path(
            "workstreams/finite_asymptotic_theory/"
            "golay_ba3_concave_majorant.json"
        ),
    )
    parser.add_argument("--maximum-depth", type=int, default=50)
    parser.add_argument("--maximum-boxes", type=int, default=100_000)
    parser.add_argument("--progress-every", type=int, default=100)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "workstreams/finite_asymptotic_theory/"
            "golay_ba3_rm2sub_joint_interval_d11.json"
        ),
    )
    args = parser.parse_args()

    segments = load_segments(args.majorant)
    envelope = load_envelope(args.selection, float(DELTA))
    classes = [int(value) for value in envelope.classes]
    counts = [int(value) for value in envelope.counts]
    association = [[int(value) for value in row] for row in envelope.association]
    state_space = 1 << envelope.state_bits

    stack = [
        (
            JointBox(
                segment,
                Fraction(1, 10_000),
                ONE,
                segment.lower,
                segment.upper,
            ),
            0,
        )
        for segment in reversed(segments)
    ]
    accepted = []
    processed = 0
    deepest = 0
    unresolved = None
    while stack:
        box, depth = stack.pop()
        processed += 1
        deepest = max(deepest, depth)
        alpha_mid = (box.alpha_lo + box.alpha_hi) * HALF
        x_mid = (box.x_lo + box.x_hi) * HALF
        optimized = fast_witness(envelope, float(alpha_mid), float(x_mid))
        p = point_fraction(float(optimized["candidate_probability"]))
        y = point_fraction(float(optimized["value_probability"]))
        z = point_fraction(float(optimized["z"]))
        actual_bit_probability = p * y

        float_transfer, _ = envelope.transfer(
            candidate_probability=2.0 * float(actual_bit_probability),
            surprisal=-math.log(float(z)),
        )
        eigenvalues, eigenvectors = np.linalg.eig(float_transfer)
        index = int(np.argmax(np.abs(eigenvalues)))
        float_vector = np.abs(np.real(eigenvectors[:, index]))
        float_vector /= np.max(float_vector)
        vector = [
            max(Fraction(1, 10**30), point_fraction(float(value)))
            for value in float_vector
        ]
        transfer = interval_three_state_transfer(
            actual_bit_probability=actual_bit_probability,
            z=z,
            classes=classes,
            counts=counts,
            association=association,
            state_space=state_space,
            step_bits=envelope.step_bits,
        )
        radius_upper = collatz_upper(transfer, vector)
        vertex_uppers = [
            exponent_upper(
                alpha=alpha,
                row_density=density,
                segment=box.segment,
                candidate_probability=p,
                value_probability=y,
                z=z,
                radius_upper=radius_upper,
            )
            for alpha, density in box.vertices()
        ]
        if max(vertex_uppers) < 0:
            accepted.append(
                {
                    "segment": box.segment.index,
                    "alpha": [str(box.alpha_lo), str(box.alpha_hi)],
                    "row_density": [str(box.x_lo), str(box.x_hi)],
                    "candidate_probability": str(p),
                    "value_probability": str(y),
                    "z": str(z),
                    "actual_bit_probability": str(actual_bit_probability),
                    "collatz_vector": [str(value) for value in vector],
                    "radius_upper": str(radius_upper),
                    "vertex_exponent_uppers_natural": vertex_uppers,
                }
            )
            continue
        if depth >= args.maximum_depth:
            unresolved = {
                "segment": box.segment.index,
                "alpha": [str(box.alpha_lo), str(box.alpha_hi)],
                "row_density": [str(box.x_lo), str(box.x_hi)],
                "depth": depth,
                "vertex_exponent_uppers_natural": vertex_uppers,
            }
            break
        left, right = box.split()
        stack.append((right, depth + 1))
        stack.append((left, depth + 1))
        if args.progress_every and processed % args.progress_every == 0:
            print(
                json.dumps(
                    {
                        "event": "progress",
                        "processed": processed,
                        "accepted": len(accepted),
                        "queued": len(stack),
                        "deepest": deepest,
                    },
                    sort_keys=True,
                ),
                flush=True,
            )
        if processed >= args.maximum_boxes:
            unresolved = {
                "reason": "maximum box count reached",
                "processed": processed,
                "queued": len(stack),
            }
            break

    accepted.sort(
        key=lambda cell: (
            int(cell["segment"]),
            Fraction(cell["alpha"][0]),
            Fraction(cell["row_density"][0]),
        )
    )
    payload = {
        "schema": "golay-ba3-rm2sub-joint-interval-v1",
        "status": "proved" if unresolved is None else "unresolved",
        "claim": (
            "For delta=0.11, the three-state RM2Sub exponent with the proved "
            "concave BA majorant is negative for every 10^-4<=alpha<=1 "
            "and 0.104<=x<=0.896."
        ),
        "interval_dps": iv.dps,
        "convexity_rule": (
            "For fixed p,y,z, an affine spectrum support, and a Collatz bound, "
            "the exponent is convex in (alpha,q=alpha*x).  Its maximum on "
            "each trapezoid occurs at one of the four listed vertices."
        ),
        "processed_boxes": processed,
        "accepted_boxes": len(accepted),
        "deepest_box": deepest,
        "largest_accepted_upper_natural": max(
            value
            for cell in accepted
            for value in cell["vertex_exponent_uppers_natural"]
        )
        if accepted
        else None,
        "boxes": accepted,
        "unresolved": unresolved,
    }
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({key: value for key, value in payload.items() if key != "boxes"}, indent=2, sort_keys=True))
    if unresolved is not None:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
