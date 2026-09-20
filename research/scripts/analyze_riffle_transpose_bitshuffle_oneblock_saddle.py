#!/usr/bin/env python3
"""Sharpen the dominant one-block tail for transpose bit shuffle.

The existing certificate uses ``inf_t E[exp(t W-t d)]`` for ``t<0``.  This
script differentiates the exact one-active-block probability generating
function, adds the one-dimensional lattice saddle prefactor, and validates
the approximation against exact dynamic programs on smaller instances.

The saddle values are diagnostics, not upper bounds.  They decide whether a
more expensive rigorous tilted-tail calculation is worth implementing.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import log_ndtr

from analyze_riffle_striped_random_outer import LOG2, log_two_power_minus_one


DEFAULT_OUTPUT = Path(
    "constructions/riffle_transpose_bitshuffle_randomstepconv/"
    "receipts/oneblock_saddle_g8.json"
)


Jet = tuple[np.ndarray, np.ndarray, np.ndarray]


def jet_identity(size: int) -> Jet:
    return (
        np.eye(size, dtype=np.float64),
        np.zeros((size, size), dtype=np.float64),
        np.zeros((size, size), dtype=np.float64),
    )


def jet_multiply(left: Jet, right: Jet) -> Jet:
    a0, a1, a2 = left
    b0, b1, b2 = right
    return (
        a0 @ b0,
        a1 @ b0 + a0 @ b1,
        a2 @ b0 + 2.0 * (a1 @ b1) + a0 @ b2,
    )


def jet_power(base: Jet, exponent: int) -> Jet:
    result = jet_identity(base[0].shape[0])
    power = base
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = jet_multiply(result, power)
        remaining >>= 1
        if remaining:
            power = jet_multiply(power, power)
    return result


def transition_jets(packet_bits: int, sigma: int, t: float) -> tuple[Jet, Jet]:
    z = math.exp(t)
    logistic = z / (1.0 + z)
    q = ((1.0 + z) * 0.5) ** packet_bits
    log_q_first = packet_bits * logistic
    log_q_second = packet_bits * logistic * (1.0 - logistic)
    q_first = q * log_q_first
    q_second = q * (log_q_second + log_q_first * log_q_first)

    state_zero = math.ldexp(1.0, -sigma)
    state_live = 1.0 - state_zero
    zero = []
    marked = []
    for derivative in (q, q_first, q_second):
        a = state_zero * derivative
        b = state_live * derivative
        zero.append(np.asarray(((0.0, 0.0), (a, b))))
        marked.append(np.asarray(((0.5 * a, 0.5 * b), (a, b))))
    zero[0][0, 0] = 1.0
    marked[0][0, 0] += 0.5
    return tuple(zero), tuple(marked)  # type: ignore[return-value]


def one_row_jet(packet_bits: int, sigma: int, positions: int, t: float) -> Jet:
    zero, marked = transition_jets(packet_bits, sigma, t)
    block = []
    for derivative in range(3):
        matrix = np.zeros((4, 4), dtype=np.float64)
        matrix[:2, :2] = zero[derivative]
        matrix[:2, 2:] = marked[derivative]
        matrix[2:, 2:] = zero[derivative]
        block.append(matrix)
    powered = jet_power(tuple(block), positions)  # type: ignore[arg-type]
    return tuple(value[:2, 2:] / positions for value in powered)  # type: ignore[return-value]


def cumulants(
    *, outer_bits: int, packet_bits: int, sigma: int, positions: int, t: float
) -> tuple[float, float, float]:
    row = one_row_jet(packet_bits, sigma, positions, t)
    total = jet_power(row, outer_bits)
    values = [float(value[0, 0] + value[0, 1]) for value in total]
    if not values[0] > 0.0:
        raise ArithmeticError("one-block generating function underflowed")
    mean = values[1] / values[0]
    variance = values[2] / values[0] - mean * mean
    return math.log(values[0]), mean, variance


def log_normal_cdf_plus_correction(w: float, u: float) -> float:
    """Return log(Phi(w)+phi(w)(1/w-1/u)) stably enough here."""
    log_phi = -0.5 * w * w - 0.5 * math.log(2.0 * math.pi)
    log_cdf = float(log_ndtr(w))
    correction_factor = 1.0 / w - 1.0 / u
    cdf = math.exp(log_cdf)
    correction = math.exp(log_phi) * correction_factor
    estimate = cdf + correction
    if not estimate > 0.0:
        return math.nan
    return math.log(estimate)


def oneblock_saddle(
    *,
    message_bits: int,
    outer_bits: int,
    packet_bits: int,
    sigma: int,
    relative_distance: float,
) -> dict[str, float | int | bool]:
    outer_dimension = outer_bits // 2
    blocks = message_bits // outer_dimension
    if blocks % packet_bits:
        raise ValueError("packet width must divide the outer-block count")
    positions = blocks // packet_bits
    output_bits = outer_bits * positions * packet_bits
    distance = math.floor(relative_distance * output_bits)

    def centered_mean(t: float) -> float:
        return cumulants(
            outer_bits=outer_bits,
            packet_bits=packet_bits,
            sigma=sigma,
            positions=positions,
            t=t,
        )[1] - distance

    lower = -1.0
    while centered_mean(lower) > 0.0:
        lower *= 2.0
    if centered_mean(0.0) <= 0.0:
        raise ArithmeticError("target is not a lower-tail event")
    t = brentq(centered_mean, lower, 0.0, xtol=1e-14, rtol=1e-13)
    log_mgf, tilted_mean, tilted_variance = cumulants(
        outer_bits=outer_bits,
        packet_bits=packet_bits,
        sigma=sigma,
        positions=positions,
        t=t,
    )
    chernoff_log = log_mgf - t * distance
    point_log = chernoff_log - 0.5 * math.log(2.0 * math.pi * tilted_variance)
    lattice_log = point_log - math.log1p(-math.exp(t))

    rate = t * distance - log_mgf
    w = -math.sqrt(2.0 * rate)
    lattice_u = 2.0 * math.sinh(0.5 * t) * math.sqrt(tilted_variance)
    lugannani_rice_log = log_normal_cdf_plus_correction(w, lattice_u)

    conditioning_penalty = -math.log1p(-math.ldexp(1.0, -outer_bits))
    outer_log2 = (
        math.log(blocks)
        + log_two_power_minus_one(outer_dimension)
        + conditioning_penalty
    ) / LOG2

    def margin(log_probability: float) -> float:
        return -(outer_log2 + log_probability / LOG2)

    return {
        "outer_bits": outer_bits,
        "outer_blocks": blocks,
        "packet_bits": packet_bits,
        "packet_positions_per_row": positions,
        "sigma": sigma,
        "output_bits": output_bits,
        "distance": distance,
        "saddle_t": t,
        "saddle_log_surprisal": math.log(-t),
        "tilted_mean": tilted_mean,
        "tilted_variance": tilted_variance,
        "outer_log2": outer_log2,
        "chernoff_inner_log2": chernoff_log / LOG2,
        "chernoff_lambda_bits": margin(chernoff_log),
        "lattice_prefactor_bits": (lattice_log - chernoff_log) / LOG2,
        "lattice_saddle_inner_log2": lattice_log / LOG2,
        "lattice_saddle_lambda_bits": margin(lattice_log),
        "lugannani_rice_valid": math.isfinite(lugannani_rice_log),
        "lugannani_rice_inner_log2": lugannani_rice_log / LOG2,
        "lugannani_rice_lambda_bits": margin(lugannani_rice_log),
    }


def one_row_matrix_batch(
    *, packet_bits: int, sigma: int, positions: int, z: np.ndarray
) -> np.ndarray:
    """Evaluate the exact one-marked-packet row PGF at complex points."""
    q = ((1.0 + z) * 0.5) ** packet_bits
    state_zero = math.ldexp(1.0, -sigma)
    a = state_zero * q
    b = (1.0 - state_zero) * q
    one_minus_b = 1.0 - b
    h = a / one_minus_b
    b_to_p = np.power(b, positions)
    b_to_q = np.power(b, positions - 1)
    geometric = (1.0 - b_to_p) / one_minus_b

    c00 = 0.5 * (1.0 + a)
    c01 = 0.5 * (1.0 - state_zero) * q
    p_minus_h = positions - geometric

    result = np.empty(z.shape + (2, 2), dtype=np.complex128)
    result[..., 0, 0] = (
        positions * c00 + c01 * h * p_minus_h
    ) / positions
    result[..., 0, 1] = c01 * geometric / positions
    result[..., 1, 0] = (
        c00 * h * p_minus_h
        + c01 * h * h * (positions - 2.0 * geometric + positions * b_to_q)
        + a * geometric
        + b * h * (geometric - positions * b_to_q)
    ) / positions
    result[..., 1, 1] = (
        c01 * h * (geometric - positions * b_to_q)
        + positions * b_to_p
    ) / positions
    return result


def matrix_power_batch(matrix: np.ndarray, exponent: int) -> np.ndarray:
    result = np.zeros_like(matrix)
    result[..., 0, 0] = 1.0
    result[..., 1, 1] = 1.0
    power = matrix
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = np.matmul(result, power)
        remaining >>= 1
        if remaining:
            power = np.matmul(power, power)
    return result


def oneblock_mgf_batch(
    *, outer_bits: int, packet_bits: int, sigma: int, positions: int, z: np.ndarray
) -> np.ndarray:
    row = one_row_matrix_batch(
        packet_bits=packet_bits,
        sigma=sigma,
        positions=positions,
        z=z,
    )
    total = matrix_power_batch(row, outer_bits)
    return total[..., 0, 0] + total[..., 0, 1]


def tilted_fourier_tail(
    *,
    outer_bits: int,
    packet_bits: int,
    sigma: int,
    positions: int,
    distance: int,
    saddle_t: float,
    transform_size: int,
    batch_size: int,
) -> dict[str, float | int]:
    """Numerically invert the optimally tilted PGF on a cyclic lattice."""
    if transform_size & (transform_size - 1):
        raise ValueError("Fourier transform size must be a power of two")
    if distance >= transform_size:
        raise ValueError("Fourier transform must exceed the target distance")
    z0 = math.exp(saddle_t)
    log_mgf, _mean, _variance = cumulants(
        outer_bits=outer_bits,
        packet_bits=packet_bits,
        sigma=sigma,
        positions=positions,
        t=saddle_t,
    )
    mgf0 = math.exp(log_mgf)
    characteristic = np.empty(transform_size, dtype=np.complex128)
    for start in range(0, transform_size, batch_size):
        stop = min(transform_size, start + batch_size)
        angles = (2.0 * math.pi / transform_size) * np.arange(start, stop)
        points = z0 * np.exp(1j * angles)
        characteristic[start:stop] = oneblock_mgf_batch(
            outer_bits=outer_bits,
            packet_bits=packet_bits,
            sigma=sigma,
            positions=positions,
            z=points,
        ) / mgf0
    imaginary_symmetry_error = float(
        np.max(
            np.abs(
                characteristic[1:] - np.conj(characteristic[:0:-1])
            )
        )
    )
    tilted_mass = np.fft.fft(characteristic).real / transform_size
    weights = np.exp(saddle_t * (distance - np.arange(distance + 1)))
    weighted_tail = float(np.dot(tilted_mass[: distance + 1], weights))
    if not weighted_tail > 0.0:
        raise ArithmeticError("Fourier tilted tail is not positive")
    chernoff_log = log_mgf - saddle_t * distance
    tail_log = chernoff_log + math.log(weighted_tail)
    return {
        "transform_size": transform_size,
        "batch_size": batch_size,
        "cyclic_alias_distance_above_target": transform_size,
        "characteristic_symmetry_error": imaginary_symmetry_error,
        "tilted_mass_sum": float(np.sum(tilted_mass)),
        "tilted_mass_minimum": float(np.min(tilted_mass)),
        "tilted_negative_l1": float(-np.sum(tilted_mass[tilted_mass < 0.0])),
        "weighted_tilted_tail": weighted_tail,
        "tail_log2": tail_log / LOG2,
        "prefactor_bits": math.log2(weighted_tail),
    }


def transition_polynomials(packet_bits: int, sigma: int) -> tuple[np.ndarray, np.ndarray]:
    coefficients = np.asarray(
        [math.comb(packet_bits, weight) for weight in range(packet_bits + 1)],
        dtype=np.float64,
    ) / (1 << packet_bits)
    state_zero = math.ldexp(1.0, -sigma)
    state_live = 1.0 - state_zero
    zero = np.zeros((packet_bits + 1, 2, 2), dtype=np.float64)
    marked = np.zeros_like(zero)
    zero[:, 1, 0] = state_zero * coefficients
    zero[:, 1, 1] = state_live * coefficients
    zero[0, 0, 0] = 1.0
    marked[:, 0, 0] = 0.5 * state_zero * coefficients
    marked[:, 0, 1] = 0.5 * state_live * coefficients
    marked[0, 0, 0] += 0.5
    marked[:, 1, :] = zero[:, 1, :]
    return zero, marked


def apply_polynomial_transition(distribution: np.ndarray, transition: np.ndarray) -> np.ndarray:
    packet_bits = transition.shape[0] - 1
    result = np.zeros((2, distribution.shape[1] + packet_bits), dtype=np.float64)
    for old_state in range(2):
        for new_state in range(2):
            for weight in range(packet_bits + 1):
                probability = transition[weight, old_state, new_state]
                if probability:
                    result[new_state, weight : weight + distribution.shape[1]] += (
                        probability * distribution[old_state]
                    )
    return result


def exact_small_tail(
    *, outer_bits: int, packet_bits: int, sigma: int, positions: int, distance: int
) -> float:
    zero, marked = transition_polynomials(packet_bits, sigma)
    current = np.asarray(((1.0,), (0.0,)), dtype=np.float64)
    for _row in range(outer_bits):
        unmarked = current
        used = np.zeros_like(current)
        for _position in range(positions):
            next_unmarked = apply_polynomial_transition(unmarked, zero)
            next_used = apply_polynomial_transition(used, zero)
            marked_now = apply_polynomial_transition(unmarked, marked)
            if next_used.shape[1] < marked_now.shape[1]:
                next_used = np.pad(next_used, ((0, 0), (0, marked_now.shape[1] - next_used.shape[1])))
            next_used[:, : marked_now.shape[1]] += marked_now
            unmarked = next_unmarked
            used = next_used
        current = used / positions
    return float(np.sum(current[:, : distance + 1]))


def validation_rows() -> list[dict[str, float | int]]:
    cases = (
        (8, 4, 7, 8),
        (12, 4, 8, 12),
        (16, 8, 9, 12),
    )
    rows = []
    for outer_bits, packet_bits, sigma, positions in cases:
        output_bits = outer_bits * packet_bits * positions
        distance = math.floor(0.09 * output_bits)
        exact = exact_small_tail(
            outer_bits=outer_bits,
            packet_bits=packet_bits,
            sigma=sigma,
            positions=positions,
            distance=distance,
        )
        # Set message_bits so the generic formula obtains the chosen number
        # of packet positions in each row.
        message_bits = outer_bits * positions * packet_bits // 2
        saddle = oneblock_saddle(
            message_bits=message_bits,
            outer_bits=outer_bits,
            packet_bits=packet_bits,
            sigma=sigma,
            relative_distance=distance / output_bits,
        )
        exact_log2 = math.log2(exact)
        rows.append(
            {
                "outer_bits": outer_bits,
                "packet_bits": packet_bits,
                "sigma": sigma,
                "positions": positions,
                "output_bits": output_bits,
                "distance": distance,
                "exact_tail_log2": exact_log2,
                "chernoff_tail_log2": float(saddle["chernoff_inner_log2"]),
                "lattice_saddle_tail_log2": float(saddle["lattice_saddle_inner_log2"]),
                "lattice_error_bits": float(saddle["lattice_saddle_inner_log2"]) - exact_log2,
                "lugannani_rice_tail_log2": float(saddle["lugannani_rice_inner_log2"]),
                "lugannani_rice_error_bits": float(saddle["lugannani_rice_inner_log2"]) - exact_log2,
            }
        )
    return rows


def parse_ints(text: str) -> list[int]:
    return [int(item) for item in text.split(",") if item]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=1 << 20)
    parser.add_argument("--g", type=int, default=8)
    parser.add_argument("--cases", default="256:16,256:17,512:13,512:14,1024:11,1024:12")
    parser.add_argument("--relative-distance", type=float, default=0.09)
    parser.add_argument("--fourier-size", type=int, default=0)
    parser.add_argument("--fourier-batch-size", type=int, default=65536)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    checks = validation_rows()
    rows = []
    for item in args.cases.split(","):
        outer_bits_text, sigma_text = item.split(":")
        row = oneblock_saddle(
            message_bits=args.message_bits,
            outer_bits=int(outer_bits_text),
            packet_bits=args.g,
            sigma=int(sigma_text),
            relative_distance=args.relative_distance,
        )
        if args.fourier_size:
            fourier = tilted_fourier_tail(
                outer_bits=int(row["outer_bits"]),
                packet_bits=args.g,
                sigma=int(row["sigma"]),
                positions=int(row["packet_positions_per_row"]),
                distance=int(row["distance"]),
                saddle_t=float(row["saddle_t"]),
                transform_size=args.fourier_size,
                batch_size=args.fourier_batch_size,
            )
            fourier_lambda = -(
                float(row["outer_log2"]) + float(fourier["tail_log2"])
            )
            fourier["lambda_bits"] = fourier_lambda
            row["tilted_fourier"] = fourier
        rows.append(row)
        print(
            f"B,{row['outer_bits']},sigma,{row['sigma']},"
            f"chernoff,{row['chernoff_lambda_bits']:.6f},"
            f"lattice,{row['lattice_saddle_lambda_bits']:.6f},"
            f"LR,{row['lugannani_rice_lambda_bits']:.6f}"
            + (
                f",fourier,{row['tilted_fourier']['lambda_bits']:.6f}"
                if "tilted_fourier" in row
                else ""
            ),
            flush=True,
        )

    payload = {
        "schema": "riffle-transpose-bitshuffle-oneblock-saddle-v1",
        "evidence_label": "EXACT_GENERATING_FUNCTION_NUMERICAL_SADDLE_DIAGNOSTIC",
        "message_bits": args.message_bits,
        "packet_bits": args.g,
        "relative_distance": args.relative_distance,
        "validation_rows": checks,
        "rows": rows,
        "scope": "The Chernoff column is an upper bound. The lattice and Lugannani-Rice columns are numerical approximations, not certificates.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({"validation_rows": checks}, indent=2))
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
