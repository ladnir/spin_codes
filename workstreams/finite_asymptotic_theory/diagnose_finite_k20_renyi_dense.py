#!/usr/bin/env python3
"""Finite Renyi/conditioning diagnostic for dense B=240 occupations.

The reduction is exact in form.  The optimizer and arithmetic are nearest
binary64, so this program does not emit a certificate.
"""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import sys

import numpy as np
import mpmath as mp
from scipy.optimize import minimize
from scipy.special import logsumexp


WORKSTREAM = Path(__file__).resolve().parent
REPO_ROOT = WORKSTREAM.parents[1]
sys.path.insert(0, str(WORKSTREAM))

from analyze_golay_ba_rm2sub_joint import (  # noqa: E402
    DEFAULT_SELECTION,
    load_envelope,
)


B = 240
L = 8832
N = B * L
EPOCHS = N // 128
DISTANCE = 233_164
SPECTRUM = WORKSTREAM / os.environ.get(
    "SPIN_SPECTRUM_RECEIPT",
    "golay_ba3_rm2sub_finite_B240_q1_outward_k20_d11.json",
)
OUTPUT = WORKSTREAM / "golay_ba3_rm2sub_finite_B240_renyi_dense_diagnostic.json"


def log_choose(total: int, selected: int) -> float:
    if not 0 <= selected <= total:
        return -math.inf
    return (
        math.lgamma(total + 1.0)
        - math.lgamma(selected + 1.0)
        - math.lgamma(total - selected + 1.0)
    )


def log_binomial_point_mass(total: int, selected: int, probability: float) -> float:
    """Return log Pr[Bin(total, probability) = selected], including endpoints."""
    if not 0 <= selected <= total or not 0.0 <= probability <= 1.0:
        return -math.inf
    if probability == 0.0:
        return 0.0 if selected == 0 else -math.inf
    if probability == 1.0:
        return 0.0 if selected == total else -math.inf
    return (
        log_choose(total, selected)
        + selected * math.log(probability)
        + (total - selected) * math.log1p(-probability)
    )


def logistic(value: float) -> float:
    if value >= 0.0:
        inverse = math.exp(-min(value, 745.0))
        return 1.0 / (1.0 + inverse)
    direct = math.exp(max(value, -745.0))
    return direct / (1.0 + direct)


def logit(value: float) -> float:
    value = min(1.0 - 1e-12, max(1e-12, value))
    return math.log(value / (1.0 - value))


def load_spectrum_logs() -> np.ndarray:
    payload = json.loads(SPECTRUM.read_text(encoding="utf-8"))
    result = np.full(B + 1, -math.inf, dtype=np.float64)
    for row in payload["weight_witnesses"]:
        weight = int(row["outer_weight"])
        value = float.fromhex(
            row["conditioned_expected_multiplicity_upper_hex"]
        )
        result[weight] = math.log(value)
    return result


def log_matrix_power_row_sum(matrix: np.ndarray, exponent: int) -> float:
    """Return log(e_0^T matrix^exponent 1) with scale normalization."""
    result = np.eye(matrix.shape[0], dtype=np.float64)
    result_scale = 0.0
    factor = matrix.copy()
    maximum = float(np.max(factor))
    if maximum <= 0.0:
        return -math.inf
    factor /= maximum
    factor_scale = math.log(maximum)
    power = exponent
    while power:
        if power & 1:
            result = result @ factor
            result_scale += factor_scale
            maximum = float(np.max(result))
            if maximum <= 0.0:
                return -math.inf
            result /= maximum
            result_scale += math.log(maximum)
        power >>= 1
        if power:
            factor = factor @ factor
            factor_scale *= 2.0
            maximum = float(np.max(factor))
            if maximum <= 0.0:
                return -math.inf
            factor /= maximum
            factor_scale += math.log(maximum)
    row_sum = float(np.sum(result[0]))
    if row_sum <= 0.0:
        return -math.inf
    return result_scale + math.log(row_sum)


def high_precision_transfer_log_moment(
    envelope,
    actual_bit_probability: float,
    surprisal: float,
    epochs: int = EPOCHS,
) -> float:
    """Evaluate the same three-state transfer with cancellation-safe arithmetic."""
    with mp.workdps(100):
        beta = mp.mpf(actual_bit_probability)
        z = mp.exp(-mp.mpf(surprisal))
        u = 1 - beta + beta * z
        v = beta + (1 - beta) * z
        signed = 1 - beta - beta * z
        classes = [int(value) for value in envelope.classes]
        counts = [int(value) for value in envelope.counts]
        association = [
            [int(envelope.association[i, j]) for j in range(len(classes))]
            for i in range(len(classes))
        ]
        live = [u ** (envelope.step_bits - weight) * v**weight for weight in classes]
        fourier = [
            u ** (envelope.step_bits - weight) * signed**weight
            for weight in classes
        ]
        state_space = mp.mpf(envelope.state_space)
        zero_to_zero = mp.fsum(
            mp.mpf(count) * value for count, value in zip(counts, fourier)
        ) / state_space
        zero_total = u**envelope.step_bits
        zero_to_deterministic = zero_total - zero_to_zero
        paired_total = mp.fsum(
            fourier[i] * mp.mpf(association[i][j]) * live[j]
            for i in range(len(classes))
            for j in range(len(classes))
        ) / state_space
        paired_nonzero = paired_total - zero_to_zero * live[0]
        deterministic_moment = paired_nonzero / zero_to_deterministic
        uniform_live_moment = (
            mp.fsum(mp.mpf(count) * value for count, value in zip(counts, live))
            - live[0]
        ) / mp.mpf(envelope.live_states)
        punctured_live_moment = (
            mp.mpf(envelope.live_states)
            / mp.mpf(envelope.live_states - 1)
            * uniform_live_moment
        )
        transfer = mp.matrix(
            [
                [zero_to_zero, zero_to_deterministic, 0],
                [
                    deterministic_moment / envelope.live_states,
                    0,
                    deterministic_moment,
                ],
                [
                    punctured_live_moment / envelope.live_states,
                    0,
                    punctured_live_moment,
                ],
            ]
        )
        if min(transfer) < 0:
            raise ArithmeticError("high-precision transfer has a negative entry")
        powered = transfer**epochs
        row_sum = mp.fsum(powered[0, column] for column in range(3))
        return float(mp.log(row_sum))


def outer_renyi_log_moment(
    spectrum_logs: np.ndarray, value_probability: float, order: float
) -> float:
    terms = []
    for weight in range(1, B + 1):
        log_a = float(spectrum_logs[weight])
        if not math.isfinite(log_a):
            continue
        log_reference_shell = (
            log_choose(B, weight)
            + weight * math.log(value_probability)
            + (B - weight) * math.log1p(-value_probability)
        )
        terms.append(order * log_a + (1.0 - order) * log_reference_shell)
    return float(logsumexp(np.asarray(terms, dtype=np.float64)))


def reference_bad_log_upper(
    *,
    envelope,
    occupation: int,
    candidate_probability: float,
    value_probability: float,
    surprisal: float,
    distance: int = DISTANCE,
) -> tuple[float, dict[str, float]]:
    actual_bit_probability = candidate_probability * value_probability
    if not 0.0 < actual_bit_probability < 1.0:
        return math.inf, {}
    transfer, details = envelope.transfer(
        candidate_probability=2.0 * actual_bit_probability,
        surprisal=surprisal,
    )
    if actual_bit_probability <= 0.5:
        log_moment = log_matrix_power_row_sum(transfer, EPOCHS)
        arithmetic_method = "binary64"
    else:
        log_moment = high_precision_transfer_log_moment(
            envelope, actual_bit_probability, surprisal
        )
        arithmetic_method = "mpmath-100-digit"
    log_conditioning_probability = log_binomial_point_mass(
        L, occupation, candidate_probability
    )
    raw = (
        log_moment
        + distance * surprisal
        - B * log_conditioning_probability
    )
    return min(0.0, raw), {
        "actual_bit_probability": actual_bit_probability,
        "log_inner_moment": log_moment,
        "log_region_conditioning_probability": log_conditioning_probability,
        "raw_reference_bad_log_upper": raw,
        "transfer_arithmetic_method": arithmetic_method,
        "perron_radius_display": float(
            np.max(np.abs(np.linalg.eigvals(transfer)))
        ),
        **details,
    }


def objective(
    point: np.ndarray,
    *,
    envelope,
    spectrum_logs: np.ndarray,
    occupation: int,
) -> tuple[float, dict[str, float]]:
    candidate_probability = logistic(float(point[0]))
    value_probability = logistic(float(point[1]))
    surprisal = math.exp(min(4.0, max(-16.0, float(point[2]))))
    order = 1.0 + math.exp(min(14.0, max(-14.0, float(point[3]))))
    outer_moment = outer_renyi_log_moment(
        spectrum_logs, value_probability, order
    )
    reference_bad, details = reference_bad_log_upper(
        envelope=envelope,
        occupation=occupation,
        candidate_probability=candidate_probability,
        value_probability=value_probability,
        surprisal=surprisal,
    )
    if not math.isfinite(reference_bad):
        return math.inf, {}
    conjugate_weight = (order - 1.0) / order
    value = (
        log_choose(L, occupation)
        + occupation * outer_moment / order
        + conjugate_weight * reference_bad
    )
    return value, {
        "candidate_probability": candidate_probability,
        "value_probability": value_probability,
        "surprisal": surprisal,
        "holder_order": order,
        "outer_renyi_log_moment": outer_moment,
        "conjugate_weight": conjugate_weight,
        **details,
    }


def optimize_occupation(envelope, spectrum_logs: np.ndarray, occupation: int) -> dict[str, object]:
    alpha = occupation / L
    starts = []
    for order in (2.0, 8.0, 32.0, 128.0):
        for value_probability in (0.4, 0.5, 0.6):
            candidate_probability = min(0.999, max(0.001, alpha))
            surprisal = max(1e-5, 1.6 * alpha * value_probability)
            starts.append(
                np.asarray(
                    [
                        logit(candidate_probability),
                        logit(value_probability),
                        math.log(surprisal),
                        math.log(order - 1.0),
                    ],
                    dtype=np.float64,
                )
            )
    best = None
    for start in starts:
        result = minimize(
            lambda point: objective(
                point,
                envelope=envelope,
                spectrum_logs=spectrum_logs,
                occupation=occupation,
            )[0],
            start,
            method="Nelder-Mead",
            options={"maxiter": 1200, "xatol": 1e-8, "fatol": 1e-8},
        )
        if best is None or float(result.fun) < float(best.fun):
            best = result
    assert best is not None
    value, details = objective(
        best.x,
        envelope=envelope,
        spectrum_logs=spectrum_logs,
        occupation=occupation,
    )
    return {
        "occupation": occupation,
        "occupation_density": alpha,
        "log2_upper": value / math.log(2.0),
        "margin_bits": -value / math.log(2.0),
        "optimizer_success": bool(best.success),
        "optimizer_message": str(best.message),
        **details,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--occupations",
        type=int,
        nargs="+",
        default=(65, 128, 256, 512, 1024, 2048, 3072, 4096, 6144, 8192, 8704, 8832),
    )
    parser.add_argument("--selection", type=Path, default=DEFAULT_SELECTION)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    if any(not 65 <= value <= L for value in args.occupations):
        raise ValueError("occupations must lie in [65,8832]")
    envelope = load_envelope(args.selection, DISTANCE / N)
    spectrum_logs = load_spectrum_logs()
    rows = []
    for occupation in args.occupations:
        row = optimize_occupation(envelope, spectrum_logs, occupation)
        rows.append(row)
        print(
            f"occupation={occupation},margin_bits={row['margin_bits']:.6f},"
            f"order={row['holder_order']:.6f}",
            flush=True,
        )
    payload = {
        "schema": "golay-ba3-rm2sub-finite-b240-renyi-dense-diagnostic-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "parameters": {
            "outer_bits": B,
            "outer_rows": L,
            "output_bits": N,
            "distance": DISTANCE,
            "inner_epochs": EPOCHS,
            "occupations": list(args.occupations),
        },
        "method": {
            "outer": "exact conditioned expected BA spectrum upper",
            "change_of_measure": "Renyi/Hölder against Bernoulli-y row values",
            "region_relaxation": (
                "iid Bernoulli-p candidate positions, divided by the exact "
                "Binomial(8832,p) mass at Q in each of 240 regions"
            ),
            "inner": "finite three-state RM2Sub transfer powered through 16560 epochs",
            "arithmetic": "nearest binary64",
        },
        "rows": rows,
        "all_sampled_below_2^-40": all(
            float(row["margin_bits"]) > 40.0 for row in rows
        ),
        "limitations": [
            "The occupation set is sampled, not a cover of every Q.",
            "The optimization and arithmetic are not outward rounded.",
            "The finite conditioning/Renyi reduction still requires a written proof audit.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
