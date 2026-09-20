#!/usr/bin/env python3
"""Witness optimizer for EBCH32--PF31x33--BA three-band compositions."""

from __future__ import annotations

import math
import os

import numpy as np
from scipy.optimize import minimize

from diagnose_finite_k20_band_compositions_qL import log_multinomial, softmax
from diagnose_finite_k20_renyi_dense import (
    high_precision_transfer_log_moment,
    logistic,
    logit,
)
from diagnose_ebch32_ba_k20_three_band_qL import (
    B,
    BANDS,
    EPOCHS,
    L,
    band_log_moment,
    fast_inner_log,
)


DISTANCE = int(os.environ.get("SPIN_EBCH_DISTANCE", "230686"))


def optimize_composition(
    envelope,
    spectrum: np.ndarray,
    occupation: float,
    active_counts: tuple[float, float, float],
    *,
    cover_mode: bool = False,
    seed_witnesses: tuple[dict[str, object], ...] = (),
) -> dict[str, object]:
    if not math.isclose(sum(active_counts), occupation, rel_tol=0.0, abs_tol=1e-8):
        raise ValueError("active counts must sum to occupation")
    counts = (L - occupation, *active_counts)
    frequencies = (np.asarray(counts, dtype=np.float64) + 0.5) / (L + 2.0)
    log_type_count = log_multinomial(counts)

    def decode(point: np.ndarray):
        probabilities = softmax(point[:4])
        values = np.asarray(
            [min(1.0 - 1e-12, max(1e-12, logistic(float(x)))) for x in point[4:7]]
        )
        surprisal = math.exp(min(4.0, max(-12.0, float(point[7]))))
        order = 1.0 + math.exp(min(10.0, max(-8.0, float(point[8]))))
        return probabilities, values, surprisal, order

    def evaluate(point: np.ndarray, high_precision: bool = False):
        probabilities, values, surprisal, order = decode(point)
        if np.any(probabilities <= 0.0):
            return math.inf, None
        moments = np.asarray(
            [
                band_log_moment(spectrum, band, value, order)
                for band, value in zip(BANDS, values)
            ]
        )
        bit_probability = float(np.dot(probabilities[1:], values))
        inner = (
            high_precision_transfer_log_moment(
                envelope, bit_probability, surprisal, epochs=EPOCHS
            )
            if high_precision
            else fast_inner_log(envelope, bit_probability, surprisal)
        )
        if not math.isfinite(inner):
            return math.inf, None
        log_conditioning = log_type_count + float(
            np.dot(np.asarray(counts), np.log(probabilities))
        )
        raw = inner + DISTANCE * surprisal - B * log_conditioning
        conjugate = (order - 1.0) / order
        value = (
            log_type_count
            + float(np.dot(np.asarray(active_counts), moments)) / order
            + conjugate * min(0.0, raw)
        )
        details = {
            "counts": list(counts),
            "active_counts": list(active_counts),
            "reference_type_probabilities": probabilities.tolist(),
            "reference_value_probabilities": values.tolist(),
            "reference_bit_probability": bit_probability,
            "surprisal": surprisal,
            "holder_order": order,
            "band_log2_renyi_moments": (moments / math.log(2.0)).tolist(),
            "log2_region_conditioning_probability": log_conditioning / math.log(2.0),
            "raw_reference_bad_log2_upper": raw / math.log(2.0),
        }
        return value, details

    initial_values = (0.32, 0.5, 0.68)
    alpha = occupation / L
    surprisals = (
        (max(1e-4, 0.8 * alpha), max(0.02, 2.1 * alpha))
        if cover_mode
        else (max(1e-4, 0.8 * alpha), max(0.02, 2.1 * alpha), 2.1)
    )
    orders = (32.0,) if cover_mode else (2.0, 8.0, 32.0)
    starts = []
    for witness in seed_witnesses:
        probabilities = np.asarray(
            witness["reference_type_probabilities"], dtype=np.float64
        )
        values = np.asarray(
            witness["reference_value_probabilities"], dtype=np.float64
        )
        surprisal = float(witness["surprisal"])
        order = float(witness["holder_order"])
        if (
            probabilities.shape != (4,)
            or values.shape != (3,)
            or np.any(probabilities <= 0.0)
            or np.any(values <= 0.0)
            or np.any(values >= 1.0)
            or surprisal <= 0.0
            or order <= 1.0
        ):
            continue
        probabilities /= float(np.sum(probabilities))
        starts.append(
            (
                order,
                surprisal,
                np.asarray(
                    [
                        *np.log(probabilities),
                        *(logit(float(value)) for value in values),
                        math.log(surprisal),
                        math.log(order - 1.0),
                    ]
                ),
            )
        )
    for order in orders:
        for surprisal in surprisals:
            starts.append(
                (
                    order,
                    surprisal,
                    np.asarray(
                        [
                            *np.log(frequencies),
                            *(logit(value) for value in initial_values),
                            math.log(surprisal),
                            math.log(order - 1.0),
                        ]
                    ),
                )
            )
    best = None
    best_start = None
    for start_order, start_surprisal, start in starts:
        result = minimize(
            lambda point: evaluate(point)[0],
            start,
            method="Nelder-Mead",
            options={
                "maxiter": 600 if cover_mode else 900,
                "xatol": 1e-8,
                "fatol": 1e-7,
            },
        )
        if best is None or float(result.fun) < float(best.fun):
            best = result
            best_start = (start_order, start_surprisal)
    assert best is not None and best_start is not None
    value, details = evaluate(best.x, high_precision=True)
    assert details is not None
    return {
        "occupation": occupation,
        "occupation_density": occupation / L,
        "log2_upper": value / math.log(2.0),
        "margin_bits": -value / math.log(2.0),
        "optimizer_success": bool(best.success),
        "optimizer_start_order": best_start[0],
        "optimizer_start_surprisal": best_start[1],
        **details,
    }


def optimize_tetrahedron_witness(
    envelope,
    spectrum: np.ndarray,
    vertices: np.ndarray,
    *,
    seed_witnesses: tuple[dict[str, object], ...] = (),
) -> dict[str, object]:
    """Search for one witness that minimizes the worst vertex exponent.

    This routine only discovers a binary64 witness.  The cover checker still
    tests the returned fixed witness at every vertex, and the Arb verifier
    independently validates every accepted cell.
    """
    vertices = np.asarray(vertices, dtype=np.float64)
    if vertices.shape != (4, 3):
        raise ValueError("vertices must have shape (4, 3)")
    occupations = np.sum(vertices, axis=1)
    if np.any(vertices < 0.0) or np.any(occupations > L):
        raise ValueError("vertices must be admissible active-count triples")
    center = np.mean(vertices, axis=0)
    center_occupation = float(np.sum(center))
    center_counts = np.concatenate(([L - center_occupation], center))
    frequencies = (center_counts + 0.5) / (L + 2.0)

    def decode(point: np.ndarray):
        probabilities = softmax(point[:4])
        values = np.asarray(
            [min(1.0 - 1e-12, max(1e-12, logistic(float(x)))) for x in point[4:7]]
        )
        surprisal = math.exp(min(4.0, max(-12.0, float(point[7]))))
        order = 1.0 + math.exp(min(10.0, max(-8.0, float(point[8]))))
        return probabilities, values, surprisal, order

    def evaluate(point: np.ndarray):
        probabilities, values, surprisal, order = decode(point)
        if np.any(probabilities <= 0.0):
            return math.inf, None
        moments = np.asarray(
            [
                band_log_moment(spectrum, band, value, order)
                for band, value in zip(BANDS, values)
            ]
        )
        bit_probability = float(np.dot(probabilities[1:], values))
        inner = fast_inner_log(envelope, bit_probability, surprisal)
        if not math.isfinite(inner):
            return math.inf, None
        conjugate = (order - 1.0) / order
        vertex_values = []
        vertex_raws = []
        for active_counts, occupation in zip(vertices, occupations):
            counts = np.concatenate(([L - occupation], active_counts))
            log_type_count = log_multinomial(tuple(float(value) for value in counts))
            log_conditioning = log_type_count + float(
                np.dot(counts, np.log(probabilities))
            )
            raw = inner + DISTANCE * surprisal - B * log_conditioning
            value = (
                log_type_count
                + float(np.dot(active_counts, moments)) / order
                + conjugate * min(0.0, raw)
            )
            vertex_values.append(value)
            vertex_raws.append(raw)
        maximum_value = max(vertex_values)
        maximum_raw = max(vertex_raws)
        # The cancellation step requires raw < 0 at every vertex.  Penalize
        # infeasible candidates while retaining a finite objective from an
        # infeasible starting point.
        objective = maximum_value + 10.0 * max(0.0, maximum_raw)
        return objective, {
            "reference_type_probabilities": probabilities.tolist(),
            "reference_value_probabilities": values.tolist(),
            "reference_bit_probability": bit_probability,
            "surprisal": surprisal,
            "holder_order": order,
            "band_log2_renyi_moments": (moments / math.log(2.0)).tolist(),
            "diagnostic_maximum_vertex_log2": maximum_value / math.log(2.0),
            "diagnostic_maximum_raw_log2": maximum_raw / math.log(2.0),
        }

    starts = []
    for witness in seed_witnesses:
        probabilities = np.asarray(
            witness["reference_type_probabilities"], dtype=np.float64
        )
        values = np.asarray(
            witness["reference_value_probabilities"], dtype=np.float64
        )
        surprisal = float(witness["surprisal"])
        order = float(witness["holder_order"])
        if (
            probabilities.shape != (4,)
            or values.shape != (3,)
            or np.any(probabilities <= 0.0)
            or np.any(values <= 0.0)
            or np.any(values >= 1.0)
            or surprisal <= 0.0
            or order <= 1.0
        ):
            continue
        probabilities /= float(np.sum(probabilities))
        starts.append(
            np.asarray(
                [
                    *np.log(probabilities),
                    *(logit(float(value)) for value in values),
                    math.log(surprisal),
                    math.log(order - 1.0),
                ]
            )
        )
    if not starts:
        alpha = center_occupation / L
        starts.append(
            np.asarray(
                [
                    *np.log(frequencies),
                    *(logit(value) for value in (0.32, 0.5, 0.68)),
                    math.log(max(1e-4, 0.8 * alpha)),
                    math.log(31.0),
                ]
            )
        )

    best = None
    for start in starts[:3]:
        result = minimize(
            lambda point: evaluate(point)[0],
            start,
            method="Nelder-Mead",
            options={"maxiter": 1200, "xatol": 1e-8, "fatol": 1e-7},
        )
        if best is None or float(result.fun) < float(best.fun):
            best = result
    assert best is not None
    _, details = evaluate(best.x)
    assert details is not None
    return {
        "optimizer_success": bool(best.success),
        "optimizer_mode": "tetrahedron-minimax",
        **details,
    }
