#!/usr/bin/env python3
"""OA-15 upper bound for the BCH intermediate code under RandomStepConv.

This is a binary64 diagnostic, not yet a rigorous certificate.  It imports the
same RandomStepConv transfer calculation used by the finite-theory replay,
forms the shell functional at the requested finite target, and majorizes that
functional on the rigorously allowed BCH support by an even Krawtchouk
polynomial of degree at most 14.

For a binary [256,128] code of orthogonal-array strength at least 15,

    sum_w A_w K_j(w) = 0,  1 <= j <= 15.

Consequently, if p(w) majorizes the symmetrized per-codeword failure
coefficient on every allowed weight, then

    sum_w A_w f(w) <= 2^128 * p_0,

where p_0 is the constant Krawtchouk coefficient.  The support used here is
the proven support for the intermediate BCH code: 0, the even weights 38..218,
and 256.
"""

from __future__ import annotations

import argparse
import importlib
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import linprog
from scipy.special import logsumexp


N = 256
K = 128
OA_STRENGTH = 15
ALLOWED_WEIGHTS = [0, *range(38, 219, 2), 256]
EVEN_DEGREES = list(range(0, OA_STRENGTH + 1, 2))


def normalized_krawtchouk(n: int, degree: int, weight: int) -> float:
    """Return K_degree(weight) / binom(n, degree), evaluated exactly first."""

    total = 0
    lo = max(0, degree - (n - weight))
    hi = min(degree, weight)
    for intersection in range(lo, hi + 1):
        term = math.comb(weight, intersection) * math.comb(
            n - weight, degree - intersection
        )
        total += -term if intersection & 1 else term
    return total / math.comb(n, degree)


def load_phase_module(finite_theory_root: Path):
    root = finite_theory_root.resolve()
    if not (root / "small_k_replay" / "evaluate_exact_spectra_q1_phase.py").is_file():
        raise FileNotFoundError(
            "Could not find small_k_replay/evaluate_exact_spectra_q1_phase.py "
            f"under {root}"
        )
    sys.path.insert(0, str(root))
    sys.path.insert(0, str(root / "small_k_replay"))
    return importlib.import_module("small_k_replay.evaluate_exact_spectra_q1_phase")


def randomstepconv_log_coefficients(
    phase,
    exponent: int,
    m_value: int,
    distance_fraction: float,
    tilt_start: float,
    tilt_stop: float,
    tilt_step: float,
) -> tuple[np.ndarray, dict]:
    """Return the best log coefficient for one codeword in each shell."""

    message_bits = 1 << exponent
    if message_bits % K:
        raise ValueError("message length must be divisible by 128")
    outer_rows = message_bits // K
    output_bits = 2 * message_bits
    distance = math.ceil(distance_fraction * output_bits)

    # Match the inclusive decimal grid used by the finite-theory replay.
    count = int(round((tilt_stop - tilt_start) / tilt_step)) + 1
    log_surprisals = [tilt_start + index * tilt_step for index in range(count)]

    best = np.full(N + 1, np.inf, dtype=np.float64)
    best_tilt = np.full(N + 1, np.nan, dtype=np.float64)
    log_rows = math.log(outer_rows)

    for log_surprisal in log_surprisals:
        surprisal = math.exp(log_surprisal)
        z_value = math.exp(-surprisal)
        zero, active = phase.inner.step_matrices(z_value, m_value)
        zero_region, one_region = phase.one_active_region_coefficients(
            zero, active, outer_rows
        )
        coordinates = phase.uniform_coefficients(
            zero_region,
            one_region,
            N,
            N,
        )
        moments = np.logaddexp(coordinates[:, 0, 0], coordinates[:, 0, 1])
        candidate = log_rows + np.minimum(0.0, moments + distance * surprisal)
        improved = candidate < best
        best[improved] = candidate[improved]
        best_tilt[improved] = log_surprisal

    # The union bound ranges over nonzero messages; the zero codeword is absent.
    best[0] = -np.inf
    best_tilt[0] = np.nan
    metadata = {
        "message_exponent": exponent,
        "message_bits": message_bits,
        "outer_rows": outer_rows,
        "output_bits": output_bits,
        "distance_fraction": distance_fraction,
        "distance": distance,
        "m_value": m_value,
        "tilt_grid": {
            "start": tilt_start,
            "stop": tilt_stop,
            "step": tilt_step,
            "count": count,
        },
        "best_tilt_by_weight": [None if math.isnan(x) else float(x) for x in best_tilt],
    }
    return best, metadata


def solve_oa_majorant(log_coefficients: np.ndarray) -> dict:
    """Solve the even degree-14 Krawtchouk majorant LP."""

    # Complement symmetry of the code lets us replace f(w) by
    # h(w)=(f(w)+f(256-w))/2 without changing sum A_w f(w).
    log_h = np.logaddexp(log_coefficients, log_coefficients[::-1]) - math.log(2.0)
    allowed_log_h = np.asarray([log_h[w] for w in ALLOWED_WEIGHTS])
    log_scale = float(np.max(allowed_log_h))
    target = np.exp(allowed_log_h - log_scale)

    basis = np.asarray(
        [
            [normalized_krawtchouk(N, degree, weight) for degree in EVEN_DEGREES]
            for weight in ALLOWED_WEIGHTS
        ],
        dtype=np.float64,
    )

    # OA-15 annihilates every nonconstant Krawtchouk moment in this basis.
    # The code-average of a majorant is therefore just its constant coefficient.
    objective = np.zeros(len(EVEN_DEGREES), dtype=np.float64)
    objective[0] = 1.0
    result = linprog(
        c=objective,
        A_ub=-basis,
        b_ub=-target,
        bounds=[(None, None)] * len(EVEN_DEGREES),
        method="highs",
    )
    if not result.success:
        raise RuntimeError(f"majorant LP failed: {result.message}")

    polynomial = basis @ result.x
    slack = polynomial - target
    max_violation = float(max(0.0, -np.min(slack)))
    constant = float(result.x[0])
    if constant <= 0.0:
        raise RuntimeError(f"nonpositive majorant constant: {constant}")

    log_bound = log_scale + K * math.log(2.0) + math.log(constant)
    active = [
        {
            "weight": weight,
            "target_scaled": float(target[index]),
            "majorant_scaled": float(polynomial[index]),
            "slack_scaled": float(slack[index]),
        }
        for index, weight in enumerate(ALLOWED_WEIGHTS)
        if slack[index] <= 1e-7
    ]
    return {
        "status": result.message,
        "degrees": EVEN_DEGREES,
        "coefficients_scaled_basis": [float(x) for x in result.x],
        "log_scale_natural": log_scale,
        "constant_coefficient_scaled": constant,
        "log2_functional_upper_bound": log_bound / math.log(2.0),
        "q1_margin_bits_lower_bound": -log_bound / math.log(2.0),
        "max_constraint_violation_scaled": max_violation,
        "active_constraints": active,
        "allowed_weights": ALLOWED_WEIGHTS,
    }


def logadd_counts(log_coefficients: np.ndarray, counts: dict[int, int]) -> float:
    terms = [
        math.log(count) + float(log_coefficients[weight])
        for weight, count in counts.items()
        if count > 0 and math.isfinite(float(log_coefficients[weight]))
    ]
    return float(logsumexp(terms)) if terms else -math.inf


def heuristic_comparators(log_coefficients: np.ndarray) -> dict:
    # Random linear [256,128] ensemble expectation, identical to the neutral
    # family-curve comparator.  This deliberately does not impose evenness.
    log_denominator = math.log((1 << N) - 1)
    log_numerator = math.log((1 << K) - 1)
    random_terms = []
    for weight in range(1, N + 1):
        log_count = (
            log_numerator
            + math.log(math.comb(N, weight))
            - log_denominator
        )
        random_terms.append(log_count + float(log_coefficients[weight]))
    random_log_value = float(logsumexp(random_terms))

    # A deliberately simple code-specific null model: distribute all codewords
    # other than 0 and 1^256 according to binomial mass on the allowed interior
    # even shells.  It is a heuristic comparator, not a claimed BCH spectrum.
    interior = list(range(38, 219, 2))
    normalization = sum(math.comb(N, weight) for weight in interior)
    remaining = (1 << K) - 2
    even_terms = [float(log_coefficients[N])]  # the all-one word
    for weight in interior:
        log_count = (
            math.log(remaining)
            + math.log(math.comb(N, weight))
            - math.log(normalization)
        )
        even_terms.append(log_count + float(log_coefficients[weight]))
    even_log_value = float(logsumexp(even_terms))

    # Rigorous currently known minimum multiplicities and their complements.
    known_floor_counts = {
        38: 912640,
        40: 18045952,
        42: 42289664,
        44: 267281408,
        46: 1247821056,
        48: 5525363968,
        50: 1574477312,
        206: 1574477312,
        208: 5525363968,
        210: 1247821056,
        212: 267281408,
        214: 42289664,
        216: 18045952,
        218: 912640,
        256: 1,
    }
    floor_log_value = logadd_counts(log_coefficients, known_floor_counts)

    return {
        "random_linear_ensemble": {
            "log2_functional": random_log_value / math.log(2.0),
            "q1_margin_bits": -random_log_value / math.log(2.0),
        },
        "truncated_even_binomial_null_model": {
            "description": "Binomial mass on even weights 38..218, plus 0 and 1^256",
            "log2_functional": even_log_value / math.log(2.0),
            "q1_margin_bits": -even_log_value / math.log(2.0),
        },
        "known_spectrum_floor_only": {
            "counts": {str(k): v for k, v in known_floor_counts.items()},
            "log2_functional": floor_log_value / math.log(2.0),
            "q1_margin_bits": -floor_log_value / math.log(2.0),
            "note": "This is a lower contribution, not an upper bound.",
        },
    }


def low_shell_stress_test(log_coefficients: np.ndarray) -> dict:
    """Stress the two shells that dominate the code-specific null model."""

    interior = list(range(38, 219, 2))
    normalization = sum(math.comb(N, weight) for weight in interior)
    remaining = (1 << K) - 2
    shell_logs = {}
    for weight in range(38, 129, 2):
        log_count = (
            math.log(remaining)
            + math.log(math.comb(N, weight))
            - math.log(normalization)
        )
        log_pair_coefficient = (
            float(log_coefficients[weight])
            if weight == 128
            else float(
                np.logaddexp(
                    log_coefficients[weight], log_coefficients[N - weight]
                )
            )
        )
        shell_logs[weight] = log_count + log_pair_coefficient

    stressed = [38, 40]
    stressed_log = float(logsumexp([shell_logs[w] for w in stressed]))
    rest_logs = [value for weight, value in shell_logs.items() if weight not in stressed]
    rest_logs.append(float(log_coefficients[N]))
    rest_log = float(logsumexp(rest_logs))
    target_log = -40.0 * math.log(2.0)
    break_multiplier = (
        math.exp(target_log) - math.exp(rest_log)
    ) / math.exp(stressed_log)

    scenarios = []
    for multiplier_log2 in [0, 1, 4, 8, 12, 16]:
        combined = float(
            np.logaddexp(
                rest_log,
                stressed_log + multiplier_log2 * math.log(2.0),
            )
        )
        scenarios.append(
            {
                "a38_a40_multiplier": 1 << multiplier_log2,
                "q1_margin_bits": -combined / math.log(2.0),
            }
        )

    thresholds = []
    for weight in range(38, 61, 2):
        pair_log = float(
            np.logaddexp(
                log_coefficients[weight], log_coefficients[N - weight]
            )
        )
        max_log2_count = -40.0 - pair_log / math.log(2.0)
        null_log2_count = (
            math.log(remaining)
            + math.log(math.comb(N, weight))
            - math.log(normalization)
        ) / math.log(2.0)
        thresholds.append(
            {
                "weight": weight,
                "max_log2_multiplicity_if_alone_for_40_bits": max_log2_count,
                "null_model_log2_multiplicity": null_log2_count,
                "headroom_bits": max_log2_count - null_log2_count,
                "null_model_term_log2": shell_logs[weight] / math.log(2.0),
            }
        )

    return {
        "dominant_shells": stressed,
        "joint_break_multiplier_log2": math.log2(break_multiplier),
        "joint_break_multiplier": break_multiplier,
        "scenarios": scenarios,
        "per_shell_thresholds": thresholds,
    }


def known_rate_half_bch_minimum_shells(phase) -> dict:
    """Compare known exact BCH minimum shells with the same even null model."""

    names = ["ebch8", "ebch32", "ebch128"]
    rows = []
    for name in names:
        constituent = phase.CONSTITUENTS[name]
        spectrum = phase.load_spectrum(constituent)
        n = constituent.block_bits
        k = constituent.dimension
        distance = constituent.minimum_distance
        normalization = sum(
            math.comb(n, weight)
            for weight in range(distance, n - distance + 1)
            if weight % 2 == 0
        )
        null_count = (
            ((1 << k) - 2) * math.comb(n, distance) / normalization
        )
        exact_count = spectrum[distance]
        rows.append(
            {
                "block_bits": n,
                "dimension": k,
                "minimum_distance": distance,
                "exact_minimum_shell": exact_count,
                "truncated_even_binomial_minimum_shell": null_count,
                "exact_to_null_ratio": exact_count / null_count,
                "exact_to_null_log2": math.log2(exact_count / null_count),
            }
        )
    return {
        "description": (
            "Exact minimum-shell multiplicity divided by an even binomial model "
            "conditioned on the same minimum distance and total code size"
        ),
        "cases": rows,
        "largest_observed_ratio": max(row["exact_to_null_ratio"] for row in rows),
    }


def proofish_low_prefix_reduction(
    log_coefficients: np.ndarray, johnson_certificate_paths: dict[int, Path]
) -> dict:
    """Reduce a certificate to seven BCH-specific low-weight shell bounds."""

    johnson: dict[int, dict] = {}
    for weight, path in johnson_certificate_paths.items():
        certificate = json.loads(path.read_text(encoding="utf-8"))
        if (
            certificate.get("n") != N
            or certificate.get("weight") != weight
            or certificate.get("minimum_distance") != 38
            or certificate.get("qsopt_status") != "OPTIMAL"
        ):
            raise ValueError(
                f"unexpected or incomplete weight-{weight} Johnson certificate"
            )
        johnson[weight] = certificate

    # If two weight-w supports shared a (w-18)-subset, their distance would be
    # at most 36.  Hence each such subset occurs in at most one support.
    packing_terms = []
    for weight in range(52, 129, 2):
        subset_size = weight - 18
        packing_upper = math.comb(N, subset_size) // math.comb(
            weight, subset_size
        )
        upper = min(1 << K, packing_upper)
        source = "constant-weight packing"
        if weight in johnson:
            upper = min(upper, int(johnson[weight]["integer_shell_upper"]))
            source = "exact rational Johnson-scheme Delsarte LP"
        pair_log = (
            float(log_coefficients[weight])
            if weight == 128
            else float(
                np.logaddexp(
                    log_coefficients[weight], log_coefficients[N - weight]
                )
            )
        )
        term_log = math.log(upper) + pair_log
        packing_terms.append(
            {
                "weight": weight,
                "subset_size": subset_size,
                "multiplicity_upper_bound": upper,
                "multiplicity_upper_bound_log2": math.log2(upper),
                "bound_source": source,
                "functional_term_upper_bound_log2": term_log / math.log(2.0),
                "term_log_natural": term_log,
            }
        )

    high_log = float(
        logsumexp([item["term_log_natural"] for item in packing_terms])
    )

    interior = list(range(38, 219, 2))
    normalization = sum(math.comb(N, weight) for weight in interior)
    remaining = (1 << K) - 2
    low_null_logs = []
    low_null_count_logs = {}
    for weight in range(38, 52, 2):
        log_count = (
            math.log(remaining)
            + math.log(math.comb(N, weight))
            - math.log(normalization)
        )
        low_null_count_logs[weight] = log_count
        pair_log = float(
            np.logaddexp(
                log_coefficients[weight], log_coefficients[N - weight]
            )
        )
        low_null_logs.append(log_count + pair_log)
    low_null_log = float(logsumexp(low_null_logs))
    target_log = -40.0 * math.log(2.0)
    residual = math.exp(target_log) - math.exp(high_log)
    joint_multiplier = residual / math.exp(low_null_log)
    common_multiplier_caps = [
        {
            "weight": weight,
            "null_model_log2_multiplicity": log_count / math.log(2.0),
            "sufficient_log2_multiplicity_cap": (
                log_count / math.log(2.0) + math.log2(joint_multiplier)
            ),
        }
        for weight, log_count in low_null_count_logs.items()
    ]

    return {
        "bch_specific_shells_required": list(range(38, 52, 2)),
        "generic_packing_shells": packing_terms,
        "johnson_certificates": {
            str(weight): {
                "path": str(johnson_certificate_paths[weight].resolve()),
                "integer_shell_upper": certificate["integer_shell_upper"],
                "exact_shell_upper": certificate["exact_shell_upper"],
                "improvement_over_packing_bits": certificate[
                    "improvement_over_simple_packing_bits"
                ],
            }
            for weight, certificate in johnson.items()
        },
        "weights_52_and_above_margin_bits": -high_log / math.log(2.0),
        "low_prefix_null_margin_bits": -low_null_log / math.log(2.0),
        "joint_low_prefix_multiplier_before_40_bit_failure": joint_multiplier,
        "joint_low_prefix_multiplier_log2": math.log2(joint_multiplier),
        "sufficient_common_multiplier_shell_caps": common_multiplier_caps,
        "scope": (
            "The combinatorial bounds are exact.  The reported margins still "
            "use binary64 RandomStepConv transfer coefficients."
        ),
    }


def solve_bch_sandwich_relaxation(
    log_coefficients: np.ndarray, model_path: Path, time_limit: float
) -> dict:
    """Maximize the failure functional over the coupled BCH LP relaxation."""

    from solve_float_diagnostic import (  # imported from this script's directory
        build_relaxation,
        use_krawtchouk_oa_basis,
    )

    model = json.loads(model_path.read_text(encoding="utf-8"))
    forced_zero_variables = {
        next(iter(row["coeffs"]))
        for row in model["constraints"]
        if "_support_" in row["name"]
        and row["sense"] == "eq"
        and int(row["rhs"]) == 0
        and len(row["coeffs"]) == 1
    }
    # Bounds already impose nonnegativity.  Dropping duplicate explicit rows
    # makes the diagnostic substantially easier for HiGHS to presolve.
    model["constraints"] = [
        row
        for row in model["constraints"]
        if not row["name"].startswith(("q_nonneg_", "h_nonneg_"))
    ]
    model = use_krawtchouk_oa_basis(model)
    variables, scales, a_eq, b_eq, a_ub, b_ub = build_relaxation(model)

    half_log_coefficient: dict[int, float] = {}
    for weight in range(0, 129, 2):
        if weight == 128:
            half_log_coefficient[weight] = float(log_coefficients[weight])
        else:
            half_log_coefficient[weight] = float(
                np.logaddexp(log_coefficients[weight], log_coefficients[N - weight])
            )

    log_objective = np.full(len(variables), -np.inf, dtype=np.float64)
    for index, name in enumerate(variables):
        if name in forced_zero_variables:
            continue
        prefix, weight_text = name.split("_", 1)
        multiplier = 1 if prefix == "q" else 31
        log_objective[index] = (
            half_log_coefficient[int(weight_text)]
            + math.log(multiplier)
            + math.log(float(scales[index]))
        )
    objective_log_scale = float(np.max(log_objective))
    objective = -np.exp(log_objective - objective_log_scale)

    # Explicit mass-derived upper bounds prevent HiGHS from mistaking tiny
    # coefficients in the scaled t=0 rows for exact zeros.  They are redundant
    # over the reals and therefore do not strengthen the mathematical model.
    variable_bounds = []
    coset_mass = 1 << 123
    for index, name in enumerate(variables):
        prefix, weight_text = name.split("_", 1)
        weight = int(weight_text)
        if name in forced_zero_variables or name == "h_0":
            variable_bounds.append((0.0, 0.0))
        elif name == "q_0":
            fixed = 1.0 / float(scales[index])
            variable_bounds.append((fixed, fixed))
        else:
            symmetry_factor = 1 if weight == 128 else 2
            physical_upper = coset_mass / symmetry_factor
            variable_bounds.append((0.0, physical_upper / float(scales[index])))

    result = linprog(
        objective,
        A_ub=a_ub,
        b_ub=b_ub,
        A_eq=a_eq,
        b_eq=b_eq,
        bounds=variable_bounds,
        method="highs-ds",
        options={"presolve": True, "time_limit": time_limit},
    )
    if result.status == 3:
        # The scaled model is mathematically bounded.  A presolve-unbounded
        # result can arise when a tiny tail coefficient is dropped relative to
        # a central-shell coefficient; retry the untouched model.
        result = linprog(
            objective,
            A_ub=a_ub,
            b_ub=b_ub,
            A_eq=a_eq,
            b_eq=b_eq,
            bounds=variable_bounds,
            method="highs-ipm",
            options={"presolve": False, "time_limit": time_limit},
        )
    payload = {
        "classification": "scaled binary64 LP relaxation; not a rigorous certificate",
        "model": str(model_path.resolve()),
        "solver_status": int(result.status),
        "solver_message": result.message,
        "success": bool(result.success),
        "objective_log_scale_natural": objective_log_scale,
    }
    if not result.success:
        return payload

    normalized_value = float(-result.fun)
    if normalized_value <= 0.0:
        payload["success"] = False
        payload["solver_message"] = "solver returned a nonpositive maximum"
        return payload
    log_value = objective_log_scale + math.log(normalized_value)
    payload.update(
        {
            "log2_functional_upper_bound": log_value / math.log(2.0),
            "q1_margin_bits_lower_bound": -log_value / math.log(2.0),
            "max_normalized_equality_residual": float(
                np.max(np.abs(result.eqlin.residual))
            ),
            "min_normalized_inequality_slack": float(
                np.min(result.ineqlin.residual)
            ),
        }
    )
    contributions = -objective * result.x
    largest = np.argsort(contributions)[::-1][:12]
    payload["largest_objective_terms"] = [
        {
            "variable": variables[index],
            "normalized_contribution": float(contributions[index]),
            "physical_value": float(result.x[index] * scales[index]),
        }
        for index in largest
        if contributions[index] > 0.0
    ]
    return payload


def solve_normalized_sandwich(
    log_coefficients: np.ndarray, model_path: Path
) -> dict:
    """Stable LP using the exact primal and dual BCH sandwich envelopes."""

    model = json.loads(model_path.read_text(encoding="utf-8"))
    weights = list(range(0, 129, 2))
    count = len(weights)
    weight_index = {weight: index for index, weight in enumerate(weights)}
    mass = 1 << 123

    # X variables are full symmetric-shell probabilities for Q and one nonzero
    # P/Q coset H.  Thus their entries are O(1), and each block sums to one.
    equality_rows = []
    equality_rhs = []
    for block in range(2):
        row = np.zeros(2 * count)
        row[block * count : (block + 1) * count] = 1.0
        equality_rows.append(row)
        equality_rhs.append(1.0)
        for degree in EVEN_DEGREES[1:]:
            row = np.zeros(2 * count)
            row[block * count : (block + 1) * count] = [
                normalized_krawtchouk(N, degree, weight) for weight in weights
            ]
            equality_rows.append(row)
            equality_rhs.append(0.0)

    bounds = [(0.0, 1.0)] * (2 * count)
    bounds = list(bounds)
    # Rigorous support and the unique zero/all-one pair.
    for weight in weights:
        index = weight_index[weight]
        if 0 < weight < 40:
            bounds[index] = (0.0, 0.0)  # Q has minimum distance 40
        if weight < 38:
            bounds[count + index] = (0.0, 0.0)  # nonzero cosets start at 38
    q0_probability = 2.0 / mass
    bounds[weight_index[0]] = (q0_probability, q0_probability)
    bounds[count + weight_index[0]] = (0.0, 0.0)

    inequality_rows = []
    inequality_rhs = []
    for constraint in model["constraints"]:
        name = constraint["name"]
        if name.startswith("Q_ge_L_"):
            weight = int(name.rsplit("_", 1)[1])
            symmetry = 1 if weight == 128 else 2
            lower = symmetry * int(constraint["rhs"]) / mass
            index = weight_index[weight]
            lo, hi = bounds[index]
            bounds[index] = (max(lo, lower), hi)
        elif name.startswith("P_le_U_"):
            weight = int(name.rsplit("_", 1)[1])
            symmetry = 1 if weight == 128 else 2
            upper = symmetry * int(constraint["rhs"]) / mass
            row = np.zeros(2 * count)
            index = weight_index[weight]
            row[index] = 1.0
            row[count + index] = 255.0
            inequality_rows.append(row)
            inequality_rhs.append(upper)

    # Add the dual sandwich and dual-code inclusion constraints in normalized
    # Krawtchouk coordinates.  This is algebraically the same information as
    # the enormous integer rows in the exact export, but its coefficients stay
    # in [-1,1].
    for degree in weights:
        kraw = np.asarray(
            [normalized_krawtchouk(N, degree, weight) for weight in weights]
        )

        bp_row = next(
            row
            for row in model["constraints"]
            if row["name"] == f"Bp_ge_Udual_{degree}"
        )
        bp_lower = int(bp_row["rhs"]) / (
            float(1 << 131) * float(math.comb(N, degree))
        )
        row = np.zeros(2 * count)
        row[:count] = -kraw / 256.0
        row[count:] = -255.0 * kraw / 256.0
        inequality_rows.append(row)
        inequality_rhs.append(-bp_lower)

        bq_row = next(
            row
            for row in model["constraints"]
            if row["name"] == f"Bq_le_Ldual_{degree}"
        )
        bq_upper = int(bq_row["rhs"]) / (
            float(1 << 123) * float(math.comb(N, degree))
        )
        row = np.zeros(2 * count)
        row[:count] = kraw
        inequality_rows.append(row)
        inequality_rhs.append(bq_upper)

        # P^perp is a subcode of Q^perp, so B_Q(degree)>=B_P(degree).
        row = np.zeros(2 * count)
        row[:count] = -kraw
        row[count:] = kraw
        inequality_rows.append(row)
        inequality_rhs.append(0.0)

    # Add the rigorous radius-five affine-orbit floors.
    h38 = count + weight_index[38]
    lo, hi = bounds[h38]
    bounds[h38] = (max(lo, 2.0 * 29440.0 / mass), hi)
    q40 = weight_index[40]
    lo, hi = bounds[q40]
    bounds[q40] = (max(lo, 2.0 * 12990720.0 / mass), hi)
    h40 = count + weight_index[40]
    lo, hi = bounds[h40]
    bounds[h40] = (max(lo, 2.0 * 163072.0 / mass), hi)
    for weight, q_floor, h_floor in (
        (42, 9204480.0, 1067264.0),
        (44, 81991680.0, 5977088.0),
        (46, 381757440.0, 27937536.0),
        (48, 1980725760.0, 114343168.0),
        (50, 679042560.0, 28884992.0),
    ):
        q_index = weight_index[weight]
        lo, hi = bounds[q_index]
        bounds[q_index] = (max(lo, 2.0 * q_floor / mass), hi)
        h_index = count + weight_index[weight]
        lo, hi = bounds[h_index]
        bounds[h_index] = (max(lo, 2.0 * h_floor / mass), hi)

    log_h = np.logaddexp(log_coefficients, log_coefficients[::-1]) - math.log(2.0)
    half_log_h = np.asarray([log_h[weight] for weight in weights])
    log_scale = float(np.max(half_log_h))
    scaled_h = np.exp(half_log_h - log_scale)
    objective = -np.concatenate((scaled_h, 31.0 * scaled_h))

    result = linprog(
        objective,
        A_ub=np.asarray(inequality_rows),
        b_ub=np.asarray(inequality_rhs),
        A_eq=np.asarray(equality_rows),
        b_eq=np.asarray(equality_rhs),
        bounds=bounds,
        method="highs",
    )
    payload = {
        "classification": (
            "normalized binary64 LP relaxation using OA-15, support, Q>=L, "
            "P<=U, the dual sandwich, dual inclusion, and affine-orbit floors; "
            "not a rigorous certificate"
        ),
        "success": bool(result.success),
        "solver_status": int(result.status),
        "solver_message": result.message,
    }
    if not result.success:
        return payload

    normalized_value = float(-result.fun)
    log_value = 123 * math.log(2.0) + log_scale + math.log(normalized_value)
    payload.update(
        {
            "log2_functional_upper_bound": log_value / math.log(2.0),
            "q1_margin_bits_lower_bound": -log_value / math.log(2.0),
            "max_equality_residual": float(np.max(np.abs(result.eqlin.residual))),
            "min_inequality_slack": float(np.min(result.ineqlin.residual)),
        }
    )
    contributions = -objective * result.x
    labels = [f"q_{weight}" for weight in weights] + [
        f"h_{weight}" for weight in weights
    ]
    largest = np.argsort(contributions)[::-1][:12]
    payload["largest_objective_terms"] = [
        {
            "variable": labels[index],
            "normalized_shell_probability": float(result.x[index]),
            "scaled_objective_contribution": float(contributions[index]),
        }
        for index in largest
        if contributions[index] > 0.0
    ]
    return payload


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--finite-theory-root",
        type=Path,
        default=Path(
            r"C:\Users\peter\.codex\worktrees\3061\permute_conv\workstreams\finite_asymptotic_theory"
        ),
    )
    parser.add_argument("--exponent", type=int, default=20)
    parser.add_argument("--m-value", type=int, default=22)
    parser.add_argument("--distance-fraction", type=float, default=0.10)
    parser.add_argument("--tilt-start", type=float, default=-14.0)
    parser.add_argument("--tilt-stop", type=float, default=4.0)
    parser.add_argument("--tilt-step", type=float, default=0.2)
    parser.add_argument(
        "--exact-model",
        type=Path,
        default=Path(__file__).resolve().parents[1]
        / "generated"
        / "coupled_lp_exact.json",
    )
    parser.add_argument(
        "--johnson-certificate",
        type=Path,
        default=Path(__file__).resolve().parents[1]
        / "generated"
        / "johnson_n256_w52_d38.json",
    )
    parser.add_argument(
        "--johnson-certificate-54",
        type=Path,
        default=Path(__file__).resolve().parents[1]
        / "generated"
        / "johnson_n256_w54_d38.json",
    )
    parser.add_argument(
        "--lp-time-limit",
        type=float,
        default=0.0,
        help="seconds for the legacy ill-conditioned coupled LP; zero skips it",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parents[1]
        / "generated"
        / "random_inner_oa15_majorant_diagnostic.json",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    phase = load_phase_module(args.finite_theory_root)
    log_coefficients, metadata = randomstepconv_log_coefficients(
        phase=phase,
        exponent=args.exponent,
        m_value=args.m_value,
        distance_fraction=args.distance_fraction,
        tilt_start=args.tilt_start,
        tilt_stop=args.tilt_stop,
        tilt_step=args.tilt_step,
    )
    majorant = solve_oa_majorant(log_coefficients)
    comparators = heuristic_comparators(log_coefficients)
    stress = low_shell_stress_test(log_coefficients)
    family_minimum_shells = known_rate_half_bch_minimum_shells(phase)
    prefix_reduction = proofish_low_prefix_reduction(
        log_coefficients,
        {
            52: args.johnson_certificate,
            54: args.johnson_certificate_54,
        },
    )
    sandwich = (
        solve_bch_sandwich_relaxation(
            log_coefficients, args.exact_model, args.lp_time_limit
        )
        if args.lp_time_limit > 0.0
        else {
            "classification": "legacy scaled binary64 coupled LP",
            "success": False,
            "solver_message": "skipped (pass --lp-time-limit SECONDS to run)",
        }
    )
    normalized_sandwich = solve_normalized_sandwich(
        log_coefficients, args.exact_model
    )

    receipt = {
        "classification": "binary64 diagnostic; not a rigorous certificate",
        "code_facts_used": {
            "length": N,
            "dimension": K,
            "minimum_distance": 38,
            "even": True,
            "contains_all_one": True,
            "complement_symmetric_spectrum": True,
            "orthogonal_array_strength_at_least": OA_STRENGTH,
        },
        "randomstepconv": metadata,
        "oa_majorant": majorant,
        "bch_sandwich_relaxation": sandwich,
        "normalized_bch_sandwich_relaxation": normalized_sandwich,
        "comparators": comparators,
        "low_shell_stress_test": stress,
        "known_rate_half_bch_minimum_shells": family_minimum_shells,
        "proofish_low_prefix_reduction": prefix_reduction,
        "log_shell_coefficients_natural": [
            None if not math.isfinite(float(value)) else float(value)
            for value in log_coefficients
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")

    print(f"wrote {args.output.resolve()}")
    print(
        "OA-15 majorant Q1 margin lower bound (diagnostic): "
        f"{majorant['q1_margin_bits_lower_bound']:.6f} bits"
    )
    print(
        "random-linear expected Q1 margin: "
        f"{comparators['random_linear_ensemble']['q1_margin_bits']:.6f} bits"
    )
    print(
        "truncated-even-binomial Q1 margin: "
        f"{comparators['truncated_even_binomial_null_model']['q1_margin_bits']:.6f} bits"
    )
    print(
        "known spectrum floor contribution margin: "
        f"{comparators['known_spectrum_floor_only']['q1_margin_bits']:.6f} bits"
    )
    if sandwich["success"]:
        print(
            "BCH-sandwich LP margin lower bound (diagnostic): "
            f"{sandwich['q1_margin_bits_lower_bound']:.6f} bits"
        )
    else:
        print("BCH-sandwich LP did not solve: " + sandwich["solver_message"])
    if normalized_sandwich["success"]:
        print(
            "normalized BCH-sandwich margin lower bound (diagnostic): "
            f"{normalized_sandwich['q1_margin_bits_lower_bound']:.6f} bits"
        )
    else:
        print(
            "normalized BCH-sandwich LP did not solve: "
            + normalized_sandwich["solver_message"]
        )
    print(
        "majorant active weights: "
        + ", ".join(str(item["weight"]) for item in majorant["active_constraints"])
    )
    print(
        "joint A38/A40 multiplier needed to erase 40-bit margin: "
        f"2^{stress['joint_break_multiplier_log2']:.3f}"
    )
    print(
        "largest exact/null minimum-shell ratio at B=8,32,128: "
        f"{family_minimum_shells['largest_observed_ratio']:.3f}x"
    )
    print(
        "generic bound for weights 52 and above: "
        f"{prefix_reduction['weights_52_and_above_margin_bits']:.6f} bits"
    )
    print(
        "BCH-specific counting is needed only at weights: "
        + ", ".join(map(str, prefix_reduction["bch_specific_shells_required"]))
    )


if __name__ == "__main__":
    main()
