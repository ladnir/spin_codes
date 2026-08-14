#!/usr/bin/env python3
"""Bounded pair-capacity diagnostic for the fixed g=8 three-band layout.

This script reconstructs the fixed-monomial local cube moments with the
surjective band-projection formula and applies the Chernoff relaxation in
Corollary 2 of ``explorations/g8_three_band_nonclumping_outer_lemma.md``.

The formula is exact for a fixed, unaveraged tile monomial.  Jensen's
inequality makes it an upper envelope for the committed permutation-averaged
tile factor.  The resulting numbers remain binary/high-precision discovery
diagnostics; they are not outward certificates.

An active packet-profile barycenter does not determine the number ``s`` of
nonzero outer-message blocks.  Consequently, ``--support-sizes`` produces a
restricted-s diagnostic only.  A complete Corollary-2 price is reported only
when ``--all-support-sizes`` explicitly sums every s from 0 through 16384.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import tempfile
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
from scipy.optimize import brentq
from scipy.special import logsumexp


ROOT = Path(__file__).resolve().parents[1]
BLOCKS = 16_384
TILES_PER_BAND = 256
BLOCKS_PER_TILE = 64
MESSAGE_BITS = 64
BAND_DIMENSIONS = (42, 43, 43)
PACKET_CLASSES = tuple(math.comb(8, j) for j in range(9))

DEFAULT_MANIFEST = ROOT / "G8_SUPPORT_MANIFEST.json"
DEFAULT_GEOMETRY = ROOT / "out/g8_factorized_geometry_persistence_wave.json"
DEFAULT_REPLAY = ROOT / "out/g8_fresh_outer_family_round_replay_197.json"
DEFAULT_FUGACITIES = ROOT / "out/g8_outer_family_price_probe_top8.json"
FRESH_ROW_FUGACITIES = ROOT / "out/g8_fresh_round_bl2_outer_components.json"
DEFAULT_SPECTRA = ROOT / "scripts/ebch128_fixed_band_projection_spectra.csv"
DEFAULT_PROJECTION_CERT = ROOT / "scripts/certify_bch_band_projections.py"
DEFAULT_TILE_CERT = ROOT / "scripts/certify_three_band_tile_map.py"
DEFAULT_COMPONENT_CATALOGS = (
    ROOT / "out/g8_factorized_inner_components_v1.json",
    ROOT / "out/g8_factorized_outer_components_v1.json",
    ROOT / "out/g8_cell_aware_factorized_columns.json",
    ROOT / "out/g8_corrected_bl2_outer_components.json",
    ROOT / "out/g8_corrected_total_spectrum_outer_components.json",
    ROOT / "out/g8_fresh_round_total_spectrum_outer_components.json",
    ROOT / "out/g8_fresh_round_bl2_outer_components.json",
)
DEFAULT_OUTPUT = ROOT / "out/g8_pair_capacity_nonclumping_diagnostic.json"

PINNED_SHA256 = {
    "manifest": "cc446d0c6a0c986a8712ef78c4aa7c9d6073687665b7e0f315164339a5b81616",
    "replay": "a58c97af1311a4911c8b81359ddd8f9488c44b03911cb563410ee186d4e8dacf",
    "spectra": "43ae548f4bd9b7a683989f1e859347c792cf6ffa64698a3824510e3f79342053",
    "projection_certificate": "240ecc9c17f7a85d6c4aa34f7c3a56bd8802e2ca8441b8f833c545e52939866b",
    "tile_map_certificate": "45b69f70d4a5cc459a881b3321004f942f8ac83cd7f1f45a3616f87e1bd1d2e3",
}

OUTER_FAMILY_TARGET_SCHEMA = "permute-conv.packet-group-g8-outer-family-price-probe.v1"
FRESH_ROW_TARGET_SCHEMA = "permute-conv.packet-group-g8-fresh-outer-family-components.v1"
TARGET_SOURCE_SHA256 = {
    OUTER_FAMILY_TARGET_SCHEMA: "797e097f39fdfea4b30bd62a7c8c7efb11450beb0446776b5d3ccc07da85c6e5",
    FRESH_ROW_TARGET_SCHEMA: "23f8f43b4c83793d56630bf7b5da4d66efe95e0749169bbe48f45a0e628f1af4",
}
FINNER3_VS_BL2_GAP_BITS = 59_524.0


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json_decimal(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle, parse_float=Decimal)


def json_scalar(value):
    if isinstance(value, Decimal):
        return str(value)
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {key: json_scalar(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_scalar(item) for item in value]
    return value


def atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as handle:
            json.dump(json_scalar(payload), handle, sort_keys=True, separators=(",", ":"))
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, path)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def validate_pinned_inputs(paths: dict[str, Path]) -> dict[str, object]:
    bindings: dict[str, object] = {}
    for label, path in paths.items():
        if not path.is_file():
            raise ValueError(f"missing input {label}: {path}")
        digest = sha256_file(path)
        expected = PINNED_SHA256.get(label)
        if expected is not None and digest != expected:
            raise ValueError(
                f"{label} digest mismatch: expected {expected}, observed {digest}"
            )
        bindings[label] = {"path": str(path), "sha256": digest}
    return bindings


def validate_projection_spectra(path: Path) -> dict:
    totals = {band: 0 for band in range(3)}
    zero_counts = {band: None for band in range(3)}
    weights: dict[int, set[int]] = {band: set() for band in range(3)}
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            if row["code"] != "primal":
                continue
            band = int(row["band"])
            weight = int(row["weight"])
            count = int(row["count"])
            if band not in totals or weight < 0 or count < 0:
                raise ValueError("invalid fixed-band projection spectrum row")
            totals[band] += count
            weights[band].add(weight)
            if weight == 0:
                zero_counts[band] = count
    expected = 1 << MESSAGE_BITS
    if any(totals[band] != expected for band in totals):
        raise ValueError(f"projection spectrum totals are not 2^64: {totals}")
    if any(zero_counts[band] != 1 for band in zero_counts):
        raise ValueError(f"projection spectrum zero counts are not one: {zero_counts}")
    return {
        "primal_totals": totals,
        "primal_zero_counts": zero_counts,
        "primal_weight_ranges": {
            band: [min(weights[band]), max(weights[band])] for band in weights
        },
        "role": (
            "fixed-projection closure check; local inclusion-exclusion uses the "
            "certified inside ranks 42/43/43"
        ),
    }


def group_moments(log_t: Sequence[Decimal], precision: int) -> list[Decimal]:
    """Return R_w(t), w=0..64, with Decimal arithmetic."""

    if len(log_t) != 9 or log_t[0] != 0:
        raise ValueError("log fugacities must have length nine and log_t[0]=0")
    with localcontext() as context:
        context.prec = precision
        polynomial = [Decimal(classes) * value.exp() for classes, value in zip(PACKET_CLASSES, log_t)]
        power = [Decimal(1)]
        for _ in range(8):
            next_power = [Decimal(0)] * (len(power) + 8)
            for left, left_value in enumerate(power):
                for right, right_value in enumerate(polynomial):
                    next_power[left + right] += left_value * right_value
            power = next_power
        return [
            +(power[weight] / Decimal(math.comb(64, weight)))
            for weight in range(65)
        ]


def local_fixed_monomial_moments(
    log_t: Sequence[Decimal],
    dimension: int,
    precision: int,
    max_r: int = 64,
) -> list[Decimal]:
    """Reconstruct E[G^3] for r independent uniform nonzero blocks.

    For F=E_perm[G], Jensen gives E[F^3] <= E_perm E[G^3].  Thus the returned
    fixed-monomial values are local envelopes, not exact averaged moments.
    """

    if dimension not in (42, 43):
        raise ValueError("only the certified dimensions 42 and 43 are supported")
    if not 0 <= max_r <= 64:
        raise ValueError("max_r must lie in [0,64]")
    with localcontext() as context:
        context.prec = precision
        moments = group_moments(log_t, precision)
        cubes = [value * value * value for value in moments]
        a_values = [
            sum(
                (
                    Decimal(math.comb(d, weight)) * cubes[weight]
                    for weight in range(d + 1)
                ),
                Decimal(0),
            )
            for d in range(max_r + 1)
        ]
        denominator_base = Decimal((1 << MESSAGE_BITS) - 1)
        kernel_base = Decimal(1 << (MESSAGE_BITS - dimension))
        result: list[Decimal] = []
        for r in range(max_r + 1):
            positive = Decimal(0)
            negative = Decimal(0)
            for d in range(r + 1):
                term = (
                    Decimal(math.comb(r, d))
                    * (kernel_base**d)
                    * (a_values[d] ** dimension)
                )
                if (r - d) & 1:
                    negative += term
                else:
                    positive += term
            numerator = positive - negative
            if numerator <= 0:
                raise ArithmeticError(
                    f"nonpositive inclusion-exclusion result at c={dimension}, r={r}; "
                    "increase --precision"
                )
            result.append(+(numerator / (denominator_base**r)))
        return result


def stable_local_log_table(
    log_t: Sequence[Decimal],
    dimension: int,
    precision: int,
) -> tuple[list[Decimal], dict]:
    """Compute local logs twice and require guard-digit agreement."""

    first = local_fixed_monomial_moments(log_t, dimension, precision)
    second = local_fixed_monomial_moments(log_t, dimension, precision + 32)
    with localcontext() as context:
        context.prec = precision
        first_logs = [value.ln() for value in first]
        second_logs = [value.ln() for value in second]
        maximum_error = max(abs(a - b) for a, b in zip(first_logs, second_logs))
        tolerance = Decimal(10) ** (-(precision // 2))
        if maximum_error > tolerance:
            raise ArithmeticError(
                f"local moment guard-digit disagreement {maximum_error} > {tolerance}"
            )
        return second_logs, {
            "working_decimal_digits": precision + 32,
            "crosscheck_decimal_digits": precision,
            "max_log_absolute_difference": str(maximum_error),
            "required_tolerance": str(tolerance),
        }


def fit_singleton_matched_pair_majorant(log_moments: Sequence[Decimal]) -> dict:
    """Fit the least nonnegative pair interaction after matching r=1.

    This is a discovery statistic for equation (7) in the note.  It does not
    establish a polymer convergence condition or a motif remainder bound.
    """

    if len(log_moments) != 65:
        raise ValueError("pair-majorant fit requires local moments r=0,...,64")
    with localcontext() as context:
        context.prec = 96
        log_a = +(log_moments[1] - log_moments[0]) / Decimal(3)
        candidates = []
        for r in range(2, 65):
            numerator = (
                (log_moments[r] - log_moments[0]) / Decimal(3)
                - Decimal(r) * log_a
            )
            candidates.append((+(numerator / Decimal(math.comb(r, 2))), r))
        raw_log1p_rho, worst_r = max(candidates, key=lambda row: (row[0], -row[1]))
        log1p_rho = max(Decimal(0), raw_log1p_rho)
        rho = +log1p_rho.exp() - Decimal(1)
        a = +log_a.exp()
    return {
        "log_a_natural": log_a,
        "a": a,
        "raw_max_log1p_rho_natural": raw_log1p_rho,
        "log1p_rho_natural": log1p_rho,
        "rho": rho,
        "worst_occupancy_r": worst_r,
        "clipped_at_zero": raw_log1p_rho < 0,
        "status": "DISCOVERY_ONLY_NO_POLYMER_CONVERGENCE_CLAIM",
    }


def exact_t1_formula_audit() -> dict:
    """Check all 130 integer inclusion-exclusion normalization identities."""

    base = (1 << MESSAGE_BITS) - 1
    checks = 0
    for dimension in (42, 43):
        for r in range(65):
            numerator = sum(
                (-1 if (r - d) & 1 else 1)
                * math.comb(r, d)
                * (1 << ((MESSAGE_BITS - dimension) * d))
                * ((1 << d) ** dimension)
                for d in range(r + 1)
            )
            if numerator != base**r:
                raise AssertionError(f"t=1 identity failed at c={dimension}, r={r}")
            checks += 1
    return {
        "checks": checks,
        "result": "m_c_r(1)=1 for c in {42,43}, r in {0,...,64}",
    }


def tilted_stats(log_coefficients: np.ndarray, x: float) -> tuple[float, float]:
    degrees = np.arange(65, dtype=np.float64)
    terms = log_coefficients + x * degrees
    normalizer = float(logsumexp(terms))
    probabilities = np.exp(terms - normalizer)
    return normalizer, float(np.dot(probabilities, degrees))


def tilted_stats_batch(
    log_coefficients: np.ndarray, x: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return log P(exp(x)), tilted means, and variances in one dense batch."""

    degrees = np.arange(65, dtype=np.float64)
    terms = x[:, None] * degrees[None, :] + log_coefficients[None, :]
    maxima = np.max(terms, axis=1)
    terms -= maxima[:, None]
    np.exp(terms, out=terms)
    normalizers = np.sum(terms, axis=1)
    means = (terms @ degrees) / normalizers
    second_moments = (terms @ (degrees * degrees)) / normalizers
    variances = np.maximum(second_moments - means * means, np.finfo(float).tiny)
    log_polynomials = maxima + np.log(normalizers)
    return log_polynomials, means, variances


def coefficient_chernoff_log(log_coefficients: np.ndarray, support_size: int) -> float:
    """Bound log [y^s] P(y)^256 by minimizing over positive y."""

    if support_size == 0:
        return TILES_PER_BAND * float(log_coefficients[0])
    if support_size == BLOCKS:
        return TILES_PER_BAND * float(log_coefficients[64])
    target = support_size / TILES_PER_BAND

    def equation(x: float) -> float:
        return tilted_stats(log_coefficients, x)[1] - target

    lower, upper = -64.0, 64.0
    while equation(lower) >= 0.0 and lower > -4096.0:
        lower *= 2.0
    while equation(upper) <= 0.0 and upper < 4096.0:
        upper *= 2.0
    if equation(lower) >= 0.0 or equation(upper) <= 0.0:
        raise ArithmeticError("failed to bracket the Chernoff coefficient tilt")
    x = float(brentq(equation, lower, upper, xtol=1e-12, rtol=1e-13, maxiter=100))
    log_polynomial = tilted_stats(log_coefficients, x)[0]
    return TILES_PER_BAND * log_polynomial - support_size * x


def coefficient_chernoff_logs_batch(
    log_coefficients: np.ndarray,
    support_sizes: Sequence[int] | np.ndarray,
    *,
    saddlepoint_table_size: int = 2049,
    newton_iterations: int = 8,
) -> np.ndarray:
    """Evaluate all requested coefficient bounds with one saddlepoint table.

    The two scalar Brent solves locate the extreme requested interior means.
    A monotone table initializes every other saddlepoint.  Dense Newton steps
    then invert all means together.  Thus an all-s evaluation uses two scalar
    solves, rather than one solve for each of the 16,383 interior values.
    """

    sizes = np.asarray(support_sizes, dtype=np.int64)
    if sizes.ndim != 1 or sizes.size == 0:
        raise ValueError("support_sizes must be a nonempty one-dimensional array")
    if np.any(sizes < 0) or np.any(sizes > BLOCKS):
        raise ValueError("support size lies outside [0,16384]")
    if saddlepoint_table_size < 65:
        raise ValueError("saddlepoint table must contain at least 65 points")
    if newton_iterations < 1:
        raise ValueError("at least one Newton iteration is required")

    result = np.empty(sizes.size, dtype=np.float64)
    left = sizes == 0
    right = sizes == BLOCKS
    result[left] = TILES_PER_BAND * float(log_coefficients[0])
    result[right] = TILES_PER_BAND * float(log_coefficients[64])
    interior = ~(left | right)
    if not np.any(interior):
        return result

    interior_sizes = sizes[interior]
    targets = interior_sizes.astype(np.float64) / TILES_PER_BAND

    def mean_minus(x: float, target: float) -> float:
        return tilted_stats(log_coefficients, x)[1] - target

    minimum_target = float(np.min(targets))
    maximum_target = float(np.max(targets))
    lower_bracket, upper_bracket = -64.0, 64.0
    while mean_minus(lower_bracket, minimum_target) >= 0.0:
        lower_bracket *= 2.0
        if lower_bracket < -4096.0:
            raise ArithmeticError("failed to bracket the lowest batched saddlepoint")
    while mean_minus(upper_bracket, maximum_target) <= 0.0:
        upper_bracket *= 2.0
        if upper_bracket > 4096.0:
            raise ArithmeticError("failed to bracket the highest batched saddlepoint")
    lower = float(
        brentq(
            mean_minus,
            lower_bracket,
            upper_bracket,
            args=(minimum_target,),
            xtol=1e-13,
            rtol=1e-14,
            maxiter=100,
        )
    )
    upper = float(
        brentq(
            mean_minus,
            lower_bracket,
            upper_bracket,
            args=(maximum_target,),
            xtol=1e-13,
            rtol=1e-14,
            maxiter=100,
        )
    )

    if lower == upper:
        saddlepoints = np.full(targets.shape, lower, dtype=np.float64)
    else:
        table_x = np.linspace(lower, upper, saddlepoint_table_size)
        _table_logs, table_means, _table_variances = tilted_stats_batch(
            log_coefficients, table_x
        )
        if np.any(np.diff(table_means) <= 0.0):
            raise ArithmeticError("saddlepoint mean table is not strictly monotone")
        saddlepoints = np.interp(targets, table_means, table_x)

    for _ in range(newton_iterations):
        _logs, means, variances = tilted_stats_batch(log_coefficients, saddlepoints)
        step = (means - targets) / variances
        proposed = saddlepoints - step
        saddlepoints = np.clip(proposed, lower, upper)
    log_polynomials, means, _variances = tilted_stats_batch(
        log_coefficients, saddlepoints
    )
    maximum_mean_residual = float(np.max(np.abs(means - targets)))
    if maximum_mean_residual > 2e-10:
        raise ArithmeticError(
            f"batched saddlepoint residual too large: {maximum_mean_residual}"
        )
    result[interior] = (
        TILES_PER_BAND * log_polynomials - interior_sizes * saddlepoints
    )
    return result


def corollary2_chernoff_log(
    local_logs_by_dimension: dict[int, np.ndarray], support_size: int, eta: float
) -> float:
    if not 0 <= support_size <= BLOCKS:
        raise ValueError("support size lies outside [0,16384]")
    degrees = np.arange(65, dtype=np.float64)
    pairs = degrees * (degrees - 1.0) / 2.0
    log_binomial = np.array([math.log(math.comb(64, r)) for r in range(65)])
    coefficient_logs: dict[int, float] = {}
    for dimension in (42, 43):
        log_coefficients = (
            log_binomial + local_logs_by_dimension[dimension] - 3.0 * eta * pairs
        )
        coefficient_logs[dimension] = coefficient_chernoff_log(
            log_coefficients, support_size
        )
    active_assignments = support_size * math.log((1 << MESSAGE_BITS) - 1)
    pair_penalty = eta * support_size * (support_size - 1) / 2.0
    band_term = (
        coefficient_logs[42] + 2.0 * coefficient_logs[43]
    ) / 3.0
    return active_assignments + pair_penalty + band_term


def corollary2_chernoff_logs_batch(
    local_logs_by_dimension: dict[int, np.ndarray],
    support_sizes: Sequence[int] | np.ndarray,
    eta: float,
) -> np.ndarray:
    sizes = np.asarray(support_sizes, dtype=np.int64)
    degrees = np.arange(65, dtype=np.float64)
    pairs = degrees * (degrees - 1.0) / 2.0
    log_binomial = np.array([math.log(math.comb(64, r)) for r in range(65)])
    coefficient_logs: dict[int, np.ndarray] = {}
    for dimension in (42, 43):
        log_coefficients = (
            log_binomial + local_logs_by_dimension[dimension] - 3.0 * eta * pairs
        )
        coefficient_logs[dimension] = coefficient_chernoff_logs_batch(
            log_coefficients, sizes
        )
    active_assignments = sizes.astype(np.float64) * math.log(
        (1 << MESSAGE_BITS) - 1
    )
    pair_penalty = eta * sizes.astype(np.float64) * (sizes - 1.0) / 2.0
    band_term = (coefficient_logs[42] + 2.0 * coefficient_logs[43]) / 3.0
    return active_assignments + pair_penalty + band_term


def parse_support_sizes(text: str) -> list[int]:
    values: set[int] = set()
    for token in text.split(","):
        token = token.strip()
        if not token:
            continue
        if "-" in token:
            left, right = token.split("-", 1)
            begin, end = int(left), int(right)
            if end < begin:
                raise ValueError(f"descending support range {token}")
            values.update(range(begin, end + 1))
        else:
            values.add(int(token))
    if not values or min(values) < 0 or max(values) > BLOCKS:
        raise ValueError("support sizes must be a nonempty subset of [0,16384]")
    return sorted(values)


def eta_grid(count: int, minimum_exponent: float, maximum_exponent: float) -> list[float]:
    if not 2 <= count <= 64:
        raise ValueError("eta-count must lie in [2,64], including eta=0")
    if minimum_exponent > maximum_exponent:
        raise ValueError("eta-min-exp must not exceed eta-max-exp")
    positive = np.linspace(minimum_exponent, maximum_exponent, count - 1)
    return [0.0] + [2.0 ** float(exponent) for exponent in positive]


def validate_barycenter_exact(value, label: str) -> list[str]:
    if not isinstance(value, list) or len(value) != 9:
        raise ValueError(f"{label} must contain nine exact rational coordinates")
    result = []
    for coordinate in value:
        if not isinstance(coordinate, str):
            raise ValueError(f"{label} coordinates must be rational strings")
        try:
            parsed = Fraction(coordinate)
        except (ValueError, ZeroDivisionError) as error:
            raise ValueError(f"invalid rational coordinate in {label}: {coordinate}") from error
        if parsed < 0:
            raise ValueError(f"negative rational coordinate in {label}: {coordinate}")
        result.append(coordinate)
    return result


def validate_log_variables(value, label: str) -> list[Decimal]:
    if not isinstance(value, list) or len(value) != 9:
        raise ValueError(f"{label} must contain nine log variables")
    result = []
    for entry in value:
        if isinstance(entry, bool) or not isinstance(entry, (int, Decimal)):
            raise ValueError(f"{label} entries must be JSON numbers")
        converted = entry if isinstance(entry, Decimal) else Decimal(entry)
        if not converted.is_finite():
            raise ValueError(f"{label} entries must be finite")
        result.append(converted)
    if result[0] != 0:
        raise ValueError(f"{label} must have log_t[0]=0")
    return result


def adapt_outer_family_targets(rows: list) -> list[dict]:
    targets = []
    for index, row in enumerate(rows):
        if not isinstance(row, dict):
            raise ValueError(f"outer-family row {index} is not an object")
        node_id = row.get("node_id")
        h2_cell = row.get("h2_cell")
        branch = row.get("three_band_linear_bl")
        if not isinstance(node_id, str) or not isinstance(h2_cell, int):
            raise ValueError(f"outer-family row {index} lacks node_id/h2_cell")
        if node_id not in (f"h2:{h2_cell:03d}", f"h2:{h2_cell:03d}/R"):
            raise ValueError(f"outer-family row {index} has inconsistent h2 identity")
        if not isinstance(branch, dict) or not isinstance(branch.get("parameters"), dict):
            raise ValueError(f"outer-family row {index} lacks three_band_linear_bl.parameters")
        optimizer = branch.get("optimizer")
        if not isinstance(optimizer, dict) or optimizer.get("success") is not True:
            raise ValueError(f"outer-family row {index} does not have a successful frozen optimizer")
        targets.append(
            {
                "source_node_id": node_id,
                "source_h2_cell": h2_cell,
                "source_barycenter_exact": validate_barycenter_exact(
                    row.get("barycenter_exact"), f"outer-family row {index} barycenter_exact"
                ),
                "log_variables": validate_log_variables(
                    branch["parameters"].get("log_variables"),
                    f"outer-family row {index} three_band_linear_bl log_variables",
                ),
                "target_source_schema": OUTER_FAMILY_TARGET_SCHEMA,
            }
        )
    return targets


def adapt_fresh_row_targets(rows: list) -> list[dict]:
    targets = []
    for index, row in enumerate(rows):
        if not isinstance(row, dict) or row.get("accepted") is not True:
            raise ValueError(f"fresh-row target {index} is not an accepted row")
        node_id = row.get("source_node_id")
        h2_cell = row.get("source_h2_cell")
        parameters = row.get("parameters")
        if not isinstance(node_id, str) or not isinstance(h2_cell, int):
            raise ValueError(f"fresh-row target {index} lacks source node/cell")
        if node_id not in (f"h2:{h2_cell:03d}", f"h2:{h2_cell:03d}/R"):
            raise ValueError(f"fresh-row target {index} has inconsistent h2 identity")
        if not isinstance(parameters, dict):
            raise ValueError(f"fresh-row target {index} lacks parameters")
        targets.append(
            {
                "source_node_id": node_id,
                "source_h2_cell": h2_cell,
                "source_barycenter_exact": validate_barycenter_exact(
                    row.get("source_barycenter_exact"),
                    f"fresh-row target {index} source_barycenter_exact",
                ),
                "log_variables": validate_log_variables(
                    parameters.get("log_variables"),
                    f"fresh-row target {index} log_variables",
                ),
                "target_source_schema": FRESH_ROW_TARGET_SCHEMA,
            }
        )
    return targets


def load_targets(path: Path, replay_path: Path) -> list[dict]:
    source = load_json_decimal(path)
    schema = source.get("schema")
    expected_digest = TARGET_SOURCE_SHA256.get(schema)
    if expected_digest is None:
        raise ValueError(f"unsupported fugacity target schema: {schema!r}")
    observed_digest = sha256_file(path)
    if observed_digest != expected_digest:
        raise ValueError(
            f"fugacity target digest mismatch for {schema}: "
            f"expected {expected_digest}, observed {observed_digest}"
        )
    rows = source.get("rows")
    if not isinstance(rows, list) or len(rows) != 8:
        raise ValueError("fugacity source must contain exactly eight frozen rows")
    if schema == OUTER_FAMILY_TARGET_SCHEMA:
        targets = adapt_outer_family_targets(rows)
    elif schema == FRESH_ROW_TARGET_SCHEMA:
        targets = adapt_fresh_row_targets(rows)
    else:  # The digest table and adapters must remain synchronized.
        raise AssertionError(f"missing strict adapter for schema {schema}")
    replay = load_json_decimal(replay_path)
    replay_nodes = {leaf["node_id"] for leaf in replay.get("leaves", [])}
    for target in targets:
        if target["source_node_id"] not in replay_nodes:
            raise ValueError(
                f"fugacity source node {target['source_node_id']} is absent from replay"
            )
    return targets


def evaluate_target(
    target: dict,
    support_sizes: Sequence[int],
    etas: Sequence[float],
    precision: int,
) -> dict:
    log_t = target["log_variables"]
    local_decimal: dict[int, list[Decimal]] = {}
    precision_audits = {}
    for dimension in (42, 43):
        logs, audit = stable_local_log_table(log_t, dimension, precision)
        local_decimal[dimension] = logs
        precision_audits[str(dimension)] = audit
    local_float = {
        dimension: np.array([float(value) for value in values], dtype=np.float64)
        for dimension, values in local_decimal.items()
    }
    sizes = np.asarray(support_sizes, dtype=np.int64)
    candidate_matrix = np.vstack(
        [corollary2_chernoff_logs_batch(local_float, sizes, eta) for eta in etas]
    )
    best_indices = np.argmin(candidate_matrix, axis=0)
    best_logs = candidate_matrix[best_indices, np.arange(sizes.size)]
    eta_zero_logs = candidate_matrix[0]
    support_rows = []
    for column, support_size in enumerate(support_sizes):
        best_log = float(best_logs[column])
        best_eta = etas[int(best_indices[column])]
        support_rows.append(
            {
                "active_outer_message_blocks": support_size,
                "best_eta": best_eta,
                "chernoff_corollary2_log2_upper": best_log / math.log(2.0),
                "eta_zero_log2_upper": float(eta_zero_logs[column]) / math.log(2.0),
                "pair_capacity_saving_bits_on_grid": (
                    float(eta_zero_logs[column]) - best_log
                )
                / math.log(2.0),
            }
        )
    result = {
        **target,
        "local_fixed_monomial_log_moments_natural": {
            "band_0_dimension_42": local_decimal[42],
            "band_1_dimension_43": local_decimal[43],
            "band_2_dimension_43": local_decimal[43],
        },
        "local_moment_interpretation": (
            "exact for fixed monomial G; Jensen upper envelope for F=E_perm[G]"
        ),
        "precision_audit": precision_audits,
        "singleton_matched_pair_majorant": {
            "dimension_42": fit_singleton_matched_pair_majorant(local_decimal[42]),
            "dimension_43": fit_singleton_matched_pair_majorant(local_decimal[43]),
            "scope": (
                "Equation-7 discovery fit only; no polymer convergence or motif "
                "remainder theorem is asserted"
            ),
        },
        "conditional_finner_parameters": {
            "band_p": ["1/3", "1/3", "1/3"],
            "moment_order": 3,
            "status": "SOUND_SYMMETRIC_CONDITIONAL_FINNER_ONLY",
        },
        "support_rows": support_rows,
    }
    if list(support_sizes) == list(range(BLOCKS + 1)):
        eta_zero_logsum = log2sumexp_scalars(
            row["eta_zero_log2_upper"] for row in support_rows
        )
        optimized_logsum = log2sumexp_scalars(
            row["chernoff_corollary2_log2_upper"] for row in support_rows
        )
        saving = eta_zero_logsum - optimized_logsum
        result.update(
            {
                "full_s_eta_zero_logsum_log2_upper": eta_zero_logsum,
                "full_s_optimized_grid_logsum_log2_upper": optimized_logsum,
                "full_s_pair_capacity_saving_bits": saving,
                "finner3_vs_bl2_gap_bits": FINNER3_VS_BL2_GAP_BITS,
                "pair_capacity_erases_finner3_vs_bl2_gap": (
                    saving >= FINNER3_VS_BL2_GAP_BITS
                ),
            }
        )
    return result


def log2sumexp_scalars(values: Iterable[float]) -> float:
    materialized = list(values)
    maximum = max(materialized)
    return maximum + math.log2(sum(2.0 ** (value - maximum) for value in materialized))


def static_self_test(args: argparse.Namespace) -> None:
    input_paths = {
        "manifest": args.manifest,
        "geometry": args.geometry,
        "replay": args.replay,
        "fugacities": args.fugacities,
        "spectra": args.spectra,
        "projection_certificate": args.projection_certificate,
        "tile_map_certificate": args.tile_map_certificate,
    }
    bindings = validate_pinned_inputs(input_paths)
    catalog_bindings = []
    for path in args.component_catalogs:
        if not path.is_file():
            raise ValueError(f"missing component catalogue: {path}")
        catalog_bindings.append({"path": str(path), "sha256": sha256_file(path)})
    bindings["component_catalogs"] = catalog_bindings
    spectrum_audit = validate_projection_spectra(args.spectra)
    targets = load_targets(args.fugacities, args.replay)
    fresh_targets = load_targets(FRESH_ROW_FUGACITIES, args.replay)
    if targets[0]["source_node_id"] != "h2:073":
        raise AssertionError("default outer-family target source does not start at h2:073")
    if any(target["target_source_schema"] != OUTER_FAMILY_TARGET_SCHEMA for target in targets):
        raise AssertionError("outer-family adapter emitted the wrong schema label")
    if any(
        target["target_source_schema"] != FRESH_ROW_TARGET_SCHEMA
        for target in fresh_targets
    ):
        raise AssertionError("fresh-row adapter emitted the wrong schema label")
    t1_audit = exact_t1_formula_audit()

    # Exercise the Decimal reconstruction away from t=1 without launching a
    # frozen-target diagnostic.
    synthetic_log_t = [Decimal(0)] + [Decimal(-j) / Decimal(20) for j in range(1, 9)]
    low = local_fixed_monomial_moments(synthetic_log_t, 42, 72, max_r=4)
    high = local_fixed_monomial_moments(synthetic_log_t, 42, 96, max_r=4)
    with localcontext() as context:
        context.prec = 60
        relative = max(abs(a / b - 1) for a, b in zip(low, high))
    if relative > Decimal("1e-55"):
        raise AssertionError(f"Decimal reconstruction instability: {relative}")

    synthetic_pair_logs = [
        Decimal("0.7")
        + Decimal(3)
        * (
            Decimal(r) * Decimal("0.2")
            + Decimal(math.comb(r, 2)) * Decimal("0.05")
        )
        for r in range(65)
    ]
    pair_fit = fit_singleton_matched_pair_majorant(synthetic_pair_logs)
    if pair_fit["log_a_natural"] != Decimal("0.2"):
        raise AssertionError("singleton-matched pair fit did not recover log(a)")
    if pair_fit["log1p_rho_natural"] != Decimal("0.05"):
        raise AssertionError("singleton-matched pair fit did not recover log(1+rho)")
    unit_pair_fit = fit_singleton_matched_pair_majorant([Decimal(0)] * 65)
    if unit_pair_fit["rho"] != 0 or unit_pair_fit["worst_occupancy_r"] != 2:
        raise AssertionError("unit pair-majorant fit failed")

    unit_logs = {42: np.zeros(65), 43: np.zeros(65)}
    for support_size in (0, 1, 8192, 16384):
        observed = corollary2_chernoff_log(unit_logs, support_size, 0.0)
        exact = (
            support_size * math.log((1 << MESSAGE_BITS) - 1)
            + math.lgamma(BLOCKS + 1)
            - math.lgamma(support_size + 1)
            - math.lgamma(BLOCKS - support_size + 1)
        )
        if observed + 1e-7 < exact:
            raise AssertionError("Chernoff coefficient fell below the exact coefficient")

    representative_sizes = np.array(
        [0, 1, 2, 63, 64, 8191, 8192, 8193, 16320, 16383, 16384],
        dtype=np.int64,
    )
    synthetic_log_coefficients = np.array(
        [
            math.log(math.comb(64, r))
            + 0.017 * r
            - 0.00091 * r * (r - 1) / 2.0
            for r in range(65)
        ],
        dtype=np.float64,
    )
    synthetic_all_s = coefficient_chernoff_logs_batch(
        synthetic_log_coefficients, np.arange(BLOCKS + 1, dtype=np.int64)
    )
    batched = synthetic_all_s[representative_sizes]
    scalar = np.array(
        [
            coefficient_chernoff_log(synthetic_log_coefficients, int(size))
            for size in representative_sizes
        ]
    )
    maximum_batch_error = float(np.max(np.abs(batched - scalar)))
    if maximum_batch_error > 2e-8:
        raise AssertionError(
            f"batched/scalar saddlepoint mismatch: {maximum_batch_error}"
        )
    if synthetic_all_s[0] != TILES_PER_BAND * synthetic_log_coefficients[0]:
        raise AssertionError("batched s=0 endpoint is not exact")
    if synthetic_all_s[-1] != TILES_PER_BAND * synthetic_log_coefficients[64]:
        raise AssertionError("batched s=16384 endpoint is not exact")

    strong_tilt_coefficients = np.array(
        [
            math.log(math.comb(64, r))
            + 0.017 * r
            - 1.5 * r * (r - 1) / 2.0
            for r in range(65)
        ],
        dtype=np.float64,
    )
    strong_batch = coefficient_chernoff_logs_batch(
        strong_tilt_coefficients, representative_sizes
    )
    strong_scalar = np.array(
        [
            coefficient_chernoff_log(strong_tilt_coefficients, int(size))
            for size in representative_sizes
        ]
    )
    maximum_strong_tilt_error = float(np.max(np.abs(strong_batch - strong_scalar)))
    if maximum_strong_tilt_error > 2e-8:
        raise AssertionError(
            f"strong-tilt batched/scalar mismatch: {maximum_strong_tilt_error}"
        )

    batch_corollary = corollary2_chernoff_logs_batch(
        unit_logs, representative_sizes, 2.0**-12
    )
    scalar_corollary = np.array(
        [
            corollary2_chernoff_log(unit_logs, int(size), 2.0**-12)
            for size in representative_sizes
        ]
    )
    maximum_corollary_error = float(
        np.max(np.abs(batch_corollary - scalar_corollary))
    )
    if maximum_corollary_error > 2e-8:
        raise AssertionError(
            f"batched/scalar Corollary-2 mismatch: {maximum_corollary_error}"
        )

    restricted = parse_support_sizes("0,2-4,16384")
    if restricted != [0, 2, 3, 4, 16384]:
        raise AssertionError("support parser failed")
    if len(range(BLOCKS + 1)) != 16_385:
        raise AssertionError("full support range omitted an endpoint")
    logsum_test = log2sumexp_scalars([1000.0, 1000.0])
    if abs(logsum_test - 1001.0) > 1e-12:
        raise AssertionError("base-2 full-s logsum helper failed")

    # At t=1 and eta=0 the unrelaxed coefficient is C(16384,s), so summing
    # over all s gives (1+(2^64-1))^16384 = 2^(64*16384).
    if (1 + ((1 << MESSAGE_BITS) - 1)) ** BLOCKS != 1 << (MESSAGE_BITS * BLOCKS):
        raise AssertionError("global t=1 normalization identity failed")

    print("pair-capacity static self-test: PASS")
    print(
        f"pinned_inputs={len(input_paths)} component_catalogs={len(catalog_bindings)} "
        f"projection_totals={spectrum_audit['primal_totals']}"
    )
    print(
        f"default_fugacity_targets={len(targets)} first={targets[0]['source_node_id']} "
        f"fresh_schema_targets={len(fresh_targets)} exact_t1_checks={t1_audit['checks']}"
    )
    print("support_contract=restricted-s-is-diagnostic; full-price-requires-0..16384")
    print("averaging_contract=fixed-G-exact; permutation-averaged-F-via-Jensen-upper-envelope")
    print(
        "batch_equivalence="
        f"coefficient_abs_error<={maximum_batch_error:.3g},"
        f"strong_tilt_abs_error<={maximum_strong_tilt_error:.3g},"
        f"corollary_abs_error<={maximum_corollary_error:.3g}"
    )
    print("conditional_bl=symmetric-p=1/3-only; asymmetric-conditioned-formula-rejected")
    print("pair_majorant=singleton-matched-static-pass; no-polymer-convergence-claim")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    support = parser.add_mutually_exclusive_group()
    support.add_argument(
        "--support-sizes",
        help="comma-separated sizes/ranges; explicitly a restricted-s diagnostic",
    )
    support.add_argument(
        "--all-support-sizes",
        action="store_true",
        help="scan and sum every s=0,...,16384 (the only full-price mode)",
    )
    parser.add_argument("--eta-count", type=int, default=33)
    parser.add_argument("--eta-min-exp", type=float, default=-32.0)
    parser.add_argument("--eta-max-exp", type=float, default=-1.0)
    parser.add_argument("--precision", type=int, default=128)
    parser.add_argument(
        "--target-limit",
        type=int,
        default=8,
        help="evaluate only the first N frozen targets (first full probe should use 1)",
    )
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--geometry", type=Path, default=DEFAULT_GEOMETRY)
    parser.add_argument("--replay", type=Path, default=DEFAULT_REPLAY)
    parser.add_argument("--fugacities", type=Path, default=DEFAULT_FUGACITIES)
    parser.add_argument("--spectra", type=Path, default=DEFAULT_SPECTRA)
    parser.add_argument(
        "--projection-certificate", type=Path, default=DEFAULT_PROJECTION_CERT
    )
    parser.add_argument("--tile-map-certificate", type=Path, default=DEFAULT_TILE_CERT)
    parser.add_argument(
        "--component-catalog",
        action="append",
        type=Path,
        dest="component_catalogs",
        help="component catalogue to bind; repeatable (defaults to current bank)",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    if args.precision < 64:
        parser.error("--precision must be at least 64 decimal digits")
    if not 1 <= args.target_limit <= 8:
        parser.error("--target-limit must lie in [1,8]")
    if args.component_catalogs is None:
        args.component_catalogs = list(DEFAULT_COMPONENT_CATALOGS)
    if args.self_test:
        static_self_test(args)
        return
    if args.support_sizes is None and not args.all_support_sizes:
        parser.error("choose --support-sizes or --all-support-sizes")

    paths = {
        "manifest": args.manifest,
        "geometry": args.geometry,
        "replay": args.replay,
        "fugacities": args.fugacities,
        "spectra": args.spectra,
        "projection_certificate": args.projection_certificate,
        "tile_map_certificate": args.tile_map_certificate,
    }
    bindings = validate_pinned_inputs(paths)
    bindings["component_catalogs"] = [
        {"path": str(path), "sha256": sha256_file(path)}
        for path in args.component_catalogs
    ]
    projection_audit = validate_projection_spectra(args.spectra)
    all_targets = load_targets(args.fugacities, args.replay)
    targets = all_targets[: args.target_limit]
    support_sizes = (
        list(range(BLOCKS + 1))
        if args.all_support_sizes
        else parse_support_sizes(args.support_sizes)
    )
    etas = eta_grid(args.eta_count, args.eta_min_exp, args.eta_max_exp)
    target_results = [
        evaluate_target(target, support_sizes, etas, args.precision)
        for target in targets
    ]
    complete = support_sizes == list(range(BLOCKS + 1))
    artifact = {
        "schema": "permute-conv.g8-pair-capacity-nonclumping-diagnostic.v1",
        "status": (
            "DIAGNOSTIC_FULL_S_SCAN_NO_OUTWARD_CLAIM"
            if complete
            else "RESTRICTED_S_SCAN_DIAGNOSTIC_NOT_A_FULL_OUTER_BOUND"
        ),
        "scope_limit": (
            "High-precision local fixed-G moments and binary64 Chernoff prices. "
            "No motif remainder, graph/puncture correction, or outward replay."
        ),
        "local_moment_semantics": (
            "The inclusion-exclusion formula is exact for fixed G. Jensen bounds "
            "the cube moment of F=E_perm[G] by the permutation average of fixed-G moments."
        ),
        "support_coverage": {
            "active_outer_message_blocks": support_sizes,
            "complete_0_through_16384": complete,
            "profile_barycenter_does_not_determine_s": True,
        },
        "configuration": {
            "band_dimensions": BAND_DIMENSIONS,
            "tiles_per_band": TILES_PER_BAND,
            "blocks_per_tile": BLOCKS_PER_TILE,
            "eta_grid": etas,
            "decimal_precision": args.precision,
            "chernoff_coefficient": True,
            "target_limit": args.target_limit,
            "available_frozen_targets": len(all_targets),
            "selected_frozen_targets": len(targets),
            "finner3_vs_bl2_gap_gate_bits": FINNER3_VS_BL2_GAP_BITS,
            "coefficient_solver": (
                "two endpoint Brent solves plus monotone table and batched Newton per "
                "(target,eta,dimension)"
            ),
        },
        "source_bindings": bindings,
        "projection_spectrum_audit": projection_audit,
        "t1_local_formula_audit": exact_t1_formula_audit(),
        "t1_global_normalization": "2^(64*16384)=2^1048576",
        "targets": target_results,
    }
    if complete:
        artifact["full_s_target_summaries"] = [
            {
                "source_node_id": target["source_node_id"],
                "full_s_eta_zero_logsum_log2_upper": target[
                    "full_s_eta_zero_logsum_log2_upper"
                ],
                "full_s_optimized_grid_logsum_log2_upper": target[
                    "full_s_optimized_grid_logsum_log2_upper"
                ],
                "full_s_pair_capacity_saving_bits": target[
                    "full_s_pair_capacity_saving_bits"
                ],
                "pair_capacity_erases_finner3_vs_bl2_gap": target[
                    "pair_capacity_erases_finner3_vs_bl2_gap"
                ],
            }
            for target in target_results
        ]
    atomic_json(args.output, artifact)
    print(f"wrote {args.output}")
    print(f"sha256={sha256_file(args.output)}")
    print(f"status={artifact['status']}")


if __name__ == "__main__":
    main()
