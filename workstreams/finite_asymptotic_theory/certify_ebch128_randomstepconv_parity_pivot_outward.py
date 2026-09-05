#!/usr/bin/env python3
"""Outward verifier for repeated EBCH128 plus RandomStepConv.

The verifier treats every stored floating-point witness as an exact dyadic
number.  Arb encloses transcendental operations.  Positive binary64 matrix
recurrences advance one ULP after every addition and multiplication.

For Q>=100, the verifier replaces the pivot polynomial by
(Q+1) times its largest coefficient term.  This is looser than the
diagnostic log-sum-exp evaluation but makes the coefficient check small and
independent of binary64 summation accuracy.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path

from flint import arb, ctx
import numpy as np

import evaluate_ebch128_randomstepconv_g1 as base


WORKSTREAM = Path(__file__).resolve().parent
DIAGNOSTIC = WORKSTREAM / "ebch128_randomstepconv_g1_s30_parity_pivot_d11_diagnostic.json"
OUTPUT = WORKSTREAM / "ebch128_randomstepconv_g1_s30_parity_pivot_outward_d11.json"
COMPLEMENT_DIAGNOSTIC = None
SPECTRUM = base.SPECTRUM
SPARSE_MAX_Q = 99
MIN_EXPONENT = np.int64(-(1 << 60))


def up(value: float) -> float:
    if math.isnan(value):
        raise ArithmeticError("NaN in outward calculation")
    return value if value == math.inf else math.nextafter(value, math.inf)


def up_array(values: np.ndarray) -> np.ndarray:
    if np.isnan(values).any():
        raise ArithmeticError("NaN in outward calculation")
    return np.nextafter(values, math.inf)


def mul_up(left: np.ndarray, right: np.ndarray | float) -> np.ndarray:
    return up_array(left * right)


def add_up(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    return up_array(left + right)


def matmul_up(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    result = np.empty(np.broadcast_shapes(left.shape, right.shape), dtype=np.float64)
    result[..., 0, 0] = add_up(
        mul_up(left[..., 0, 0], right[..., 0, 0]),
        mul_up(left[..., 0, 1], right[..., 1, 0]),
    )
    result[..., 0, 1] = add_up(
        mul_up(left[..., 0, 0], right[..., 0, 1]),
        mul_up(left[..., 0, 1], right[..., 1, 1]),
    )
    result[..., 1, 0] = add_up(
        mul_up(left[..., 1, 0], right[..., 0, 0]),
        mul_up(left[..., 1, 1], right[..., 1, 0]),
    )
    result[..., 1, 1] = add_up(
        mul_up(left[..., 1, 0], right[..., 0, 1]),
        mul_up(left[..., 1, 1], right[..., 1, 1]),
    )
    return result


def arb_of_float(value: float) -> arb:
    numerator, denominator = value.as_integer_ratio()
    return arb(numerator) / denominator


def upper_float(value: arb) -> float:
    if not value.is_finite():
        raise ArithmeticError("nonfinite Arb enclosure")
    return math.nextafter(float(value.upper()), math.inf)


def transition_upper(
    surprisal: float, memory_bits: int
) -> tuple[np.ndarray, np.ndarray]:
    ctx.prec = 192
    s = arb_of_float(surprisal)
    z = (-s).exp()
    collision = arb(1) / (1 << memory_bits)
    bit_moment = (1 + z) / 2
    terminate = collision * bit_moment
    survive = (1 - collision) * bit_moment
    t = upper_float(terminate)
    r = upper_float(survive)
    zero = np.asarray(((1.0, 0.0), (t, r)), dtype=np.float64)
    active = np.asarray(((t, r), (t, r)), dtype=np.float64)
    return zero, active


def normalize(
    matrices: np.ndarray, exponents: np.ndarray, valid: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    maxima = np.max(matrices, axis=(-2, -1))
    nonzero = valid & (maxima > 0.0)
    shifts = np.zeros(len(matrices), dtype=np.int64)
    if np.any(nonzero):
        _, raw = np.frexp(maxima[nonzero])
        shifts[nonzero] = raw.astype(np.int64)
        matrices[nonzero] = np.ldexp(
            matrices[nonzero], -shifts[nonzero, None, None]
        )
        matrices[nonzero] = up_array(matrices[nonzero])
        exponents[nonzero] += shifts[nonzero]
    return matrices, exponents


def align_up(
    mantissas: np.ndarray, exponents: np.ndarray, common: np.ndarray
) -> np.ndarray:
    valid = exponents != MIN_EXPONENT
    differences = np.zeros(len(exponents), dtype=np.int64)
    differences[valid] = exponents[valid] - common[valid]
    aligned = np.zeros_like(mantissas)
    if np.any(valid):
        aligned[valid] = np.ldexp(
            mantissas[valid], differences[valid, None, None]
        )
        aligned[valid] = up_array(aligned[valid])
    return aligned


def uniform_coefficients_upper(
    zero: np.ndarray,
    candidate: np.ndarray,
    positions: int,
    maximum_degree: int,
) -> tuple[np.ndarray, np.ndarray]:
    current = np.zeros((maximum_degree + 1, 2, 2), dtype=np.float64)
    current[0] = np.eye(2)
    exponents = np.full(maximum_degree + 1, MIN_EXPONENT, dtype=np.int64)
    exponents[0] = 0
    old_maximum = 0
    for completed in range(positions):
        new_maximum = min(completed + 1, maximum_degree)
        size = new_maximum + 1
        old_size = old_maximum + 1
        zero_terms = np.zeros((size, 2, 2), dtype=np.float64)
        zero_exponents = np.full(size, MIN_EXPONENT, dtype=np.int64)
        zero_terms[:old_size] = matmul_up(current[:old_size], zero)
        zero_weights = up_array(
            np.arange(completed + 1, completed + 1 - old_size, -1, dtype=np.float64)
            / float(completed + 1)
        )
        zero_terms[:old_size] = mul_up(
            zero_terms[:old_size], zero_weights[:, None, None]
        )
        zero_exponents[:old_size] = exponents[:old_size]

        selected_terms = np.zeros_like(zero_terms)
        selected_exponents = np.full(size, MIN_EXPONENT, dtype=np.int64)
        selected = np.arange(1, new_maximum + 1, dtype=np.float64)
        selected_weights = up_array(selected / float(completed + 1))
        selected_terms[1:] = mul_up(
            matmul_up(current[:new_maximum], candidate),
            selected_weights[:, None, None],
        )
        selected_exponents[1:] = exponents[:new_maximum]

        common = np.maximum(zero_exponents, selected_exponents)
        updated = add_up(
            align_up(zero_terms, zero_exponents, common),
            align_up(selected_terms, selected_exponents, common),
        )
        valid = common != MIN_EXPONENT
        updated, common = normalize(updated, common, valid)
        current.fill(0.0)
        current[:size] = updated
        exponents.fill(MIN_EXPONENT)
        exponents[:size] = common
        old_maximum = new_maximum
    return current, exponents


def scaled_power(
    matrices: np.ndarray, matrix_exponents: np.ndarray, exponent: int
) -> tuple[np.ndarray, np.ndarray]:
    count = len(matrices)
    result = np.zeros_like(matrices)
    result[:, 0, 0] = 1.0
    result[:, 1, 1] = 1.0
    result_exponents = np.zeros(count, dtype=np.int64)
    factor = matrices.copy()
    factor_exponents = matrix_exponents.copy()
    remaining = exponent
    valid = np.ones(count, dtype=bool)
    while remaining:
        if remaining & 1:
            result = matmul_up(result, factor)
            result_exponents += factor_exponents
            result, result_exponents = normalize(result, result_exponents, valid)
        remaining >>= 1
        if remaining:
            factor = matmul_up(factor, factor)
            factor_exponents *= 2
            factor, factor_exponents = normalize(factor, factor_exponents, valid)
    return result, result_exponents


def row_sum_log2_upper(matrices: np.ndarray, exponents: np.ndarray) -> np.ndarray:
    totals = add_up(matrices[:, 0, 0], matrices[:, 0, 1])
    result = np.empty(len(matrices), dtype=np.float64)
    ctx.prec = 160
    log_two = arb(2).log()
    for index, total in enumerate(totals):
        result[index] = up(
            upper_float(arb_of_float(float(total)).log() / log_two)
            + float(exponents[index])
        )
    return result


def sum_positive_up(values: np.ndarray) -> np.ndarray:
    """Bound a positive binary64 reduction by the standard gamma factor."""

    count = len(values)
    if count == 0:
        return np.zeros(values.shape[1:], dtype=np.float64)
    rounded = np.sum(values, axis=0)
    epsilon = np.finfo(np.float64).eps / 2.0
    denominator = 1.0 - (count - 1) * epsilon
    if denominator <= 0.0:
        raise ArithmeticError("positive reduction is too long")
    factor = math.nextafter(1.0 / denominator, math.inf)
    return mul_up(rounded, factor)


def fugacity_weights_upper(r: float, maximum: int) -> tuple[np.ndarray, np.ndarray]:
    """Enclose r^m/m! as scaled positive binary64 values."""

    mantissas = np.empty(maximum + 1, dtype=np.float64)
    exponents = np.empty(maximum + 1, dtype=np.int64)
    mantissas[0] = 1.0
    exponents[0] = 0
    for m in range(1, maximum + 1):
        value = up(up(mantissas[m - 1] * r) / float(m))
        mantissa, shift = math.frexp(value)
        mantissas[m] = up(mantissa)
        exponents[m] = exponents[m - 1] + shift
    return mantissas, exponents


def complement_pivot_matrix_upper(
    region_by_zeros: np.ndarray,
    region_exponents: np.ndarray,
    *,
    inactive_rows: int,
    active_rows: int,
    r: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Enclose the exact complement-indexed pivot polynomial."""

    weights, weight_exponents = fugacity_weights_upper(r, active_rows)
    indices = inactive_rows + np.arange(active_rows + 1)
    terms = mul_up(
        region_by_zeros[indices],
        weights[:, None, None],
    )
    term_exponents = region_exponents[indices] + weight_exponents
    common = int(np.max(term_exponents))
    differences = term_exponents - common
    aligned = np.ldexp(terms, differences[:, None, None])
    aligned = up_array(aligned)
    result = sum_positive_up(aligned)[None, :, :]
    exponents = np.asarray([common], dtype=np.int64)
    return normalize(result, exponents, np.ones(1, dtype=bool))


def sparse_inner_log2_upper(
    zero: np.ndarray, candidate: np.ndarray, maximum_q: int
) -> np.ndarray:
    regions, region_exponents = uniform_coefficients_upper(
        zero, candidate, base.L, maximum_q
    )
    fair, fair_exponents = scaled_power(
        regions, region_exponents, base.FAIR_REGIONS
    )
    zero_regions, zero_exponents = scaled_power(
        zero[None, :, :], np.zeros(1, dtype=np.int64), base.L
    )
    complete = matmul_up(fair, zero_regions[0])
    complete_exponents = fair_exponents + zero_exponents[0]
    complete, complete_exponents = normalize(
        complete, complete_exponents, np.ones(len(complete), dtype=bool)
    )
    return row_sum_log2_upper(complete, complete_exponents)


def q1_inner_sum_upper(
    zero: np.ndarray, active: np.ndarray, spectrum: dict[int, int]
) -> arb:
    bit_regions, bit_exponents = uniform_coefficients_upper(
        zero, active, base.L, 1
    )
    if np.any(bit_exponents < -900) or np.any(bit_exponents > 900):
        raise ArithmeticError("unexpected q1 region exponent")
    materialized = np.ldexp(bit_regions, bit_exponents[:, None, None])
    coordinates, coordinate_exponents = uniform_coefficients_upper(
        materialized[0], materialized[1], base.B, base.B
    )
    ctx.prec = 192
    total = arb(0)
    for weight, count in spectrum.items():
        if weight == 0 or count == 0:
            continue
        row_sum = up(coordinates[weight, 0, 0] + coordinates[weight, 0, 1])
        value = arb_of_float(row_sum) * (arb(2) ** int(coordinate_exponents[weight]))
        total += count * value
    return total


def maximum_pivot_term_indices(q: int, r: float) -> list[int]:
    """Return a proved superset of the maximizers of the pivot terms."""
    numerator, denominator = r.as_integer_ratio()

    def difference(m: int) -> int:
        # Positive iff term[m+1] < term[m].
        return (
            (m + 1) * (q - m) * denominator
            - numerator * (base.L - q + m + 1)
        )

    candidates = {0, q}
    if q == 0 or difference(0) > 0:
        return sorted(candidates)
    vertex = (q - 1.0 - r) / 2.0
    vertex_candidates = {
        max(0, min(q - 1, int(math.floor(vertex)) + offset))
        for offset in (-2, -1, 0, 1, 2, 3)
    }
    peak = max(vertex_candidates, key=difference)
    if difference(peak) <= 0:
        return sorted(candidates)
    low = 0
    high = peak
    while high - low > 1:
        middle = (low + high) // 2
        if difference(middle) > 0:
            high = middle
        else:
            low = middle
    candidates.update((max(0, high - 1), high, min(q, high + 1)))
    return sorted(candidates)


def pivot_cost_log2_upper(q: int, r: float) -> float:
    ctx.prec = 192
    log_r = arb_of_float(r).log()
    candidates = maximum_pivot_term_indices(q, r)
    term_upper: list[float] = []
    for m in candidates:
        selected = q - m
        log_binomial = (
            arb(base.L + 1).lgamma()
            - arb(selected + 1).lgamma()
            - arb(base.L - selected + 1).lgamma()
        )
        term = m * log_r - arb(m + 1).lgamma() - log_binomial
        term_upper.append(upper_float(term))
    maximum_term = arb_of_float(max(term_upper))
    cost = (
        base.B * (maximum_term + arb(q + 1).log()) - q * log_r
    ) / arb(2).log()
    return upper_float(cost)


def dense_matrix_log2_upper(
    zero: np.ndarray, candidate: np.ndarray, x_values: np.ndarray
) -> np.ndarray:
    matrices = add_up(
        zero[None, :, :], mul_up(candidate[None, :, :], x_values[:, None, None])
    )
    exponents = np.zeros(len(matrices), dtype=np.int64)
    matrices, exponents = normalize(
        matrices, exponents, np.ones(len(matrices), dtype=bool)
    )
    powered, powered_exponents = scaled_power(matrices, exponents, base.N)
    return row_sum_log2_upper(powered, powered_exponents)


def exact_log2_upper(integer: int) -> float:
    ctx.prec = 192
    return upper_float(arb(integer).log() / arb(2).log())


def exact_fraction_log2_upper(value: Fraction) -> float:
    ctx.prec = 192
    return upper_float(
        (arb(value.numerator) / value.denominator).log() / arb(2).log()
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--diagnostic", type=Path, default=DIAGNOSTIC)
    parser.add_argument("--complement-diagnostic", type=Path, default=COMPLEMENT_DIAGNOSTIC)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--outer-rows", type=int, default=base.L)
    parser.add_argument("--memory-bits", type=int, default=base.MEMORY)
    parser.add_argument("--distance-numerator", type=int, default=11)
    parser.add_argument("--distance-denominator", type=int, default=100)
    parser.add_argument("--certified-margin-bits", type=int, default=20)
    args = parser.parse_args()
    if args.memory_bits < 1:
        raise ValueError("memory bits must be positive")
    if args.certified_margin_bits < 0:
        raise ValueError("certified margin must be nonnegative")
    base.configure_outer_rows(
        args.outer_rows,
        args.distance_numerator,
        args.distance_denominator,
    )

    diagnostic = json.loads(args.diagnostic.read_text(encoding="utf-8"))
    if diagnostic.get("schema") != "ebch128-randomstepconv-parity-pivot-diagnostic-v1":
        raise ValueError("unexpected diagnostic schema")
    parameters = diagnostic["parameters"]
    expected = {
        "outer_rows": base.L,
        "active_outer_rows": base.ACTIVE_ROWS,
        "message_bits": base.K * base.ACTIVE_ROWS,
        "output_bits": base.N,
        "distance_cutoff": base.D,
        "memory_bits": args.memory_bits,
    }
    for key, value in expected.items():
        if parameters.get(key) != value:
            raise ValueError(f"diagnostic parameter mismatch: {key}")

    rows = diagnostic["occupation_rows"]
    if [int(row["active_outer_rows"]) for row in rows] != list(
        range(1, base.ACTIVE_ROWS + 1)
    ):
        raise ValueError("occupation rows are incomplete or unordered")
    spectrum = base.load_spectrum(SPECTRUM)
    envelope, _ = base.body_envelope(spectrum)
    effective_mass = envelope + 1
    per_row_log2 = exact_fraction_log2_upper(effective_mass)
    log_two = arb(2).log()
    pointwise = np.full(base.ACTIVE_ROWS + 1, math.inf)

    # Q=1 uses the exact BCH spectrum.
    row = rows[0]
    s = float.fromhex(row["surprisal_hex"])
    zero, active = transition_upper(s, args.memory_bits)
    q1_moment = q1_inner_sum_upper(zero, active, spectrum)
    ctx.prec = 192
    q1 = (
        arb(base.ACTIVE_ROWS).log()
        + q1_moment.log()
        + base.D * arb_of_float(s)
    ) / log_two
    pointwise[1] = upper_float(q1)

    # Q=2..99 uses the aligned-pivot deletion bound and exact region counts.
    sparse_groups: dict[str, list[int]] = {}
    for q in range(2, SPARSE_MAX_Q + 1):
        sparse_groups.setdefault(rows[q - 1]["surprisal_hex"], []).append(q)
    for group_index, (s_hex, occupations) in enumerate(sorted(sparse_groups.items())):
        s = float.fromhex(s_hex)
        zero, active = transition_upper(s, args.memory_bits)
        candidate = mul_up(add_up(zero, active), 0.5)
        inner = sparse_inner_log2_upper(zero, candidate, max(occupations))
        correction = upper_float(base.D * arb_of_float(s) / log_two)
        for q in occupations:
            inner_bound = min(0.0, up(inner[q] + correction))
            outer = up(
                exact_log2_upper(math.comb(base.ACTIVE_ROWS, q))
                + up(q * per_row_log2)
            )
            pointwise[q] = up(outer + inner_bound)
        print(
            f"sparse_group,{group_index + 1},{len(sparse_groups)},"
            f"occupations,{len(occupations)},surprisal,{s_hex}",
            flush=True,
        )

    # Q>=100 uses the dispersed-pivot coefficient bound.
    dense_groups: dict[str, list[int]] = {}
    for q in range(SPARSE_MAX_Q + 1, base.ACTIVE_ROWS + 1):
        dense_groups.setdefault(rows[q - 1]["surprisal_hex"], []).append(q)
    for group_index, (s_hex, occupations) in enumerate(sorted(dense_groups.items())):
        s = float.fromhex(s_hex)
        zero, active = transition_upper(s, args.memory_bits)
        candidate = mul_up(add_up(zero, active), 0.5)
        x_values = np.asarray(
            [float.fromhex(rows[q - 1]["region_fugacity_hex"]) for q in occupations]
        )
        matrix_logs = dense_matrix_log2_upper(zero, candidate, x_values)
        correction = upper_float(base.D * arb_of_float(s) / log_two)
        for local, q in enumerate(occupations):
            r = float.fromhex(rows[q - 1]["pivot_fugacity_hex"])
            pivot = pivot_cost_log2_upper(q, r)
            x = float(x_values[local])
            ctx.prec = 192
            scalar = (
                arb(q + 1).lgamma()
                - q * arb(base.B).log()
                - q * (base.B - 1) * arb_of_float(x).log()
            ) / log_two
            conditional = up(
                up(upper_float(scalar) + pivot) + matrix_logs[local]
            )
            inner_bound = min(0.0, up(conditional + correction))
            outer = up(
                exact_log2_upper(math.comb(base.ACTIVE_ROWS, q))
                + up(q * per_row_log2)
            )
            pointwise[q] = up(outer + inner_bound)
        print(
            f"dense_group,{group_index + 1},{len(dense_groups)},"
            f"occupations,{len(occupations)},surprisal,{s_hex}",
            flush=True,
        )

    complement_rows: list[dict[str, object]] = []
    if args.complement_diagnostic is not None:
        complement = json.loads(
            args.complement_diagnostic.read_text(encoding="utf-8")
        )
        if complement.get("schema") != "ebch128-randomstepconv-dense-complement-diagnostic-v1":
            raise ValueError("unexpected complement diagnostic schema")
        complement_parameters = complement["parameters"]
        complement_expected = {
            "outer_rows": base.L,
            "active_message_rows": base.ACTIVE_ROWS,
            "output_bits": base.N,
            "distance_cutoff": base.D,
            "memory_bits": args.memory_bits,
        }
        for key, value in complement_expected.items():
            if complement_parameters.get(key) != value:
                raise ValueError(f"complement diagnostic parameter mismatch: {key}")
        complement_rows = complement["occupation_rows"]
        inactive_values = [int(row["inactive_outer_rows"]) for row in complement_rows]
        if inactive_values != list(range(len(complement_rows))):
            raise ValueError("complement rows are incomplete or unordered")

        complement_groups: dict[str, list[dict[str, object]]] = {}
        for row in complement_rows:
            complement_groups.setdefault(str(row["surprisal_hex"]), []).append(row)
        for group_index, (s_hex, group_rows) in enumerate(
            sorted(complement_groups.items())
        ):
            s = float.fromhex(s_hex)
            zero, active = transition_upper(s, args.memory_bits)
            candidate = mul_up(add_up(zero, active), 0.5)
            region_by_zeros, region_exponents = uniform_coefficients_upper(
                candidate,
                zero,
                base.L,
                base.L,
            )
            correction = upper_float(base.D * arb_of_float(s) / log_two)
            for row in group_rows:
                inactive = int(row["inactive_outer_rows"])
                q = int(row["active_outer_rows"])
                if q != base.L - inactive:
                    raise ValueError("invalid complement occupation")
                r = float.fromhex(str(row["pivot_fugacity_hex"]))
                pivot_matrix, pivot_exponents = complement_pivot_matrix_upper(
                    region_by_zeros,
                    region_exponents,
                    inactive_rows=inactive,
                    active_rows=q,
                    r=r,
                )
                powered, powered_exponents = scaled_power(
                    pivot_matrix,
                    pivot_exponents,
                    base.B,
                )
                matrix_log2 = row_sum_log2_upper(powered, powered_exponents)[0]
                scalar = (
                    arb(q + 1).lgamma()
                    - q * arb(base.B).log()
                    - q * arb_of_float(r).log()
                ) / log_two
                conditional = up(upper_float(scalar) + matrix_log2)
                inner_bound = min(0.0, up(conditional + correction))
                outer = up(
                    exact_log2_upper(math.comb(base.ACTIVE_ROWS, q))
                    + up(q * per_row_log2)
                )
                pointwise[q] = up(outer + inner_bound)
            print(
                f"complement_group,{group_index + 1},{len(complement_groups)},"
                f"occupations,{len(group_rows)},surprisal,{s_hex}",
                flush=True,
            )

    if not np.isfinite(pointwise[1:]).all():
        raise ArithmeticError("some occupations lack outward witnesses")
    ctx.prec = 256
    aggregate = arb(0)
    for value in pointwise[1:]:
        aggregate += arb(2) ** arb_of_float(float(value))
    aggregate_log2 = upper_float(aggregate.log() / arb(2).log())
    if not aggregate_log2 < -args.certified_margin_bits:
        raise AssertionError(
            f"{args.certified_margin_bits}-bit outward certificate did not close"
        )
    dominant = int(np.argmax(pointwise[1:])) + 1
    payload = {
        "schema": (
            "ebch128-randomstepconv-parity-pivot-outward-v2"
            if complement_rows
            else "ebch128-randomstepconv-parity-pivot-outward-v1"
        ),
        "status": "OUTWARD_DISTANCE_CERTIFICATE",
        "claim": {
            "minimum_distance_at_least": base.D,
            "relative_distance_lower": base.D / base.N,
            "bad_probability_log2_upper": aggregate_log2,
            "margin_bits_lower": -aggregate_log2,
            "certified_margin_bits": args.certified_margin_bits,
            "meets_certified_margin": True,
            "dominant_active_outer_rows": dominant,
            "dominant_pointwise_log2_upper": float(pointwise[dominant]),
        },
        "parameters": {
            **expected,
            "outer_code": "extended BCH [128,64,22]",
            "outer_bits": base.B,
            "outer_dimension": base.K,
            "shortened_outer_rows": base.L - base.ACTIVE_ROWS,
            "sparse_exact_max_occupation": SPARSE_MAX_Q,
        },
        "probability_space": {
            "outer": "one fixed EBCH constituent repeated in every outer row",
            "routing": "independent uniform row-coordinate and region permutations",
            "inner": (
                f"independent RandomStepConv-M{args.memory_bits} matrices "
                "sampled once and shared by every message"
            ),
        },
        "proof_accounting": {
            "occupation_one": "exact BCH spectrum",
            "occupations_2_through_99": "exact uniform-subset regions after aligned parity-input deletion",
            "occupations_100_through_16384": (
                "dispersed parity pivots and two positive coefficient bounds; "
                "the supplied complement receipt replaces its final occupations"
                if complement_rows
                else "dispersed parity pivots and two positive coefficient bounds"
            ),
            "complement_occupations": (
                "exact region coefficients indexed by noncandidate positions, "
                "followed by one positive pivot fugacity"
                if complement_rows
                else None
            ),
            "pivot_polynomial": "bounded by (Q+1) times its largest exactly located term",
            "arithmetic": "Arb transcendental enclosures plus one-ULP outward positive binary64 matrix recurrences",
            "diagnostic_witness_source": args.diagnostic.name,
            "complement_witness_source": (
                args.complement_diagnostic.name
                if args.complement_diagnostic is not None
                else None
            ),
        },
        "occupation_rows": [
            {
                "active_outer_rows": q,
                "pointwise_log2_upper": float(pointwise[q]),
                "surprisal_hex": rows[q - 1]["surprisal_hex"],
                "region_fugacity_hex": rows[q - 1]["region_fugacity_hex"],
                "pivot_fugacity_hex": rows[q - 1]["pivot_fugacity_hex"],
            }
            for q in range(1, base.ACTIVE_ROWS + 1)
        ],
        "limitations": [
            "The theorem concerns RandomStepConv, not RM2Sub.",
            "The random inner is a proof comparator; no efficient representation or implementation is claimed.",
            "The fixed BCH spectrum is imported from scripts/EBCH128_64.wd and must be hash-bound in the final manifest.",
        ],
    }
    if not complement_rows:
        del payload["proof_accounting"]["complement_occupations"]
        del payload["proof_accounting"]["complement_witness_source"]
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
