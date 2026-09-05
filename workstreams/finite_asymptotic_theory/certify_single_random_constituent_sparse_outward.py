#!/usr/bin/env python3
"""Outward sparse certificate for one repeated random [512,256] code.

The script verifies exact spectrum-shell transfers for occupations one and
two.  It then verifies the shared-category two-band transfer for occupations
3 through 159.  Stored binary64 tilts are witnesses only.  Arb encloses all
transcendental values, and each positive binary64 recurrence operation is
advanced by one ULP.
"""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path

from flint import arb, ctx
import numpy as np

from certify_single_random_constituent_dense_outward import (
    B,
    D,
    K,
    L,
    MEMORY,
    N,
    P_CENTRAL,
    P_DEFECT,
    a,
    exact_band_majorant,
    exact_caps,
    fraction_of_float,
    upper_float,
)


WORKSTREAM = Path(__file__).resolve().parent
Q1_DIAGNOSTIC = WORKSTREAM / "single_random_constituent_B512_highprob_q1_s22.json"
Q2_DIAGNOSTIC = WORKSTREAM / "single_random_constituent_B512_highprob_q2_s22.json"
SHARED_DIAGNOSTICS = (
    WORKSTREAM / "single_random_constituent_B512_shared_two_band_q3_16_s22.json",
    WORKSTREAM / "single_random_constituent_B512_shared_two_band_q17_32_s22.json",
    WORKSTREAM / "single_random_constituent_B512_shared_two_band_q33_64_s22.json",
    WORKSTREAM / "single_random_constituent_B512_shared_two_band_q65_128_s22.json",
    WORKSTREAM / "single_random_constituent_B512_shared_two_band_q129_159_s22.json",
)
OUTPUT = WORKSTREAM / "single_random_constituent_B512_sparse_outward_s22.json"
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


def matmul2_up(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    result = np.empty_like(left)
    result[..., 0, 0] = add_up(
        mul_up(left[..., 0, 0], right[0, 0]),
        mul_up(left[..., 0, 1], right[1, 0]),
    )
    result[..., 0, 1] = add_up(
        mul_up(left[..., 0, 0], right[0, 1]),
        mul_up(left[..., 0, 1], right[1, 1]),
    )
    result[..., 1, 0] = add_up(
        mul_up(left[..., 1, 0], right[0, 0]),
        mul_up(left[..., 1, 1], right[1, 0]),
    )
    result[..., 1, 1] = add_up(
        mul_up(left[..., 1, 0], right[0, 1]),
        mul_up(left[..., 1, 1], right[1, 1]),
    )
    return result


def normalize(matrices: np.ndarray, exponents: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    maxima = np.max(matrices, axis=(-2, -1))
    valid = (exponents != MIN_EXPONENT) & (maxima > 0.0)
    shifts = np.zeros(exponents.shape, dtype=np.int64)
    if np.any(valid):
        _, raw = np.frexp(maxima[valid])
        shifts[valid] = raw.astype(np.int64)
        matrices[valid] = np.ldexp(matrices[valid], -shifts[valid, None, None])
        matrices[valid] = up_array(matrices[valid])
        exponents[valid] += shifts[valid]
    return matrices, exponents


def align(mantissas: np.ndarray, exponents: np.ndarray, common: np.ndarray) -> np.ndarray:
    valid = exponents != MIN_EXPONENT
    result = np.zeros_like(mantissas)
    if np.any(valid):
        differences = exponents[valid] - common[valid]
        result[valid] = np.ldexp(mantissas[valid], differences[:, None, None])
        result[valid] = up_array(result[valid])
    return result


def combine_terms(terms: list[tuple[np.ndarray, np.ndarray]]) -> tuple[np.ndarray, np.ndarray]:
    common = terms[0][1].copy()
    for _, exponents in terms[1:]:
        common = np.maximum(common, exponents)
    result = align(terms[0][0], terms[0][1], common)
    for matrices, exponents in terms[1:]:
        result = add_up(result, align(matrices, exponents, common))
    return normalize(result, common)


def nonnegative_upper(value: arb) -> float:
    if value == 0:
        return 0.0
    result = upper_float(value)
    if result < 0:
        raise ArithmeticError("negative transition entry")
    return result


def transition_matrices(log_surprisal: float, probabilities: tuple[Fraction, ...]) -> tuple[np.ndarray, arb]:
    """Return zero and Bernoulli-mixture matrices with outward entries."""
    ctx.prec = 192
    u = a(fraction_of_float(log_surprisal))
    surprisal = u.exp()
    z = (-surprisal).exp()
    b = (1 + z) / 2
    q = arb(1) / (1 << MEMORY)
    terminate = q * b
    survive = (1 - q) * b
    zero = (arb(1), arb(0), terminate, survive)
    active = (terminate, survive, terminate, survive)
    rows = [np.asarray(tuple(nonnegative_upper(value) for value in zero), dtype=np.float64).reshape(2, 2)]
    for probability in probabilities:
        p = a(probability)
        mixed = tuple((1 - p) * zero[index] + p * active[index] for index in range(4))
        rows.append(np.asarray(tuple(nonnegative_upper(value) for value in mixed), dtype=np.float64).reshape(2, 2))
    return np.asarray(rows), surprisal


def coefficients_one_dimension(
    zero: np.ndarray,
    zero_exponent: int,
    candidate: np.ndarray,
    candidate_exponent: int,
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
        zero_matrices = np.zeros((size, 2, 2), dtype=np.float64)
        zero_exponents = np.full(size, MIN_EXPONENT, dtype=np.int64)
        zero_matrices[:old_size] = matmul2_up(current[:old_size], zero)
        zero_weights = up_array(
            np.arange(completed + 1, completed + 1 - old_size, -1, dtype=np.float64)
            / float(completed + 1)
        )
        zero_matrices[:old_size] = mul_up(zero_matrices[:old_size], zero_weights[:, None, None])
        zero_exponents[:old_size] = exponents[:old_size] + zero_exponent

        selected_matrices = np.zeros_like(zero_matrices)
        selected_exponents = np.full(size, MIN_EXPONENT, dtype=np.int64)
        if new_maximum:
            selected_matrices[1:] = matmul2_up(current[:new_maximum], candidate)
            selected_weights = up_array(np.arange(1, new_maximum + 1, dtype=np.float64) / float(completed + 1))
            selected_matrices[1:] = mul_up(selected_matrices[1:], selected_weights[:, None, None])
            selected_exponents[1:] = exponents[:new_maximum] + candidate_exponent
        updated, updated_exponents = combine_terms(
            [(zero_matrices, zero_exponents), (selected_matrices, selected_exponents)]
        )
        current.fill(0.0)
        current[:size] = updated
        exponents.fill(MIN_EXPONENT)
        exponents[:size] = updated_exponents
        old_maximum = new_maximum
    return current, exponents


def row_sum_log2_upper(matrices: np.ndarray, exponents: np.ndarray) -> np.ndarray:
    totals = add_up(matrices[..., 0, 0], matrices[..., 0, 1])
    result = np.full(exponents.shape, -math.inf)
    ctx.prec = 160
    log_two = arb(2).log()
    valid = (exponents != MIN_EXPONENT) & (totals > 0)
    for index in zip(*np.nonzero(valid)):
        total = float(totals[index])
        value = a(fraction_of_float(total)).log() / log_two + int(exponents[index])
        result[index] = upper_float(value)
    return result


def powered_row_sum_log2_upper(
    matrices: np.ndarray, exponents: np.ndarray, power: int
) -> np.ndarray:
    count = len(matrices)
    vectors = np.zeros((count, 2), dtype=np.float64)
    vectors[:, 0] = 1.0
    vector_exponents = np.zeros(count, dtype=np.int64)
    for _ in range(power):
        next_zero = add_up(
            mul_up(vectors[:, 0], matrices[:, 0, 0]),
            mul_up(vectors[:, 1], matrices[:, 1, 0]),
        )
        next_live = add_up(
            mul_up(vectors[:, 0], matrices[:, 0, 1]),
            mul_up(vectors[:, 1], matrices[:, 1, 1]),
        )
        vectors[:, 0] = next_zero
        vectors[:, 1] = next_live
        vector_exponents += exponents
        maxima = np.max(vectors, axis=1)
        _, shifts = np.frexp(maxima)
        shifts = shifts.astype(np.int64)
        vectors = np.ldexp(vectors, -shifts[:, None])
        vectors = up_array(vectors)
        vector_exponents += shifts
    totals = add_up(vectors[:, 0], vectors[:, 1])
    result = np.empty(count)
    for index, total in enumerate(totals):
        value = a(fraction_of_float(float(total))).log() / arb(2).log() + int(vector_exponents[index])
        result[index] = upper_float(value)
    return result


def pair_support_coefficients_upper(
    regions: np.ndarray, region_exponents: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    current = np.zeros((B + 1, B + 1, 2, 2), dtype=np.float64)
    current[0, 0] = np.eye(2)
    exponents = np.full((B + 1, B + 1), MIN_EXPONENT, dtype=np.int64)
    exponents[0, 0] = 0
    for completed in range(B):
        size = completed + 1
        old = current[:size, :size]
        old_exponents = exponents[:size, :size]
        next_size = size + 1
        unselected = up_array((size - np.arange(size, dtype=np.float64)) / float(size))
        selected = up_array(np.arange(1, next_size, dtype=np.float64) / float(size))
        terms: list[tuple[np.ndarray, np.ndarray]] = []
        for first_selected, second_selected, region_index in (
            (False, False, 0), (True, False, 1), (False, True, 1), (True, True, 2)
        ):
            matrices = np.zeros((next_size, next_size, 2, 2), dtype=np.float64)
            term_exponents = np.full((next_size, next_size), MIN_EXPONENT, dtype=np.int64)
            row_slice = slice(1, next_size) if first_selected else slice(0, size)
            column_slice = slice(1, next_size) if second_selected else slice(0, size)
            weights_first = selected if first_selected else unselected
            weights_second = selected if second_selected else unselected
            product = matmul2_up(old, regions[region_index])
            product = mul_up(product, weights_first[:, None, None, None])
            product = mul_up(product, weights_second[None, :, None, None])
            matrices[row_slice, column_slice] = product
            term_exponents[row_slice, column_slice] = old_exponents + region_exponents[region_index]
            terms.append((matrices, term_exponents))
        updated, updated_exponents = combine_terms(terms)
        current[:next_size, :next_size] = updated
        exponents[:next_size, :next_size] = updated_exponents
        if (completed + 1) % 128 == 0:
            print(f"q2_coordinates,{completed + 1},{B}", flush=True)
    return current, exponents


def shared_two_band_regions_upper(
    zero: np.ndarray,
    defect: np.ndarray,
    central: np.ndarray,
    maximum_defects: int,
    maximum_centrals: int,
    maximum_occupation: int,
) -> tuple[np.ndarray, np.ndarray]:
    shape = (maximum_defects + 1, maximum_centrals + 1)
    current = np.zeros((*shape, 2, 2), dtype=np.float64)
    current[0, 0] = np.eye(2)
    exponents = np.full(shape, MIN_EXPONENT, dtype=np.int64)
    exponents[0, 0] = 0
    defects, centrals = np.indices(shape)
    totals = defects + centrals
    for completed in range(L):
        terms: list[tuple[np.ndarray, np.ndarray]] = []
        inactive_matrices = matmul2_up(current, zero)
        inactive_exponents = exponents.copy()
        inactive_valid = (totals <= completed) & (totals <= maximum_occupation)
        inactive_weights = np.zeros(shape)
        inactive_weights[inactive_valid] = up_array(
            (completed + 1 - totals[inactive_valid]).astype(np.float64)
            / float(completed + 1)
        )
        inactive_matrices = mul_up(inactive_matrices, inactive_weights[..., None, None])
        inactive_exponents[~inactive_valid] = MIN_EXPONENT
        terms.append((inactive_matrices, inactive_exponents))

        defect_matrices = np.zeros_like(current)
        defect_exponents = np.full(shape, MIN_EXPONENT, dtype=np.int64)
        if maximum_defects:
            product = matmul2_up(current[:-1], defect)
            weights = up_array(defects[1:].astype(np.float64) / float(completed + 1))
            valid = (totals[:-1] <= completed) & (totals[1:] <= maximum_occupation)
            product = mul_up(product, weights[..., None, None])
            defect_matrices[1:] = product
            defect_exponents[1:] = exponents[:-1]
            defect_exponents[1:][~valid] = MIN_EXPONENT
        terms.append((defect_matrices, defect_exponents))

        central_matrices = np.zeros_like(current)
        central_exponents = np.full(shape, MIN_EXPONENT, dtype=np.int64)
        if maximum_centrals:
            product = matmul2_up(current[:, :-1], central)
            weights = up_array(centrals[:, 1:].astype(np.float64) / float(completed + 1))
            valid = (totals[:, :-1] <= completed) & (totals[:, 1:] <= maximum_occupation)
            product = mul_up(product, weights[..., None, None])
            central_matrices[:, 1:] = product
            central_exponents[:, 1:] = exponents[:, :-1]
            central_exponents[:, 1:][~valid] = MIN_EXPONENT
        terms.append((central_matrices, central_exponents))
        current, exponents = combine_terms(terms)
    return current, exponents


def log2_integer_upper(value: int) -> float:
    return upper_float(arb(value).log() / arb(2).log())


def correction_log2_upper(surprisal: arb) -> float:
    return upper_float(D * surprisal / arb(2).log())


def q1_bound(caps: list[int]) -> tuple[float, dict[str, object]]:
    diagnostic = json.loads(Q1_DIAGNOSTIC.read_text(encoding="utf-8"))
    groups: dict[float, list[int]] = {}
    for row in diagnostic["weight_rows"]:
        groups.setdefault(float(row["log_surprisal"]), []).append(int(row["weight"]))
    pointwise: dict[int, float] = {}
    for group_index, (u, weights) in enumerate(sorted(groups.items())):
        transitions, surprisal = transition_matrices(u, (Fraction(1),))
        regions, region_exponents = coefficients_one_dimension(
            transitions[0], 0, transitions[1], 0, L, 1
        )
        coordinates, coordinate_exponents = coefficients_one_dimension(
            regions[0], int(region_exponents[0]), regions[1], int(region_exponents[1]), B, B
        )
        moments = row_sum_log2_upper(coordinates, coordinate_exponents)
        correction = correction_log2_upper(surprisal)
        for weight in weights:
            inner = min(0.0, up(moments[weight] + correction))
            outer = log2_integer_upper(L * caps[weight])
            pointwise[weight] = up(outer + inner)
        print(f"q1_witness_group,{group_index + 1},{len(groups)}", flush=True)
    maximum = max(pointwise.values())
    aggregate = up(maximum + log2_integer_upper(len(pointwise)))
    return aggregate, {"shells": len(pointwise), "maximum_pointwise_log2_upper": maximum}


def q2_bound(caps: list[int]) -> tuple[float, dict[str, object]]:
    support = [weight for weight in range(1, B + 1) if caps[weight]]
    best = np.full((B + 1, B + 1), math.inf)
    # A single common witness suffices after the conservative shell-pair union.
    # The denser diagnostic grid was useful for optimization, not soundness.
    u_values = np.asarray([-7.5])
    for index, u in enumerate(u_values):
        transitions, surprisal = transition_matrices(float(u), (Fraction(1),))
        regions, region_exponents = coefficients_one_dimension(
            transitions[0], 0, transitions[1], 0, L, 2
        )
        pairs, pair_exponents = pair_support_coefficients_upper(regions, region_exponents)
        moments = row_sum_log2_upper(pairs, pair_exponents)
        inner = up_array(moments + correction_log2_upper(surprisal))
        inner = np.minimum(0.0, inner)
        best = np.minimum(best, inner)
        print(f"q2_tilt,{index + 1},{len(u_values)},u,{float(u):.6f}", flush=True)
    maximum = -math.inf
    row_pairs = math.comb(L, 2)
    for first in support:
        first_outer = log2_integer_upper(row_pairs * caps[first])
        for second in support:
            outer = up(first_outer + log2_integer_upper(caps[second]))
            maximum = max(maximum, up(outer + float(best[first, second])))
    term_count = len(support) ** 2
    aggregate = up(maximum + log2_integer_upper(term_count))
    return aggregate, {"shell_pairs": term_count, "maximum_pointwise_log2_upper": maximum}


def sparse_shared_bound(caps: list[int]) -> tuple[float, dict[str, object]]:
    low = exact_band_majorant(caps, 42, 79, P_DEFECT)
    high = exact_band_majorant(caps, 433, 470, 1 - P_DEFECT)
    central_majorant = exact_band_majorant(caps, 80, 432, P_CENTRAL)
    defect_majorant = 2 * max(low, high)
    defect_log2 = upper_float(a(defect_majorant).log() / arb(2).log())
    central_log2 = upper_float(a(central_majorant).log() / arb(2).log())

    groups: dict[float, list[tuple[int, int]]] = {}
    sources = []
    seen: set[tuple[int, int]] = set()
    for path in SHARED_DIAGNOSTICS:
        receipt = json.loads(path.read_text(encoding="utf-8"))
        sources.append({"file": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
        for row in receipt["rows"]:
            composition = (int(row["defect_rows"]), int(row["central_rows"]))
            if composition in seen:
                raise AssertionError("duplicate sparse composition")
            seen.add(composition)
            groups.setdefault(float(row["log_surprisal"]), []).append(composition)
    expected = {(defects, occupation - defects) for occupation in range(3, 160) for defects in range(occupation + 1)}
    if seen != expected:
        raise AssertionError("sparse shared-category receipts do not cover Q=3,...,159")

    pointwise: list[float] = []
    worst_composition = None
    worst_value = -math.inf
    for group_index, (u, compositions) in enumerate(sorted(groups.items())):
        maximum_defects = max(value[0] for value in compositions)
        maximum_centrals = max(value[1] for value in compositions)
        maximum_occupation = max(sum(value) for value in compositions)
        transitions, surprisal = transition_matrices(u, (P_DEFECT, P_CENTRAL))
        regions, region_exponents = shared_two_band_regions_upper(
            transitions[0], transitions[1], transitions[2], maximum_defects, maximum_centrals, maximum_occupation
        )
        selected_matrices = np.asarray([regions[composition] for composition in compositions])
        selected_exponents = np.asarray([region_exponents[composition] for composition in compositions])
        moments = powered_row_sum_log2_upper(selected_matrices, selected_exponents, B)
        correction = correction_log2_upper(surprisal)
        for composition, moment in zip(compositions, moments):
            defects, centrals = composition
            inactive = L - defects - centrals
            multinomial = math.comb(L, defects) * math.comb(L - defects, centrals)
            outer = log2_integer_upper(multinomial)
            outer = up(outer + up(defects * defect_log2))
            outer = up(outer + up(centrals * central_log2))
            value = up(outer + min(0.0, up(float(moment) + correction)))
            pointwise.append(value)
            if value > worst_value:
                worst_value = value
                worst_composition = (defects, centrals)
        print(
            f"shared_witness_group,{group_index + 1},{len(groups)},u,{u:.6f},compositions,{len(compositions)}",
            flush=True,
        )
    ctx.prec = 192
    total = sum(((arb(2) ** a(fraction_of_float(value))) for value in pointwise), arb(0))
    aggregate = upper_float(total.log() / arb(2).log())
    return aggregate, {
        "compositions": len(pointwise),
        "witness_groups": len(groups),
        "maximum_pointwise_log2_upper": worst_value,
        "worst_composition": list(worst_composition) if worst_composition else None,
        "witness_sources": sources,
    }


def main() -> None:
    ctx.prec = 192
    caps, setup_failure = exact_caps()
    q1, q1_details = q1_bound(caps)
    q2, q2_details = q2_bound(caps)
    shared, shared_details = sparse_shared_bound(caps)
    conditional_terms = [q1, q2, shared]
    conditional_sum = sum(
        ((arb(2) ** a(fraction_of_float(value))) for value in conditional_terms),
        arb(0),
    )
    conditional = upper_float(conditional_sum.log() / arb(2).log())
    if conditional >= -40:
        raise AssertionError("sparse conditional bound does not close at 40 bits")
    payload = {
        "schema": "single-random-constituent-sparse-outward-v1",
        "status": "OUTWARD_SPARSE_CERTIFICATE",
        "claim": {
            "minimum_occupation": 1,
            "maximum_occupation": 159,
            "q1_log2_upper": q1,
            "q2_log2_upper": q2,
            "q3_159_log2_upper": shared,
            "aggregate_log2_upper": conditional,
            "aggregate_margin_bits_lower": -conditional,
        },
        "parameters": {
            "outer_code": "one uniform binary [512,256] generator, repeated in all 4096 rows",
            "outer_rows": L,
            "output_bits": N,
            "distance_cutoff": D,
            "memory_bits": MEMORY,
        },
        "details": {"q1": q1_details, "q2": q2_details, "q3_159": shared_details},
        "constituent_event": {
            "total_failure_upper": upper_float(a(setup_failure)),
            "margin_bits_lower": -upper_float(a(setup_failure).log() / arb(2).log()),
        },
        "arithmetic": (
            "192-bit Arb transcendental endpoints and exact rational spectrum data; positive binary64 recurrences advanced one ULP after every addition, multiplication, division, and power-of-two scaling"
        ),
        "limitations": [
            "Occupations 160 through 4096 use the separate outward dense cover.",
            "The inner ensemble is RandomStepConv-M22, not RM2Sub.",
        ],
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(f"output={OUTPUT}")


if __name__ == "__main__":
    main()
