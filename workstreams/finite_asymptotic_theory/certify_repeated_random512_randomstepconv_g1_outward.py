#!/usr/bin/env python3
"""Outward certificate for repeated random [512,256] + RandomStepConv.

The witness tilts come from the factor-4 binary64 diagnostic. This verifier
takes the certified spectrum factor to be 13/2, treats
each stored tilt as an exact dyadic number, obtains its exponential through
Arb, and evaluates the remaining positive transfer recurrences with
one-binary64-ULP outward rounding after every arithmetic operation.  Per-
degree powers of two prevent underflow from losing positive mass.

The Bernoulli reference's globally zero tuple is retained.  This is a looser
upper bound than the diagnostic, but avoids every subtraction.
"""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path

from flint import arb, ctx
import numpy as np


WORKSTREAM = Path(__file__).resolve().parent
DIAGNOSTIC = (
    WORKSTREAM / "repeated_random512_randomstepconv_g1_s30_allq_d11.json"
)
OUTPUT = (
    WORKSTREAM
    / "repeated_random512_randomstepconv_g1_s30_allq_outward_d11.json"
)

B = 512
K = 256
L = 4140
ACTIVE_ROWS = 4096
N = B * L
DIAGNOSTIC_DISTANCE = (11 * N) // 100
D = (11 * N + 99) // 100
MEMORY = 30
DIAGNOSTIC_SPECTRUM_FACTOR = 4
SPECTRUM_FACTOR = Fraction(13, 2)
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
    """Upper matrix products for left[...,2,2] and fixed right[2,2]."""
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


def arb_of_float(value: float) -> arb:
    numerator, denominator = value.as_integer_ratio()
    return arb(numerator) / denominator


def upper_float(value: arb) -> float:
    if not value.is_finite():
        raise ArithmeticError("nonfinite Arb enclosure")
    return math.nextafter(float(value.upper()), math.inf)


def transition_upper(surprisal: float) -> tuple[np.ndarray, np.ndarray]:
    """Return outward upper entries for T_0 and T_1 at z=exp(-s)."""
    ctx.prec = 192
    s = arb_of_float(surprisal)
    z = (-s).exp()
    q = arb(1) / (1 << MEMORY)
    b = (1 + z) / 2
    terminate = q * b
    survive = (1 - q) * b
    t = upper_float(terminate)
    r = upper_float(survive)
    zero = np.asarray(((1.0, 0.0), (t, r)), dtype=np.float64)
    active = np.asarray(((t, r), (t, r)), dtype=np.float64)
    return zero, active


def normalize_matrices(
    matrices: np.ndarray, exponents: np.ndarray, valid: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    maxima = np.max(matrices, axis=(1, 2))
    nonzero = valid & (maxima > 0.0)
    shifts = np.zeros(len(matrices), dtype=np.int64)
    if np.any(nonzero):
        _, raw_shifts = np.frexp(maxima[nonzero])
        shifts[nonzero] = raw_shifts.astype(np.int64)
        matrices[nonzero] = np.ldexp(
            matrices[nonzero], -shifts[nonzero, None, None]
        )
        # ldexp is exact except at underflow; advancing covers that boundary.
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


def region_coefficients_upper(
    zero: np.ndarray, candidate: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Uniform-subset transfer averages, as mantissa * 2**exponent."""
    current = np.zeros((ACTIVE_ROWS + 1, 2, 2), dtype=np.float64)
    current[0] = np.eye(2)
    exponents = np.full(ACTIVE_ROWS + 1, MIN_EXPONENT, dtype=np.int64)
    exponents[0] = 0
    old_maximum = 0

    for completed in range(L):
        new_maximum = min(completed + 1, ACTIVE_ROWS)
        size = new_maximum + 1
        zero_terms = np.zeros((size, 2, 2), dtype=np.float64)
        zero_exponents = np.full(size, MIN_EXPONENT, dtype=np.int64)
        old_size = old_maximum + 1
        zero_terms[:old_size] = matmul2_up(current[:old_size], zero)
        zero_weights = np.arange(
            completed + 1, completed + 1 - old_size, -1, dtype=np.float64
        )
        zero_weights = up_array(zero_weights / float(completed + 1))
        zero_terms[:old_size] = mul_up(
            zero_terms[:old_size], zero_weights[:, None, None]
        )
        zero_exponents[:old_size] = exponents[:old_size]

        candidate_terms = np.zeros_like(zero_terms)
        candidate_exponents = np.full(size, MIN_EXPONENT, dtype=np.int64)
        selected = np.arange(1, new_maximum + 1, dtype=np.float64)
        candidate_weights = up_array(selected / float(completed + 1))
        candidate_terms[1:] = mul_up(
            matmul2_up(current[:new_maximum], candidate),
            candidate_weights[:, None, None],
        )
        candidate_exponents[1:] = exponents[:new_maximum]

        common = np.maximum(zero_exponents, candidate_exponents)
        updated = add_up(
            align_up(zero_terms, zero_exponents, common),
            align_up(candidate_terms, candidate_exponents, common),
        )
        valid = common != MIN_EXPONENT
        updated, common = normalize_matrices(updated, common, valid)
        current.fill(0.0)
        current[:size] = updated
        exponents.fill(MIN_EXPONENT)
        exponents[:size] = common
        old_maximum = new_maximum
    return current, exponents


def region_powers_log2_upper(
    regions: np.ndarray, region_exponents: np.ndarray
) -> np.ndarray:
    """Upper log2 of [1,0] R_q^B [1,1]^T for every q."""
    count = len(regions)
    vectors = np.zeros((count, 2), dtype=np.float64)
    vectors[:, 0] = 1.0
    vector_exponents = np.zeros(count, dtype=np.int64)
    for _ in range(B):
        next_zero = add_up(
            mul_up(vectors[:, 0], regions[:, 0, 0]),
            mul_up(vectors[:, 1], regions[:, 1, 0]),
        )
        next_live = add_up(
            mul_up(vectors[:, 0], regions[:, 0, 1]),
            mul_up(vectors[:, 1], regions[:, 1, 1]),
        )
        vectors[:, 0] = next_zero
        vectors[:, 1] = next_live
        vector_exponents += region_exponents
        maxima = np.max(vectors, axis=1)
        _, shifts = np.frexp(maxima)
        shifts = shifts.astype(np.int64)
        vectors = np.ldexp(vectors, -shifts[:, None])
        vectors = up_array(vectors)
        vector_exponents += shifts

    totals = add_up(vectors[:, 0], vectors[:, 1])
    result = np.empty(count, dtype=np.float64)
    ctx.prec = 128
    log_two = arb(2).log()
    for index, total in enumerate(totals):
        value = arb_of_float(float(total)).log() / log_two
        result[index] = up(upper_float(value) + float(vector_exponents[index]))
    return result


def setup_event_exact() -> tuple[Fraction, dict[str, object]]:
    message_count = (1 << K) - 1
    tail = Fraction(0)
    central = Fraction(0)
    for weight in range(1, B + 1):
        mean = Fraction(message_count * math.comb(B, weight), 1 << B)
        if mean < 1 / SPECTRUM_FACTOR:
            tail += mean
        else:
            # Cantelli's inequality, using Var(A_w) <= E[A_w].
            central += 1 / (1 + (SPECTRUM_FACTOR - 1) ** 2 * mean)
    rank = sum(Fraction(1 << index, 1 << B) for index in range(K))
    failure = tail + central + rank

    def integer_hash(value: int) -> str:
        encoded = value.to_bytes((value.bit_length() + 7) // 8, "big")
        return hashlib.sha256(encoded).hexdigest()

    details = {
        "low_mean_markov_union_upper": upper_float(arb(tail.numerator) / tail.denominator),
        "central_cantelli_union_upper": upper_float(
            arb(central.numerator) / central.denominator
        ),
        "rank_failure_union_upper": upper_float(arb(rank.numerator) / rank.denominator),
        "total_failure_upper": upper_float(
            arb(failure.numerator) / failure.denominator
        ),
        "success_probability_lower": float(
            (arb(failure.denominator - failure.numerator) / failure.denominator).lower()
        ),
        "exact_total_failure_numerator_bits": failure.numerator.bit_length(),
        "exact_total_failure_denominator_bits": failure.denominator.bit_length(),
        "exact_total_failure_numerator_sha256": integer_hash(failure.numerator),
        "exact_total_failure_denominator_sha256": integer_hash(failure.denominator),
    }
    return failure, details


def exact_log2_upper(integer: int) -> float:
    ctx.prec = 160
    return upper_float(arb(integer).log() / arb(2).log())


def exact_fraction_log2_upper(value: Fraction) -> float:
    ctx.prec = 160
    quotient = arb(value.numerator) / value.denominator
    return upper_float(quotient.log() / arb(2).log())


def main() -> None:
    diagnostic = json.loads(DIAGNOSTIC.read_text(encoding="utf-8"))
    parameters = diagnostic["parameters"]
    expected = {
        "outer_bits": B,
        "outer_dimension": K,
        "outer_rows": L,
        "active_outer_rows": ACTIVE_ROWS,
        "message_bits": K * ACTIVE_ROWS,
        "output_bits": N,
        "distance_cutoff": DIAGNOSTIC_DISTANCE,
        "memory_bits": MEMORY,
        "spectrum_factor": DIAGNOSTIC_SPECTRUM_FACTOR,
    }
    for key, value in expected.items():
        if parameters.get(key) != value:
            raise ValueError(f"diagnostic parameter mismatch: {key}")

    witness_by_q = {
        int(row["active_outer_rows"]): float(row["surprisal"])
        for row in diagnostic["occupation_rows"]
    }
    groups: dict[float, list[int]] = {}
    for occupation, surprisal in witness_by_q.items():
        groups.setdefault(surprisal, []).append(occupation)

    per_row_log2 = exact_fraction_log2_upper(
        SPECTRUM_FACTOR * ((1 << K) - 1)
    )
    pointwise = np.full(ACTIVE_ROWS + 1, math.inf)
    witness_rows: list[dict[str, object]] = []
    for group_index, (surprisal, occupations) in enumerate(sorted(groups.items())):
        zero, active = transition_upper(surprisal)
        candidate = up_array(mul_up(add_up(zero, active), 0.5))
        regions, region_exponents = region_coefficients_upper(zero, candidate)
        inner_log2 = region_powers_log2_upper(regions, region_exponents)
        correction = upper_float(
            arb(D) * arb_of_float(surprisal) / arb(2).log()
        )
        for occupation in occupations:
            inner = min(0.0, up(inner_log2[occupation] + correction))
            combinatorial = exact_log2_upper(
                math.comb(ACTIVE_ROWS, occupation)
            )
            value = up(
                up(combinatorial + up(occupation * per_row_log2)) + inner
            )
            pointwise[occupation] = value
            witness_rows.append(
                {
                    "active_outer_rows": occupation,
                    "surprisal_hex": surprisal.hex(),
                    "pointwise_log2_upper": value,
                    "inner_log2_upper": inner,
                }
            )
        print(
            f"witness_group,{group_index + 1},{len(groups)},"
            f"occupations,{len(occupations)},surprisal,{surprisal.hex()}",
            flush=True,
        )

    if not np.isfinite(pointwise[1:]).all():
        raise ArithmeticError("some occupations lack outward witnesses")
    maximum = up(float(np.max(pointwise[1:])))
    coarse_union_log2 = up(maximum + math.log2(ACTIVE_ROWS))
    failure, setup = setup_event_exact()
    if failure >= 1:
        raise AssertionError("constituent event has no certified probability")
    if not coarse_union_log2 < -40.0:
        raise AssertionError("40-bit distance certificate did not close")
    setup_attempts = 28
    conservative_joint = failure**setup_attempts + Fraction(1, 1 << 108)
    if not conservative_joint < Fraction(1, 1 << 40):
        raise AssertionError("bounded setup plus distance did not close")
    ctx.prec = 192
    joint_arb = arb(conservative_joint.numerator) / conservative_joint.denominator
    unconditional_log2 = upper_float(joint_arb.log() / arb(2).log())

    dominant_q = int(np.argmax(pointwise[1:])) + 1
    witness_rows.sort(key=lambda row: int(row["active_outer_rows"]))
    payload = {
        "schema": "repeated-random512-randomstepconv-g1-allq-outward-v1",
        "status": "OUTWARD_DISTANCE_CERTIFICATE",
        "claim": {
            "minimum_distance_at_least": D,
            "relative_distance_lower": D / N,
            "conditional_bad_probability_log2_upper": coarse_union_log2,
            "conditional_margin_bits_lower": -coarse_union_log2,
            "comparison_to_2^-40": True,
            "maximum_pointwise_log2_upper": maximum,
            "dominant_active_outer_rows": dominant_q,
            "occupation_count": ACTIVE_ROWS,
            "bounded_setup_attempts": setup_attempts,
            "unconditional_abort_or_bad_log2_upper": unconditional_log2,
            "unconditional_comparison_to_2^-40": True,
        },
        "parameters": {
            **{
                key: value
                for key, value in expected.items()
                if key not in {"spectrum_factor", "distance_cutoff"}
            },
            "distance_cutoff": D,
            "spectrum_factor_numerator": SPECTRUM_FACTOR.numerator,
            "spectrum_factor_denominator": SPECTRUM_FACTOR.denominator,
            "parent_message_bits": K * L,
            "shortened_outer_rows": L - ACTIVE_ROWS,
        },
        "constituent_event": {
            "event": (
                "the one sampled [512,256] generator has full row rank and "
                "A_w <= (13/2) E[A_w] for every 1<=w<=512"
            ),
            **setup,
        },
        "probability_space": {
            "outer_selection": (
                "try at most 28 independent uniform binary 256-by-512 "
                "generators, retain the first satisfying the constituent "
                "event, and repeat only that one code"
            ),
            "routing": (
                "independent uniform row-coordinate permutations and "
                "independent uniform region permutations"
            ),
            "inner": (
                "independent uniform binary (M+1)-by-(M+1) linear maps at "
                "all N steps, sampled once and shared by all messages"
            ),
            "conditioning": (
                "the distance probability is conditional on any fixed outer "
                "constituent satisfying the constituent event"
            ),
        },
        "proof_accounting": {
            "outer_row_envelope_log2_upper": per_row_log2,
            "all_zero_reference_tuple": "retained in the positive upper bound",
            "occupation_union": (
                "sum_Q U_Q <= 4096 * max_Q U_Q"
            ),
            "arithmetic": (
                "Arb transcendental endpoints and exact setup rationals; "
                "positive binary64 recurrences advanced one ULP after every "
                "addition, multiplication, division, and power-of-two scaling"
            ),
            "diagnostic_witness_source": DIAGNOSTIC.name,
        },
        "occupation_rows": witness_rows,
        "limitations": [
            "The exact constituent-event test may enumerate all 2^256 row messages; no efficient test or deterministic construction is supplied.",
            "The claim concerns the RandomStepConv proof model, not RM2Sub.",
            "No implementation or running-time claim is made for RandomStepConv.",
        ],
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(json.dumps(payload["constituent_event"], indent=2))
    print(f"output={OUTPUT}")


if __name__ == "__main__":
    main()
