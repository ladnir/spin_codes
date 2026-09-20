#!/usr/bin/env python3
"""Certify the k=2^20 pairwise-systematic outer comparator.

The outer constituent maps u in GF(2^256) to

    (u, a*u + b*u^2),

for one uniformly sampled pair (a,b).  The existing sparse and dense
receipts are conditional only on their integer shell caps.  This verifier
checks that the same caps hold for the new ensemble and combines the three
failure terms with outward Arb arithmetic.
"""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path

from flint import arb, ctx

from certify_single_random_constituent_dense_outward import (
    B,
    D,
    K,
    L,
    MEMORY,
    N,
    a,
    exact_caps,
    fraction_of_float,
    upper_float,
)


WORKSTREAM = Path(__file__).resolve().parent
SPARSE = WORKSTREAM / "single_random_constituent_B512_sparse_outward_s22.json"
DENSE = (
    WORKSTREAM
    / "single_random_constituent_B512_shared_two_band_dense_cover_outward_s22.json"
)
OUTPUT = (
    WORKSTREAM
    / "pairwise_systematic_constituent_B512_combined_outward_s22.json"
)

FIELD_DEGREE = 256
FIELD_MODULUS = (1 << 256) | (1 << 10) | (1 << 5) | (1 << 2) | 1


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial_mod(value: int, modulus: int) -> int:
    modulus_degree = modulus.bit_length() - 1
    while value.bit_length() - 1 >= modulus_degree:
        value ^= modulus << (value.bit_length() - 1 - modulus_degree)
    return value


def polynomial_gcd(first: int, second: int) -> int:
    while second:
        first, second = second, polynomial_mod(first, second)
    return first


def field_multiply(first: int, second: int) -> int:
    result = 0
    while second:
        if second & 1:
            result ^= first
        second >>= 1
        first <<= 1
        if first >> FIELD_DEGREE:
            first ^= FIELD_MODULUS
    return result


def verify_field_modulus() -> dict[str, object]:
    """Apply Rabin's irreducibility criterion for degree 2^8."""
    x = 2
    frobenius = x
    halfway_gcd = None
    for exponent in range(1, FIELD_DEGREE + 1):
        frobenius = field_multiply(frobenius, frobenius)
        if exponent == FIELD_DEGREE // 2:
            halfway_gcd = polynomial_gcd(frobenius ^ x, FIELD_MODULUS)
    if halfway_gcd != 1 or frobenius != x:
        raise AssertionError("the selected degree-256 polynomial is reducible")
    return {
        "modulus_hex": hex(FIELD_MODULUS),
        "degree": FIELD_DEGREE,
        "gcd_x_2^128_minus_x": halfway_gcd,
        "x_2^256_equals_x_mod_modulus": True,
    }


def systematic_shell_mean(weight: int) -> Fraction:
    """Return E[A_weight] for u -> (u, a*u+b*u^2), u != 0."""
    removed_zero_message_term = math.comb(K, weight) if weight <= K else 0
    return Fraction(
        math.comb(B, weight) - removed_zero_message_term,
        1 << K,
    )


def verify_vandermonde_means() -> None:
    for weight in range(1, B + 1):
        direct = sum(
            math.comb(K, input_weight) * math.comb(K, weight - input_weight)
            for input_weight in range(1, K + 1)
            if 0 <= weight - input_weight <= K
        )
        if Fraction(direct, 1 << K) != systematic_shell_mean(weight):
            raise AssertionError(f"Vandermonde check failed at weight {weight}")


def systematic_event_failure(caps: list[int]) -> Fraction:
    """Verify Markov/Cantelli with Var(A_w) <= E[A_w]."""
    failure = Fraction(0)
    for weight in range(1, B + 1):
        mean = systematic_shell_mean(weight)
        cap = caps[weight]
        if cap == 0:
            failure += mean
            continue
        deviation = Fraction(cap + 1) - mean
        if deviation <= 0:
            raise AssertionError(f"cap does not exceed the mean at weight {weight}")
        failure += mean / (mean + deviation * deviation)
    return failure


def exp2_of_upper_float(value: float) -> arb:
    return arb(2) ** a(fraction_of_float(value))


def main() -> None:
    ctx.prec = 192
    field_check = verify_field_modulus()
    verify_vandermonde_means()

    sparse = json.loads(SPARSE.read_text(encoding="utf-8"))
    dense = json.loads(DENSE.read_text(encoding="utf-8"))
    if sparse.get("status") != "OUTWARD_SPARSE_CERTIFICATE":
        raise AssertionError("sparse conditional receipt is not accepting")
    if dense.get("status") != "OUTWARD_DENSE_CERTIFICATE":
        raise AssertionError("dense conditional receipt is not accepting")
    if sparse["claim"]["maximum_occupation"] + 1 != dense["claim"]["minimum_occupation"]:
        raise AssertionError("occupation cover has a gap or overlap")
    if dense["claim"]["maximum_occupation"] != L:
        raise AssertionError("occupation cover does not reach all outer rows")

    caps, _ = exact_caps()
    support = [weight for weight in range(1, B + 1) if caps[weight]]
    if support[0] != 42 or support[-1] != 470:
        raise AssertionError("unexpected cap support")
    event_failure = systematic_event_failure(caps)

    sparse_log2 = float(sparse["claim"]["aggregate_log2_upper"])
    dense_log2 = float(dense["claim"]["aggregate_log2_upper"])
    conditional_mass = exp2_of_upper_float(sparse_log2) + exp2_of_upper_float(
        dense_log2
    )
    conditional_log2 = upper_float(conditional_mass.log() / arb(2).log())
    total_failure = a(event_failure) + conditional_mass
    total_log2 = upper_float(total_failure.log() / arb(2).log())
    if total_log2 >= -40:
        raise AssertionError("pairwise-systematic failure does not close at 40 bits")

    payload = {
        "schema": "pairwise-systematic-constituent-combined-outward-v1",
        "status": "OUTWARD_FINITE_DISTANCE_CERTIFICATE",
        "claim": {
            "minimum_distance_at_least": D,
            "relative_distance_lower": D / N,
            "conditional_bad_probability_log2_upper": conditional_log2,
            "conditional_margin_bits_lower": -conditional_log2,
            "constituent_event_failure_upper": upper_float(a(event_failure)),
            "constituent_event_margin_bits_lower": -upper_float(
                a(event_failure).log() / arb(2).log()
            ),
            "unconditional_setup_or_bad_log2_upper": total_log2,
            "unconditional_margin_bits_lower": -total_log2,
            "comparison_to_2^-40": True,
        },
        "parameters": {
            "message_bits": K * L,
            "output_bits": N,
            "outer_block_bits": B,
            "outer_dimension": K,
            "outer_rows": L,
            "distance_cutoff": D,
            "memory_bits": MEMORY,
            "outer_sampler_bits": 2 * K,
            "outer_field_multiplications_per_active_row": 2,
        },
        "outer_constituent": {
            "field": "GF(2^256)",
            "field_check": field_check,
            "sample": "independent uniform a,b in GF(2^256)",
            "map": "u -> (u, a*u + b*u^2) = (u, u*(a+b*u))",
            "linearity": "Frobenius squaring and multiplication by fixed field elements are GF(2)-linear",
            "injectivity": "the first 256 output bits equal u",
            "pairwise_independence": "for distinct nonzero u,v, the coefficient determinant is u*v*(u+v), which is nonzero",
            "reuse": "one sampled pair (a,b) defines one constituent repeated in all 4096 outer rows",
        },
        "probability_space": {
            "outer": "sample a,b independently and uniformly in GF(2^256); no rejection test is performed",
            "routing": "independently sample one uniform 512-coordinate permutation per outer row and one uniform 4096-position permutation per transposed region",
            "inner": "independently sample one uniform binary 23-by-23 linear map at each of the 2^21 RandomStepConv positions; all maps are fixed and shared by every message",
        },
        "proof_accounting": {
            "shell_caps": "the exact integer cap vector used by the uniform-matrix certificate",
            "shell_mean": "(binom(512,w)-binom(256,w))/2^256 for w<=256 and binom(512,w)/2^256 otherwise",
            "shell_variance": "pairwise independence gives Var(A_w)=sum_u p_u(1-p_u)<=E[A_w]",
            "outer_event": "all nonzero shell counts are at most their integer caps; systematic form makes rank failure impossible",
            "failure_union": "Pr[outer event fails] plus the conditional expected number of nonzero messages below the distance cutoff",
        },
        "component_receipts": [
            {
                "file": SPARSE.name,
                "sha256": sha256(SPARSE),
                "role": "cap-conditional occupations 1 through 159",
            },
            {
                "file": DENSE.name,
                "sha256": sha256(DENSE),
                "role": "cap-conditional occupations 160 through 4096",
            },
        ],
        "limitations": [
            "The reused sparse receipt labels its original uniform-matrix source; its numerical proof depends only on the identical integer shell caps.",
            "The theorem concerns RandomStepConv-M22, not RM2Sub.",
            "The inner samples a fresh 23-by-23 binary matrix at every position; no compact inner implementation claim is made.",
            "The two-field-multiplication outer cost is an operation count, not a benchmark result.",
            "The theorem is a finite k=2^20 result. An asymptotic family requires a separate parameter schedule.",
        ],
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(f"output={OUTPUT}")


if __name__ == "__main__":
    main()
