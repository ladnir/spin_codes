#!/usr/bin/env python3
"""Certify an arbitrary fixed [512,256] outer plus a linearized mixer."""

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
OUTPUT = WORKSTREAM / "fixed_outer_linearized_mixer_B512_combined_outward_s22.json"

FIELD_DEGREE = 512
FIELD_MODULUS = (1 << 512) | (1 << 8) | (1 << 5) | (1 << 2) | 1


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial_mod(value: int, modulus: int) -> int:
    degree = modulus.bit_length() - 1
    while value.bit_length() - 1 >= degree:
        value ^= modulus << (value.bit_length() - 1 - degree)
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
    x = 2
    frobenius = x
    halfway_gcd = None
    for exponent in range(1, FIELD_DEGREE + 1):
        frobenius = field_multiply(frobenius, frobenius)
        if exponent == FIELD_DEGREE // 2:
            halfway_gcd = polynomial_gcd(frobenius ^ x, FIELD_MODULUS)
    if halfway_gcd != 1 or frobenius != x:
        raise AssertionError("the selected degree-512 polynomial is reducible")
    return {
        "modulus_hex": hex(FIELD_MODULUS),
        "degree": FIELD_DEGREE,
        "gcd_x_2^256_minus_x": halfway_gcd,
        "x_2^512_equals_x_mod_modulus": True,
    }


def exact_event_failure(caps: list[int]) -> tuple[Fraction, Fraction, Fraction]:
    """Return shell failure, rank failure, and their union bound."""
    code_size = 1 << K
    field_size = 1 << B
    shell_failure = Fraction(0)
    for weight in range(1, B + 1):
        mean = Fraction((code_size - 1) * math.comb(B, weight), field_size)
        cap = caps[weight]
        if cap == 0:
            shell_failure += mean
            continue
        deviation = Fraction(cap + 1) - mean
        if deviation <= 0:
            raise AssertionError(f"cap does not exceed the mean at weight {weight}")
        shell_failure += mean / (mean + deviation * deviation)

    # For b != 0 and a != 0, ker(x -> ax+bx^2)={0,a/b}. The
    # restriction fails exactly when a/b is a nonzero base-code word.
    # The zero map (a,b)=(0,0) also fails.
    rank_failure = Fraction(
        1 + (field_size - 1) * (code_size - 1),
        field_size * field_size,
    )
    return shell_failure, rank_failure, shell_failure + rank_failure


def exp2_of_upper_float(value: float) -> arb:
    return arb(2) ** a(fraction_of_float(value))


def main() -> None:
    ctx.prec = 192
    field_check = verify_field_modulus()
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
    shell_failure, rank_failure, event_failure = exact_event_failure(caps)
    sparse_log2 = float(sparse["claim"]["aggregate_log2_upper"])
    dense_log2 = float(dense["claim"]["aggregate_log2_upper"])
    conditional_mass = exp2_of_upper_float(sparse_log2) + exp2_of_upper_float(
        dense_log2
    )
    conditional_log2 = upper_float(conditional_mass.log() / arb(2).log())
    total_failure = a(event_failure) + conditional_mass
    total_log2 = upper_float(total_failure.log() / arb(2).log())
    if total_log2 >= -40:
        raise AssertionError("linearized-mixer failure does not close at 40 bits")

    payload = {
        "schema": "fixed-outer-linearized-mixer-combined-outward-v1",
        "status": "OUTWARD_FINITE_DISTANCE_CERTIFICATE",
        "claim": {
            "minimum_distance_at_least": D,
            "relative_distance_lower": D / N,
            "conditional_bad_probability_log2_upper": conditional_log2,
            "conditional_margin_bits_lower": -conditional_log2,
            "shell_event_failure_upper": upper_float(a(shell_failure)),
            "rank_failure_upper": upper_float(a(rank_failure)),
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
            "outer_sampler_bits": 2 * B,
            "mixer_field_multiplications_per_active_row": 2,
        },
        "base_constituent": {
            "theorem_scope": "any fixed binary [512,256] linear code",
            "concrete_instance": "direct sum of four fixed EBCH [128,64,22] codes",
            "accumulator_stages_required": 0,
        },
        "linearized_mixer": {
            "field": "GF(2^512)",
            "field_check": field_check,
            "sample": "independent uniform a,b in GF(2^512)",
            "map": "x -> a*x+b*x^2 = x*(a+b*x)",
            "linearity": "Frobenius squaring and multiplication by fixed field elements are GF(2)-linear",
            "pairwise_independence": "for distinct nonzero x,y, the coefficient determinant is x*y*(x+y), which is nonzero",
            "reuse": "one sampled pair (a,b) is repeated in all 4096 outer rows",
        },
        "probability_space": {
            "base": "the [512,256] base constituent is fixed before setup",
            "mixer": "sample a,b independently and uniformly in GF(2^512); no rejection test is performed",
            "routing": "independently sample one uniform 512-coordinate permutation per outer row and one uniform 4096-position permutation per transposed region",
            "inner": "independently sample one uniform binary 23-by-23 linear map at each of the 2^21 RandomStepConv positions; all maps are fixed and shared by every message",
        },
        "proof_accounting": {
            "shell_mean": "(2^256-1)*binom(512,w)/2^512",
            "shell_variance": "pairwise independence gives Var(A_w)<=E[A_w]",
            "rank_failure": "the one-dimensional mixer kernel intersects the fixed base code",
            "failure_union": "rank failure, shell-cap failure, and the conditional distance failure",
        },
        "component_receipts": [
            {"file": SPARSE.name, "sha256": sha256(SPARSE)},
            {"file": DENSE.name, "sha256": sha256(DENSE)},
        ],
        "limitations": [
            "The algebraic mixer changes the BCH-plus-BA construction; it is a certified fallback, not a proof of BA-only concentration.",
            "The theorem concerns RandomStepConv-M22, not RM2Sub.",
            "The two-field-multiplication mixer cost is an operation count, not a benchmark result.",
            "The theorem is finite at k=2^20; an asymptotic family needs a separate parameter schedule.",
        ],
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(f"output={OUTPUT}")


if __name__ == "__main__":
    main()
