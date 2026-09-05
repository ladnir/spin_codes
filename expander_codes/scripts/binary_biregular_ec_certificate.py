#!/usr/bin/env python3
"""Generate and verify binary two-sided regular wrapped-EC certificates.

The rate-half profiles have odd right degree.  Exact regional transfers cover
the smallest supports.  Positive coefficient blocks cover the outer ranges,
and a Poisson-binomial local bound covers the central range.  Generation uses
SciPy only to select decimal markers.  Verification reads fixed markers and
uses outward-rounded Arb arithmetic for every reported probability bound.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from flint import arb, arb_mat, ctx
from scipy.optimize import minimize_scalar

from binary_biregular_diagnostic import (
    biregular_parity_shell_counts,
    exact_biregular_logterm,
    geometric_blocks,
    saddle_biregular_block_logterm,
)
from check_binary_biregular_ec_certificate import (
    DEGREE_FIVE_SCHEMA,
    ODD_DEGREE_SCHEMA,
    validate_certificate_structure,
)
from ea_certificate import decimal_marker, hamming_ball_bound_arb
from regular_ec_certificate import (
    binomial_shell_bound_arb,
    uniform_slice_transfer_matrices_arb,
    zero_matrix,
)
from regular_ec_certificate import wrapping_input_matrices_arb
from regular_ec_diagnostic import wrapping_input_matrices
from expander_bounds import log_weight_mgf


@dataclass(frozen=True)
class BinaryBiregularECResult:
    success: bool
    total_bound: arb
    security_bits: arb
    exact_bound: arb
    intermediate_bound: arb
    dense_bound: arb
    full_support_bound: arb
    largest_exact_r: int
    largest_exact_term: arb
    intermediate_blocks: int
    dense_blocks: int


def endpoint_mass(matrix: arb_mat, exponent: int, start_state: int) -> arb:
    power = matrix**exponent
    return sum((power[start_state, column] for column in range(power.ncols())), arb(0))


def exact_term_arb(
    *, k: int, left_degree: int, right_degree: int, cutoff: int,
    memory: int, message_weight: int, output_marker: str,
) -> arb:
    region_length = k // right_degree
    denominator, counts = biregular_parity_shell_counts(
        left_vertices=k,
        right_degree=right_degree,
        support_size=message_weight,
    )
    z = arb(output_marker)
    slices = uniform_slice_transfer_matrices_arb(
        length=region_length,
        max_weight=message_weight,
        output_marker=z,
        memory=memory,
    )
    region = zero_matrix(memory + 1)
    for weight, count in counts.items():
        region += slices[weight] * (arb(count) / denominator)
    tail = endpoint_mass(region, left_degree, memory) * z ** (-cutoff)
    return arb(math.comb(k, message_weight)) * tail


def parity_probability_arb(right_degree: int, x: arb) -> arb:
    return (1 - ((1 - x) / (1 + x)) ** right_degree) / 2


def coefficient_block_arb(
    *, k: int, left_degree: int, right_degree: int, cutoff: int,
    memory: int, lo: int, hi: int, input_marker: str, output_marker: str,
) -> arb:
    region_length = k // right_degree
    n = left_degree * region_length
    x = arb(input_marker)
    z = arb(output_marker)
    q = parity_probability_arb(right_degree, x)
    zero, one = wrapping_input_matrices_arb(z, memory)
    transfer = zero * (1 - q) + one * q
    transform = endpoint_mass(transfer, n, memory) * z ** (-cutoff)

    def endpoint(r: int) -> arb:
        choose = binomial_shell_bound_arb(k, r)
        return choose ** (1 - left_degree) * x ** (-left_degree * r)

    # For fixed markers, the logarithm of the shell bound is convex in r.
    # The sum of both endpoint bounds safely dominates their maximum.
    endpoints = endpoint(lo) + endpoint(hi)
    return (
        arb(hi - lo + 1)
        * endpoints
        * (1 + x) ** (left_degree * k)
        * transform
    )


def odd_half_count_variance_arb(p: arb) -> arb:
    """Variance of (J-1)/2 for Bin(5,p) conditioned on odd J."""
    one_minus_p = 1 - p
    odd_probability = (1 - (1 - 2 * p) ** 5) / 2
    mass_one = 10 * p**3 * one_minus_p**2
    mass_two = p**5
    mean = (mass_one + 2 * mass_two) / odd_probability
    return (mass_one + 4 * mass_two) / odd_probability - mean**2


def odd_degree_variance_lower_arb(p: arb, right_degree: int) -> arb:
    """Lower-bound both conditional variances for ``p <= s <= 1/2``.

    For right degree ``2h+1``, each conditional probability-generating
    polynomial has only negative real roots.  It is therefore a product of
    ``h`` Bernoulli generating polynomials.  Each factor variance is unimodal
    on the parameter interval, so its minimum occurs at an endpoint.
    """
    if right_degree < 3 or right_degree % 2 != 1:
        raise ValueError("right_degree must be odd and at least three")
    h = (right_degree - 1) // 2
    odds_ratio = (1 - p) / p

    def factor_variance(rho: arb) -> arb:
        return rho / (1 + rho) ** 2

    def parity_sum(angles: list[arb]) -> arb:
        total = arb(0)
        for angle in angles:
            tangent_squared = angle.tan() ** 2
            at_lo = factor_variance(odds_ratio**2 * tangent_squared)
            at_half = factor_variance(tangent_squared)
            lo_value = at_lo.lower()
            half_value = at_half.lower()
            total += lo_value if lo_value < half_value else half_value
        return total.lower()

    even = parity_sum([
        arb(2 * j + 1) * arb.pi() / (2 * right_degree)
        for j in range(h)
    ])
    odd = parity_sum([
        arb(j + 1) * arb.pi() / right_degree
        for j in range(h)
    ])
    return even if even < odd else odd


def full_support_marker(*, n: int, cutoff: int, memory: int) -> str:
    """Select an output marker for the deterministic all-one EC input."""
    def objective(log_z: float) -> float:
        _, one = wrapping_input_matrices(math.exp(log_z), memory)
        return log_weight_mgf(one, n, memory) - cutoff * log_z

    optimum = minimize_scalar(
        objective,
        bounds=(-40.0, 0.0),
        method="bounded",
        options={"xatol": 2e-10, "maxiter": 160},
    )
    return decimal_marker(math.exp(float(optimum.x)))


def full_support_term_arb(
    *, n: int, cutoff: int, memory: int, output_marker: str
) -> arb:
    z = arb(output_marker)
    _, one = wrapping_input_matrices_arb(z, memory)
    return endpoint_mass(one, n, memory) * z ** (-cutoff)


def dense_low_half_block_arb(
    *, k: int, left_degree: int, right_degree: int = 5,
    cutoff: int, lo: int, hi: int,
) -> arb:
    """Bound one block in k/d_R <= lo <= hi <= (k-1)/2."""
    region_length = k // right_degree
    n = left_degree * region_length
    p_lo = arb(lo) / k
    p_hi = arb(hi) / k
    q_zero_lo = (1 + (1 - 2 * p_lo) ** right_degree) / 2
    variance_lo = (
        odd_half_count_variance_arb(p_lo)
        if right_degree == 5
        else odd_degree_variance_lower_arb(p_lo, right_degree)
    )
    local_mass = (arb.pi() / (8 * region_length * variance_lo)).sqrt()
    if not local_mass < 1:
        raise ValueError("local-limit branch is not below the trivial bound")
    choose_hi = binomial_shell_bound_arb(k, hi)
    binomial_mass_hi = (
        choose_hi * p_hi**hi * (1 - p_hi) ** (k - hi)
    )
    regional_point_mass = (
        q_zero_lo**region_length * local_mass / binomial_mass_hi
    )
    return (
        arb(hi - lo + 1)
        * choose_hi
        * hamming_ball_bound_arb(n, cutoff)
        * regional_point_mass**left_degree
    )


def generate_certificate(
    *, k: int = 1_048_575, left_degree: int = 10,
    right_degree: int = 5, cutoff: int = 230_729, memory: int = 15,
    exact_limit: int = 16, intermediate_relative_width: float = 0.2,
    dense_relative_width: float = 0.02, target_bits: int = 20,
    precision_bits: int = 192,
) -> dict[str, Any]:
    if right_degree < 3 or right_degree % 2 == 0:
        raise ValueError("right_degree must be odd and at least three")
    if k % right_degree or left_degree * k != 2 * k * right_degree:
        raise ValueError("parameters must define a rate-half two-sided regular graph")
    central_start = k // right_degree
    if not 1 <= exact_limit < central_start:
        raise ValueError("invalid exact limit")
    n = left_degree * (k // right_degree)

    exact_weights = []
    for r in range(1, exact_limit + 1):
        selected = exact_biregular_logterm(
            code="ec",
            k=k,
            left_degree=left_degree,
            right_degree=right_degree,
            cutoff=cutoff,
            memory=memory,
            message_weight=r,
        )
        exact_weights.append({
            "r": r,
            "z": decimal_marker(selected.output_marker),
        })

    low_blocks = geometric_blocks(
        exact_limit + 1, central_start - 1, intermediate_relative_width
    )
    complement_blocks = geometric_blocks(
        1, central_start - 1, intermediate_relative_width
    )
    high_blocks = [
        (k - complement_hi, k - complement_lo)
        for complement_lo, complement_hi in reversed(complement_blocks)
    ]
    outer_blocks = []
    for lo, hi in low_blocks + high_blocks:
        selected = saddle_biregular_block_logterm(
            code="ec",
            k=k,
            left_degree=left_degree,
            right_degree=right_degree,
            cutoff=cutoff,
            memory=memory,
            support_start=lo,
            support_limit=hi,
        )
        outer_blocks.append({
            "lo": lo,
            "hi": hi,
            "x": decimal_marker(selected.input_marker),
            "z": decimal_marker(selected.output_marker),
        })

    half = (k - 1) // 2
    central_blocks = [
        {"lo": lo, "hi": hi, "include_complement": True}
        for lo, hi in geometric_blocks(central_start, half, dense_relative_width)
    ]
    certificate = {
        "schema": (
            DEGREE_FIVE_SCHEMA if right_degree == 5 else ODD_DEGREE_SCHEMA
        ),
        "parameters": {
            "k": k,
            "n": n,
            "cutoff": cutoff,
            "left_degree": left_degree,
            "right_degree": right_degree,
            "memory": memory,
            "target_bits": target_bits,
        },
        "verification": {"precision_bits": precision_bits},
        "exact_weights": exact_weights,
        "outer_blocks": outer_blocks,
        "central_blocks": central_blocks,
        "full_support": {
            "r": k,
            "z": full_support_marker(n=n, cutoff=cutoff, memory=memory),
        },
    }
    validate_certificate_structure(certificate)
    return certificate


def verify_certificate(certificate: dict[str, Any]) -> BinaryBiregularECResult:
    validate_certificate_structure(certificate)
    parameters = certificate["parameters"]
    k = int(parameters["k"])
    cutoff = int(parameters["cutoff"])
    left_degree = int(parameters["left_degree"])
    right_degree = int(parameters["right_degree"])
    memory = int(parameters["memory"])
    target_bits = int(parameters["target_bits"])
    n = int(parameters["n"])
    ctx.prec = int(certificate["verification"]["precision_bits"])

    exact_terms = []
    for item in certificate["exact_weights"]:
        exact_terms.append(exact_term_arb(
            k=k,
            left_degree=left_degree,
            right_degree=right_degree,
            cutoff=cutoff,
            memory=memory,
            message_weight=int(item["r"]),
            output_marker=str(item["z"]),
        ))
    exact = sum(exact_terms, arb(0))
    largest_index = max(range(len(exact_terms)), key=lambda index: exact_terms[index])

    intermediate = arb(0)
    for item in certificate["outer_blocks"]:
        intermediate += coefficient_block_arb(
            k=k,
            left_degree=left_degree,
            right_degree=right_degree,
            cutoff=cutoff,
            memory=memory,
            lo=int(item["lo"]),
            hi=int(item["hi"]),
            input_marker=str(item["x"]),
            output_marker=str(item["z"]),
        )

    dense_low = sum((dense_low_half_block_arb(
        k=k,
        left_degree=left_degree,
        right_degree=right_degree,
        cutoff=cutoff,
        lo=int(item["lo"]),
        hi=int(item["hi"]),
    ) for item in certificate["central_blocks"]), arb(0))
    # The structural checker requires every central block to include its
    # complement.  Odd k pairs the two ranges without a fixed midpoint.
    dense = 2 * dense_low

    full_support = full_support_term_arb(
        n=n,
        cutoff=cutoff,
        memory=memory,
        output_marker=str(certificate["full_support"]["z"]),
    )

    total = exact + intermediate + dense + full_support
    return BinaryBiregularECResult(
        success=bool(total < arb(2) ** (-target_bits)),
        total_bound=total,
        security_bits=-total.log() / arb(2).log(),
        exact_bound=exact,
        intermediate_bound=intermediate,
        dense_bound=dense,
        full_support_bound=full_support,
        largest_exact_r=int(certificate["exact_weights"][largest_index]["r"]),
        largest_exact_term=exact_terms[largest_index],
        intermediate_blocks=len(certificate["outer_blocks"]),
        dense_blocks=len(certificate["central_blocks"]),
    )


def verify_candidate(
    *, k: int = 1_048_575, left_degree: int = 10,
    right_degree: int = 5, cutoff: int = 230_729, memory: int = 15,
    exact_limit: int = 16, intermediate_relative_width: float = 0.2,
    dense_relative_width: float = 0.02, target_bits: int = 20,
    precision_bits: int = 192,
) -> BinaryBiregularECResult:
    """Generate markers in memory and verify them; retained for diagnostics."""
    return verify_certificate(generate_certificate(
        k=k,
        left_degree=left_degree,
        right_degree=right_degree,
        cutoff=cutoff,
        memory=memory,
        exact_limit=exact_limit,
        intermediate_relative_width=intermediate_relative_width,
        dense_relative_width=dense_relative_width,
        target_bits=target_bits,
        precision_bits=precision_bits,
    ))


def print_result(result: BinaryBiregularECResult) -> None:
    print(f"success: {result.success}")
    print(f"total bound: {result.total_bound}")
    print(f"security bits: {result.security_bits}")
    print(f"exact contribution: {result.exact_bound}")
    print(f"intermediate contribution: {result.intermediate_bound}")
    print(f"dense contribution: {result.dense_bound}")
    print(f"full-support contribution: {result.full_support_bound}")
    print(f"largest exact term: r={result.largest_exact_r}, {result.largest_exact_term}")
    print(f"intermediate blocks: {result.intermediate_blocks}")
    print(f"dense blocks: {result.dense_blocks}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    generate = subparsers.add_parser("generate")
    generate.add_argument("--output", type=Path, required=True)
    generate.add_argument("--k", type=int, default=1_048_575)
    generate.add_argument("--left-degree", type=int, default=10)
    generate.add_argument("--right-degree", type=int, default=5)
    generate.add_argument("--cutoff", type=int, default=230_729)
    generate.add_argument("--memory", type=int, default=15)
    generate.add_argument("--exact-limit", type=int, default=16)
    generate.add_argument("--target-bits", type=int, default=20)
    generate.add_argument("--precision-bits", type=int, default=192)
    verify = subparsers.add_parser("verify")
    verify.add_argument("certificate", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.command == "generate":
        certificate = generate_certificate(
            k=args.k,
            left_degree=args.left_degree,
            right_degree=args.right_degree,
            cutoff=args.cutoff,
            memory=args.memory,
            exact_limit=args.exact_limit,
            target_bits=args.target_bits,
            precision_bits=args.precision_bits,
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(certificate, indent=2) + "\n", encoding="utf-8"
        )
        print(f"wrote {args.output}")
        return
    certificate = json.loads(args.certificate.read_text(encoding="utf-8"))
    result = verify_certificate(certificate)
    print_result(result)
    if not result.success:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
