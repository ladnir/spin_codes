#!/usr/bin/env python3
"""Generate and verify rigorous finite EA first-moment certificates.

Certificate generation uses floating-point optimization only to select useful
values of the Chernoff marker z.  Verification treats those values as fixed
decimal rationals and recomputes every bound with Arb ball arithmetic.  The
validity of a certificate therefore does not depend on the optimizer or on
binary floating-point arithmetic.
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

from expander_bounds import activation_probability, ea_matrix, log_weight_mgf


SCHEMA = "ea-first-moment-v1"


def spectral_radius_float(q: float, z: float) -> float:
    discriminant = (1.0 - q) ** 2 * (1.0 - z) ** 2 + 4.0 * q * q * z
    return 0.5 * ((1.0 - q) * (1.0 + z) + math.sqrt(discriminant))


def optimize_exact_marker(q: float, n: int, cutoff: int) -> float:
    def objective(log_z: float) -> float:
        z = math.exp(log_z)
        return log_weight_mgf(ea_matrix(q, z), n, 0) - cutoff * log_z

    optimum = minimize_scalar(
        objective,
        bounds=(-40.0, 0.0),
        method="bounded",
        options={"xatol": 1e-13, "maxiter": 200},
    )
    return math.exp(float(optimum.x))


def optimize_spectral_marker(q: float, n: int, cutoff: int) -> float:
    def objective(log_z: float) -> float:
        z = math.exp(log_z)
        return (
            0.5 * math.log1p(z)
            + n * math.log(spectral_radius_float(q, z))
            - cutoff * log_z
        )

    optimum = minimize_scalar(
        objective,
        bounds=(-40.0, 0.0),
        method="bounded",
        options={"xatol": 1e-13, "maxiter": 200},
    )
    return math.exp(float(optimum.x))


def decimal_marker(value: float) -> str:
    # Seventeen significant digits recover the selected binary64 value.  The
    # verifier parses this string as an Arb ball containing that decimal value.
    return format(value, ".17g")


def geometric_blocks(start: int, stop: int, growth: float) -> list[tuple[int, int]]:
    blocks: list[tuple[int, int]] = []
    lo = start
    while lo <= stop:
        hi = min(stop, max(lo, math.floor(lo * growth)))
        blocks.append((lo, hi))
        lo = hi + 1
    return blocks


def generate_certificate(
    *,
    k: int,
    n: int,
    cutoff: int,
    row_weight: int,
    target_bits: int,
    exact_limit: int,
    block_growth: float,
    precision_bits: int,
) -> dict[str, Any]:
    if not (1 <= cutoff < n // 2):
        raise ValueError("cutoff must lie in [1,n/2)")
    if not (1 <= exact_limit < k // 2):
        raise ValueError("exact_limit must lie in [1,k/2)")
    if not (0 < row_weight < n // 2):
        raise ValueError("row_weight must induce a density in (0,1/2)")
    if block_growth <= 1.0:
        raise ValueError("block_growth must exceed one")

    density = row_weight / n
    exact: list[dict[str, Any]] = []
    for r in range(1, exact_limit + 1):
        q = activation_probability(density, r)
        exact.append(
            {"r": r, "z": decimal_marker(optimize_exact_marker(q, n, cutoff))}
        )

    blocks: list[dict[str, Any]] = []
    midpoint = k // 2
    for lo, hi in geometric_blocks(exact_limit + 1, midpoint, block_growth):
        q = activation_probability(density, lo)
        blocks.append(
            {
                "lo": lo,
                "hi": hi,
                "count_bound": "entropy-at-hi",
                "z": decimal_marker(optimize_spectral_marker(q, n, cutoff)),
            }
        )

    if midpoint < k:
        lo = midpoint + 1
        q = activation_probability(density, lo)
        blocks.append(
            {
                "lo": lo,
                "hi": k,
                "count_bound": "all-messages",
                "z": decimal_marker(optimize_spectral_marker(q, n, cutoff)),
            }
        )

    return {
        "schema": SCHEMA,
        "parameters": {
            "k": k,
            "n": n,
            "cutoff": cutoff,
            "row_weight": row_weight,
            "target_bits": target_bits,
        },
        "verification": {
            "precision_bits": precision_bits,
            "exact_limit": exact_limit,
            "block_growth": decimal_marker(block_growth),
        },
        "exact_weights": exact,
        "blocks": blocks,
    }


def q_arb(row_weight: int, n: int, r: int) -> arb:
    p = arb(row_weight) / n
    return (arb(1) - (arb(1) - 2 * p) ** r) / 2


def exact_tail_bound_arb(q: arb, z: arb, n: int, cutoff: int) -> arb:
    matrix = arb_mat([[1 - q, q * z], [q, (1 - q) * z]])
    powered = matrix**n
    mgf = powered[0, 0] + powered[0, 1]
    return mgf * z ** (-cutoff)


def spectral_tail_bound_arb(q: arb, z: arb, n: int, cutoff: int) -> arb:
    discriminant = (1 - q) ** 2 * (1 - z) ** 2 + 4 * q * q * z
    eigenvalue = ((1 - q) * (1 + z) + discriminant.sqrt()) / 2
    return (1 + z).sqrt() * eigenvalue**n * z ** (-cutoff)


def entropy_count_bound_arb(k: int, lo: int, hi: int) -> arb:
    if not (1 <= lo <= hi <= k // 2):
        raise ValueError("entropy block must lie in [1,k/2]")
    x = arb(hi) / k
    entropy = -x * x.log() - (1 - x) * (1 - x).log()
    return arb(hi - lo + 1) * (arb(k) * entropy).exp()


def hamming_ball_bound_arb(length: int, cutoff: int) -> arb:
    """Upper-bound a binary Hamming ball using its largest shell.

    For ``cutoff < length/2``, consecutive shells below ``cutoff`` are
    bounded by a geometric series. Arb evaluates the log-gamma expression
    for the largest binomial coefficient with outward rounding.
    """
    if not (0 <= cutoff and 2 * cutoff < length):
        raise ValueError("require 0 <= cutoff < length/2")
    if cutoff == 0:
        return arb(1)
    log_binomial = (
        arb(length + 1).lgamma()
        - arb(cutoff + 1).lgamma()
        - arb(length - cutoff + 1).lgamma()
    )
    geometric_factor = arb(length - cutoff + 1) / (length - 2 * cutoff + 1)
    return log_binomial.exp() * geometric_factor


def validate_coverage(certificate: dict[str, Any]) -> None:
    parameters = certificate["parameters"]
    k = int(parameters["k"])
    expected = 1
    for item in certificate["exact_weights"]:
        r = int(item["r"])
        if r != expected:
            raise ValueError(f"exact-weight coverage breaks at r={expected}")
        expected += 1
    for block in certificate["blocks"]:
        lo, hi = int(block["lo"]), int(block["hi"])
        if lo != expected or hi < lo:
            raise ValueError(f"block coverage breaks at r={expected}")
        expected = hi + 1
    if expected != k + 1:
        raise ValueError(f"certificate stops at r={expected - 1}, expected {k}")


@dataclass(frozen=True)
class VerificationResult:
    success: bool
    total_bound: arb
    target: arb
    security_bits: arb
    exact_bound: arb
    block_bound: arb
    largest_exact_r: int
    largest_exact_term: arb
    largest_block: tuple[int, int]
    largest_block_term: arb


def verify_certificate(certificate: dict[str, Any]) -> VerificationResult:
    if certificate.get("schema") != SCHEMA:
        raise ValueError(f"unsupported certificate schema: {certificate.get('schema')!r}")
    validate_coverage(certificate)

    parameters = certificate["parameters"]
    k = int(parameters["k"])
    n = int(parameters["n"])
    cutoff = int(parameters["cutoff"])
    row_weight = int(parameters["row_weight"])
    target_bits = int(parameters["target_bits"])
    precision_bits = int(certificate["verification"]["precision_bits"])
    if k < 1 or n < 1:
        raise ValueError("k and n must be positive")
    if not (0 <= cutoff < n):
        raise ValueError("cutoff must lie in [0,n)")
    if not (0 < 2 * row_weight < n):
        raise ValueError("row_weight must induce a density in (0,1/2)")
    if target_bits < 0:
        raise ValueError("target_bits must be nonnegative")
    if precision_bits < 64:
        raise ValueError("verification precision must be at least 64 bits")

    ctx.prec = precision_bits
    exact_total = arb(0)
    largest_exact_r = 0
    largest_exact_term = arb(0)
    for item in certificate["exact_weights"]:
        r = int(item["r"])
        z = arb(str(item["z"]))
        if not (z > 0 and z <= 1):
            raise ValueError(f"invalid marker for r={r}")
        q = q_arb(row_weight, n, r)
        term = arb(math.comb(k, r)) * exact_tail_bound_arb(q, z, n, cutoff)
        exact_total += term
        if largest_exact_r == 0 or term > largest_exact_term:
            largest_exact_r = r
            largest_exact_term = term

    block_total = arb(0)
    largest_block = (0, 0)
    largest_block_term = arb(0)
    for block in certificate["blocks"]:
        lo, hi = int(block["lo"]), int(block["hi"])
        z = arb(str(block["z"]))
        if not (z > 0 and z <= 1):
            raise ValueError(f"invalid marker for block [{lo},{hi}]")
        q = q_arb(row_weight, n, lo)
        method = block["count_bound"]
        if method == "entropy-at-hi":
            count_bound = entropy_count_bound_arb(k, lo, hi)
        elif method == "all-messages":
            count_bound = arb(2) ** k
        else:
            raise ValueError(f"unsupported count bound: {method!r}")
        term = count_bound * spectral_tail_bound_arb(q, z, n, cutoff)
        block_total += term
        if largest_block == (0, 0) or term > largest_block_term:
            largest_block = (lo, hi)
            largest_block_term = term

    total = exact_total + block_total
    target = arb(2) ** (-target_bits)
    security_bits = -total.log() / arb(2).log()
    return VerificationResult(
        success=bool(total < target),
        total_bound=total,
        target=target,
        security_bits=security_bits,
        exact_bound=exact_total,
        block_bound=block_total,
        largest_exact_r=largest_exact_r,
        largest_exact_term=largest_exact_term,
        largest_block=largest_block,
        largest_block_term=largest_block_term,
    )


def result_json(result: VerificationResult) -> dict[str, Any]:
    def text(value: arb) -> str:
        return str(value)

    return {
        "success": result.success,
        "total_bound": text(result.total_bound),
        "target": text(result.target),
        "security_bits": text(result.security_bits),
        "exact_bound": text(result.exact_bound),
        "block_bound": text(result.block_bound),
        "largest_exact": {
            "r": result.largest_exact_r,
            "bound": text(result.largest_exact_term),
        },
        "largest_block": {
            "lo": result.largest_block[0],
            "hi": result.largest_block[1],
            "bound": text(result.largest_block_term),
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    generate = subparsers.add_parser("generate")
    generate.add_argument("--output", type=Path, required=True)
    generate.add_argument("--k", type=int, default=2**20)
    generate.add_argument("--rate", type=float, default=0.2)
    generate.add_argument("--delta", type=float, default=0.05)
    generate.add_argument("--row-weight", type=int, required=True)
    generate.add_argument("--target-bits", type=int, default=20)
    generate.add_argument("--exact-limit", type=int, default=128)
    generate.add_argument("--block-growth", type=float, default=1.04)
    generate.add_argument("--precision-bits", type=int, default=192)

    verify = subparsers.add_parser("verify")
    verify.add_argument("certificate", type=Path)

    table = subparsers.add_parser("table")
    table.add_argument("certificates", type=Path, nargs="+")

    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.command == "generate":
        n = int(round(args.k / args.rate))
        cutoff = math.floor(args.delta * n)
        certificate = generate_certificate(
            k=args.k,
            n=n,
            cutoff=cutoff,
            row_weight=args.row_weight,
            target_bits=args.target_bits,
            exact_limit=args.exact_limit,
            block_growth=args.block_growth,
            precision_bits=args.precision_bits,
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(certificate, indent=2) + "\n", encoding="utf-8")
        print(f"wrote {args.output}")
        print(f"exact weights: {len(certificate['exact_weights'])}")
        print(f"spectral blocks: {len(certificate['blocks'])}")
        return

    if args.command == "verify":
        certificate = json.loads(args.certificate.read_text(encoding="utf-8"))
        result = verify_certificate(certificate)
        print(json.dumps(result_json(result), indent=2))
        if not result.success:
            raise SystemExit(1)
        return

    rows: list[dict[str, Any]] = []
    all_success = True
    for path in args.certificates:
        certificate = json.loads(path.read_text(encoding="utf-8"))
        result = verify_certificate(certificate)
        parameters = certificate["parameters"]
        rows.append(
            {
                "certificate": str(path),
                "k": parameters["k"],
                "n": parameters["n"],
                "cutoff": parameters["cutoff"],
                "row_weight": parameters["row_weight"],
                **result_json(result),
            }
        )
        all_success &= result.success
    print(json.dumps(rows, indent=2))
    if not all_success:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
