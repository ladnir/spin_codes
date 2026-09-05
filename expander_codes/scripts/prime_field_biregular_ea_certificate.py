#!/usr/bin/env python3
"""Generate and verify Arb certificates for two-sided regular prime-field EA.

Floating-point optimization selects positive coefficient markers.  Verification
parses them as fixed decimals and recomputes every inequality with outward-
rounded Arb arithmetic.  Small message supports use the exact regional
occupancy transfer.  The remaining supports use uniform three-marker blocks.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from flint import arb, arb_mat, ctx

from ea_certificate import decimal_marker
from expander_bounds import log_binom, log_weight_mgf
from prime_field_regular_ea_trace_diagnostic import (
    _optimize_three_markers,
    _optimize_two_markers,
    balanced_occupancy_size_distribution,
    balanced_slot_trace_matrix,
    trace_input_matrices,
    uniform_occupancy_trace_transfers,
)
from regular_ec_certificate import (
    combine_uniform_slice_transfers_arb,
    identity_matrix,
    zero_matrix,
)


SCHEMA = "prime-field-biregular-ea-trace-v1"


def trace_input_matrices_arb(z: arb, v: arb) -> tuple[arb_mat, arb_mat]:
    return (
        arb_mat([[arb(1), arb(0)], [arb(0), z]]),
        arb_mat([[v, z], [v, z]]),
    )


def balanced_occupancy_size_distribution_arb(
    region_length: int, group_size: int, draws: int
) -> list[arb]:
    total_slots = region_length * group_size
    if region_length < 1 or group_size < 1 or not 0 <= draws <= total_slots:
        raise ValueError("invalid balanced occupancy parameters")
    distribution = [arb(1)]
    for used in range(draws):
        following = [arb(0) for _ in range(len(distribution) + 1)]
        remaining = total_slots - used
        for occupied, probability in enumerate(distribution):
            following[occupied] += (
                probability * (group_size * occupied - used) / remaining
            )
            following[occupied + 1] += (
                probability * group_size * (region_length - occupied) / remaining
            )
        distribution = following
    return distribution


def uniform_occupancy_trace_transfers_arb(
    *, region_length: int, max_occupied: int, z: arb, v: arb
) -> list[arb_mat]:
    empty, occupied = trace_input_matrices_arb(z, v)
    result: tuple[int, list[arb_mat]] = (0, [identity_matrix(2)])
    power: tuple[int, list[arb_mat]] = (1, [empty, occupied])
    remaining = region_length
    while remaining:
        if remaining & 1:
            result = combine_uniform_slice_transfers_arb(
                result, power, max_occupied
            )
        remaining >>= 1
        if remaining:
            power = combine_uniform_slice_transfers_arb(
                power, power, max_occupied
            )
    return result[1]


def exact_region_trace_matrix_arb(
    *, region_length: int, group_size: int, message_weight: int,
    z: arb, v: arb,
) -> arb_mat:
    occupancy = balanced_occupancy_size_distribution_arb(
        region_length, group_size, message_weight
    )
    slices = uniform_occupancy_trace_transfers_arb(
        region_length=region_length,
        max_occupied=message_weight,
        z=z,
        v=v,
    )
    region = zero_matrix(2)
    for occupied, probability in enumerate(occupancy):
        if not probability.is_zero():
            region += slices[occupied] * probability
    return region


def exact_markers_float(
    *, prime: int, k: int, n: int, cutoff: int, region_count: int,
    message_weight: int,
) -> dict[str, list[str]]:
    region_length = n // region_count
    group_size = k // region_length
    occupancy = balanced_occupancy_size_distribution(
        region_length, group_size, message_weight
    )

    def region_matrix(z: float, v: float) -> np.ndarray:
        slices = uniform_occupancy_trace_transfers(
            region_length=region_length,
            max_occupied=message_weight,
            z=z,
            v=v,
        )
        region = np.zeros((2, 2), dtype=np.float64)
        for occupied, probability in enumerate(occupancy):
            if probability:
                region += probability * slices[occupied]
        return region

    def base(log_z: float, log_v: float) -> float:
        return (
            log_weight_mgf(
                region_matrix(math.exp(log_z), math.exp(log_v)),
                region_count,
                0,
            )
            - cutoff * log_z
        )

    r = message_weight
    log_s = math.log(prime - 1)

    def structural(point: np.ndarray) -> float:
        log_z, log_v = map(float, point)
        return base(log_z, log_v) - (r - 1) * log_v

    def field(point: np.ndarray) -> float:
        log_z, log_v = map(float, point)
        return -log_s + base(log_z, log_v) - r * log_v

    z_start = max(1e-9, min(0.95, cutoff / n))
    v_start = max(1e-12, min(0.8, r / n))
    lower = -max(120.0, log_s + 12.0)
    _, z_a, v_a = _optimize_two_markers(
        objective=structural,
        z_start=z_start,
        v_start=v_start,
        lower_z=lower,
        lower_v=lower,
    )
    _, z_b, v_b = _optimize_two_markers(
        objective=field,
        z_start=z_start,
        v_start=max(math.exp(-log_s), v_start),
        lower_z=lower,
        lower_v=-log_s,
    )
    v_b = max(v_b, math.nextafter(1.0 / (prime - 1), math.inf))
    return {
        "structural": [decimal_marker(z_a), decimal_marker(v_a)],
        "field": [decimal_marker(z_b), decimal_marker(v_b)],
    }


def exact_marker_schedule_float(
    *, prime: int, k: int, n: int, cutoff: int, region_count: int,
    limit: int,
) -> list[dict[str, Any]]:
    """Optimize geometric anchors and log-interpolate the remaining markers."""
    anchors = [1]
    while anchors[-1] < limit:
        anchors.append(min(limit, 2 * anchors[-1]))
    anchors = sorted(set(anchors) | set(range(10, min(limit, 16) + 1)))
    anchor_markers: dict[int, dict[str, list[str]]] = {}
    for r in anchors:
        anchor_markers[r] = exact_markers_float(
            prime=prime, k=k, n=n, cutoff=cutoff,
            region_count=region_count, message_weight=r,
        )
        print(f"optimized exact anchor r={r}", flush=True)

    schedule: list[dict[str, Any]] = []
    for r in range(1, limit + 1):
        if r in anchor_markers:
            markers = anchor_markers[r]
        else:
            upper = next(anchor for anchor in anchors if anchor > r)
            lower = anchors[anchors.index(upper) - 1]
            fraction = (r - lower) / (upper - lower)
            markers = {}
            for branch in ("structural", "field"):
                values = []
                for index in range(2):
                    low_value = float(anchor_markers[lower][branch][index])
                    high_value = float(anchor_markers[upper][branch][index])
                    interpolated = math.exp(
                        (1 - fraction) * math.log(low_value)
                        + fraction * math.log(high_value)
                    )
                    values.append(decimal_marker(interpolated))
                markers[branch] = values
            field_v = max(
                float(markers["field"][1]),
                math.nextafter(1.0 / (prime - 1), math.inf),
            )
            markers["field"][1] = decimal_marker(field_v)
        schedule.append({"r": r, "markers": markers})
    return schedule


def exact_term_arb(
    *, prime: int, k: int, n: int, cutoff: int, region_count: int,
    message_weight: int, markers: dict[str, list[str]],
) -> arb:
    region_length = n // region_count
    group_size = k // region_length
    r = message_weight

    def branch(values: list[str], field: bool) -> arb:
        z, v = map(arb, values)
        if not (z > 0 and z <= 1 and v > 0 and v <= 1):
            raise ValueError(f"invalid exact markers at r={r}")
        if field and not v >= arb(1) / (prime - 1):
            raise ValueError(f"field marker is below 1/(p-1) at r={r}")
        region = exact_region_trace_matrix_arb(
            region_length=region_length,
            group_size=group_size,
            message_weight=r,
            z=z,
            v=v,
        )
        powered = region**region_count
        mgf = powered[0, 0] + powered[0, 1]
        if field:
            return mgf * z ** (-cutoff) * v ** (-r) / (prime - 1)
        return mgf * z ** (-cutoff) * v ** (-(r - 1))

    return arb(math.comb(k, r)) * (
        branch(markers["structural"], False)
        + branch(markers["field"], True)
    )


def block_markers_float(
    *, prime: int, k: int, n: int, cutoff: int, region_count: int,
    lo: int, hi: int,
) -> tuple[dict[str, list[str]], float]:
    region_length = n // region_count
    group_size = k // region_length
    log_s = math.log(prime - 1)
    common = math.log(hi - lo + 1)

    def base(log_z: float, log_v: float, log_x: float) -> tuple[float, float]:
        z = math.exp(log_z)
        v = math.exp(log_v)
        x = math.exp(log_x)
        log_b = region_count * log_x + log_v
        endpoint_factor = max(
            (1 - region_count) * log_binom(k, lo) - lo * log_b,
            (1 - region_count) * log_binom(k, hi) - hi * log_b,
        )
        value = (
            common
            + endpoint_factor
            + log_weight_mgf(
                balanced_slot_trace_matrix(group_size, x, z, v), n, 0
            )
            - cutoff * log_z
        )
        return value, v

    def structural(point: np.ndarray) -> float:
        log_z, log_v, log_x = map(float, point)
        value, v = base(log_z, log_v, log_x)
        return value + math.log(v)

    def field(point: np.ndarray) -> float:
        log_z, log_v, log_x = map(float, point)
        value, _ = base(log_z, log_v, log_x)
        return value - log_s

    midpoint = (lo + hi) // 2
    x_start = midpoint / max(1.0, k - midpoint)
    z_start = max(1e-9, min(0.95, cutoff / n))
    v_start = max(1e-12, min(0.8, midpoint / n))
    lower = -max(120.0, log_s + 12.0)
    upper_x = math.log(max(2.0, math.expm1(650.0 / group_size)))
    a, z_a, v_a, x_a = _optimize_three_markers(
        objective=structural,
        z_start=z_start,
        v_start=v_start,
        t_start=x_start,
        lower_z=lower,
        lower_v=lower,
        upper_t=upper_x,
    )
    b, z_b, v_b, x_b = _optimize_three_markers(
        objective=field,
        z_start=z_start,
        v_start=max(math.exp(-log_s), v_start),
        t_start=x_start,
        lower_z=lower,
        lower_v=-log_s,
        upper_t=upper_x,
    )
    v_b = max(v_b, math.nextafter(1.0 / (prime - 1), math.inf))
    markers = {
        "structural": [
            decimal_marker(x_a), decimal_marker(z_a), decimal_marker(v_a)
        ],
        "field": [
            decimal_marker(x_b), decimal_marker(z_b), decimal_marker(v_b)
        ],
    }
    return markers, float(np.logaddexp(a, b))


def block_term_arb(
    *, prime: int, k: int, n: int, cutoff: int, region_count: int,
    lo: int, hi: int, markers: dict[str, list[str]],
) -> arb:
    region_length = n // region_count
    group_size = k // region_length
    common = arb(hi - lo + 1)

    def branch(values: list[str], field: bool) -> arb:
        x, z, v = map(arb, values)
        if not (x > 0 and z > 0 and z <= 1 and v > 0 and v <= 1):
            raise ValueError(f"invalid block markers on [{lo},{hi}]")
        if field and not v >= arb(1) / (prime - 1):
            raise ValueError(f"field marker is below 1/(p-1) on [{lo},{hi}]")
        empty, occupied = trace_input_matrices_arb(z, v)
        matrix = empty + ((1 + x) ** group_size - 1) * occupied
        powered = matrix**n
        mgf = powered[0, 0] + powered[0, 1]
        b = x**region_count * v
        lo_factor = arb(math.comb(k, lo)) ** (1 - region_count) * b ** (-lo)
        if lo == hi:
            endpoint_factor = lo_factor
        else:
            hi_factor = (
                arb(math.comb(k, hi)) ** (1 - region_count) * b ** (-hi)
            )
            # The logarithm of the r-dependent factor is convex, so its
            # maximum on an integer interval occurs at an endpoint.  Their
            # sum is a comparison-free Arb upper bound for that maximum.
            endpoint_factor = lo_factor + hi_factor
        result = common * endpoint_factor * mgf * z ** (-cutoff)
        return result / (prime - 1) if field else result * v

    return branch(markers["structural"], False) + branch(
        markers["field"], True
    )


def full_support_markers_float(
    *, prime: int, n: int, cutoff: int, message_weight: int
) -> dict[str, list[str]]:
    r = message_weight
    log_s = math.log(prime - 1)

    def base(log_z: float, log_v: float) -> float:
        z, v = math.exp(log_z), math.exp(log_v)
        _, occupied = trace_input_matrices(z, v)
        return log_weight_mgf(occupied, n, 0) - cutoff * log_z

    def structural(point: np.ndarray) -> float:
        log_z, log_v = map(float, point)
        return base(log_z, log_v) - (r - 1) * log_v

    def field(point: np.ndarray) -> float:
        log_z, log_v = map(float, point)
        return -log_s + base(log_z, log_v) - r * log_v

    z_start = max(1e-9, min(0.95, cutoff / n))
    v_start = max(1e-12, min(0.8, r / n))
    lower = -max(120.0, log_s + 12.0)
    _, z_a, v_a = _optimize_two_markers(
        objective=structural,
        z_start=z_start,
        v_start=v_start,
        lower_z=lower,
        lower_v=lower,
    )
    _, z_b, v_b = _optimize_two_markers(
        objective=field,
        z_start=z_start,
        v_start=max(math.exp(-log_s), v_start),
        lower_z=lower,
        lower_v=-log_s,
    )
    v_b = max(v_b, math.nextafter(1.0 / (prime - 1), math.inf))
    return {
        "structural": [decimal_marker(z_a), decimal_marker(v_a)],
        "field": [decimal_marker(z_b), decimal_marker(v_b)],
    }


def full_support_term_arb(
    *, prime: int, n: int, cutoff: int, message_weight: int,
    markers: dict[str, list[str]],
) -> arb:
    r = message_weight

    def branch(values: list[str], field: bool) -> arb:
        z, v = map(arb, values)
        if field and not v >= arb(1) / (prime - 1):
            raise ValueError("full-support field marker is below 1/(p-1)")
        _, occupied = trace_input_matrices_arb(z, v)
        powered = occupied**n
        mgf = powered[0, 0] + powered[0, 1]
        if field:
            return mgf * z ** (-cutoff) * v ** (-r) / (prime - 1)
        return mgf * z ** (-cutoff) * v ** (-(r - 1))

    return branch(markers["structural"], False) + branch(
        markers["field"], True
    )


def generate_blocks(
    *, prime: int, k: int, n: int, cutoff: int, region_count: int,
    start: int, stop: int, target_log2: float,
) -> list[dict[str, Any]]:
    target = target_log2 * math.log(2.0)
    blocks: list[dict[str, Any]] = []
    lo = start
    while lo <= stop:
        markers, value = block_markers_float(
            prime=prime, k=k, n=n, cutoff=cutoff,
            region_count=region_count, lo=lo, hi=lo,
        )
        if value > target:
            raise ValueError(
                f"singleton saddle block r={lo} has log2 bound "
                f"{value / math.log(2):.6f}"
            )
        hi = lo
        step = 1
        while hi < stop:
            candidate = min(stop, hi + step)
            candidate_markers, candidate_value = block_markers_float(
                prime=prime, k=k, n=n, cutoff=cutoff,
                region_count=region_count, lo=lo, hi=candidate,
            )
            if candidate_value > target:
                break
            hi, markers, value = candidate, candidate_markers, candidate_value
            step *= 2
        low, high = hi + 1, min(stop, hi + step - 1)
        while low <= high:
            middle = (low + high) // 2
            candidate_markers, candidate_value = block_markers_float(
                prime=prime, k=k, n=n, cutoff=cutoff,
                region_count=region_count, lo=lo, hi=middle,
            )
            if candidate_value <= target:
                hi, markers, value = middle, candidate_markers, candidate_value
                low = middle + 1
            else:
                high = middle - 1
        blocks.append({"lo": lo, "hi": hi, "markers": markers})
        print(
            f"block [{lo},{hi}] selected log2={value / math.log(2):.6f}",
            flush=True,
        )
        lo = hi + 1
    return blocks


def generate_certificate(
    *, prime: int, k: int, n: int, cutoff: int, region_count: int,
    exact_limit: int, target_bits: int, target_block_log2: float,
    precision_bits: int,
) -> dict[str, Any]:
    if n % region_count:
        raise ValueError("region_count must divide n")
    region_length = n // region_count
    if k % region_length:
        raise ValueError("region_length must divide k")
    exact = exact_marker_schedule_float(
        prime=prime, k=k, n=n, cutoff=cutoff,
        region_count=region_count, limit=exact_limit,
    )
    blocks = generate_blocks(
        prime=prime, k=k, n=n, cutoff=cutoff,
        region_count=region_count, start=exact_limit + 1, stop=k - 1,
        target_log2=target_block_log2,
    )
    return {
        "schema": SCHEMA,
        "parameters": {
            "prime": str(prime), "k": k, "n": n, "cutoff": cutoff,
            "region_count": region_count, "target_bits": target_bits,
        },
        "verification": {
            "precision_bits": precision_bits, "exact_limit": exact_limit,
            "target_block_log2": target_block_log2,
        },
        "exact_weights": exact,
        "blocks": blocks,
        "full_support": {
            "r": k,
            "markers": full_support_markers_float(
                prime=prime, n=n, cutoff=cutoff, message_weight=k
            ),
        },
    }


def validate_coverage(certificate: dict[str, Any]) -> None:
    k = int(certificate["parameters"]["k"])
    expected = 1
    for item in certificate["exact_weights"]:
        if int(item["r"]) != expected:
            raise ValueError(f"exact coverage breaks at r={expected}")
        expected += 1
    for block in certificate["blocks"]:
        lo, hi = int(block["lo"]), int(block["hi"])
        if lo != expected or hi < lo:
            raise ValueError(f"block coverage breaks at r={expected}")
        expected = hi + 1
    if expected != k:
        raise ValueError(f"block coverage stops at r={expected - 1}")
    if int(certificate["full_support"]["r"]) != k:
        raise ValueError("full-support entry is missing")


@dataclass(frozen=True)
class VerificationResult:
    success: bool
    total_bound: arb
    security_bits: arb
    exact_bound: arb
    block_bound: arb
    full_support_bound: arb
    largest_exact_r: int
    largest_exact_term: arb
    largest_block: tuple[int, int]
    largest_block_term: arb


def verify_certificate(certificate: dict[str, Any]) -> VerificationResult:
    if certificate.get("schema") != SCHEMA:
        raise ValueError("unsupported certificate schema")
    validate_coverage(certificate)
    parameters = certificate["parameters"]
    prime = int(parameters["prime"])
    k, n = int(parameters["k"]), int(parameters["n"])
    cutoff = int(parameters["cutoff"])
    region_count = int(parameters["region_count"])
    target_bits = int(parameters["target_bits"])
    ctx.prec = int(certificate["verification"]["precision_bits"])

    exact_total = arb(0)
    largest_exact_r, largest_exact_term = 0, arb(0)
    for item in certificate["exact_weights"]:
        r = int(item["r"])
        term = exact_term_arb(
            prime=prime, k=k, n=n, cutoff=cutoff,
            region_count=region_count, message_weight=r,
            markers=item["markers"],
        )
        exact_total += term
        if largest_exact_r == 0 or term > largest_exact_term:
            largest_exact_r, largest_exact_term = r, term
        print(f"verified exact r={r}", flush=True)

    block_total = arb(0)
    largest_block, largest_block_term = (0, 0), arb(0)
    for item in certificate["blocks"]:
        lo, hi = int(item["lo"]), int(item["hi"])
        term = block_term_arb(
            prime=prime, k=k, n=n, cutoff=cutoff,
            region_count=region_count, lo=lo, hi=hi,
            markers=item["markers"],
        )
        block_total += term
        if largest_block == (0, 0) or term > largest_block_term:
            largest_block, largest_block_term = (lo, hi), term
        print(f"verified block [{lo},{hi}]", flush=True)

    full = certificate["full_support"]
    full_term = full_support_term_arb(
        prime=prime, n=n, cutoff=cutoff, message_weight=k,
        markers=full["markers"],
    )
    total = exact_total + block_total + full_term
    security_bits = -total.log() / arb(2).log()
    return VerificationResult(
        success=bool(total < arb(2) ** (-target_bits)),
        total_bound=total,
        security_bits=security_bits,
        exact_bound=exact_total,
        block_bound=block_total,
        full_support_bound=full_term,
        largest_exact_r=largest_exact_r,
        largest_exact_term=largest_exact_term,
        largest_block=largest_block,
        largest_block_term=largest_block_term,
    )


def result_json(result: VerificationResult) -> dict[str, Any]:
    return {
        "success": result.success,
        "total_bound": str(result.total_bound),
        "security_bits": str(result.security_bits),
        "exact_bound": str(result.exact_bound),
        "block_bound": str(result.block_bound),
        "full_support_bound": str(result.full_support_bound),
        "largest_exact": {
            "r": result.largest_exact_r,
            "bound": str(result.largest_exact_term),
        },
        "largest_block": {
            "lo": result.largest_block[0], "hi": result.largest_block[1],
            "bound": str(result.largest_block_term),
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    generate = subparsers.add_parser("generate")
    generate.add_argument("--output", type=Path, required=True)
    generate.add_argument(
        "--prime", type=int,
        default=170141183460469231731687303715884105727,
    )
    generate.add_argument("--k", type=int, default=1_048_575)
    generate.add_argument("--n", type=int, default=2_097_150)
    generate.add_argument("--cutoff", type=int, default=1_032_064)
    generate.add_argument("--regions", type=int, default=30)
    generate.add_argument("--exact-limit", type=int, default=64)
    generate.add_argument("--target-bits", type=int, default=20)
    generate.add_argument("--target-block-log2", type=float, default=-32.0)
    generate.add_argument("--precision-bits", type=int, default=192)
    verify = subparsers.add_parser("verify")
    verify.add_argument("certificate", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.command == "generate":
        certificate = generate_certificate(
            prime=args.prime, k=args.k, n=args.n, cutoff=args.cutoff,
            region_count=args.regions, exact_limit=args.exact_limit,
            target_bits=args.target_bits,
            target_block_log2=args.target_block_log2,
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
    print(json.dumps(result_json(result), indent=2))
    if not result.success:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
