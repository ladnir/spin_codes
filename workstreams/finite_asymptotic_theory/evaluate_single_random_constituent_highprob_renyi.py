#!/usr/bin/env python3
"""High-probability spectrum and Renyi diagnostic for one random [256,128]."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

import evaluate_ebch128_randomstepconv_g1 as transfer


WORKSTREAM = Path(__file__).resolve().parent
B = 256
K = 128
L = 8192
N = 1 << 21
D = (109 * N + 999) // 1000
MEMORY = 22
LOG2 = math.log(2.0)
DEFAULT_OUTPUT = (
    WORKSTREAM / "single_random_constituent_B256_highprob_renyi_probe.json"
)


def log_choose(total: int, selected: int) -> float:
    return (
        math.lgamma(total + 1)
        - math.lgamma(selected + 1)
        - math.lgamma(total - selected + 1)
    )


def logistic(value: float) -> float:
    if value >= 0.0:
        inverse = math.exp(-min(value, 745.0))
        return 1.0 / (1.0 + inverse)
    direct = math.exp(max(value, -745.0))
    return direct / (1.0 + direct)


def logit(value: float) -> float:
    value = min(1.0 - 1e-12, max(1e-12, value))
    return math.log(value / (1.0 - value))


def spectrum_caps(
    per_shell_failure: float,
    *,
    block_bits: int = B,
    dimension: int = K,
) -> tuple[np.ndarray, dict[str, float]]:
    """Return integer shell caps and their Cantelli/Markov union bound."""

    message_count = (1 << dimension) - 1
    caps = np.zeros(block_bits + 1, dtype=object)
    failure_mp = mp.mpf(0)
    with mp.workdps(max(100, 2 * block_bits)):
        delta = mp.mpf(per_shell_failure)
        denominator = mp.mpf(1 << block_bits)
        for weight in range(1, block_bits + 1):
            shell_size = math.comb(block_bits, weight)
            mean = mp.mpf(message_count) * shell_size / denominator
            if mean <= delta:
                cap = 0
                shell_failure = mean
            else:
                deviation = mp.sqrt(mean * (1 - delta) / delta)
                cap = max(0, int(mp.ceil(mean + deviation)) - 1)
                cap = min(cap, message_count, shell_size)
                actual_deviation = mp.mpf(cap + 1) - mean
                if actual_deviation <= 0:
                    raise ArithmeticError("spectrum cap does not exceed its mean")
                shell_failure = mean / (
                    mean + actual_deviation * actual_deviation
                )
            caps[weight] = cap
            failure_mp += shell_failure
        failure = float(failure_mp)

    rank_failure = sum(
        math.ldexp(1.0, index - block_bits) for index in range(dimension)
    )
    total_mass = sum(int(value) for value in caps)
    return caps, {
        "per_shell_failure_target": per_shell_failure,
        "spectrum_failure_upper": failure,
        "rank_failure_upper": rank_failure,
        "total_event_failure_upper": failure + rank_failure,
        "total_cap_mass_log2": math.log2(total_mass),
        "true_nonzero_mass_log2": math.log2(message_count),
        "cap_mass_excess_bits": math.log2(total_mass / message_count),
    }


def spectrum_logs(caps: np.ndarray) -> np.ndarray:
    result = np.full(B + 1, -math.inf)
    for weight in range(1, B + 1):
        cap = int(caps[weight])
        if cap:
            result[weight] = math.log(cap)
    return result


def outer_renyi_log_moment(
    cap_logs: np.ndarray,
    value_probability: float,
    order: float,
) -> float:
    terms = []
    for weight in range(1, B + 1):
        log_cap = float(cap_logs[weight])
        if not math.isfinite(log_cap):
            continue
        log_reference_shell = (
            log_choose(B, weight)
            + weight * math.log(value_probability)
            + (B - weight) * math.log1p(-value_probability)
        )
        terms.append(order * log_cap + (1.0 - order) * log_reference_shell)
    return float(logsumexp(np.asarray(terms)))


def reference_bad_log_upper(
    *,
    occupation: int,
    candidate_probability: float,
    value_probability: float,
    surprisal: float,
    memory_bits: int,
) -> tuple[float, dict[str, float]]:
    bit_probability = candidate_probability * value_probability
    z = math.exp(-surprisal)
    zero, active = transfer.step_matrices(z, memory_bits)
    mixed = (1.0 - bit_probability) * zero + bit_probability * active
    powered = transfer.log_power(transfer.log_entries(mixed), N)
    log_moment = float(np.logaddexp(powered[0, 0], powered[0, 1]))
    log_conditioning = (
        log_choose(L, occupation)
        + occupation * math.log(candidate_probability)
        + (L - occupation) * math.log1p(-candidate_probability)
    )
    raw = log_moment + D * surprisal - B * log_conditioning
    return min(0.0, raw), {
        "candidate_probability": candidate_probability,
        "value_probability": value_probability,
        "reference_bit_probability": bit_probability,
        "surprisal": surprisal,
        "log_surprisal": math.log(surprisal),
        "log_inner_moment": log_moment,
        "log_region_conditioning_probability": log_conditioning,
        "raw_reference_bad_log_upper": raw,
    }


def objective(
    point: np.ndarray,
    *,
    cap_logs: np.ndarray,
    occupation: int,
    memory_bits: int,
) -> tuple[float, dict[str, float]]:
    candidate_probability = min(
        1.0 - 1e-15, max(1e-15, logistic(float(point[0])))
    )
    value_probability = min(
        1.0 - 1e-15, max(1e-15, logistic(float(point[1])))
    )
    surprisal = math.exp(min(4.0, max(-16.0, float(point[2]))))
    order = 1.0 + math.exp(min(10.0, max(-16.0, float(point[3]))))
    outer_moment = outer_renyi_log_moment(
        cap_logs, value_probability, order
    )
    reference_bad, details = reference_bad_log_upper(
        occupation=occupation,
        candidate_probability=candidate_probability,
        value_probability=value_probability,
        surprisal=surprisal,
        memory_bits=memory_bits,
    )
    conjugate = (order - 1.0) / order
    value = (
        log_choose(L, occupation)
        + occupation * outer_moment / order
        + conjugate * reference_bad
    )
    return value, {
        "holder_order": order,
        "conjugate_weight": conjugate,
        "outer_renyi_log_moment": outer_moment,
        **details,
    }


def optimize_occupation(
    cap_logs: np.ndarray,
    occupation: int,
    memory_bits: int,
) -> dict[str, object]:
    density = occupation / L
    starts = []
    for order in (1.05, 1.2, 2.0, 8.0, 32.0):
        for value_probability in (0.35, 0.5, 0.65):
            candidate_probability = min(0.999999, max(1.0 / L, density))
            bit_probability = candidate_probability * value_probability
            surprisal = max(1e-5, 2.2 * bit_probability)
            starts.append(
                np.asarray(
                    [
                        logit(candidate_probability),
                        logit(value_probability),
                        math.log(surprisal),
                        math.log(order - 1.0),
                    ]
                )
            )
    best = None
    for start in starts:
        result = minimize(
            lambda point: objective(
                point,
                cap_logs=cap_logs,
                occupation=occupation,
                memory_bits=memory_bits,
            )[0],
            start,
            method="Nelder-Mead",
            options={"maxiter": 1600, "xatol": 1e-9, "fatol": 1e-8},
        )
        if best is None or float(result.fun) < float(best.fun):
            best = result
    assert best is not None
    value, details = objective(
        best.x,
        cap_logs=cap_logs,
        occupation=occupation,
        memory_bits=memory_bits,
    )
    return {
        "occupation": occupation,
        "occupation_density": density,
        "log2_upper": value / LOG2,
        "margin_bits": -value / LOG2,
        "optimizer_success": bool(best.success),
        "optimizer_message": str(best.message),
        **details,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--occupations",
        type=int,
        nargs="+",
        default=(1, 2, 3, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 6144, 8192),
    )
    parser.add_argument("--memory-bits", type=int, default=MEMORY)
    parser.add_argument("--per-shell-failure-exponent", type=int, default=51)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if any(not 1 <= q <= L for q in args.occupations):
        parser.error(f"occupations must lie in [1,{L}]")

    caps, event = spectrum_caps(
        math.ldexp(1.0, -args.per_shell_failure_exponent)
    )
    cap_logs = spectrum_logs(caps)
    rows = []
    for occupation in args.occupations:
        row = optimize_occupation(cap_logs, occupation, args.memory_bits)
        rows.append(row)
        print(
            f"Q={occupation},margin={row['margin_bits']:.9f},"
            f"order={row['holder_order']:.6f}",
            flush=True,
        )

    payload = {
        "schema": "single-random-constituent-highprob-renyi-v1",
        "status": "BINARY64_SAMPLED_OCCUPATION_DIAGNOSTIC",
        "parameters": {
            "outer_code": "one uniform binary [256,128] generator",
            "outer_rows": L,
            "output_bits": N,
            "distance_cutoff": D,
            "memory_bits": args.memory_bits,
            "occupations": args.occupations,
        },
        "spectrum_event": event,
        "rows": rows,
        "worst_sampled": min(rows, key=lambda row: float(row["margin_bits"])),
        "limitations": [
            "The occupations are sampled and do not form an all-Q cover.",
            "The optimizer and arithmetic use nearest binary64.",
            "The conditioning and Renyi reduction require an outward verifier.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"event": event, "worst": payload["worst_sampled"]}, indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
