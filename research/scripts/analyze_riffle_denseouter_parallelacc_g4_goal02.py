#!/usr/bin/env python3
"""Optimize the Goal 02 joint-lane moment bound.

The script prints a JSON receipt. It does not modify the workspace.
"""

from __future__ import annotations

import argparse
import json
import math

from scipy.optimize import minimize
from scipy.special import gammaln


def log_binom(n: int, h: int) -> float:
    return float(gammaln(n + 1) - gammaln(h + 1) - gammaln(n - h + 1))


def log_one_minus_exp_neg(x: float) -> float:
    """Return log(1-exp(-x)) for x>0 without cancellation."""
    if x <= math.log(2.0):
        return math.log(-math.expm1(-x))
    return math.log1p(-math.exp(-x))


def feasible_final_weights(g: int, w: int, h: int) -> list[int]:
    if h == 1:
        return [w]
    return [q for q in range(min(g, w) + 1) if (q - w) % 2 == 0]


def extremal_state_counts(g: int, w: int, h: int, q: int) -> tuple[int, list[int]]:
    """Relax the path to the convex-extremal state-weight histogram."""
    edge_lower = (w + q + 1) // 2
    terminal_nonzero = int(q > 0)
    capacity_needed = max(0, edge_lower - q)
    remaining_needed = (capacity_needed + g - 1) // g
    nonzero_states = max((h + 1) // 2, terminal_nonzero + remaining_needed)
    if nonzero_states > h:
        raise ValueError("infeasible nonzero-state count")
    remaining_nonzero = nonzero_states - terminal_nonzero
    total = max(edge_lower, q + remaining_nonzero)
    remaining = total - q
    if remaining_nonzero == 0:
        if remaining != 0:
            raise ValueError("infeasible terminal state")
        full = residual = 0
    else:
        extra = remaining - remaining_nonzero
        full, residual = divmod(extra, g - 1)
        if full > remaining_nonzero or (full == remaining_nonzero and residual):
            raise ValueError("infeasible relaxed state weight")
    counts = [0] * (g + 1)
    counts[q] += 1
    counts[g] += full
    if residual:
        counts[1 + residual] += 1
    counts[1] += remaining_nonzero - full - int(bool(residual))
    counts[0] += h - sum(counts)
    return total, counts


def optimize_profile_q(
    n: int, g: int, d: int, w: int, h: int, q: int
) -> dict[str, float | int | bool | list[int]]:
    total, counts = extremal_state_counts(g, w, h, q)
    lb = log_binom(n, h)

    def objective(xy: object) -> float:
        u, scaled_v = xy  # type: ignore[misc]
        u = float(u)
        v = float(scaled_v) / n
        value = (d - total) * u + (n - h) * v - lb
        for a, count in enumerate(counts):
            value -= count * log_one_minus_exp_neg(v + a * u)
        return value

    def gradient(xy: object) -> list[float]:
        u, scaled_v = xy  # type: ignore[misc]
        u = float(u)
        v = float(scaled_v) / n
        grad_u = float(d - total)
        grad_v = float(n - h)
        for a, count in enumerate(counts):
            reciprocal = 1.0 / math.expm1(v + a * u)
            grad_u -= count * a * reciprocal
            grad_v -= count * reciprocal
        return [grad_u, grad_v / n]

    start_u = 1.0 / max(1.0, d - total)
    start_scaled_v = max(1.0, float(h))
    result = minimize(
        objective,
        (start_u, start_scaled_v),
        jac=gradient,
        method="L-BFGS-B",
        bounds=((1e-12, 100.0), (1e-9, float(n))),
        options={"ftol": 1e-13, "gtol": 1e-10, "maxiter": 1000},
    )

    u = float(result.x[0])
    v = float(result.x[1]) / n
    log_bound = min(0.0, float(result.fun))
    return {
        "w": w,
        "h": h,
        "final_weight": q,
        "state_weight_sum": total,
        "extremal_state_counts": counts,
        "z": math.exp(-u),
        "tau": math.exp(-v),
        "log2_bound": log_bound / math.log(2.0),
        "root": math.exp(log_bound / w),
        "optimizer_success": bool(result.success),
    }


def optimize_profile(n: int, g: int, d: int, w: int, h: int) -> dict[str, object]:
    candidates = [
        optimize_profile_q(n, g, d, w, h, q)
        for q in feasible_final_weights(g, w, h)
    ]
    worst = max(candidates, key=lambda item: float(item["root"]))
    result = dict(worst)
    result["final_weight_candidates"] = len(candidates)
    return result


def exact_ordered_two_packet_count(
    n: int, d: int, first_weight: int, terminal_weight: int
) -> int:
    """Count failing supports for one fixed order of two active packets."""
    count = 0
    max_gap = min(n - 1, d // first_weight)
    for gap in range(1, max_gap + 1):
        remaining = d - first_weight * gap
        if terminal_weight == 0:
            count += n - gap
        else:
            count += min(n - gap, remaining // terminal_weight)
    return count


def exact_small_support_profile(
    n: int, g: int, d: int, w: int, h: int
) -> dict[str, object]:
    """Return the exact worst fixed-word probability for H in {1,2}."""
    if h == 1:
        if not 1 <= w <= g:
            raise ValueError("infeasible one-packet weight")
        probability = min(n, d // w) / n
        return {
            "w": w,
            "h": h,
            "method": "exact-support-count",
            "packet_values": [next(x for x in range(1, 1 << g) if x.bit_count() == w)],
            "log2_bound": math.log2(probability) if probability else -math.inf,
            "root": probability ** (1.0 / w),
            "optimizer_success": True,
        }

    if h != 2:
        raise ValueError("exact small-support counter only covers H in {1,2}")

    denominator = math.comb(n, 2)
    worst: dict[str, object] | None = None
    for x in range(1, 1 << g):
        for y in range(x, 1 << g):
            if x.bit_count() + y.bit_count() != w:
                continue
            terminal_weight = (x ^ y).bit_count()
            forward = exact_ordered_two_packet_count(
                n, d, x.bit_count(), terminal_weight
            )
            if x == y:
                numerator = float(forward)
            else:
                reverse = exact_ordered_two_packet_count(
                    n, d, y.bit_count(), terminal_weight
                )
                numerator = (forward + reverse) / 2.0
            probability = numerator / denominator
            candidate = {
                "w": w,
                "h": h,
                "method": "exact-support-and-order-count",
                "packet_values": [x, y],
                "packet_weights": [x.bit_count(), y.bit_count()],
                "terminal_weight": terminal_weight,
                "log2_bound": math.log2(probability) if probability else -math.inf,
                "root": probability ** (1.0 / w),
                "optimizer_success": True,
            }
            if worst is None or candidate["root"] > worst["root"]:
                worst = candidate

    if worst is None:
        raise ValueError("infeasible two-packet weight")
    return worst


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--packets", type=int, default=524_352)
    parser.add_argument("--width", type=int, default=4)
    parser.add_argument("--distance", type=int, default=40)
    parser.add_argument("--max-weight", type=int)
    parser.add_argument("--exact-small-support", action="store_true")
    parser.add_argument("--full", action="store_true")
    args = parser.parse_args()

    profiles = []
    worst = None
    max_weight = args.max_weight or 2 * args.distance
    for w in range(1, max_weight + 1):
        h_min = (w + args.width - 1) // args.width
        for h in range(h_min, w + 1):
            if args.exact_small_support and h <= 2:
                profile = exact_small_support_profile(
                    args.packets, args.width, args.distance, w, h
                )
            else:
                profile = optimize_profile(
                    args.packets, args.width, args.distance, w, h
                )
            profiles.append(profile)
            if worst is None or profile["root"] > worst["root"]:
                worst = profile

    assert worst is not None
    receipt = {
        "schema": (
            "riffle-denseouter-parallelacc-g4-goal03-hybrid-v1"
            if args.exact_small_support
            else "riffle-denseouter-parallelacc-g4-goal02-v1"
        ),
        "parameters": {
            "packets": args.packets,
            "packet_width": args.width,
            "binary_length": args.packets * args.width,
            "distance_threshold": args.distance,
            "exact_small_support": args.exact_small_support,
        },
        "profile_count": len(profiles),
        "worst_profile": worst,
        "uniform_rho": worst["root"],
        "uniform_rho_log2": math.log2(worst["root"]),
        "failed_optimizers": sum(not p["optimizer_success"] for p in profiles),
    }
    if args.full:
        receipt["profiles"] = profiles
    else:
        receipt["top_profiles"] = sorted(
            profiles, key=lambda p: p["root"], reverse=True
        )[:20]
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
