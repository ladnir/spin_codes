#!/usr/bin/env python3
"""Audit the direct fixed-region-weight recurrence used for tail shifts."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from analyze_riffle_bitshuffle_splitstate_regular_bulk import (
    LOG2,
    binomial_transform,
    candidate_epoch_logs,
    load_nonzero_spectrum,
    load_uniform_nonactivation,
    log_choose,
    logsumexp,
    regular_region_log_matrices,
    splitstate_impulse_matrices,
)
from analyze_riffle_body_arbitrary_tail_direct import (
    arbitrary_tail_region_envelopes,
)


def max_finite_difference(left: np.ndarray, right: np.ndarray) -> float:
    finite = np.isfinite(left) & np.isfinite(right)
    if np.any(np.isfinite(left) != np.isfinite(right)):
        return math.inf
    return float(np.max(np.abs(left[finite] - right[finite])))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--activation", type=Path, required=True)
    parser.add_argument("--live-spectrum", type=Path, required=True)
    parser.add_argument("--step-bits", type=int, default=128)
    parser.add_argument("--state-bits", type=int, default=16)
    parser.add_argument("--region-bits", type=int, default=2048)
    parser.add_argument("--constituent-distance", type=int, default=48)
    parser.add_argument("--live-moment-order", type=int, default=3)
    parser.add_argument("--log-surprisal", type=float, default=-7.0)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.region_bits % args.step_bits:
        raise ValueError("step bits must divide region bits")

    nonactivation = load_uniform_nonactivation(args.activation, args.step_bits)
    live_spectrum = load_nonzero_spectrum(
        args.live_spectrum, args.step_bits, args.state_bits
    )
    surprisal = math.exp(args.log_surprisal)
    impulses = splitstate_impulse_matrices(
        z=math.exp(-surprisal),
        step_bits=args.step_bits,
        state_bits=args.state_bits,
        constituent_distance=args.constituent_distance,
        live_moment_order=args.live_moment_order,
        live_model="support-averaged-preaddmul",
        nonactivation=nonactivation,
        live_spectrum=live_spectrum,
    )
    with np.errstate(divide="ignore"):
        impulse_logs = np.log(impulses)
    epochs = args.region_bits // args.step_bits
    fixed = regular_region_log_matrices(
        impulse_logs, args.step_bits, epochs, args.region_bits
    )

    candidate_logs = candidate_epoch_logs(
        impulses, binomial_transform(args.step_bits)
    )
    regular = regular_region_log_matrices(
        candidate_logs, args.step_bits, epochs, args.region_bits
    )
    current = fixed.copy()
    fair_zero_shift_error = 0.0
    direct_binomial_error = 0.0
    checks = ((1, 0), (2, 1), (3, 4), (5, 7), (8, 13))
    direct_rows: dict[tuple[int, int], np.ndarray] = {}
    for body_count, forced_count in checks:
        values = []
        for body_weight in range(body_count + 1):
            values.append(
                fixed[forced_count + body_weight]
                + log_choose(body_count, body_weight)
                - body_count * LOG2
            )
        direct_rows[(body_count, forced_count)] = np.asarray(
            [
                logsumexp(
                    np.asarray([value.reshape(4)[entry] for value in values])
                )
                for entry in range(4)
            ]
        ).reshape(2, 2)

    for body_count in range(1, args.region_bits + 1):
        current = np.logaddexp(current[:-1], current[1:]) - LOG2
        fair_zero_shift_error = max(
            fair_zero_shift_error,
            max_finite_difference(current[0], regular[body_count]),
        )
        for checked_body, forced_count in checks:
            if checked_body == body_count:
                direct_binomial_error = max(
                    direct_binomial_error,
                    max_finite_difference(
                        current[forced_count],
                        direct_rows[(checked_body, forced_count)],
                    ),
                )

    envelopes, maximizers = arbitrary_tail_region_envelopes(
        fixed, args.region_bits - 1
    )
    maximizer_range_ok = all(
        0 <= int(maximizers[b - 1, entry]) <= args.region_bits - b
        for b in range(1, args.region_bits)
        for entry in range(4)
    )
    result = {
        "schema": "riffle-arbitrary-tail-region-recurrence-audit-v1",
        "parameters": {
            "step_bits": args.step_bits,
            "state_bits": args.state_bits,
            "region_bits": args.region_bits,
            "epochs": epochs,
            "log_surprisal": args.log_surprisal,
        },
        "max_log_error_fair_zero_shift_vs_regular": fair_zero_shift_error,
        "max_log_error_recurrence_vs_direct_binomial": direct_binomial_error,
        "maximizer_range_ok": maximizer_range_ok,
        "all_finite_envelopes": bool(np.all(np.isfinite(envelopes))),
        "passed": bool(
            fair_zero_shift_error < 2e-11
            and direct_binomial_error < 2e-11
            and maximizer_range_ok
            and np.all(np.isfinite(envelopes))
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, sort_keys=True))
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
