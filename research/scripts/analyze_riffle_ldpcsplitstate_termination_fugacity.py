#!/usr/bin/env python3
"""Measure termination occupation in selected shared-group profiles.

Multiply every U-to-Z epoch transfer by exp(y).  The logarithmic derivative
of the complete moment at y=0 is the expected number of termination
transitions under the Chernoff-tilted path measure.  The second derivative is
the corresponding variance.

The script uses centered finite differences around the optimized all-quad
occupation-128 point.  It is a floating-point diagnostic, not a certificate.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from analyze_riffle_ldpcsplitstate_occupation_ladder import DEFAULT_ACTIVATION
from analyze_riffle_ldpcsplitstate_shared_groups import evaluate_profile_point


DEFAULT_QUAD_RECEIPT = Path(
    "constructions/riffle_ldpcsplitstate_g4_t256_s64/"
    "receipts/shared_group_occupation128_quad_optimized.json"
)
DEFAULT_OUTPUT = Path(
    "constructions/riffle_ldpcsplitstate_g4_t256_s64/"
    "receipts/termination_fugacity_occupation128.json"
)


def parse_profiles(text: str) -> list[tuple[int, int, int, int]]:
    result = []
    for item in text.split(";"):
        profile = tuple(int(value) for value in item.split(":"))
        if len(profile) != 4:
            raise ValueError(f"invalid profile: {item}")
        if sum((index + 1) * count for index, count in enumerate(profile)) != 128:
            raise ValueError(f"profile does not have occupation 128: {item}")
        result.append(profile)
    return result


def quad_point(path: Path) -> tuple[float, float]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    for row in payload["occupation_rows"][0]["profiles"]:
        if row["profile_n1_n2_n3_n4"] == [0, 0, 0, 32]:
            return float(row["log_surprisal"]), float(row["placement_logit"])
    raise ValueError("optimized all-quad row not found")


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    activation = json.loads(args.activation.read_text(encoding="utf-8"))
    activation_upper = [0.0] * 129
    activation_upper[0] = 1.0
    for item in activation["by_total_weight"]:
        weight = int(item["total_weight"])
        if weight > 128:
            break
        activation_upper[weight] = float(
            item["maximum_distinct_conditioned_upper_bound"]
        )

    log_surprisal, placement_logit = quad_point(args.quad_receipt)
    target_distance = math.floor(args.relative_distance * (1 << 21))
    rows = []
    for profile in parse_profiles(args.profiles):
        values = {}
        center_row = None
        for y in (-args.epsilon, 0.0, args.epsilon):
            row = evaluate_profile_point(
                profile=profile,
                log_surprisal=log_surprisal,
                placement_logit=placement_logit,
                activation_upper=activation_upper,
                constituent_distance=args.constituent_distance,
                live_moment_order=args.live_moment_order,
                target_distance=target_distance,
                termination_scale=math.exp(y),
            )
            values[y] = float(row["log2_inner_moment"])
            if y == 0.0:
                center_row = row
        no_termination = evaluate_profile_point(
            profile=profile,
            log_surprisal=log_surprisal,
            placement_logit=placement_logit,
            activation_upper=activation_upper,
            constituent_distance=args.constituent_distance,
            live_moment_order=args.live_moment_order,
            target_distance=target_distance,
            termination_scale=0.0,
        )
        minus = values[-args.epsilon]
        center = values[0.0]
        plus = values[args.epsilon]
        expected = math.log(2.0) * (plus - minus) / (2.0 * args.epsilon)
        variance = (
            math.log(2.0)
            * (plus - 2.0 * center + minus)
            / (args.epsilon**2)
        )
        rows.append(
            {
                "profile_n1_n2_n3_n4": list(profile),
                "active_packet_groups": sum(profile),
                "full_log2_inner_moment": center,
                "no_termination_log2_inner_moment": no_termination[
                    "log2_inner_moment"
                ],
                "termination_log2_inflation": center
                - float(no_termination["log2_inner_moment"]),
                "tilted_expected_terminations": expected,
                "tilted_termination_variance": variance,
                "region_matrix": center_row["region_matrix"],
            }
        )
        print(
            f"profile,{list(profile)},groups,{sum(profile)},"
            f"inflation,{rows[-1]['termination_log2_inflation']:.6f},"
            f"expected_terminations,{expected:.6f},variance,{variance:.6f}",
            flush=True,
        )
    return {
        "schema": "riffle-ldpcsplitstate-termination-fugacity-v1",
        "candidate": "Riffle LDPCSplitState g=4 t=256 s=64",
        "parameters": {
            "occupation": 128,
            "log_surprisal": log_surprisal,
            "placement_logit": placement_logit,
            "epsilon": args.epsilon,
            "constituent_distance": args.constituent_distance,
            "live_moment_order": args.live_moment_order,
        },
        "rows": rows,
        "scope": (
            "Centered finite-difference diagnostic under the Chernoff-tilted "
            "state-path measure. Expected termination counts are not encoder "
            "probabilities under the unweighted setup distribution."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--activation", type=Path, default=DEFAULT_ACTIVATION)
    parser.add_argument("--quad-receipt", type=Path, default=DEFAULT_QUAD_RECEIPT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--profiles",
        default="128:0:0:0;2:0:42:0;0:32:0:16;0:0:0:32",
    )
    parser.add_argument("--epsilon", type=float, default=1e-3)
    parser.add_argument("--relative-distance", type=float, default=0.09)
    parser.add_argument("--constituent-distance", type=int, default=40)
    parser.add_argument(
        "--live-moment-order", type=int, choices=(1, 2, 3), default=3
    )
    args = parser.parse_args()
    payload = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
