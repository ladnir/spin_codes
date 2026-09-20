#!/usr/bin/env python3
"""All-message pivot-deletion diagnostic for fixed BA plus RandomStepConv.

The fixed outer and route are the setup from the outward Toeplitz
certificate. In each local row, a deterministic greedy algorithm selects one
full-rank column basis using an independent fixed priority order. Uniform
messages give independent fair values at these basis coordinates. Set every
other input coordinate to zero. This deletes input ones, so RandomStepConv
monotonicity makes the low-output moment no smaller.

The three-state recurrence sums all nonzero pivot assignments without
subtracting the zero word. Nearest binary64 log arithmetic is diagnostic.
"""

from __future__ import annotations

from array import array
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np


WORKSTREAM = Path(__file__).resolve().parent
sys.path.insert(0, str(WORKSTREAM))

import certify_ba240_repeated_toeplitz_prefix_outward as frozen  # noqa: E402


DEFAULT_OUTPUT = WORKSTREAM / "fixed_ba240_randomstepconv_pivot_s30_d11.json"
TOEPLITZ_RECEIPT = WORKSTREAM / "ba3_B240_repeated_toeplitz_prefix_outward.json"
LOG2 = math.log(2.0)
BASIS_SEED = 0x52414E4442415349  # ASCII bytes "RANDBASI", read big-endian.


def little_u16(values: list[int]) -> bytes:
    encoded = array("H", values)
    if sys.byteorder != "little":
        encoded.byteswap()
    return encoded.tobytes()


def rebuild_pivot_flags(basis_mode: str) -> tuple[bytearray, dict[str, object]]:
    rng = frozen.SplitMix64(frozen.SETUP_SEED)
    basis_rng = frozen.SplitMix64(BASIS_SEED)
    columns = frozen.golay_direct_sum_columns()
    accumulator_hashes = []
    for _ in range(2):
        permutation = rng.permutation(frozen.B)
        accumulator_hashes.append(
            hashlib.sha256(little_u16(permutation)).hexdigest()
        )
        columns = frozen.accumulate(columns, permutation)

    increments = np.zeros((frozen.B, frozen.L), dtype=np.uint8)
    local_hash = hashlib.sha256()
    for row in range(frozen.L):
        permutation = rng.permutation(frozen.B)
        local_hash.update(little_u16(permutation))
        routed_columns = [columns[coordinate] for coordinate in permutation]
        if basis_mode == "prefix":
            priority = list(range(frozen.B))
        elif basis_mode == "randomized":
            priority = basis_rng.permutation(frozen.B)
        else:
            raise ValueError(f"unknown basis mode: {basis_mode}")
        basis = [0] * frozen.LOCAL_DIMENSION
        rank = 0
        for region in priority:
            if frozen.add_column(basis, routed_columns[region]):
                increments[region, row] = 1
                rank += 1
        if rank != frozen.LOCAL_DIMENSION:
            raise ArithmeticError("a local routed generator lost rank")

    flags = bytearray()
    route_hash = hashlib.sha256()
    pivots_by_region = []
    for region in range(frozen.B):
        permutation = rng.permutation(frozen.L)
        route_hash.update(little_u16(permutation))
        count = 0
        for row in permutation:
            flag = int(increments[region, row])
            flags.append(flag)
            count += flag
        pivots_by_region.append(count)

    if len(flags) != frozen.N:
        raise ArithmeticError("the routed pivot stream has the wrong length")
    if sum(flags) != frozen.PARENT_DIMENSION:
        raise ArithmeticError("the routed pivot stream has the wrong rank")

    reference = json.loads(TOEPLITZ_RECEIPT.read_text(encoding="utf-8"))[
        "fixed_setup"
    ]
    observed = {
        "accumulator_permutation_sha256": accumulator_hashes,
        "local_coordinate_permutations_sha256": local_hash.hexdigest(),
        "region_row_permutations_sha256": route_hash.hexdigest(),
    }
    for key, value in observed.items():
        if value != reference[key]:
            raise ArithmeticError(f"frozen setup hash mismatch for {key}")

    return flags, {
        **observed,
        "pivot_flags_sha256": hashlib.sha256(flags).hexdigest(),
        "pivot_count": sum(flags),
        "pivots_by_region": pivots_by_region,
        "basis_mode": basis_mode,
        "basis_seed_hex": (
            None if basis_mode == "prefix" else f"0x{BASIS_SEED:016x}"
        ),
        "basis_splitmix64_words_consumed": basis_rng.words_consumed,
        "basis_splitmix64_rejected_words": basis_rng.rejected_words,
        "splitmix64_words_consumed": rng.words_consumed,
        "splitmix64_rejected_words": rng.rejected_words,
    }


def logadd(*values: float) -> float:
    maximum = max(values)
    if maximum == -math.inf:
        return maximum
    return maximum + math.log(sum(math.exp(value - maximum) for value in values))


def evaluate_tilt(flags: bytearray, memory_bits: int, u: float) -> dict[str, float]:
    surprisal = math.exp(u)
    z = math.exp(-surprisal)
    q = math.ldexp(1.0, -memory_bits)
    b = (1.0 + z) / 2.0
    alpha = q * b
    beta = (1.0 - q) * b
    log_alpha = math.log(alpha)
    log_beta = math.log(beta)
    log_one_plus_alpha = math.log1p(alpha)
    log_two_alpha = math.log(2.0 * alpha)
    log_two_beta = math.log(2.0 * beta)

    # The inactive state represents the unique all-zero pivot prefix. Its
    # logarithmic mass is always zero and is not stored. The other entries
    # sum moments over active assignments ending in zero or live inner state.
    active_zero = -math.inf
    active_live = -math.inf
    for flag in flags:
        if flag:
            next_zero = logadd(
                log_alpha,
                active_zero + log_one_plus_alpha,
                active_live + log_two_alpha,
            )
            next_live = logadd(
                log_beta,
                active_zero + log_beta,
                active_live + log_two_beta,
            )
        else:
            next_zero = logadd(active_zero, active_live + log_alpha)
            next_live = active_live + log_beta
        active_zero, active_live = next_zero, next_live

    moment = logadd(active_zero, active_live)
    log_bad_upper = moment + frozen.DISTANCE_CUTOFF * surprisal
    return {
        "log_surprisal": u,
        "surprisal": surprisal,
        "z": z,
        "log2_nonzero_moment": moment / LOG2,
        "log2_expected_bad_upper": log_bad_upper / LOG2,
        "margin_bits": -log_bad_upper / LOG2,
    }


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    flags, setup = rebuild_pivot_flags(args.basis_mode)
    grid_count = int(
        math.floor((args.grid_max - args.grid_min) / args.grid_step + 0.5)
    ) + 1
    rows = []
    for index in range(grid_count):
        u = args.grid_min + index * args.grid_step
        row = evaluate_tilt(flags, args.memory_bits, u)
        rows.append(row)
        print(
            f"tilt,{index + 1},{grid_count},u,{u:.6f},"
            f"margin,{row['margin_bits']:.6f}",
            flush=True,
        )
    best = max(rows, key=lambda row: float(row["margin_bits"]))
    return {
        "schema": "fixed-ba240-randomstepconv-pivot-deletion-v1",
        "status": "BINARY64_DIAGNOSTIC_EXACT_PIVOT_TRANSFER",
        "claim": {
            "outer": (
                "the fixed repeated Golay--BA-3 outer and fixed route from "
                "the Toeplitz certificate"
            ),
            "inner": (
                "independent uniform linear (M+1)-by-(M+1) maps at every "
                "bit position, sampled once and shared by all messages"
            ),
            "all_nonzero_parent_messages_covered": True,
            "best": best,
            "closes_40_bits": float(best["margin_bits"]) > 40.0,
        },
        "parameters": {
            "memory_bits": args.memory_bits,
            "outer_bits": frozen.B,
            "outer_dimension": frozen.LOCAL_DIMENSION,
            "outer_rows": frozen.L,
            "output_bits": frozen.N,
            "parent_dimension": frozen.PARENT_DIMENSION,
            "target_dimension_after_zero_shortening": frozen.TARGET_DIMENSION,
            "distance_cutoff": frozen.DISTANCE_CUTOFF,
        },
        "fixed_setup": setup,
        "bound": {
            "pivot_rule": (
                "retain one full-rank column basis in each routed local row; "
                "basis_mode specifies the deterministic priority order"
            ),
            "deleted_coordinates": frozen.N - frozen.PARENT_DIMENSION,
            "monotonicity": (
                "For every fixed input word and 0<z<1, changing an input one "
                "to zero cannot decrease the RandomStepConv output moment."
            ),
            "message_sum": (
                "The inactive/active-zero/active-live recurrence sums every "
                "nonzero assignment to the independent pivot values and "
                "never includes or subtracts the all-zero message."
            ),
        },
        "grid": {
            "minimum": args.grid_min,
            "maximum": args.grid_max,
            "step": args.grid_step,
            "count": grid_count,
        },
        "tilt_rows": rows,
        "limitations": [
            "The logarithmic recurrence uses nearest binary64 arithmetic.",
            "The result is not yet an outward-rounded certificate.",
            "The certificate would concern the fixed proof-model setup, not an implementation.",
        ],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--memory-bits", type=int, default=30)
    parser.add_argument(
        "--basis-mode", choices=("randomized", "prefix"), default="randomized"
    )
    parser.add_argument("--grid-min", type=float, default=-4.0)
    parser.add_argument("--grid-max", type=float, default=1.0)
    parser.add_argument("--grid-step", type=float, default=0.25)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    payload = evaluate(args)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
