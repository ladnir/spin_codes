#!/usr/bin/env python3
"""Floating-point all-one-word completion for the packet s=256 certificate.

Each outer block has one exceptional all-one word in addition to its regular
spectrum.  This checker partitions mixed configurations by the number a of
regular active blocks.

* For a<=105, delete every regular input and all but one forced all-one input.
* For 106<=a<=136, retain one forced packet and floor((a-3)/4) packed regular
  rank-four packets.  The subtraction of three allows every regular block in
  the retained forced packet to be discarded.
* For a>=137, delete every all-one input and reuse the regular-only inner
  bound.  The regular margin in this range absorbs all choices of all-one
  blocks.

Deleting input packets lowers the region's nonzero-packet count.  The
entrywise suffix envelope from the regular-word certificate converts that
count comparison into a matrix upper bound without assuming that the exact
region matrices are monotone.  Every mixed outer configuration is counted by
C(M,a) m^a (2^(M-a)-1).
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import gammaln, logsumexp

from analyze_riffle_fieldcheckpoint_regular_bulk_logdp import (
    log_matrix_power_moments_batch,
)
from analyze_riffle_packet4_fieldcheckpoint_s256_holder import log_matpow
from analyze_riffle_striped_random_outer import LOG2, log_two_power_minus_one
from certify_riffle_packet4_fieldcheckpoint_s256_support import (
    deterministic_support_region_envelope,
)


ROOT = Path(
    "constructions/riffle_bchperm_transpose_packetshuffle_fieldcheckpoint_g4_s256"
)
DEFAULT_REGULAR = ROOT / "receipts/all_profiles_regular_support_relaxation_delta09.json"
DEFAULT_OUTPUT = ROOT / "receipts/all_one_completion_delta09.json"


def mixed_middle_inner_logs(
    *,
    state_bits: int,
    packet_bits: int,
    packet_slots_per_epoch: int,
    epochs_per_region: int,
    regions: int,
    distance: int,
    log_surprisals: tuple[float, ...],
    maximum_full_groups: int,
    live_normalization_correction: float,
    groups: int,
) -> tuple[np.ndarray, np.ndarray]:
    best = np.full(maximum_full_groups + 1, math.inf)
    best_tilt = np.full(maximum_full_groups + 1, math.nan)
    full_nonzero = 1.0 - math.ldexp(1.0, -packet_bits)
    for log_surprisal in log_surprisals:
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        envelope = deterministic_support_region_envelope(
            state_bits=state_bits,
            packet_slots_per_epoch=packet_slots_per_epoch,
            epochs_per_region=epochs_per_region,
            groups=groups,
            z=z,
        )
        region = np.full((maximum_full_groups + 1, 2, 2), -math.inf)
        for full_groups in range(packet_slots_per_epoch):
            if full_groups > maximum_full_groups:
                break
            counts = np.arange(full_groups + 1, dtype=np.float64)
            log_probabilities = (
                gammaln(full_groups + 1.0)
                - gammaln(counts + 1.0)
                - gammaln(full_groups - counts + 1.0)
                + counts * math.log(full_nonzero)
                + (full_groups - counts) * math.log1p(-full_nonzero)
            )
            # One deterministic forced packet shifts the nonzero-packet count.
            region[full_groups] = logsumexp(
                log_probabilities[:, None, None]
                + envelope[1 : full_groups + 2],
                axis=0,
            )
        candidates = (
            log_matrix_power_moments_batch(region, regions)
            + live_normalization_correction
            + distance * surprisal
        )
        improved = candidates < best
        best[improved] = candidates[improved]
        best_tilt[improved] = log_surprisal
    return best, best_tilt


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    regular = json.loads(args.regular_receipt.read_text(encoding="utf-8"))
    parameters = regular["parameters"]
    outer_blocks = int(parameters["outer_blocks"])
    output_bits = int(parameters["output_bits"])
    distance = int(parameters["distance"])
    state_bits = int(parameters["state_bits"])
    packet_bits = int(parameters["packet_bits"])
    epochs_per_region = int(parameters["epochs_per_region"])
    regions = int(parameters["outer_bits"])
    packet_slots_per_epoch = state_bits // packet_bits
    groups = outer_blocks // packet_bits
    total_epochs = regions * epochs_per_region
    log_live_normalization = -math.log1p(-math.ldexp(1.0, -state_bits))
    live_normalization_correction = total_epochs * log_live_normalization
    if args.low_end >= args.high_start:
        raise ValueError("the mixed middle interval must be nonempty")

    best_forced = math.inf
    best_forced_tilt = math.nan
    for log_surprisal in args.forced_log_surprisals:
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        region = deterministic_support_region_envelope(
            state_bits=state_bits,
            packet_slots_per_epoch=packet_slots_per_epoch,
            epochs_per_region=epochs_per_region,
            groups=groups,
            z=z,
        )[1]
        all_regions = log_matpow(region, regions)
        moment = float(np.logaddexp(all_regions[0, 0], all_regions[0, 1]))
        candidate = moment + live_normalization_correction + distance * surprisal
        if candidate < best_forced:
            best_forced = candidate
            best_forced_tilt = log_surprisal

    maximum_middle_groups = (args.high_start - 1 - 3) // packet_bits
    middle_best, middle_tilt = mixed_middle_inner_logs(
        state_bits=state_bits,
        packet_bits=packet_bits,
        packet_slots_per_epoch=packet_slots_per_epoch,
        epochs_per_region=epochs_per_region,
        regions=regions,
        distance=distance,
        log_surprisals=tuple(args.middle_log_surprisals),
        maximum_full_groups=maximum_middle_groups,
        live_normalization_correction=live_normalization_correction,
        groups=groups,
    )

    rows = []
    contributions = []
    for regular_count in range(outer_blocks):
        all_one_choices = log_two_power_minus_one(outer_blocks - regular_count)
        if regular_count:
            regular_row = regular["occupation_rows"][regular_count - 1]
            regular_outer = float(regular_row["outer_log2_envelope"]) * LOG2
        else:
            regular_outer = 0.0
        mixed_outer = regular_outer + all_one_choices

        if regular_count <= args.low_end:
            inner = best_forced
            method = "one forced packet; all regular inputs deleted"
            tilt = best_forced_tilt
            retained_groups = 0
        elif regular_count < args.high_start:
            retained_groups = max(0, regular_count - 3) // packet_bits
            inner = float(middle_best[retained_groups])
            method = "one forced packet plus pessimistically packed regular packets"
            tilt = float(middle_tilt[retained_groups])
        else:
            regular_row = regular["occupation_rows"][regular_count - 1]
            inner = float(regular_row["inner_log2_upper"]) * LOG2
            method = "all all-one inputs deleted; regular bound reused"
            tilt = float(regular_row["best_log_surprisal"])
            retained_groups = None

        contribution = mixed_outer + inner
        contributions.append(contribution)
        rows.append(
            {
                "regular_active_blocks": regular_count,
                "retained_regular_full_groups": retained_groups,
                "method": method,
                "best_log_surprisal": tilt,
                "mixed_outer_log2_upper": mixed_outer / LOG2,
                "inner_log2_upper": inner / LOG2,
                "pointwise_log2_upper": contribution / LOG2,
            }
        )

    all_one_log = float(logsumexp(np.asarray(contributions)))
    regular_log = float(regular["regular_log2_upper"]) * LOG2
    total_log = float(np.logaddexp(regular_log, all_one_log))
    return {
        "schema": "riffle-packet4-fieldcheckpoint-s256-all-one-completion-v1",
        "candidate": regular["candidate"],
        "parameters": {
            "message_bits": int(parameters["message_bits"]),
            "output_bits": output_bits,
            "distance": distance,
            "relative_distance": float(parameters["relative_distance"]),
            "low_end": args.low_end,
            "high_start": args.high_start,
            "forced_log_surprisals": args.forced_log_surprisals,
            "middle_log_surprisals": args.middle_log_surprisals,
            "log2_live_normalization_correction": (
                live_normalization_correction / LOG2
            ),
        },
        "method": {
            "outer_mixed_count": "C(M,a) m^a (2^(M-a)-1)",
            "input_comparison": (
                "delete packets, then apply the entrywise packet-count "
                "suffix envelope"
            ),
            "arithmetic": "nearest binary64 log semiring",
            "outward_rounded": False,
        },
        "one_forced": {
            "best_log_surprisal": best_forced_tilt,
            "inner_log2_upper": best_forced / LOG2,
        },
        "all_one_class_log2_upper": all_one_log / LOG2,
        "all_one_class_lambda_bits_lower_float": -all_one_log / LOG2,
        "regular_class_log2_upper": regular_log / LOG2,
        "combined_log2_upper": total_log / LOG2,
        "combined_lambda_bits_lower_float": -total_log / LOG2,
        "dominant_all_one_rows": sorted(
            rows,
            key=lambda row: float(row["pointwise_log2_upper"]),
            reverse=True,
        )[:30],
        "mixed_rows": rows,
        "scope": (
            "Every message under the modeled regular spectrum, including every "
            "number and placement of exceptional all-one outer words.  Outward "
            "rounding and an outer-spectrum theorem remain open."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--regular-receipt", type=Path, default=DEFAULT_REGULAR)
    parser.add_argument("--low-end", type=int, default=105)
    parser.add_argument("--high-start", type=int, default=137)
    parser.add_argument(
        "--forced-log-surprisals",
        type=float,
        nargs="+",
        default=(-3.8, -3.775, -3.75),
    )
    parser.add_argument(
        "--middle-log-surprisals",
        type=float,
        nargs="+",
        default=(-5.0, -4.5, -4.0, -3.5, -3.0, -2.5),
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        "all_one_lambda_bits,"
        f"{payload['all_one_class_lambda_bits_lower_float']:.12f}"
    )
    print(
        "combined_lambda_bits,"
        f"{payload['combined_lambda_bits_lower_float']:.12f}"
    )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
