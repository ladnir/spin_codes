#!/usr/bin/env python3
"""Probe occupation two with one Bernoulli reference per exact RM shell.

The two active outer rows use independent coordinate permutations. For each
ordered pair of exact RM weights, the probe dominates each fixed-weight shell
by a Bernoulli product measure with matching mean. Arithmetic is binary64.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.special import logsumexp


HERE = Path(__file__).resolve().parent
WORKSTREAM = HERE.parent
REPOSITORY = WORKSTREAM.parent.parent
SCRIPTS = REPOSITORY / "scripts"
sys.path.insert(0, str(WORKSTREAM))
sys.path.insert(0, str(SCRIPTS))

from analyze_riffle_bitshuffle_splitstate_regular_bulk import (  # noqa: E402
    load_nonzero_spectrum,
    load_uniform_nonactivation,
    regular_region_log_matrices,
    splitstate_impulse_matrices,
)
from small_k_replay.evaluate_exact_spectra_q1_phase import (  # noqa: E402
    CONSTITUENTS,
    load_spectrum,
)
from small_k_replay.evaluate_rm2sub_primary_tranche import configuration  # noqa: E402
from small_k_replay.probe_rm2sub_three_group_bridge import (  # noqa: E402
    log_matrix_power_moments_binary,
)
from small_k_replay.probe_rm2sub_two_band_bridge import (  # noqa: E402
    verify_exact_nonactivation,
)


DEFAULT_OUTPUT = HERE / "rm2sub_rm49_q2_shell_reference_d100.json"
BLOCK_BITS = 512
DIMENSION = 256
MESSAGE_BITS = 1 << 16
OUTER_ROWS = MESSAGE_BITS // DIMENSION
OUTPUT_BITS = 2 * MESSAGE_BITS
BAD_WEIGHT = OUTPUT_BITS // 10
STEP_BITS = 64
STATE_BITS = 14


def shifted_probability(weight: int, logit_shift: float) -> float:
    probability = weight / BLOCK_BITS
    if probability == 1.0:
        return 1.0
    odds = probability / (1.0 - probability)
    shifted_odds = odds * math.exp(logit_shift)
    return shifted_odds / (1.0 + shifted_odds)


def shell_envelope_log(weight: int, count: int, probability: float) -> float:
    shell = math.lgamma(BLOCK_BITS + 1) - math.lgamma(weight + 1) - math.lgamma(BLOCK_BITS - weight + 1)
    if probability == 1.0:
        reference = 0.0
    else:
        reference = weight * math.log(probability) + (BLOCK_BITS - weight) * math.log1p(-probability)
    return math.log(count) - shell - reference


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tilt-minimum-tenths", type=int, default=-80)
    parser.add_argument("--tilt-maximum-tenths", type=int, default=-10)
    parser.add_argument(
        "--logit-shifts",
        type=float,
        nargs="+",
        default=[0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0],
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    constituent = CONSTITUENTS["rm49"]
    spectrum = load_spectrum(constituent)
    weights = np.asarray(
        sorted(weight for weight, count in spectrum.items() if weight and count),
        dtype=np.int64,
    )
    counts = np.asarray([spectrum[int(weight)] for weight in weights], dtype=object)
    config = configuration(STEP_BITS, STATE_BITS + 6)
    nonactivation_path = Path(config["b"])
    nonactivation = load_uniform_nonactivation(nonactivation_path, STEP_BITS)
    verify_exact_nonactivation(nonactivation_path, nonactivation)
    live_spectrum = load_nonzero_spectrum(
        Path(config["a"]), STEP_BITS, STATE_BITS
    )

    best = np.full(len(weights) * len(weights), math.inf)
    witnesses = np.zeros(best.shape, dtype=np.int64)
    shift_witnesses = np.full(best.shape, math.nan)
    location_log = math.log(math.comb(OUTER_ROWS, 2))
    for logit_shift in args.logit_shifts:
        probabilities = np.asarray(
            [shifted_probability(int(weight), logit_shift) for weight in weights]
        )
        envelope_logs = np.asarray(
            [
                shell_envelope_log(int(weight), int(count), float(probability))
                for weight, count, probability in zip(
                    weights, counts, probabilities, strict=True
                )
            ]
        )
        first_probabilities = np.repeat(probabilities, len(weights))
        second_probabilities = np.tile(probabilities, len(weights))
        first_envelopes = np.repeat(envelope_logs, len(weights))
        second_envelopes = np.tile(envelope_logs, len(weights))
        distributions = np.column_stack(
            (
                (1.0 - first_probabilities) * (1.0 - second_probabilities),
                first_probabilities * (1.0 - second_probabilities)
                + (1.0 - first_probabilities) * second_probabilities,
                first_probabilities * second_probabilities,
            )
        )
        if np.max(np.abs(np.sum(distributions, axis=1) - 1.0)) > 4e-16:
            raise ArithmeticError("two-shell reference distribution lost mass")
        for witness_tenths in range(
            args.tilt_minimum_tenths, args.tilt_maximum_tenths + 1
        ):
            log_surprisal = witness_tenths / 10.0
            surprisal = math.exp(log_surprisal)
            z = math.exp(-surprisal)
            impulses = splitstate_impulse_matrices(
                z=z,
                step_bits=STEP_BITS,
                state_bits=STATE_BITS,
                constituent_distance=int(config["a_distance"]),
                live_moment_order=3,
                live_model="support-averaged-preaddmul",
                nonactivation=nonactivation,
                live_spectrum=live_spectrum,
            )
            for input_weight in range(1, STEP_BITS + 1):
                impulses[input_weight, 0, 1] = (
                    max(0.0, 1.0 - float(nonactivation[input_weight]))
                    * z**input_weight
                )
            with np.errstate(divide="ignore"):
                impulse_logs = np.log(impulses)
            region_logs = regular_region_log_matrices(
                impulse_logs, STEP_BITS, OUTER_ROWS // STEP_BITS, 2
            )
            regions = np.exp(region_logs).reshape((3, 4))
            references = (distributions @ regions).reshape((-1, 2, 2))
            with np.errstate(divide="ignore"):
                reference_logs = np.log(references)
            moments = log_matrix_power_moments_binary(reference_logs, BLOCK_BITS)
            inner = np.minimum(0.0, moments + BAD_WEIGHT * surprisal)
            candidates = location_log + first_envelopes + second_envelopes + inner
            improved = candidates < best
            best[improved] = candidates[improved]
            witnesses[improved] = witness_tenths
            shift_witnesses[improved] = logit_shift

    aggregate = float(logsumexp(best))
    dominant_index = int(np.argmax(best))
    rows = []
    for index, value in enumerate(best):
        first_index, second_index = divmod(index, len(weights))
        rows.append(
            {
                "first_weight": int(weights[first_index]),
                "second_weight": int(weights[second_index]),
                "log_surprisal_exact": f"{int(witnesses[index])}/10",
                "logit_shift": float(shift_witnesses[index]),
                "log2_upper_diagnostic": float(value / math.log(2.0)),
                "margin_bits_diagnostic": float(-value / math.log(2.0)),
            }
        )
    dominant = rows[dominant_index]
    payload = {
        "schema": "rm2sub-rm49-q2-shell-reference-diagnostic-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "construction": "one fixed RM(4,9) constituent repeated in every row, structured routing, fixed RM2Sub t64,s14 maps, and independent nonzero epoch scalars",
        "probability_space": "row-coordinate and region permutations plus RM2Sub epoch scalars",
        "parameters": {
            "message_bits": MESSAGE_BITS,
            "output_bits": OUTPUT_BITS,
            "outer_rows": OUTER_ROWS,
            "outer_block_bits": BLOCK_BITS,
            "step_bits": STEP_BITS,
            "state_bits": STATE_BITS,
            "bad_weight": BAD_WEIGHT,
            "occupation": 2,
            "ordered_weight_pair_count": len(rows),
            "tilt_tenths": [args.tilt_minimum_tenths, args.tilt_maximum_tenths],
            "logit_shifts": args.logit_shifts,
        },
        "reduction": "For each exact outer weight w, dominate its uniformly permuted shell by Bernoulli(w/512) with the exact pointwise density ratio. Apply the positive forced-support RM2Sub transfer to the sum of the two Bernoulli bits in each region.",
        "log2_upper_diagnostic": aggregate / math.log(2.0),
        "margin_bits_diagnostic": -aggregate / math.log(2.0),
        "dominant_pair": dominant,
        "pair_rows": rows,
        "scope": "Nearest binary64 discovery only. An outward checker must authenticate the exact shell ratios and fixed witnesses.",
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"margin_bits_diagnostic,{payload['margin_bits_diagnostic']:.12f}")
    print(f"dominant_pair,{dominant['first_weight']},{dominant['second_weight']}")
    print(f"dominant_witness,{dominant['log_surprisal_exact']}")
    print(f"dominant_logit_shift,{dominant['logit_shift']}")
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
