#!/usr/bin/env python3
"""Expected spectra of four EBCH128 blocks followed by BA-t.

The four-block direct sum is deterministic.  Before each accumulator, setup
samples one independent uniform permutation of the 512 coordinates.  The
reported post-accumulator spectrum is the ensemble expectation; it is not the
spectrum of a selected realization and is not a concentration certificate.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import re
import sys

from flint import fmpz_poly
import numpy as np
from scipy.special import logsumexp


WORKSTREAM = Path(__file__).resolve().parent
REPO_ROOT = WORKSTREAM.parents[1]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from analyze_systematic_expander_accumulate_spectrum import (  # noqa: E402
    LN2,
    apply_accumulator,
    log2_comb_scalar,
)


SOURCE = REPO_ROOT / "scripts" / "EBCH128_64.wd"
OUTPUT = WORKSTREAM / "ebch128x4_ba0_6_B512_expected_spectra.json"
LOCAL_LENGTH = 128
LOCAL_DIMENSION = 64
BLOCKS = 4
B = LOCAL_LENGTH * BLOCKS
K = LOCAL_DIMENSION * BLOCKS


def read_wd(path: Path) -> list[int]:
    spectrum = [0] * (LOCAL_LENGTH + 1)
    pattern = re.compile(r"^\s*(\d+)\s+(\d+)\s*$")
    for line in path.read_text(encoding="utf-8").splitlines():
        match = pattern.match(line)
        if match is None:
            continue
        weight, count = map(int, match.groups())
        if not 0 <= weight <= LOCAL_LENGTH:
            raise ValueError("spectrum weight outside EBCH128")
        spectrum[weight] = count
    if spectrum[0] != 1 or sum(spectrum) != 1 << LOCAL_DIMENSION:
        raise AssertionError("imported EBCH128 spectrum has the wrong mass")
    return spectrum


def macwilliams_dual_spectrum(spectrum: list[int]) -> list[int]:
    """Return the exact binary MacWilliams transform."""

    code_size = sum(spectrum)
    result = []
    for dual_weight in range(LOCAL_LENGTH + 1):
        numerator = 0
        for weight, multiplicity in enumerate(spectrum):
            krawtchouk = sum(
                (-1) ** overlap
                * math.comb(weight, overlap)
                * math.comb(LOCAL_LENGTH - weight, dual_weight - overlap)
                for overlap in range(max(0, dual_weight - (LOCAL_LENGTH - weight)), min(weight, dual_weight) + 1)
            )
            numerator += multiplicity * krawtchouk
        if numerator % code_size:
            raise AssertionError("MacWilliams coefficient is not integral")
        result.append(numerator // code_size)
    return result


def log2_integer(value: int) -> float:
    if value <= 0:
        return -math.inf
    shift = max(0, value.bit_length() - 53)
    return math.log2(value >> shift) + shift


def direct_sum_spectrum(local: list[int]) -> np.ndarray:
    polynomial = fmpz_poly(local) ** BLOCKS
    result = np.full(B + 1, -math.inf)
    for weight in range(B + 1):
        coefficient = int(polynomial[weight])
        if coefficient:
            result[weight] = log2_integer(coefficient)
    return result


def random_spectrum() -> np.ndarray:
    result = np.full(B + 1, -math.inf)
    result[0] = 0.0
    log_ratio = math.log2((1 << K) - 1) - math.log2((1 << B) - 1)
    for weight in range(1, B + 1):
        result[weight] = log2_comb_scalar(B, weight) + log_ratio
    return result


def apply_dual_accumulator(log_spectrum: np.ndarray) -> np.ndarray:
    """Apply the expected spectrum kernel for the dual BA stage.

    A primal stage first permutes and then accumulates.  Its action on the
    dual first permutes and then applies the inverse transpose of the
    accumulator.  Bijection counting gives

        Pr[wt(A^{-T} P x)=w | wt(x)=h]
          = N_B(w,h) / C(B,h),

    where N_B(w,h) is the usual accumulator input-output enumerator with
    input weight w and output weight h.
    """

    block_size = len(log_spectrum) - 1
    result = np.full(block_size + 1, -math.inf)
    result[0] = log_spectrum[0]
    for output_weight in range(1, block_size + 1):
        half_down = output_weight // 2
        half_up = (output_weight + 1) // 2
        low_input = half_up
        high_input = block_size - half_down
        input_weights = np.arange(low_input, high_input + 1)
        log_pair_counts = (
            np.asarray(
                [
                    log2_comb_scalar(block_size - int(weight), half_down)
                    for weight in input_weights
                ]
            )
            + np.asarray(
                [
                    log2_comb_scalar(int(weight) - 1, half_up - 1)
                    for weight in input_weights
                ]
            )
        )
        contributions = (
            log_spectrum[input_weights]
            + log_pair_counts
            - np.asarray(
                [log2_comb_scalar(block_size, int(weight)) for weight in input_weights]
            )
        )
        result[output_weight] = float(np.logaddexp2.reduce(contributions))
    return result


def log2_sum(values: np.ndarray) -> float:
    finite = np.isfinite(values)
    if not finite.any():
        return -math.inf
    return float(logsumexp(values[finite] * LN2) / LN2)


def summarize(stage: int, spectrum: np.ndarray, random: np.ndarray) -> dict[str, object]:
    mass = log2_sum(spectrum)
    if abs(mass - K) > 1e-7:
        raise AssertionError(f"stage {stage} spectrum mass is {mass}, expected {K}")
    nonzero = np.arange(1, B + 1)
    excess = spectrum[1:] - random[1:]
    finite = np.isfinite(spectrum[1:])
    maximum_index = int(np.argmax(np.where(finite, excess, -math.inf))) + 1
    bands = ((1, 79), (42, 79), (80, 432), (433, 470), (471, 512))
    band_rows = []
    for lower, upper in bands:
        indices = np.arange(lower, upper + 1)
        stage_mass = log2_sum(spectrum[indices])
        random_mass = log2_sum(random[indices])
        pointwise = spectrum[indices] - random[indices]
        finite_band = np.isfinite(spectrum[indices])
        maximizing = int(indices[np.argmax(np.where(finite_band, pointwise, -math.inf))])
        band_rows.append(
            {
                "weights": [lower, upper],
                "expected_mass_log2": stage_mass,
                "random_mass_log2": random_mass,
                "mass_excess_bits": stage_mass - random_mass,
                "maximum_pointwise_excess_bits": float(pointwise[maximizing - lower]),
                "maximizing_weight": maximizing,
            }
        )
    low_cumulative = []
    for cutoff in (10, 20, 30, 40, 41, 42, 50, 60, 79):
        low_cumulative.append(
            {
                "cutoff": cutoff,
                "expected_nonzero_mass_log2": log2_sum(spectrum[1 : cutoff + 1]),
                "random_nonzero_mass_log2": log2_sum(random[1 : cutoff + 1]),
            }
        )
    rows = [
        {
            "weight": weight,
            "log2_expected_multiplicity": (
                None if not math.isfinite(float(spectrum[weight])) else float(spectrum[weight])
            ),
            "random_linear_log2_expected_multiplicity": float(random[weight]),
            "excess_over_random_bits": (
                None if not math.isfinite(float(spectrum[weight])) else float(spectrum[weight] - random[weight])
            ),
        }
        for weight in range(B + 1)
    ]
    return {
        "accumulators": stage,
        "total_mass_log2": mass,
        "maximum_pointwise_excess_bits": float(excess[maximum_index - 1]),
        "maximum_pointwise_excess_weight": maximum_index,
        "band_summaries": band_rows,
        "low_cumulative": low_cumulative,
        "spectrum": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maximum-accumulators", type=int, default=6)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    if args.maximum_accumulators < 0:
        parser.error("maximum accumulator count must be nonnegative")

    local = read_wd(SOURCE)
    local_dual = macwilliams_dual_spectrum(local)
    if local_dual != local:
        raise AssertionError("imported EBCH128 enumerator is not formally self-dual")
    dual_distance = next(weight for weight, count in enumerate(local_dual) if weight and count)
    if dual_distance != 22:
        raise AssertionError("imported EBCH128 dual distance is not 22")
    expected = direct_sum_spectrum(local)
    # The imported extended BCH [128,64,22] weight enumerator is formally
    # self-dual.  Hence the four-block direct sum starts with the same primal
    # and dual spectra.  The two spectra diverge after the first BA stage.
    expected_dual = expected.copy()
    random = random_spectrum()
    stages = []
    dual_stages = []
    for stage in range(args.maximum_accumulators + 1):
        stages.append(summarize(stage, expected, random))
        dual_stages.append(summarize(stage, expected_dual, random))
        summary = stages[-1]
        dual_summary = dual_stages[-1]
        print(
            f"stage,{stage},max_excess,{summary['maximum_pointwise_excess_bits']:.9f},"
            f"weight,{summary['maximum_pointwise_excess_weight']},"
            f"dual_max_excess,{dual_summary['maximum_pointwise_excess_bits']:.9f},"
            f"dual_weight,{dual_summary['maximum_pointwise_excess_weight']}",
            flush=True,
        )
        if stage < args.maximum_accumulators:
            expected = apply_accumulator(expected)
            expected_dual = apply_dual_accumulator(expected_dual)

    payload = {
        "schema": "ebch128x4-ba-expected-spectra-v2",
        "status": "BINARY64_ENSEMBLE_EXPECTATION_DIAGNOSTIC",
        "construction": {
            "base": "direct sum of four extended BCH [128,64,22] codes",
            "block_bits": B,
            "dimension": K,
            "accumulator_stages": [0, args.maximum_accumulators],
            "sampler": "one independent uniform 512-coordinate permutation before each accumulator",
            "reuse": "one sampled BA-t constituent is repeated in every SPIN outer row",
        },
        "source_spectrum": str(SOURCE.relative_to(REPO_ROOT)).replace("\\", "/"),
        "method": {
            "direct_sum": "exact integer polynomial fourth power",
            "local_dual": "exact integer MacWilliams transform",
            "accumulators": "exact accumulator IOWE composed in binary64 log space",
            "comparison": "expected spectrum of a uniform random binary [512,256] linear code",
        },
        "validation": {
            "local_enumerator_formally_self_dual": True,
            "local_dual_distance": dual_distance,
            "base_dual_distance": dual_distance,
            "base_coordinate_independence_order": dual_distance - 1,
        },
        "stages": stages,
        "dual_stages": dual_stages,
        "limitations": [
            "An expected spectrum does not imply a high-probability spectrum event for one sampled and repeated constituent.",
            "The dual spectra use the inverse-transpose accumulator kernel derived by bijection counting.",
            "Nearest binary64 arithmetic makes the numerical comparison diagnostic.",
            "The calculation does not yet include the structured SPIN routing or RandomStepConv transfer.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
