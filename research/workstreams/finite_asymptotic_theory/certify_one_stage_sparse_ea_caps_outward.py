#!/usr/bin/env python3
"""Build simultaneous shell caps for the one-stage sparse-EA constituent.

The constituent maps 256 message bits through 512 independent weight-33
parity rows and then through the zero-initialized accumulator.  The script
combines the certified defect-shell variance receipts with the cancellation-
free central receipt.  It uses exact rational Cantelli checks to construct the
integer caps and outward first moments for the zero-cap tails and kernel.
"""

from __future__ import annotations

import hashlib
import json
import math
import sys
from fractions import Fraction
from pathlib import Path

from flint import ctx
import numpy as np


WORKSTREAM = Path(__file__).resolve().parent
PURE = WORKSTREAM / "pure_expander_accumulate"
sys.path.insert(0, str(PURE))

from certify_dominant_character_deviation_outward import reference_intervals  # noqa: E402
from certify_primal_schur_diagonal_outward import (  # noqa: E402
    add_upper_scalar,
    coefficient_upper,
    positive_multiply_upper,
)


LOW = PURE / "dominant_character_likelihood_K256_B512_r33_w42_79_outward.json"
CENTRAL = PURE / "dominant_character_deviation_K256_B512_r33_w80_432_outward.json"
HIGH = PURE / "dominant_character_likelihood_K256_B512_r33_w433_470_outward.json"
OUTPUT = WORKSTREAM / "one_stage_sparse_ea_K256_B512_r33_spectrum_caps_outward.json"

MESSAGE_BITS = 256
OUTPUT_BITS = 512
RIGHT_DEGREE = 33
LOWER_SUPPORT = 42
UPPER_SUPPORT = 470
PER_SHELL_DELTA = Fraction(1, 1 << 51)
RANK_ATTEMPTS = 16


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def upper_fraction(value: float) -> Fraction:
    if not math.isfinite(value) or value < 0:
        raise ArithmeticError(f"invalid nonnegative endpoint: {value}")
    return Fraction(*value.as_integer_ratio())


def random_mean(shell: int) -> Fraction:
    return Fraction(math.comb(OUTPUT_BITS, shell), 1 << (OUTPUT_BITS - MESSAGE_BITS))


def receipt_rows(path: Path, lower: int, upper: int) -> dict[int, dict[str, object]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("status") != "OUTWARD_ARB_AND_BINARY64_CERTIFICATE":
        raise AssertionError(f"nonaccepting variance receipt: {path.name}")
    parameters = payload.get("parameters", {})
    expected = {
        "message_bits": MESSAGE_BITS,
        "output_bits": OUTPUT_BITS,
        "right_degree": RIGHT_DEGREE,
    }
    for name, value in expected.items():
        if parameters.get(name) != value:
            raise AssertionError(f"parameter mismatch in {path.name}: {name}")
    rows = {int(row["shell_weight"]): row for row in payload["shells"]}
    if set(rows) != set(range(lower, upper + 1)):
        raise AssertionError(f"incomplete shell range in {path.name}")
    return rows


def mean_ratio_upper(shell: int) -> Fraction:
    _one_low, one_high, _one_deviation, _pair_deviation = reference_intervals(
        MESSAGE_BITS, OUTPUT_BITS, RIGHT_DEGREE, shell
    )
    result = np.float64(0)
    for weight in range(1, MESSAGE_BITS + 1):
        message_fraction = np.ldexp(
            coefficient_upper(math.comb(MESSAGE_BITS, weight)), -MESSAGE_BITS
        )
        term = positive_multiply_upper(
            np.asarray([message_fraction]), np.asarray([one_high[weight]])
        )[0]
        result = add_upper_scalar(result, term)
    return upper_fraction(float(result))


def cantelli_failure(
    cap: int, mean_upper: Fraction, variance_upper: Fraction
) -> Fraction:
    deviation = Fraction(cap + 1) - mean_upper
    if deviation <= 0:
        return Fraction(1)
    return variance_upper / (variance_upper + deviation * deviation)


def cantelli_cap(mean_upper: Fraction, variance_upper: Fraction) -> tuple[int, Fraction]:
    if mean_upper <= 0 or variance_upper <= 0:
        raise AssertionError("positive-cap shell needs positive moments")
    # Use an integer upper bound on the square root.  This avoids converting
    # the central-shell integers, whose binary64 spacing is enormous, into a
    # purportedly minimal cap.
    radicand = variance_upper * (1 - PER_SHELL_DELTA) / PER_SHELL_DELTA
    radicand_ceiling = (
        radicand.numerator + radicand.denominator - 1
    ) // radicand.denominator
    deviation = math.isqrt(radicand_ceiling)
    if deviation * deviation < radicand_ceiling:
        deviation += 1
    shifted_mean = mean_upper - 1
    shifted_mean_ceiling = (
        shifted_mean.numerator + shifted_mean.denominator - 1
    ) // shifted_mean.denominator
    cap = max(0, shifted_mean_ceiling + deviation)

    failure = cantelli_failure(cap, mean_upper, variance_upper)
    if failure > PER_SHELL_DELTA:
        raise AssertionError("integer square-root cap failed its exact check")
    return cap, failure


def band_majorant(caps: list[int], lower: int, upper: int, probability: Fraction) -> Fraction:
    result = Fraction(0)
    for shell in range(lower, upper + 1):
        if not caps[shell]:
            continue
        shell_probability = (
            Fraction(math.comb(OUTPUT_BITS, shell))
            * probability**shell
            * (1 - probability) ** (OUTPUT_BITS - shell)
        )
        result = max(result, Fraction(caps[shell], 1) / shell_probability)
    if result <= 0:
        raise AssertionError("empty band")
    return result


def log2_upper(value: Fraction) -> float:
    estimate = math.log2(value.numerator) - math.log2(value.denominator)
    return math.nextafter(estimate, math.inf)


def main() -> None:
    ctx.prec = 1536
    ctx.threads = 1
    low = receipt_rows(LOW, 42, 79)
    central = receipt_rows(CENTRAL, 80, 432)
    high = receipt_rows(HIGH, 433, 470)
    rows = {**low, **central, **high}
    if set(rows) != set(range(LOWER_SUPPORT, UPPER_SUPPORT + 1)):
        raise AssertionError("variance receipts do not cover the cap support")

    sys.path.insert(0, str(WORKSTREAM))
    from certify_single_random_constituent_dense_outward import (  # noqa: E402
        exact_band_majorant,
        exact_caps,
    )

    p_defect = Fraction(79, OUTPUT_BITS)
    p_central = Fraction(1, 2)
    frozen_caps, _ = exact_caps()
    frozen_low = exact_band_majorant(frozen_caps, 42, 79, p_defect)
    frozen_high = exact_band_majorant(frozen_caps, 433, 470, 1 - p_defect)
    frozen_central = exact_band_majorant(frozen_caps, 80, 432, p_central)
    frozen_majorants = {
        "low": frozen_low,
        "central": frozen_central,
        "high": frozen_high,
    }

    moments: dict[int, tuple[Fraction, Fraction]] = {}
    independent_caps: dict[int, tuple[int, Fraction]] = {}
    for shell in range(LOWER_SUPPORT, UPPER_SUPPORT + 1):
        row = rows[shell]
        mean_upper = random_mean(shell) * upper_fraction(
            float(row["mean_ratio_upper"])
        )
        variance_upper = mean_upper * upper_fraction(
            float(row["variance_to_mean_upper"])
        )
        moments[shell] = (mean_upper, variance_upper)
        independent_caps[shell] = cantelli_cap(mean_upper, variance_upper)

    # Intersect two sufficient constraints.  The independent 2^-51 cap keeps
    # Q=1 and Q=2 sharp.  The frozen band envelope preserves every Q>=3
    # majorant hypothesis in the existing transfer.
    caps = [0 for _ in range(OUTPUT_BITS + 1)]
    cap_rows = []
    cap_failure = Fraction(0)
    for shell in range(LOWER_SUPPORT, UPPER_SUPPORT + 1):
        if shell <= 79:
            band = "low"
            probability = p_defect
        elif shell <= 432:
            band = "central"
            probability = p_central
        else:
            band = "high"
            probability = 1 - p_defect
        shell_probability = (
            Fraction(math.comb(OUTPUT_BITS, shell))
            * probability**shell
            * (1 - probability) ** (OUTPUT_BITS - shell)
        )
        envelope_value = frozen_majorants[band] * shell_probability
        envelope_cap = envelope_value.numerator // envelope_value.denominator
        mean_upper, variance_upper = moments[shell]
        independent_cap, independent_failure = independent_caps[shell]
        cap = min(independent_cap, envelope_cap)
        failure = cantelli_failure(cap, mean_upper, variance_upper)
        if failure >= 1:
            raise AssertionError(f"frozen envelope lies below the mean at shell {shell}")
        caps[shell] = cap
        cap_failure += failure
        cap_rows.append(
            {
                "shell_weight": shell,
                "mean_upper_log2": log2_upper(mean_upper),
                "variance_upper_log2": log2_upper(variance_upper),
                "integer_cap": cap,
                "cantelli_failure_upper_log2": log2_upper(failure),
                "independent_2^-51_cap": independent_cap,
                "frozen_band_envelope_cap": envelope_cap,
                "independent_2^-51_failure_upper_log2": log2_upper(
                    independent_failure
                ),
                "band": band,
                "variance_source": (
                    LOW.name
                    if shell <= 79
                    else CENTRAL.name
                    if shell <= 432
                    else HIGH.name
                ),
            }
        )

    zero_shells = [
        *range(1, LOWER_SUPPORT),
        *range(UPPER_SUPPORT + 1, OUTPUT_BITS + 1),
    ]
    mean_ratio_cache: dict[int, Fraction] = {}
    tail_failure = Fraction(0)
    tail_rows = []
    for index, shell in enumerate([0, *zero_shells]):
        ratio = mean_ratio_upper(shell)
        mean = random_mean(shell) * ratio
        mean_ratio_cache[shell] = ratio
        if shell:
            tail_failure += mean
        tail_rows.append(
            {
                "shell_weight": shell,
                "mean_upper_log2": log2_upper(mean),
                "role": "kernel" if shell == 0 else "zero-cap tail",
            }
        )
        if (index + 1) % 16 == 0:
            print(f"first_moment_shells,{index + 1},{len(zero_shells)+1}", flush=True)

    kernel_failure = random_mean(0) * mean_ratio_cache[0]
    if kernel_failure >= 1:
        raise AssertionError("rank lower bound is nonpositive")
    accepted_cap_failure = (cap_failure + tail_failure) / (1 - kernel_failure)
    abort_failure = kernel_failure**RANK_ATTEMPTS
    setup_failure = accepted_cap_failure + abort_failure

    low_majorant = band_majorant(caps, 42, 79, p_defect)
    high_majorant = band_majorant(caps, 433, 470, 1 - p_defect)
    central_majorant = band_majorant(caps, 80, 432, p_central)

    domination = {
        "low": low_majorant <= frozen_low,
        "high": high_majorant <= frozen_high,
        "central": central_majorant <= frozen_central,
    }

    payload = {
        "schema": "one-stage-sparse-ea-spectrum-caps-outward-v1",
        "status": "OUTWARD_CAP_CERTIFICATE",
        "parameters": {
            "message_bits": MESSAGE_BITS,
            "output_bits": OUTPUT_BITS,
            "right_degree": RIGHT_DEGREE,
            "per_shell_delta": f"{PER_SHELL_DELTA.numerator}/{PER_SHELL_DELTA.denominator}",
            "rank_test_attempts": RANK_ATTEMPTS,
        },
        "probability_space": {
            "one_attempt": "sample 512 independent uniform weight-33 row vectors in F_2^256; apply the zero-initialized length-512 accumulator",
            "bounded_setup": "sample independent attempts and accept the first full-rank constituent; abort after 16 rank failures",
            "reuse": "the accepted constituent is fixed and reused in every SPIN outer position",
        },
        "claim": {
            "positive_cap_support": [LOWER_SUPPORT, UPPER_SUPPORT],
            "all_positive_shells_capped": True,
            "all_other_nonzero_shell_caps_zero": True,
            "kernel_failure_upper_log2": log2_upper(kernel_failure),
            "zero_cap_tail_failure_upper_log2": log2_upper(tail_failure),
            "positive_cap_failure_upper_log2": log2_upper(cap_failure),
            "rank_conditioned_cap_failure_upper_log2": log2_upper(accepted_cap_failure),
            "bounded_setup_abort_upper_log2": log2_upper(abort_failure),
            "total_setup_failure_upper_log2": log2_upper(setup_failure),
            "total_setup_margin_bits_lower": -log2_upper(setup_failure),
        },
        "band_majorants_log2_upper": {
            "low": log2_upper(low_majorant),
            "central": log2_upper(central_majorant),
            "high": log2_upper(high_majorant),
        },
        "frozen_random_band_majorants_log2_upper": {
            "low": log2_upper(frozen_low),
            "central": log2_upper(frozen_central),
            "high": log2_upper(frozen_high),
        },
        "band_majorant_dominated_by_frozen_random": domination,
        "caps": caps,
        "cap_rows": cap_rows,
        "zero_cap_and_kernel_rows": tail_rows,
        "sources": [
            {"file": path.name, "sha256": sha256(path)}
            for path in (LOW, CENTRAL, HIGH)
        ],
        "proof_scope": [
            "Cantelli caps use exact rational checks after outward moment conversion.",
            "Markov's inequality controls the zero-cap tails and rank failure.",
            "Rank conditioning divides the cap-event failure by the certified full-rank probability; sixteen independent failures bound setup abort.",
            "Band domination alone transfers only those SPIN bounds whose routing and inner-ensemble hypotheses match the frozen verifier.",
        ],
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"claim": payload["claim"], "domination": domination}, indent=2))
    print(f"output={OUTPUT}")


if __name__ == "__main__":
    main()
