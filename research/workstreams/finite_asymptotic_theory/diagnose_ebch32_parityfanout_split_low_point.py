"""Diagnostic four-band refinement of the rejected three-band lattice point."""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import minimize


WORKSTREAM = Path(__file__).resolve().parent
sys.path.insert(0, str(WORKSTREAM))

from analyze_golay_ba_rm2sub_joint import DEFAULT_SELECTION, load_envelope  # noqa: E402
from diagnose_finite_k20_band_compositions_qL import log_multinomial, softmax  # noqa: E402
from diagnose_finite_k20_renyi_dense import (  # noqa: E402
    high_precision_transfer_log_moment,
    logistic,
    logit,
)
from diagnose_ebch32_ba_k20_three_band_qL import (  # noqa: E402
    B,
    EPOCHS,
    L,
    N,
    band_log_moment,
    conditioned_spectrum,
    fast_inner_log,
)
from diagnose_ebch32_parityfanout_ba_three_band_q import DISTANCE  # noqa: E402


if os.environ.get("SPIN_EBCH_PARITY_FANOUT", "0") != "1":
    raise RuntimeError("set SPIN_EBCH_PARITY_FANOUT=1")
if os.environ.get("SPIN_EBCH_LOWER_WEIGHT", "22") != "24":
    raise RuntimeError("set SPIN_EBCH_LOWER_WEIGHT=24")

BANDS = ((24, 54), (55, 100), (101, 155), (156, 232))
EXTERNAL_COUNTS = (3344, 1, 21)


def optimize_split(envelope, spectrum: np.ndarray, first_low_count: int):
    active_counts = np.asarray(
        (
            first_low_count,
            EXTERNAL_COUNTS[0] - first_low_count,
            EXTERNAL_COUNTS[1],
            EXTERNAL_COUNTS[2],
        ),
        dtype=np.float64,
    )
    occupation = float(np.sum(active_counts))
    counts = np.concatenate(([L - occupation], active_counts))
    frequencies = (counts + 0.5) / (L + 2.5)
    log_type_count = log_multinomial(tuple(float(value) for value in counts))

    def decode(point: np.ndarray):
        probabilities = softmax(point[:5])
        values = np.asarray(
            [min(1.0 - 1e-12, max(1e-12, logistic(float(x)))) for x in point[5:9]]
        )
        surprisal = math.exp(min(4.0, max(-12.0, float(point[9]))))
        order = 1.0 + math.exp(min(10.0, max(-8.0, float(point[10]))))
        return probabilities, values, surprisal, order

    def evaluate(point: np.ndarray, high_precision: bool = False):
        probabilities, values, surprisal, order = decode(point)
        moments = np.asarray(
            [
                band_log_moment(spectrum, band, value, order)
                for band, value in zip(BANDS, values)
            ]
        )
        bit_probability = float(np.dot(probabilities[1:], values))
        inner = (
            high_precision_transfer_log_moment(
                envelope, bit_probability, surprisal, epochs=EPOCHS
            )
            if high_precision
            else fast_inner_log(envelope, bit_probability, surprisal)
        )
        if not math.isfinite(inner):
            return math.inf, None
        log_conditioning = log_type_count + float(np.dot(counts, np.log(probabilities)))
        raw = inner + DISTANCE * surprisal - B * log_conditioning
        conjugate = (order - 1.0) / order
        value = (
            log_type_count
            + float(np.dot(active_counts, moments)) / order
            + conjugate * min(0.0, raw)
        )
        return value, {
            "split_active_counts": active_counts.astype(int).tolist(),
            "reference_type_probabilities": probabilities.tolist(),
            "reference_value_probabilities": values.tolist(),
            "reference_bit_probability": bit_probability,
            "surprisal": surprisal,
            "holder_order": order,
            "band_log2_renyi_moments": (moments / math.log(2.0)).tolist(),
            "raw_reference_bad_log2_upper": raw / math.log(2.0),
        }

    starts = []
    alpha = occupation / L
    for order in (2.0, 8.0, 32.0):
        for surprisal in (max(1e-4, 0.8 * alpha), max(0.02, 2.1 * alpha), 2.1):
            starts.append(
                np.asarray(
                    [
                        *np.log(frequencies),
                        *(logit(value) for value in (0.25, 0.42, 0.5, 0.68)),
                        math.log(surprisal),
                        math.log(order - 1.0),
                    ]
                )
            )
    best = None
    for start in starts:
        result = minimize(
            lambda point: evaluate(point)[0],
            start,
            method="Nelder-Mead",
            options={"maxiter": 1000, "xatol": 1e-8, "fatol": 1e-7},
        )
        if best is None or float(result.fun) < float(best.fun):
            best = result
    assert best is not None
    value, details = evaluate(best.x, high_precision=True)
    assert details is not None
    return {
        "first_low_count": first_low_count,
        "log2_upper": value / math.log(2.0),
        "margin_bits": -value / math.log(2.0),
        "optimizer_success": bool(best.success),
        **details,
    }


def evaluate_fixed_witness(envelope, first_low_count: int, witness):
    active_counts = np.asarray(
        (
            first_low_count,
            EXTERNAL_COUNTS[0] - first_low_count,
            EXTERNAL_COUNTS[1],
            EXTERNAL_COUNTS[2],
        ),
        dtype=np.float64,
    )
    occupation = float(np.sum(active_counts))
    counts = np.concatenate(([L - occupation], active_counts))
    log_type_count = log_multinomial(tuple(float(value) for value in counts))
    probabilities = np.asarray(witness["reference_type_probabilities"], dtype=np.float64)
    values = np.asarray(witness["reference_value_probabilities"], dtype=np.float64)
    moments = np.asarray(witness["band_log2_renyi_moments"], dtype=np.float64) * math.log(2.0)
    surprisal = float(witness["surprisal"])
    order = float(witness["holder_order"])
    bit_probability = float(np.dot(probabilities[1:], values))
    inner = high_precision_transfer_log_moment(
        envelope, bit_probability, surprisal, epochs=EPOCHS
    )
    log_conditioning = log_type_count + float(np.dot(counts, np.log(probabilities)))
    raw = inner + DISTANCE * surprisal - B * log_conditioning
    conjugate = (order - 1.0) / order
    value = (
        log_type_count
        + float(np.dot(active_counts, moments)) / order
        + conjugate * min(0.0, raw)
    )
    return {
        "first_low_count": first_low_count,
        "margin_bits": -value / math.log(2.0),
        "raw_reference_bad_log2_upper": raw / math.log(2.0),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts", default="0,1672,3344")
    args = parser.parse_args()
    counts = tuple(int(value) for value in args.counts.split(","))
    if any(value < 0 or value > EXTERNAL_COUNTS[0] for value in counts):
        raise ValueError("split counts must lie in [0,3344]")
    spectrum, _ = conditioned_spectrum()
    envelope = load_envelope(DEFAULT_SELECTION, DISTANCE / N)
    rows = []
    for value in counts:
        row = optimize_split(envelope, spectrum, value)
        row["fixed_witness_endpoint_checks"] = [
            evaluate_fixed_witness(envelope, endpoint, row)
            for endpoint in (0, EXTERNAL_COUNTS[0])
        ]
        rows.append(row)
        print(
            f"first_low_count={value},margin_bits={row['margin_bits']:.9f},"
            f"raw_log2={row['raw_reference_bad_log2_upper']:.9f}",
            flush=True,
        )
    print(
        json.dumps(
            {
                "status": "BINARY64_DIAGNOSTIC",
                "external_three_band_counts": list(EXTERNAL_COUNTS),
                "split_bands": [list(band) for band in BANDS],
                "rows": rows,
                "limitations": [
                    "The split-count set is sampled, not a gap-free cover.",
                    "Witness optimization uses nearest binary64 arithmetic.",
                    "No outward verifier is implemented for this refinement.",
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
