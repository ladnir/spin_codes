#!/usr/bin/env python3
"""Sample banded all-active compositions for layered BCH250 fanout.

The input spectrum is an expected, mass-coupled envelope after layered
ParityFanout.  Therefore this script is diagnostic even when its numerical
optimization closes.  It tests whether a later one-sample concentration
theorem would have enough RM2Sub margin to be useful.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize, minimize_scalar


WORKSTREAM = Path(__file__).resolve().parent
sys.path.insert(0, str(WORKSTREAM))

from analyze_golay_ba_rm2sub_joint import DEFAULT_SELECTION, load_envelope  # noqa: E402
from diagnose_finite_k20_renyi_dense import (  # noqa: E402
    high_precision_transfer_log_moment,
    log_matrix_power_row_sum,
)


B = 250
L = 8448
N = B * L
EPOCHS = N // 128
DISTANCE = math.floor(0.11 * N)
THREE_BANDS = ((1, 60), (61, 189), (190, 249))
SEVEN_BANDS = (
    (1, 20),
    (21, 40),
    (41, 60),
    (61, 189),
    (190, 209),
    (210, 229),
    (230, 249),
)
DEFAULT_SPECTRUM = (
    WORKSTREAM / "bch250_parityfanout31x33_l64_packing_expected_envelope.json"
)


def load_spectrum(path: Path) -> np.ndarray:
    payload = json.loads(path.read_text(encoding="utf-8"))
    result = np.full(B + 1, -math.inf, dtype=np.float64)
    for row in payload["spectrum"]:
        weight = int(row["weight"])
        if weight == 0:
            continue
        result[weight] = float(row["log2_expected_multiplicity"]) * math.log(2.0)
    return result


def log_choose(length: int, weight: int) -> float:
    return math.lgamma(length + 1) - math.lgamma(weight + 1) - math.lgamma(
        length - weight + 1
    )


def band_majorant(
    spectrum: np.ndarray, lower: int, upper: int
) -> dict[str, float | int]:
    def objective(value_probability: float) -> float:
        return max(
            float(spectrum[weight])
            - log_choose(B, weight)
            - weight * math.log(value_probability)
            - (B - weight) * math.log1p(-value_probability)
            for weight in range(lower, upper + 1)
            if math.isfinite(float(spectrum[weight]))
        )

    result = minimize_scalar(
        objective,
        bounds=(max(1e-8, lower / B - 0.25), min(1 - 1e-8, upper / B + 0.25)),
        method="bounded",
        options={"xatol": 1e-13},
    )
    return {
        "lower_weight": lower,
        "upper_weight": upper,
        "value_probability": float(result.x),
        "log_majorant": objective(float(result.x)),
    }


def log_multinomial(counts: tuple[int, ...]) -> float:
    return math.lgamma(sum(counts) + 1) - sum(
        math.lgamma(count + 1) for count in counts
    )


def softmax(values: np.ndarray) -> np.ndarray:
    shifted = values - np.max(values)
    result = np.exp(shifted)
    return result / np.sum(result)


def integer_composition(fractions: np.ndarray) -> tuple[int, ...]:
    raw = fractions * L
    counts = np.floor(raw).astype(int)
    remainder = L - int(np.sum(counts))
    order = np.argsort(-(raw - counts))
    counts[order[:remainder]] += 1
    return tuple(int(value) for value in counts)


def sampled_compositions(
    band_count: int,
    random_per_concentration: int,
    pure_only: bool,
    edge_only: bool,
    central_spokes_only: bool,
    central_index: int | None,
) -> list[tuple[int, ...]]:
    count = band_count
    result: set[tuple[int, ...]] = set()
    for index in range(count):
        fractions = np.zeros(count)
        fractions[index] = 1.0
        result.add(integer_composition(fractions))
    if pure_only:
        return sorted(result)
    if central_spokes_only:
        if central_index is None:
            raise ValueError("central-spokes mode requires a central band")
        result = {tuple(L if index == central_index else 0 for index in range(count))}
        for tail_index in range(count):
            if tail_index == central_index:
                continue
            for tail_count in (1, 2, 4, 8, 16):
                counts = [0] * count
                counts[tail_index] = tail_count
                counts[central_index] = L - tail_count
                result.add(tuple(counts))
        return sorted(result)
    if edge_only:
        edge_counts = (1, 2, 4, 8, 16, 32, 64, 128, 256, 512)
        for left in range(count):
            for right in range(count):
                if left == right:
                    continue
                for left_count in edge_counts:
                    counts = [0] * count
                    counts[left] = left_count
                    counts[right] = L - left_count
                    result.add(tuple(counts))
        return sorted(result)
    for left in range(count):
        for right in range(left + 1, count):
            for share in (0.1, 0.25, 0.5, 0.75, 0.9):
                fractions = np.zeros(count)
                fractions[left] = share
                fractions[right] = 1.0 - share
                result.add(integer_composition(fractions))
    result.add(integer_composition(np.full(count, 1.0 / count)))
    generator = np.random.default_rng(0xB250F064)
    for concentration in (0.2, 0.5, 1.0, 2.0, 5.0):
        for _ in range(random_per_concentration):
            result.add(
                integer_composition(
                    generator.dirichlet(np.full(count, concentration))
                )
            )
    return sorted(result)


def optimize_composition(
    envelope,
    bands: list[dict[str, float | int]],
    counts: tuple[int, ...],
) -> dict[str, object]:
    active = [index for index, count in enumerate(counts) if count]
    active_counts = tuple(counts[index] for index in active)
    frequencies = np.asarray(active_counts, dtype=np.float64) / L
    values = np.asarray(
        [float(bands[index]["value_probability"]) for index in active]
    )
    outer_log = sum(
        count * float(bands[index]["log_majorant"])
        for index, count in enumerate(counts)
    )
    log_type_count = log_multinomial(counts)

    def evaluate(point: np.ndarray, high_precision: bool) -> float:
        probabilities = softmax(point[:-1])
        surprisal = math.exp(min(4.0, max(-8.0, float(point[-1]))))
        bit_probability = float(np.dot(probabilities, values))
        if high_precision:
            inner_log = high_precision_transfer_log_moment(
                envelope, bit_probability, surprisal, epochs=EPOCHS
            )
        else:
            transfer, _ = envelope.transfer(
                candidate_probability=2.0 * bit_probability,
                surprisal=surprisal,
            )
            if np.min(transfer) < -1e-11:
                return math.inf
            inner_log = log_matrix_power_row_sum(
                np.maximum(transfer, 0.0), EPOCHS
            )
        log_conditioning = log_type_count + sum(
            count * math.log(probability)
            for count, probability in zip(active_counts, probabilities)
        )
        reference_bad = min(
            0.0,
            inner_log + DISTANCE * surprisal - B * log_conditioning,
        )
        return log_type_count + outer_log + reference_bad

    start_logits = np.log(np.maximum(frequencies, 1e-12))
    starts = [
        np.concatenate((start_logits, [math.log(surprisal)]))
        for surprisal in (0.7, 1.3, 2.0, 2.7)
    ]
    best = None
    for start in starts:
        result = minimize(
            lambda point: evaluate(point, False),
            start,
            method="Nelder-Mead",
            options={"maxiter": 1000, "xatol": 1e-8, "fatol": 1e-7},
        )
        if best is None or float(result.fun) < float(best.fun):
            best = result
    assert best is not None
    value = evaluate(best.x, True)
    probabilities = softmax(best.x[:-1])
    full_probabilities = [0.0] * len(counts)
    for index, probability in zip(active, probabilities):
        full_probabilities[index] = float(probability)
    return {
        "counts": list(counts),
        "fractions": [count / L for count in counts],
        "log2_upper": value / math.log(2.0),
        "margin_bits": -value / math.log(2.0),
        "reference_type_probabilities": full_probabilities,
        "reference_bit_probability": float(np.dot(probabilities, values)),
        "surprisal": math.exp(min(4.0, max(-8.0, float(best.x[-1])))),
        "optimizer_success": bool(best.success),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectrum", type=Path, default=DEFAULT_SPECTRUM)
    parser.add_argument("--selection", type=Path, default=DEFAULT_SELECTION)
    parser.add_argument("--random-per-concentration", type=int, default=8)
    parser.add_argument(
        "--fine-seven-band",
        action="store_true",
        help="split each endpoint tail into three twenty-weight bands",
    )
    parser.add_argument(
        "--pure-only",
        action="store_true",
        help="check only compositions supported on a single band",
    )
    parser.add_argument(
        "--edge-only",
        action="store_true",
        help="check pure cases and powers-of-two distances from pairwise vertices",
    )
    parser.add_argument(
        "--central-spokes-only",
        action="store_true",
        help="check one, two, four, eight, or sixteen singleton-tail rows against the central band",
    )
    parser.add_argument(
        "--singleton-range",
        type=int,
        nargs=2,
        metavar=("LOWER", "UPPER"),
        help="replace the standard bands by singleton weights in this range",
    )
    parser.add_argument(
        "--symmetric-cutoff",
        type=int,
        help="use endpoint bands 1..T and (250-T)..249 with one central band",
    )
    parser.add_argument(
        "--sparse-singletons-through",
        type=int,
        help="exclude endpoints 1..16, use singleton bands through T, and one central band",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    spectrum = load_spectrum(args.spectrum)
    central_index = None
    if args.singleton_range is not None:
        lower, upper = args.singleton_range
        if not (1 <= lower <= upper < B):
            parser.error("singleton range must satisfy 1 <= LOWER <= UPPER < 250")
        band_limits = tuple((weight, weight) for weight in range(lower, upper + 1))
    elif args.sparse_singletons_through is not None:
        cutoff = args.sparse_singletons_through
        if not (17 <= cutoff < B // 2):
            parser.error("sparse singleton cutoff must satisfy 17 <= T < 125")
        low = tuple((weight, weight) for weight in range(17, cutoff + 1))
        central = ((cutoff + 1, B - cutoff - 1),)
        high = tuple(
            (weight, weight) for weight in range(B - cutoff, B - 16)
        )
        band_limits = low + central + high
        central_index = len(low)
    elif args.symmetric_cutoff is not None:
        cutoff = args.symmetric_cutoff
        if not (1 <= cutoff < B // 2):
            parser.error("symmetric cutoff must satisfy 1 <= T < 125")
        band_limits = (
            (1, cutoff),
            (cutoff + 1, B - cutoff - 1),
            (B - cutoff, B - 1),
        )
        central_index = 1
    else:
        band_limits = SEVEN_BANDS if args.fine_seven_band else THREE_BANDS
    bands = [band_majorant(spectrum, *limits) for limits in band_limits]
    envelope = load_envelope(args.selection, DISTANCE / N)
    compositions = sampled_compositions(
        len(bands),
        args.random_per_concentration,
        args.pure_only,
        args.edge_only,
        args.central_spokes_only,
        central_index,
    )
    rows = []
    for index, counts in enumerate(compositions, start=1):
        row = optimize_composition(envelope, bands, counts)
        rows.append(row)
        print(
            f"composition,{index},{len(compositions)},"
            f"margin,{row['margin_bits']:.9f},counts,{counts}",
            flush=True,
        )
    worst = min(rows, key=lambda row: float(row["margin_bits"]))
    payload = {
        "schema": "bch250-parityfanout-band-qL-diagnostic-v2",
        "status": "BINARY64_OPTIMIZATION_HIGH_PRECISION_FINAL_EVALUATION",
        "parameters": {
            "outer_bits": B,
            "outer_rows": L,
            "output_bits": N,
            "distance": DISTANCE,
            "relative_distance": DISTANCE / N,
            "spectrum": str(args.spectrum),
            "sampled_compositions": len(compositions),
            "pure_only": args.pure_only,
        },
        "bands": [
            {
                **band,
                "log2_majorant": float(band["log_majorant"]) / math.log(2.0),
            }
            for band in bands
        ],
        "worst_sampled": worst,
        "all_sampled_below_2^-40": all(
            float(row["margin_bits"]) > 40.0 for row in rows
        ),
        "rows": rows,
        "limitations": [
            "The spectrum is an expected mass-coupled envelope, not a fixed-code spectrum.",
            "The composition set is sampled and does not cover the integer simplex.",
            "Optimization uses binary64; final transfer evaluation uses 100 digits but is not outward rounded.",
            "No concentration theorem for one sampled fanout composition is supplied.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"wrote,{args.output}")
    print(f"worst_margin_bits,{worst['margin_bits']:.9f}")


if __name__ == "__main__":
    main()
