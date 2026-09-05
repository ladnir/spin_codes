#!/usr/bin/env python3
"""Use an exact Poisson-binomial comparison on thin mixed faces.

The fixed BCH250 outer is repeated in every row.  Each row then receives an
independent draw of the layered fanout wrapper.  The expected row counting
measure is therefore the input spectrum in this diagnostic, and products of
row expectations are valid.  This probability space differs from one
fanout composition reused in every row.

For a fixed band composition, the row references are independent Bernoulli
variables with band-dependent probabilities.  A region permutation makes
their law exchangeable.  The pointwise comparison with iid Bernoulli bits of
the same mean is exactly the maximum ratio between a Poisson-binomial mass
and the matching binomial mass.  This avoids the earlier categorical-type
conditioning loss.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.signal import fftconvolve
from scipy.special import gammaln, logsumexp
from scipy.stats import binom


WORKSTREAM = Path(__file__).resolve().parent
sys.path.insert(0, str(WORKSTREAM))

from analyze_golay_ba_rm2sub_joint import DEFAULT_SELECTION, load_envelope  # noqa: E402
from diagnose_bch250_parityfanout_three_band_qL import (  # noqa: E402
    B,
    DISTANCE,
    EPOCHS,
    L,
    N,
    band_majorant,
    load_spectrum,
    log_multinomial,
)
from diagnose_finite_k20_renyi_dense import (  # noqa: E402
    high_precision_transfer_log_moment,
    log_matrix_power_row_sum,
)


BANDS = ((1, 16), (17, 233), (234, 249))
DEFAULT_SPECTRUM = (
    WORKSTREAM / "bch250_parityfanout31x33_l128_packing_expected_envelope.json"
)


def log_choose_vector(length: int) -> np.ndarray:
    values = np.arange(length + 1, dtype=np.float64)
    return gammaln(length + 1.0) - gammaln(values + 1.0) - gammaln(length - values + 1.0)


def poisson_binomial_log_ratio(
    defect_count: int, defect_probability: float
) -> tuple[float, int, float]:
    if defect_count == 0:
        return 0.0, 0, 0.5
    central_count = L - defect_count
    mean_probability = (
        0.5 * central_count + defect_probability * defect_count
    ) / L
    log_choose_central = log_choose_vector(central_count)
    defect_indices = np.arange(defect_count + 1, dtype=np.float64)
    log_defect = (
        log_choose_vector(defect_count)
        + defect_indices * math.log(defect_probability)
        + (defect_count - defect_indices) * math.log1p(-defect_probability)
    )
    log_choose_total = log_choose_vector(L)
    log_binomial = (
        log_choose_total
        + np.arange(L + 1) * math.log(mean_probability)
        + (L - np.arange(L + 1)) * math.log1p(-mean_probability)
    )
    ratios = np.full(L + 1, -math.inf, dtype=np.float64)
    central_normalization = central_count * math.log(0.5)
    for total in range(L + 1):
        lower = max(0, total - central_count)
        upper = min(defect_count, total)
        if lower > upper:
            continue
        indices = np.arange(lower, upper + 1, dtype=int)
        central_weights = total - indices
        log_pb = logsumexp(
            log_defect[indices]
            + log_choose_central[central_weights]
            + central_normalization
        )
        ratios[total] = log_pb - log_binomial[total]
    maximizing_weight = int(np.argmax(ratios))
    return float(ratios[maximizing_weight]), maximizing_weight, mean_probability


def grouped_poisson_binomial_log_ratio(
    counts: tuple[int, int, int], probabilities: tuple[float, float, float]
) -> tuple[float, int, float]:
    occupation = sum(counts)
    if not 1 <= occupation <= L:
        raise ValueError("composition occupation must lie in [1,L]")
    mean_probability = sum(
        count * probability for count, probability in zip(counts, probabilities)
    ) / L
    mass = np.asarray([1.0], dtype=np.float64)
    for count, probability in zip(counts, probabilities):
        if count == 0:
            continue
        component = binom.pmf(np.arange(count + 1), count, probability)
        mass = fftconvolve(mass, component)
        mass = np.maximum(mass, 0.0)
        mass /= float(np.sum(mass))
    padded_mass = np.zeros(L + 1, dtype=np.float64)
    padded_mass[: occupation + 1] = mass
    mass = padded_mass
    weights = np.arange(L + 1, dtype=np.float64)
    log_reference = (
        log_choose_vector(L)
        + weights * math.log(mean_probability)
        + (L - weights) * math.log1p(-mean_probability)
    )
    log_ratio = np.full(L + 1, -math.inf, dtype=np.float64)
    # FFT convolution leaves roundoff-sized positive noise in remote tails.
    # The threshold is diagnostic only; the final theorem needs a rigorous
    # ratio bound or an outward recurrence.
    positive = mass > float(np.max(mass)) * 1e-12
    log_ratio[positive] = np.log(mass[positive]) - log_reference[positive]
    endpoint_zero = sum(
        count * math.log1p(-probability)
        for count, probability in zip(counts, probabilities)
    ) - L * math.log1p(-mean_probability)
    endpoint_one = -math.inf
    if occupation == L:
        endpoint_one = sum(
            count * math.log(probability)
            for count, probability in zip(counts, probabilities)
        ) - L * math.log(mean_probability)
    log_ratio[0] = endpoint_zero
    log_ratio[-1] = endpoint_one
    maximizing_weight = int(np.argmax(log_ratio))
    return float(log_ratio[maximizing_weight]), maximizing_weight, mean_probability


def fast_inner_log(envelope, bit_probability: float, surprisal: float) -> float:
    transfer, _ = envelope.transfer(
        candidate_probability=2.0 * bit_probability,
        surprisal=surprisal,
    )
    if np.min(transfer) < -1e-12:
        return math.inf
    return log_matrix_power_row_sum(np.maximum(transfer, 0.0), EPOCHS)


def evaluate_edge(
    envelope,
    bands: list[dict[str, float | int]],
    defect_band: int,
    defect_count: int,
) -> dict[str, object]:
    central_count = L - defect_count
    defect_probability = float(bands[defect_band]["value_probability"])
    central_probability = float(bands[1]["value_probability"])
    if abs(central_probability - 0.5) > 1e-8:
        raise ValueError("central reference is not half Bernoulli")
    log_ratio, maximizing_weight, bit_probability = poisson_binomial_log_ratio(
        defect_count, defect_probability
    )
    counts = [0, central_count, 0]
    counts[defect_band] = defect_count
    outer_log = (
        central_count * float(bands[1]["log_majorant"])
        + defect_count * float(bands[defect_band]["log_majorant"])
    )
    type_log = log_multinomial(tuple(counts))

    def objective(surprisal: float) -> float:
        inner = fast_inner_log(envelope, bit_probability, surprisal)
        return (
            type_log
            + outer_log
            + B * log_ratio
            + min(0.0, inner + DISTANCE * surprisal)
        )

    result = minimize_scalar(
        objective,
        bounds=(0.05, 20.0),
        method="bounded",
        options={"xatol": 1e-12},
    )
    surprisal = float(result.x)
    high_precision_inner = high_precision_transfer_log_moment(
        envelope,
        bit_probability,
        surprisal,
        epochs=EPOCHS,
    )
    value = (
        type_log
        + outer_log
        + B * log_ratio
        + min(0.0, high_precision_inner + DISTANCE * surprisal)
    )
    return {
        "counts": counts,
        "margin_bits": -value / math.log(2.0),
        "log2_upper": value / math.log(2.0),
        "defect_probability": defect_probability,
        "reference_bit_probability": bit_probability,
        "poisson_binomial_log2_ratio_per_region": log_ratio / math.log(2.0),
        "poisson_binomial_maximizing_weight": maximizing_weight,
        "surprisal": surprisal,
        "optimizer_success": bool(result.success),
    }


def evaluate_composition(
    envelope,
    bands: list[dict[str, float | int]],
    counts: tuple[int, int, int],
) -> dict[str, object]:
    probabilities = tuple(float(band["value_probability"]) for band in bands)
    log_ratio, maximizing_weight, bit_probability = (
        grouped_poisson_binomial_log_ratio(counts, probabilities)
    )
    outer_log = sum(
        count * float(band["log_majorant"])
        for count, band in zip(counts, bands)
    )
    occupation = sum(counts)
    type_log = (
        math.lgamma(L + 1)
        - math.lgamma(L - occupation + 1)
        - sum(math.lgamma(count + 1) for count in counts)
    )

    def objective(surprisal: float) -> float:
        inner = fast_inner_log(envelope, bit_probability, surprisal)
        return (
            type_log
            + outer_log
            + B * log_ratio
            + min(0.0, inner + DISTANCE * surprisal)
        )

    result = minimize_scalar(
        objective,
        bounds=(0.05, 20.0),
        method="bounded",
        options={"xatol": 1e-12},
    )
    surprisal = float(result.x)
    high_precision_inner = high_precision_transfer_log_moment(
        envelope, bit_probability, surprisal, epochs=EPOCHS
    )
    value = (
        type_log
        + outer_log
        + B * log_ratio
        + min(0.0, high_precision_inner + DISTANCE * surprisal)
    )
    return {
        "counts": list(counts),
        "margin_bits": -value / math.log(2.0),
        "log2_upper": value / math.log(2.0),
        "reference_bit_probability": bit_probability,
        "poisson_binomial_log2_ratio_per_region": log_ratio / math.log(2.0),
        "poisson_binomial_maximizing_weight": maximizing_weight,
        "surprisal": surprisal,
        "optimizer_success": bool(result.success),
    }


def sampled_compositions() -> list[tuple[int, int, int]]:
    result: set[tuple[int, int, int]] = {
        (L, 0, 0),
        (0, L, 0),
        (0, 0, L),
    }
    for left in range(3):
        for right in range(left + 1, 3):
            for share in (0.1, 0.25, 0.5, 0.75, 0.9):
                left_count = int(round(L * share))
                counts = [0, 0, 0]
                counts[left] = left_count
                counts[right] = L - left_count
                result.add(tuple(counts))
    for defect_count in (1, 2, 4, 8, 16, 32, 64, 128, 256, 512):
        for defect_band in (0, 2):
            counts = [0, L - defect_count, 0]
            counts[defect_band] = defect_count
            result.add(tuple(counts))
    generator = np.random.default_rng(0xB250128)
    for concentration in (0.2, 0.5, 1.0, 2.0, 5.0):
        for _ in range(8):
            raw = generator.dirichlet(np.full(3, concentration)) * L
            counts = np.floor(raw).astype(int)
            remainder = L - int(np.sum(counts))
            order = np.argsort(-(raw - counts))
            counts[order[:remainder]] += 1
            result.add(tuple(int(value) for value in counts))
    return sorted(result)


def sampled_compositions_at_occupation(occupation: int) -> list[tuple[int, int, int]]:
    result: set[tuple[int, int, int]] = {
        (occupation, 0, 0),
        (0, occupation, 0),
        (0, 0, occupation),
    }
    for left in range(3):
        for right in range(left + 1, 3):
            for share in (0.1, 0.25, 0.5, 0.75, 0.9):
                left_count = int(round(occupation * share))
                counts = [0, 0, 0]
                counts[left] = left_count
                counts[right] = occupation - left_count
                result.add(tuple(counts))
    for defect_count in (1, 2, 4, 8, 16):
        if defect_count >= occupation:
            continue
        for defect_band in (0, 2):
            counts = [0, occupation - defect_count, 0]
            counts[defect_band] = defect_count
            result.add(tuple(counts))
    return sorted(result)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectrum", type=Path, default=DEFAULT_SPECTRUM)
    parser.add_argument("--selection", type=Path, default=DEFAULT_SELECTION)
    parser.add_argument("--sample-compositions", action="store_true")
    parser.add_argument("--sample-occupations", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    spectrum = load_spectrum(args.spectrum)
    bands = [band_majorant(spectrum, *limits) for limits in BANDS]
    envelope = load_envelope(args.selection, DISTANCE / N)
    rows = []
    if args.sample_occupations:
        occupations = (
            101, 128, 192, 256, 384, 512, 768, 1024,
            1536, 2048, 3072, 4096, 6144, 8192, L,
        )
        compositions = [
            counts
            for occupation in occupations
            for counts in sampled_compositions_at_occupation(occupation)
        ]
        for index, counts in enumerate(compositions, start=1):
            row = evaluate_composition(envelope, bands, counts)
            rows.append(row)
            print(
                f"composition,{index},{len(compositions)},counts,{counts},"
                f"margin,{row['margin_bits']:.9f},"
                f"ratio,{row['poisson_binomial_log2_ratio_per_region']:.12f}",
                flush=True,
            )
    elif args.sample_compositions:
        compositions = sampled_compositions()
        for index, counts in enumerate(compositions, start=1):
            row = evaluate_composition(envelope, bands, counts)
            rows.append(row)
            print(
                f"composition,{index},{len(compositions)},counts,{counts},"
                f"margin,{row['margin_bits']:.9f},"
                f"ratio,{row['poisson_binomial_log2_ratio_per_region']:.12f}",
                flush=True,
            )
    else:
        for defect_count in (1, 2, 4, 8, 16):
            for defect_band in (0, 2):
                row = evaluate_edge(envelope, bands, defect_band, defect_count)
                rows.append(row)
                print(
                    f"counts,{tuple(row['counts'])},margin,{row['margin_bits']:.9f},"
                    f"ratio,{row['poisson_binomial_log2_ratio_per_region']:.12f}",
                    flush=True,
                )
    worst = min(rows, key=lambda row: float(row["margin_bits"]))
    payload = {
        "schema": "bch250-rowlocal-fanout-poisson-binomial-edges-v1",
        "status": "BINARY64_DIAGNOSTIC_WITH_HIGH_PRECISION_FINAL_TRANSFER",
        "parameters": {
            "outer_bits": B,
            "outer_rows": L,
            "output_bits": N,
            "distance": DISTANCE,
            "relative_distance": DISTANCE / N,
            "bands": [list(band) for band in BANDS],
            "spectrum": str(args.spectrum),
            "selection": str(args.selection),
            "fanout_law": "independent layered fanout composition in every row",
            "sampled_compositions": len(rows),
        },
        "band_majorants": [
            {
                **band,
                "log2_majorant": float(band["log_majorant"]) / math.log(2.0),
            }
            for band in bands
        ],
        "worst_sampled": worst,
        "rows": rows,
        "limitations": [
            "Only edge compositions with one, two, four, eight, or sixteen defect rows are evaluated.",
            "The layered fanout transition and Poisson-binomial ratios use nearest binary64 arithmetic.",
            "General-composition FFT ratios discard masses below 1e-12 of the peak to avoid roundoff tails; this is diagnostic, not a certificate rule.",
            "The final RM2Sub transfer uses 100 decimal digits but is not outward rounded.",
            "A complete proof must cover every composition and every occupation.",
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
