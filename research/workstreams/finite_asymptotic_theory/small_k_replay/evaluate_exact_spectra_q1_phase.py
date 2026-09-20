#!/usr/bin/env python3
"""Sweep exact constituent spectra over total dimension and inner memory.

Each case repeats one fixed rate-half constituent in every outer row.  The
row-coordinate permutations, region permutations, and RandomStepConv maps
are sampled independently once and fixed for the resulting code.

The calculation covers occupation Q=1.  It uses the exact constituent weight
spectrum and the exact uniform-routing coefficient in the proof model.  A
degree-one matrix-polynomial power computes a region with exactly one active
position in O(log L) matrix operations.

Nearest binary64 log arithmetic makes every result diagnostic.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.special import logsumexp


HERE = Path(__file__).resolve().parent
WORKSTREAM = HERE.parent
REPOSITORY = WORKSTREAM.parent.parent
sys.path.insert(0, str(WORKSTREAM))

import evaluate_ebch128_randomstepconv_g1 as inner  # noqa: E402
from evaluate_ebch128_randomstepconv_q1_exact import uniform_coefficients  # noqa: E402


LOG2 = math.log(2.0)
DEFAULT_OUTPUT = HERE / "exact_spectra_q1_small_k_memory_phase.json"


@dataclass(frozen=True)
class Constituent:
    name: str
    block_bits: int
    dimension: int
    minimum_distance: int
    spectrum_path: Path


CONSTITUENTS = {
    "ebch8": Constituent(
        "extended BCH [8,4,4]",
        8,
        4,
        4,
        HERE / "spectra" / "ebch8_4_spectrum.csv",
    ),
    "ebch32": Constituent(
        "extended BCH [32,16,8]",
        32,
        16,
        8,
        WORKSTREAM / "ebch32_16_delta8_spectrum.csv",
    ),
    "ebch128": Constituent(
        "extended BCH [128,64,22]",
        128,
        64,
        22,
        REPOSITORY / "scripts" / "EBCH128_64.wd",
    ),
    "xbch64": Constituent(
        "Philips shortened-XBCH [64,32,12]",
        64,
        32,
        12,
        HERE / "spectra" / "xbch64_32_philips_spectrum.csv",
    ),
    "rm13": Constituent(
        "RM(1,3) [8,4,4]",
        8,
        4,
        4,
        HERE / "spectra" / "ebch8_4_spectrum.csv",
    ),
    "rm25": Constituent(
        "RM(2,5) [32,16,8]",
        32,
        16,
        8,
        WORKSTREAM / "ebch32_16_delta8_spectrum.csv",
    ),
    "rm37": Constituent(
        "RM(3,7) [128,64,16]",
        128,
        64,
        16,
        HERE / "spectra" / "rm37_128_64_spectrum.csv",
    ),
    "rm49": Constituent(
        "RM(4,9) [512,256,32]",
        512,
        256,
        32,
        REPOSITORY / "scripts" / "rm512_256_spectrum.csv",
    ),
}


def load_spectrum(constituent: Constituent) -> dict[int, int]:
    path = constituent.spectrum_path
    first = path.read_text(encoding="utf-8").lstrip().splitlines()[0]
    if "," in first:
        with path.open(newline="", encoding="utf-8") as source:
            spectrum = {
                int(row["weight"]): int(row["count"])
                for row in csv.DictReader(source)
            }
    else:
        spectrum = {}
        for line in path.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            weight, count = line.split()
            spectrum[int(weight)] = int(count)
    if sum(spectrum.values()) != 1 << constituent.dimension:
        raise AssertionError(f"{constituent.name}: spectrum mass mismatch")
    if spectrum.get(0) != 1:
        raise AssertionError(f"{constituent.name}: missing unique zero word")
    observed_distance = min(weight for weight, count in spectrum.items() if weight and count)
    if observed_distance != constituent.minimum_distance:
        raise AssertionError(f"{constituent.name}: minimum distance mismatch")
    return spectrum


def log_poly_mul_degree_one(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    result = np.full((2, 2, 2), -math.inf)
    result[0] = inner.log_matmul(left[0], right[0])
    result[1] = np.logaddexp(
        inner.log_matmul(left[0], right[1]),
        inner.log_matmul(left[1], right[0]),
    )
    return result


def one_active_region_coefficients(
    zero: np.ndarray, active: np.ndarray, region_length: int
) -> tuple[np.ndarray, np.ndarray]:
    """Return average transfer matrices for zero or one active position."""
    result = np.full((2, 2, 2), -math.inf)
    result[0] = inner.log_identity()
    power = np.stack((inner.log_entries(zero), inner.log_entries(active)))
    exponent = region_length
    while exponent:
        if exponent & 1:
            result = log_poly_mul_degree_one(result, power)
        exponent >>= 1
        if exponent:
            power = log_poly_mul_degree_one(power, power)
    # Coefficient one is the sum over the L possible active positions.
    # Uniform routing averages those positions.
    result[1] -= math.log(region_length)
    return result[0], result[1]


def tilt_grid(args: argparse.Namespace) -> list[float]:
    coarse = np.arange(
        args.grid_min, args.grid_max + args.grid_step / 2, args.grid_step
    )
    fine = np.arange(
        args.fine_min, args.fine_max + args.fine_step / 2, args.fine_step
    )
    return sorted(set(map(float, coarse)) | set(map(float, fine)))


def evaluate_case(
    *,
    constituent: Constituent,
    spectrum: dict[int, int],
    message_exponent: int,
    memory_bits: int,
    distance_numerator: int,
    distance_denominator: int,
    tilts: list[float],
) -> dict[str, object]:
    message_bits = 1 << message_exponent
    if message_bits % constituent.dimension:
        raise ValueError("constituent dimension must divide total message bits")
    outer_rows = message_bits // constituent.dimension
    output_bits = 2 * message_bits
    if output_bits != outer_rows * constituent.block_bits:
        raise AssertionError("rate-half length identity failed")
    distance_cutoff = (
        distance_numerator * output_bits + distance_denominator - 1
    ) // distance_denominator
    weights = sorted(weight for weight, count in spectrum.items() if weight and count)
    best = {weight: math.inf for weight in weights}
    witnesses = {weight: math.nan for weight in weights}

    for log_surprisal in tilts:
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        zero, active = inner.step_matrices(z, memory_bits)
        zero_region, one_region = one_active_region_coefficients(
            zero, active, outer_rows
        )
        coordinates = uniform_coefficients(
            zero_region,
            one_region,
            constituent.block_bits,
            constituent.block_bits,
        )
        moments = np.logaddexp(
            coordinates[:, 0, 0], coordinates[:, 0, 1]
        )
        correction = distance_cutoff * surprisal
        for weight in weights:
            value = (
                math.log(outer_rows)
                + math.log(spectrum[weight])
                + min(0.0, float(moments[weight]) + correction)
            )
            if value < best[weight]:
                best[weight] = value
                witnesses[weight] = log_surprisal

    aggregate_log2 = float(logsumexp(list(best.values())) / LOG2)
    dominant_weight = max(weights, key=lambda weight: best[weight])
    return {
        "constituent": constituent.name,
        "message_exponent": message_exponent,
        "message_bits": message_bits,
        "output_bits": output_bits,
        "outer_rows": outer_rows,
        "memory_bits": memory_bits,
        "distance_cutoff": distance_cutoff,
        "relative_distance_lower": distance_cutoff / output_bits,
        "q1_log2_upper_diagnostic": aggregate_log2,
        "q1_margin_bits_diagnostic": -aggregate_log2,
        "q1_closes_40_bits_diagnostic": aggregate_log2 < -40,
        "dominant_weight": dominant_weight,
        "dominant_multiplicity": str(spectrum[dominant_weight]),
        "dominant_pointwise_log2": best[dominant_weight] / LOG2,
        "dominant_log_surprisal": witnesses[dominant_weight],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--constituents",
        nargs="+",
        choices=sorted(CONSTITUENTS),
        default=["ebch32", "ebch128", "rm49"],
    )
    parser.add_argument(
        "--message-exponents", type=int, nargs="+", default=list(range(8, 21))
    )
    parser.add_argument(
        "--memory-bits", type=int, nargs="+", default=[4, 6, 8, 10, 12, 16, 22]
    )
    parser.add_argument("--distance-numerator", type=int, default=109)
    parser.add_argument("--distance-denominator", type=int, default=1000)
    parser.add_argument("--grid-min", type=float, default=-12.0)
    parser.add_argument("--grid-max", type=float, default=1.0)
    parser.add_argument("--grid-step", type=float, default=0.5)
    parser.add_argument("--fine-min", type=float, default=-10.0)
    parser.add_argument("--fine-max", type=float, default=-5.0)
    parser.add_argument("--fine-step", type=float, default=0.1)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if not 0 < args.distance_numerator < args.distance_denominator:
        parser.error("distance fraction must lie in (0,1)")
    tilts = tilt_grid(args)

    selected = [CONSTITUENTS[name] for name in args.constituents]
    spectra = {item.name: load_spectrum(item) for item in selected}
    cases = []
    for constituent in selected:
        minimum_exponent = constituent.dimension.bit_length() - 1
        for message_exponent in args.message_exponents:
            if message_exponent < minimum_exponent:
                continue
            for memory_bits in args.memory_bits:
                print(
                    f"case,{constituent.name},kexp,{message_exponent},memory,{memory_bits}",
                    flush=True,
                )
                cases.append(
                    evaluate_case(
                        constituent=constituent,
                        spectrum=spectra[constituent.name],
                        message_exponent=message_exponent,
                        memory_bits=memory_bits,
                        distance_numerator=args.distance_numerator,
                        distance_denominator=args.distance_denominator,
                        tilts=tilts,
                    )
                )

    payload = {
        "schema": "exact-spectrum-small-k-randomstepconv-q1-phase-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "probability_space": {
            "outer": "one fixed deterministic constituent, repeated in every outer row",
            "routing": "independent uniform coordinate permutation per row and independent uniform position permutation per transposed region",
            "inner": "independent RandomStepConv map at every position, sampled once and shared by all messages",
        },
        "parameters": {
            "message_exponents": args.message_exponents,
            "memory_bits": args.memory_bits,
            "distance_numerator": args.distance_numerator,
            "distance_denominator": args.distance_denominator,
            "tilt_count": len(tilts),
        },
        "spectrum_sources": [
            {
                "constituent": item.name,
                "path": str(item.spectrum_path),
                "sha256": hashlib.sha256(item.spectrum_path.read_bytes()).hexdigest(),
            }
            for item in selected
        ],
        "cases": cases,
        "limitations": [
            "Every value uses nearest binary64 log arithmetic and is diagnostic.",
            "The finite tilt grid supplies valid candidate witnesses but is not outward rounded.",
            "Only occupation Q=1 is covered; a full distance certificate must cover every occupation.",
            "The inner is RandomStepConv, not RM2Sub.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    summary = [
        {
            "constituent": row["constituent"],
            "kexp": row["message_exponent"],
            "memory": row["memory_bits"],
            "margin": row["q1_margin_bits_diagnostic"],
            "dominant_weight": row["dominant_weight"],
        }
        for row in cases
        if row["q1_closes_40_bits_diagnostic"]
    ]
    print(json.dumps({"q1_closing_cases": summary}, indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
