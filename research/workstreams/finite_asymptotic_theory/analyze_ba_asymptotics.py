#!/usr/bin/env python3
"""Compare finite BA weighted spectra with their sparse boundary paths.

This script is diagnostic.  It uses the existing exact direct-sum polynomial
and binary64 accumulator recursion from
``scripts/analyze_designed_block_accumulate_spectrum.py``.  It writes no
files and runs configurations sequentially.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY / "scripts"))

from analyze_designed_block_accumulate_spectrum import analyze, read_spectrum  # noqa: E402


@dataclass(frozen=True)
class Configuration:
    name: str
    spectrum: Path
    local_length: int
    local_dimension: int


CONFIGURATIONS = {
    "ebch32": Configuration(
        "EBCH [32,16,8]",
        REPOSITORY / "scripts" / "ebch32_16_delta8_spectrum.csv",
        32,
        16,
    ),
    "xbch64": Configuration(
        "shortened-XBCH [64,32,12]",
        REPOSITORY / "scripts" / "xbch64_32_philips_spectrum.csv",
        64,
        32,
    ),
    "ebch128": Configuration(
        "EBCH [128,64,22]",
        REPOSITORY / "scripts" / "ebch128_64_spectrum.csv",
        128,
        64,
    ),
}


def log2_comb(n: int, k: int) -> float:
    if k < 0 or k > n:
        return -math.inf
    return (
        math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
    ) / math.log(2.0)


def accumulator_log2_probability(block_size: int, input_weight: int, output_weight: int) -> float:
    half_down = input_weight // 2
    half_up = (input_weight + 1) // 2
    return (
        log2_comb(output_weight - 1, half_up - 1)
        + log2_comb(block_size - output_weight, half_down)
        - log2_comb(block_size, input_weight)
    )


def minimum_weight_data(configuration: Configuration) -> tuple[int, int]:
    spectrum = read_spectrum(configuration.spectrum, configuration.local_length)
    for weight, count in enumerate(spectrum[1:], start=1):
        if count:
            return weight, count
    raise ValueError(f"{configuration.name} has no nonzero words")


def boundary_path(
    configuration: Configuration,
    block_size: int,
    accumulators: int,
    rho: float,
) -> dict[str, object]:
    minimum_weight, local_count = minimum_weight_data(configuration)
    blocks = block_size // configuration.local_length
    log2_mass = math.log2(blocks * local_count)
    path = [minimum_weight]
    current = minimum_weight
    for _ in range(accumulators):
        following = (current + 1) // 2
        log2_mass += accumulator_log2_probability(block_size, current, following)
        path.append(following)
        current = following
    log2_mass += current * math.log2(rho)
    polynomial_exponent = sum(path[1:]) - 1
    return {
        "minimum_weight": minimum_weight,
        "local_minimum_multiplicity": local_count,
        "path": path,
        "predicted_polynomial_exponent": polynomial_exponent,
        "log2_weighted_lower_bound": log2_mass,
    }


def full_weighted_mass(result: dict[str, object], rho: float) -> float:
    terms = []
    for weight, row in enumerate(result["spectrum"]):
        if weight == 0 or row["log2_expected_multiplicity"] is None:
            continue
        terms.append(row["log2_expected_multiplicity"] + weight * math.log2(rho))
    if not terms:
        return -math.inf
    return float(np.logaddexp2.reduce(np.asarray(terms, dtype=np.float64)))


def slope(rows: list[dict[str, object]], key: str) -> float | None:
    finite = [row for row in rows if math.isfinite(float(row[key]))]
    if len(finite) < 2:
        return None
    tail = finite[-min(3, len(finite)) :]
    x = np.log2(np.asarray([row["block_size"] for row in tail], dtype=np.float64))
    y = np.asarray([row[key] for row in tail], dtype=np.float64)
    return float(np.polyfit(x, y, 1)[0])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--configuration",
        choices=sorted(CONFIGURATIONS),
        default="ebch128",
    )
    parser.add_argument("--blocks", nargs="+", type=int, default=[4, 8, 16, 32])
    parser.add_argument("--accumulators", type=int, default=2)
    parser.add_argument("--rho", nargs="+", type=float, default=[0.05, 0.1, 0.2])
    args = parser.parse_args()

    configuration = CONFIGURATIONS[args.configuration]
    if not configuration.spectrum.is_file():
        raise FileNotFoundError(
            f"authenticated spectrum is unavailable: {configuration.spectrum}"
        )
    if any(blocks <= 0 for blocks in args.blocks):
        raise ValueError("block counts must be positive")
    if any(not 0.0 < rho < 1.0 for rho in args.rho):
        raise ValueError("rho values must lie in (0,1)")

    rows_by_rho: dict[str, list[dict[str, object]]] = {
        str(rho): [] for rho in args.rho
    }
    for blocks in args.blocks:
        block_size = blocks * configuration.local_length
        result = analyze(
            configuration.spectrum,
            configuration.name,
            configuration.local_length,
            configuration.local_dimension,
            blocks,
            args.accumulators,
        )
        for rho in args.rho:
            boundary = boundary_path(
                configuration, block_size, args.accumulators, rho
            )
            full = full_weighted_mass(result, rho)
            rows_by_rho[str(rho)].append(
                {
                    "block_size": block_size,
                    "full_log2_weighted_mass": full,
                    "boundary_log2_lower_bound": boundary[
                        "log2_weighted_lower_bound"
                    ],
                    "full_minus_boundary_bits": full
                    - float(boundary["log2_weighted_lower_bound"]),
                    "first_weight_expected_at_least_1": result["summary"][
                        "first_weight_expected_at_least_1"
                    ],
                    "boundary": boundary,
                }
            )

    output = {
        "status": "binary64 diagnostic",
        "configuration": {
            "key": args.configuration,
            "name": configuration.name,
            "local_length": configuration.local_length,
            "local_dimension": configuration.local_dimension,
            "spectrum": str(configuration.spectrum.relative_to(REPOSITORY)),
            "accumulators": args.accumulators,
        },
        "results": {
            rho: {
                "rows": rows,
                "tail_loglog_slope_full": slope(rows, "full_log2_weighted_mass"),
                "tail_loglog_slope_boundary": slope(
                    rows, "boundary_log2_lower_bound"
                ),
            }
            for rho, rows in rows_by_rho.items()
        },
    }
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
