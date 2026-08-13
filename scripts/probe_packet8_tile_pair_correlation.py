#!/usr/bin/env python3
"""Diagnose exact one-block correlations between two outer tile factors.

For fixed packet-profile fugacities, a tile factor is the product of the
64-row lane moment ``R_w`` over the coordinates in one EBCH band.  Conditional
on the one message block shared by two cross-band tiles, the other 63 rows are
independent.  Hence the conditional tile expectation depends only on the
shared block's projected band weight.  The exact [85,64] and [86,64] split
spectra then give the exact binary64 pair moment (up to final floating log
evaluation).

This is a diagnostic for the size and sign of the correlation discarded by
the global exponent-three Finner bound; it is not itself a global inequality.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

from probe_packet8_three_band_profile_enumerator import (
    BINOMIAL_LOG_PROBABILITIES,
    group_log_moments,
    optimize_outer,
)
from probe_packet8_weight_profile_scalar import parse_profile


ROOT = Path(__file__).resolve().parent.parent
PAIR01 = ROOT / "out" / "ebch85_band01_split_spectrum.csv"
PAIR12 = ROOT / "out" / "ebch86_band12_split_spectrum.csv"


def load_spectrum(path: Path) -> list[tuple[int, int, int]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = [
            (
                int(row["band0_weight"]),
                int(row["band1_weight"]),
                int(row["count"]),
            )
            for row in csv.DictReader(handle)
        ]
    if sum(count for _a, _b, count in rows) != 1 << 64:
        raise SystemExit(f"tile pair: spectrum mass mismatch in {path}")
    return rows


def logadd(log_x: float, log_y: float) -> float:
    return float(np.logaddexp(log_x, log_y))


def conditional_logs(group_logs: np.ndarray, power: float) -> tuple[float, float]:
    """Logs of E[R_(Bin(63)+bit)^power] for bit zero and one."""

    log_probabilities = np.array(
        [math.log(math.comb(63, weight)) - 63 * math.log(2.0) for weight in range(64)]
    )
    zero = logsumexp(log_probabilities + power * group_logs[:64])
    one = logsumexp(log_probabilities + power * group_logs[1:])
    return float(zero), float(one)


def tile_log_mean(group_logs: np.ndarray, band_size: int, power: float) -> float:
    one_column = logsumexp(BINOMIAL_LOG_PROBABILITIES + power * group_logs)
    return float(band_size * one_column)


def pair_log_moment(
    spectrum: list[tuple[int, int, int]],
    sizes: tuple[int, int],
    conditional: tuple[float, float],
) -> float:
    log_zero, log_one = conditional
    terms = []
    for weight0, weight1, count in spectrum:
        terms.append(
            math.log(count)
            - 64 * math.log(2.0)
            + (sizes[0] - weight0) * log_zero
            + weight0 * log_one
            + (sizes[1] - weight1) * log_zero
            + weight1 * log_one
        )
    return float(logsumexp(np.array(terms)))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", required=True)
    parser.add_argument(
        "--powers", default="1,1.5,2,3", help="comma-separated tile-factor powers"
    )
    args = parser.parse_args()
    try:
        profile = parse_profile(args.profile)
        powers = [float(value) for value in args.powers.split(",")]
    except ValueError as error:
        raise SystemExit(f"tile pair: {error}") from error

    coefficient_log, result, components = optimize_outer(profile)
    _value, _log_s3, group_logs = components
    pair01 = load_spectrum(PAIR01)
    pair12 = load_spectrum(PAIR12)
    print("packet-8 exact one-block tile-pair correlation probe")
    print(f"profile={','.join(map(str, profile))}")
    print(f"optimizer_success={result.success} message={result.message}")
    print(f"outer_profile_count_log2_upper={coefficient_log/math.log(2.0):.12f}")
    for power in powers:
        conditional = conditional_logs(group_logs, power)
        for label, sizes, spectrum in (
            ("01", (42, 43), pair01),
            ("12", (43, 43), pair12),
        ):
            joint = pair_log_moment(spectrum, sizes, conditional)
            independent = sum(
                tile_log_mean(group_logs, size, power) for size in sizes
            )
            print(
                f"power={power:.6g} pair={label} "
                f"log2_joint={joint/math.log(2.0):.12f} "
                f"log2_independent={independent/math.log(2.0):.12f} "
                f"log2_correlation={(joint-independent)/math.log(2.0):.12f}"
            )
    print("status=DIAGNOSTIC_EXACT_LOCAL_SPLIT_FLOATING_LOGS")


if __name__ == "__main__":
    main()
