#!/usr/bin/env python3
"""Contract target Schur diagonals with exact accumulator run counts.

The input receipt must contain sectors zero, one, and two through every
active primal level.  The diagonal insertion lemma bounds all higher sectors
by sector two.  Positive semidefiniteness then bounds every cross-level entry
by the geometric mean of its two diagonal entries.

The run counts are exact integers.  Input diagonals and means are currently
nondirected floating-point diagnostics, so the result is not a certificate.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from analyze_dual_walk_and_accumulator_energy import accumulator_joint_counts
from probe_one_stage_variance_scaling import accumulator_shell_law
from verify_pair_kernel_small import krawtchouk


HERE = Path(__file__).resolve().parent
DEFAULT_DIAGONALS = HERE / "primal_schur_diagonal_K256_B512_r33_j0_2_p84_probe.json"
DEFAULT_OUTPUT = HERE / "primal_schur_shell_bound_K256_B512_r33_w42_probe.json"


def expected_spectrum(
    message_bits: int, output_bits: int, right_degree: int
) -> np.ndarray:
    denominator = math.comb(message_bits, right_degree)
    result = np.zeros(output_bits + 1, dtype=np.longdouble)
    for weight in range(message_bits + 1):
        bias = krawtchouk(message_bits, right_degree, weight) / denominator
        law = accumulator_shell_law(
            np.longdouble((1 - bias) / 2), output_bits
        )
        result += np.longdouble(math.comb(message_bits, weight)) * law
    return result


def evaluate(
    path: Path, overlays: list[Path], shell_weights: list[int]
) -> dict[str, object]:
    source = json.loads(path.read_text(encoding="utf-8"))
    parameters = source["result"]
    message_bits = int(parameters["message_bits"])
    output_bits = int(parameters["output_bits"])
    right_degree = int(parameters["right_degree"])
    if output_bits != 2 * message_bits:
        raise ValueError("the present contraction expects rate one half")

    entry_rows = list(parameters["entries"])
    for overlay in overlays:
        extra = json.loads(overlay.read_text(encoding="utf-8"))["result"]
        if any(
            int(extra[key]) != expected
            for key, expected in (
                ("message_bits", message_bits),
                ("output_bits", output_bits),
                ("right_degree", right_degree),
            )
        ):
            raise ValueError("overlay parameters do not match the base receipt")
        by_key = {
            (int(row["sector"]), int(row["level"])): row for row in entry_rows
        }
        for row in extra["entries"]:
            key = (int(row["sector"]), int(row["level"]))
            if key not in by_key:
                by_key[key] = row
                continue
            merged = dict(by_key[key])
            for field in (
                "scaled_uncentered_pair_sum",
                "scaled_positive_part",
                "scaled_centered_positive_part",
                "scaled_centered_rational_log_upper",
                "scaled_locally_centered_sum",
            ):
                if row.get(field) is not None:
                    merged[field] = row[field]
            by_key[key] = merged
        entry_rows = list(by_key.values())
    diagonals = {
        (int(row["sector"]), int(row["level"])): np.longdouble(
            row["scaled_diagonal"]
        )
        for row in entry_rows
    }
    monotone_upper_diagonals = {}
    for row in entry_rows:
        sector = int(row["sector"])
        level = int(row["level"])
        if sector == 0:
            value = row.get("scaled_centered_rational_log_upper")
            if value is None:
                value = row.get("scaled_centered_positive_part")
        elif sector == 1:
            value = row.get("scaled_positive_part")
        else:
            value = row["scaled_diagonal"]
        if value is not None:
            monotone_upper_diagonals[(sector, level)] = np.longdouble(value)
    run_counts = accumulator_joint_counts(output_bits)
    spectrum = expected_spectrum(message_bits, output_bits, right_degree)
    rows = []
    for shell_weight in shell_weights:
        counts = run_counts[shell_weight]
        terms = []
        monotone_terms = []
        maximizing_sectors = {}
        for level, count in sorted(counts.items()):
            available = range(min(2, level) + 1)
            values = []
            for sector in available:
                key = (sector, level)
                if key not in diagonals:
                    raise ValueError(f"missing diagonal sector={sector}, level={level}")
                values.append((diagonals[key], sector))
            diagonal, sector = max(values)
            if diagonal < 0:
                raise ArithmeticError("negative covariance diagonal")
            maximizing_sectors[str(level)] = sector
            terms.append(np.sqrt(np.longdouble(count) * diagonal))
            monotone_values = []
            for candidate_sector in available:
                key = (candidate_sector, level)
                if key not in monotone_upper_diagonals:
                    monotone_values = []
                    break
                monotone_values.append(monotone_upper_diagonals[key])
            if monotone_values:
                monotone_terms.append(
                    np.sqrt(np.longdouble(count) * max(monotone_values))
                )
        scaled_quadratic = sum(terms, np.longdouble(0)) ** 2
        mean = spectrum[shell_weight]
        ratio = np.ldexp(scaled_quadratic, -message_bits) / mean
        monotone_ratio = None
        if len(monotone_terms) == len(terms):
            monotone_quadratic = sum(monotone_terms, np.longdouble(0)) ** 2
            monotone_ratio = np.ldexp(
                monotone_quadratic, -message_bits
            ) / mean
        rows.append(
            {
                "shell_weight": shell_weight,
                "active_level_minimum": min(counts),
                "active_level_maximum": max(counts),
                "mean_log2": float(np.log2(mean)),
                "raw_signed_diagonal_diagnostic": float(ratio),
                "raw_signed_diagonal_diagnostic_log2": float(np.log2(ratio)),
                "monotone_variance_to_mean_upper": (
                    None if monotone_ratio is None else float(monotone_ratio)
                ),
                "monotone_variance_to_mean_upper_log2": (
                    None if monotone_ratio is None else float(np.log2(monotone_ratio))
                ),
                "maximizing_sectors": maximizing_sectors,
            }
        )
    return {
        "message_bits": message_bits,
        "output_bits": output_bits,
        "right_degree": right_degree,
        "source_diagonals": str(path.resolve()),
        "overlay_diagonals": [str(overlay.resolve()) for overlay in overlays],
        "shells": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--diagonals", type=Path, default=DEFAULT_DIAGONALS)
    parser.add_argument("--overlay", type=Path, action="append", default=[])
    parser.add_argument("--shell-weight", type=int, action="append", default=[])
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    shell_weights = args.shell_weight or [42]
    result = evaluate(args.diagonals, args.overlay, shell_weights)
    payload = {
        "schema": "pure-ea-primal-schur-diagonal-shell-bound-v1",
        "status": "EXACT_REDUCTION_WITH_NONDIRECTED_NUMERICAL_INPUTS",
        "result": result,
        "proved_reductions": [
            "At levels p <= n/2, every valid sector j >= 2 has diagonal at most the sector-2 diagonal at the same level.",
            "Positive semidefiniteness bounds each cross-level block entry by its diagonal geometric mean.",
            "Accumulator level multiplicities are exact run counts.",
        ],
        "scope": [
            "The source block diagonals use nondirected binary64 with long-double reduction.",
            "The expected spectrum uses nondirected long-double arithmetic.",
            "The displayed variance bound is diagnostic until both inputs use outward rounding.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
