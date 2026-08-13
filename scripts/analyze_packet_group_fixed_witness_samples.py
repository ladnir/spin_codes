#!/usr/bin/env python3
"""Stress-test one optimized generalized-g witness on adversarial profiles."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import gammaln

from packet_group_drive_stratified import profile_classes, profile_count
from packet_group_outer_profile import (
    D,
    K,
    N,
    atom_count,
    normalization_log2,
    total_weight_outer,
)


HAMMING_BALL_LOG2_UPPER = (
    gammaln(N + 1)
    - gammaln(D + 1)
    - gammaln(N - D + 1)
    + math.log(N - D + 1)
    - math.log(N - 2 * D + 1)
) / math.log(2.0)


def largest_remainder(rows: np.ndarray, total: int) -> np.ndarray:
    raw = rows * total
    result = np.floor(raw).astype(np.int64)
    missing = total - np.sum(result, axis=1)
    fractions = raw - result
    for index, count in enumerate(missing):
        if count:
            order = np.argpartition(fractions[index], -count)[-count:]
            result[index, order] += 1
    return result


def sample_profiles(group_bits: int, seed: int, per_family: int) -> np.ndarray:
    rng = np.random.default_rng(seed + group_bits)
    dimension = group_bits + 1
    rows = [np.eye(dimension)]
    weights = np.arange(dimension, dtype=np.float64)
    for probability in np.linspace(0.02, 0.98, 25):
        distribution = np.asarray(
            [
                math.comb(group_bits, weight)
                * probability**weight
                * (1.0 - probability) ** (group_bits - weight)
                for weight in range(dimension)
            ]
        )
        rows.append(distribution[None, :] / np.sum(distribution))
    for concentration in (0.03, 0.1, 0.3, 1.0):
        rows.append(
            rng.dirichlet(np.full(dimension, concentration), per_family)
        )
    return largest_remainder(np.vstack(rows), atom_count(group_bits))


def load_rows(paths: list[Path]) -> list[dict]:
    merged: dict[int, dict] = {}
    for path in paths:
        data = json.loads(path.read_text(encoding="utf-8"))
        for row in data if isinstance(data, list) else [data]:
            group_bits = int(row["group_bits"])
            merged[group_bits] = {**merged.get(group_bits, {}), **row}
    return list(merged.values())


def load_all_rows(paths: list[Path]) -> list[dict]:
    rows = []
    for path in paths:
        data = json.loads(path.read_text(encoding="utf-8"))
        rows.extend(data if isinstance(data, list) else [data])
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, action="append", required=True)
    parser.add_argument("--extra-witness", type=Path, action="append", default=[])
    parser.add_argument("--pure-curve", type=Path, action="append", default=[])
    parser.add_argument("--groups", default="1,2,4,8,16,32,64")
    parser.add_argument("--per-family", type=int, default=500)
    parser.add_argument("--seed", type=int, default=20260811)
    parser.add_argument("--show-worst", type=int, default=5)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    selected = {int(value) for value in args.groups.split(",")}
    source = load_rows(args.input)
    extra_rows = load_all_rows(args.extra_witness) if args.extra_witness else []
    pure_curves = []
    for path in args.pure_curve:
        data = json.loads(path.read_text(encoding="utf-8"))
        pure_curves.extend(data if isinstance(data, list) else [data])
    reports = []
    for row in sorted(source, key=lambda item: item["group_bits"]):
        group_bits = int(row["group_bits"])
        if group_bits not in selected:
            continue
        anchor = np.asarray(row["profile"], dtype=np.float64)
        log_variables = np.asarray(
            row["outer_details"]["log_variables"], dtype=np.float64
        ) / math.log(2.0)
        log_fugacities = np.log2(np.asarray(row["fugacities"], dtype=np.float64))
        charge = log_variables + log_fugacities
        outer_base = float(row["outer_profile_log2"]) + float(anchor @ log_variables)
        inner_base = float(row["inner_details"]["inner_mgf_log2"]) - D * math.log2(
            float(row["pole"])
        )
        constant = outer_base + inner_base
        profiles = sample_profiles(
            group_bits, args.seed, args.per_family
        )
        normalizations = normalization_log2(group_bits, profiles)
        fixed_values = constant - profiles @ charge - normalizations
        fixed_columns = [fixed_values]
        for extra in extra_rows:
            if int(extra["group_bits"]) != group_bits:
                continue
            if extra.get("outer_type") == "total_weight":
                extra_log_fugacities = np.log2(
                    np.asarray(extra["fugacities"], dtype=np.float64)
                )
                physical_weights = profiles @ np.arange(
                    group_bits + 1, dtype=np.int64
                )
                cache = {
                    int(weight): total_weight_outer(int(weight))[0]
                    for weight in np.unique(physical_weights)
                }
                outer_values = np.asarray(
                    [cache[int(weight)] for weight in physical_weights]
                )
                inner_constant = float(
                    extra["inner_details"]["inner_mgf_log2"]
                ) - D * math.log2(float(extra["pole"]))
                fixed_columns.append(
                    outer_values
                    + inner_constant
                    - profiles @ extra_log_fugacities
                    - normalizations
                )
                continue
            extra_anchor = np.asarray(extra["profile"], dtype=np.float64)
            if "spectrum_outer_log2" in extra:
                outer_point = np.asarray(extra["outer_point"], dtype=np.float64)
                extra_log_variables_natural = np.concatenate(
                    ([0.0], outer_point[:group_bits])
                )
                extra_outer_profile = float(extra["spectrum_outer_log2"])
                replacement_is_included = False
            else:
                extra_log_variables_natural = np.asarray(
                    extra["outer_details"]["log_variables"], dtype=np.float64
                )
                extra_outer_profile = float(extra["outer_profile_log2"])
                replacement_is_included = extra.get("outer_type") == "linear_bl"
            replacement = (
                0.0
                if replacement_is_included
                else 128.0 * max(
                    extra_log_variables_natural[new]
                    - extra_log_variables_natural[old]
                    for old in range(group_bits + 1)
                    for new in range(group_bits + 1)
                    if abs(new - old) <= 1
                )
                / math.log(2.0)
            )
            extra_log_variables = extra_log_variables_natural / math.log(2.0)
            extra_log_fugacities = np.log2(
                np.asarray(extra["fugacities"], dtype=np.float64)
            )
            extra_constant = (
                extra_outer_profile
                + float(extra_anchor @ extra_log_variables)
                + replacement
                + float(extra["inner_details"]["inner_mgf_log2"])
                - D * math.log2(float(extra["pole"]))
            )
            fixed_columns.append(
                extra_constant
                - profiles @ (extra_log_variables + extra_log_fugacities)
                - normalizations
            )
        fixed_values = np.min(np.vstack(fixed_columns), axis=0)
        full_bijection = (
            N * math.log2(11) - (N - D) * math.log2(10) + K
        ) - normalizations
        classes = np.asarray(profile_classes(group_bits), dtype=np.float64)
        stars_bars = np.sum(
            gammaln(profiles + classes)
            - gammaln(profiles + 1.0)
            - gammaln(classes),
            axis=1,
        ) / math.log(2.0)
        inverse_orbit = K + HAMMING_BALL_LOG2_UPPER - normalizations + stars_bars
        values = np.minimum(fixed_values, np.minimum(full_bijection, inverse_orbit))
        leaders = np.argmin(
            np.vstack((fixed_values, full_bijection, inverse_orbit)), axis=0
        )
        # The all-zero physical profile contains only the zero outer word,
        # which is excluded from the nonzero distance union.
        pure_zero = np.all(profiles[:, 1:] == 0, axis=1)
        values[pure_zero] = -np.inf
        total_weight = profiles @ np.arange(group_bits + 1, dtype=np.int64)
        values[total_weight < 21] = -np.inf
        for curve in pure_curves:
            if int(curve["group_bits"]) != group_bits:
                continue
            pure_rows = {int(item["class_weight"]): item for item in curve["rows"]}
            for index, profile in enumerate(profiles):
                active = np.flatnonzero(profile)
                if len(active) == 1 and int(active[0]) in pure_rows:
                    values[index] = min(
                        values[index],
                        float(pure_rows[int(active[0])]["combined_log2"]),
                    )
        target = -40.0 - math.log2(profile_count(group_bits, N))
        worst_indices = np.argsort(values)[-args.show_worst :][::-1]
        report = {
            "group_bits": group_bits,
            "samples": len(profiles),
            "target_log2": target,
            "covered": int(np.sum(values <= target)),
            "uncovered": int(np.sum(values > target)),
            "worst": [
                {
                    "value_log2": float(values[index]),
                    "profile": profiles[index].tolist(),
                    "support": int(np.sum(profiles[index] > 0)),
                    "leader": ("shared_drive", "full_bijection", "inverse_orbit")[
                        int(leaders[index])
                    ],
                }
                for index in worst_indices
            ],
        }
        reports.append(report)
        print(
            f"g={group_bits} covered={report['covered']}/{report['samples']} "
            f"worst={report['worst'][0]['value_log2']:.6f} "
            f"target={target:.6f}",
            flush=True,
        )
    rendered = json.dumps(reports, indent=2)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    print("status=DIAGNOSTIC_FIXED_WITNESS_ADVERSARIAL_PROFILE_SAMPLES")


if __name__ == "__main__":
    main()
