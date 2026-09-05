#!/usr/bin/env python3
"""Directly audit an optimized equal-pair outer-support scan."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np


SCRIPT_DIRECTORY = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIRECTORY))

import analyze_riffle_double_parity as bch  # noqa: E402


def ratio_lookup(records: np.ndarray, ratio: int) -> int:
    keys = records[:, 0]
    index = int(np.searchsorted(keys, np.uint64(ratio)))
    if index < records.shape[0] and int(records[index, 0]) == ratio:
        return int(records[index, 1])
    return 0


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--histogram", type=Path, required=True)
    parser.add_argument("--optimized", type=Path, required=True)
    args = parser.parse_args()

    optimized = json.loads(args.optimized.read_text())
    blocks = int(optimized["data_blocks"])
    profile = str(optimized.get("profile", "22_22_24"))
    if profile not in ("22_22_24", "22_22_26", "24_24_22"):
        raise ValueError("unsupported equal-pair profile")
    unique_weight = {"22_22_24": 24, "22_22_26": 26, "24_24_22": 22}[profile]
    exact_key = {
        "22_22_24": "exact_finite_weight222224_outer_words",
        "22_22_26": "exact_finite_weight222226_outer_words",
        "24_24_22": "exact_finite_weight242422_outer_words",
    }[profile]
    if not 3 <= blocks <= 16:
        raise ValueError("direct audit is restricted to at most 16 data blocks")
    raw = np.memmap(args.histogram, dtype="<u8", mode="r")
    if raw.size % 2:
        raise RuntimeError("ratio histogram has a partial record")
    records = raw.reshape(-1, 2)
    if not np.all(records[1:, 0] > records[:-1, 0]):
        raise RuntimeError("ratio histogram is not sorted and unique")

    gamma = [1]
    for _ in range(1, blocks):
        gamma.append(bch.field_multiply_x(gamma[-1]))

    data_profile_supports = [0, 0, 0]
    data_outer_words = [0, 0, 0]
    for first in range(blocks - 2):
        for second in range(first + 1, blocks - 1):
            for third in range(second + 1, blocks):
                coefficients = [gamma[first], gamma[second], gamma[third]]
                for heavy in range(3):
                    lights = [index for index in range(3) if index != heavy]
                    numerator = coefficients[lights[0]] ^ coefficients[heavy]
                    denominator = coefficients[lights[1]] ^ coefficients[heavy]
                    ratio = bch.field_multiply(
                        numerator, bch.field_inverse(denominator)
                    )
                    count = ratio_lookup(records, ratio)
                    if count:
                        data_profile_supports[heavy] += 1
                        data_outer_words[heavy] += count

    p0_profile_supports = [0, 0, 0]
    p0_outer_words = [0, 0, 0]
    for first in range(blocks - 1):
        for second in range(first + 1, blocks):
            coefficients = [0, gamma[first], gamma[second]]
            for heavy in range(3):
                lights = [index for index in range(3) if index != heavy]
                numerator = coefficients[lights[0]] ^ coefficients[heavy]
                denominator = coefficients[lights[1]] ^ coefficients[heavy]
                ratio = bch.field_multiply(
                    numerator, bch.field_inverse(denominator)
                )
                count = ratio_lookup(records, ratio)
                if count:
                    p0_profile_supports[heavy] += 1
                    p0_outer_words[heavy] += count

    expected = {
        f"data_profile_supports_by_weight{unique_weight}_role": data_profile_supports,
        f"data_outer_words_by_weight{unique_weight}_role": data_outer_words,
        f"p0_profile_supports_by_weight{unique_weight}_role": p0_profile_supports,
        f"p0_outer_words_by_weight{unique_weight}_role": p0_outer_words,
    }
    for key, value in expected.items():
        if optimized[key] != value:
            raise RuntimeError(
                f"optimized/direct mismatch for {key}: {optimized[key]} != {value}"
            )
    total = sum(data_outer_words) + sum(p0_outer_words)
    if optimized[exact_key] != total:
        raise RuntimeError("optimized/direct total outer-word mismatch")

    payload = {
        "schema": f"riffle-shiftalpha64-weight{profile.replace('_', '')}-outer-audit-v1",
        "profile": profile,
        "data_blocks": blocks,
        **expected,
        exact_key: total,
        "validation": {
            "directly_enumerated_every_small_instance_support": True,
            f"all_three_weight{unique_weight}_roles_checked": True,
            "optimized_scan_matches_direct_enumeration": True,
        },
        "scope": (
            "Independent direct support audit of the optimized gap-normalized "
            "scan on the declared small data-block instance."
        ),
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
