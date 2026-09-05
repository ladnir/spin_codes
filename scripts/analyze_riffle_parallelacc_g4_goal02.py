#!/usr/bin/env python3
"""Evaluate the full-size support-33 moment bound for ParallelAcc g=4."""

from __future__ import annotations

import hashlib
import json
import math
import subprocess
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "explorations" / "riffle_dp_g4_g1_support33_refutation.json"
BINARY = ROOT / "scripts" / "analyze_riffle_parallelacc_g4_goal02_moment.exe"
CPP_SOURCE = ROOT / "scripts" / "analyze_riffle_parallelacc_g4_goal02_moment.cpp"
OUTPUT = ROOT / "constructions" / "riffle_parallelacc_g4" / "receipts" / "goal02_support33_moment.json"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def packet_values(codeword_hex: str) -> list[int]:
    word = int(codeword_hex, 16)
    values = []
    for packet in range(32):
        value = (word >> (4 * packet)) & 15
        if value:
            values.append(value)
    return values


def main() -> None:
    source = json.loads(SOURCE.read_text())
    profiles: dict[tuple[int, ...], dict[str, object]] = {}
    family_profiles = []
    for receipt in source["receipts"]:
        values = []
        for local_word in receipt["local_words"]:
            values.extend(packet_values(local_word["codeword_hex"]))
        counts = Counter(values)
        profile = tuple(counts[value] for value in range(1, 16))
        if sum(profile) != 33:
            raise RuntimeError("support-33 source changed")
        if profile not in profiles:
            profiles[profile] = {
                "profile_id": f"p{len(profiles):02d}",
                "counts": profile,
                "families": [],
            }
        profiles[profile]["families"].append(receipt["family_id"])
        family_profiles.append((receipt["family_id"], profiles[profile]["profile_id"]))

    lines = []
    for profile in profiles.values():
        lines.append(profile["profile_id"] + " " + " ".join(map(str, profile["counts"])))
    completed = subprocess.run(
        [str(BINARY)],
        input="\n".join(lines) + "\n",
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode:
        raise RuntimeError(f"moment executable failed\n{completed.stdout}{completed.stderr}")

    results = {}
    for line in completed.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) != 9:
            raise RuntimeError(f"unexpected moment output: {line}")
        profile_id = fields[0]
        results[profile_id] = {
            "support": int(fields[1]),
            "mixed_radix_states": int(fields[2]),
            "objective_evaluations": int(fields[3]),
            "log2_probability_upper_bound": float(fields[4]),
            "a_equals_minus_log_z": float(fields[5]),
            "b_equals_minus_log_tau": float(fields[6]),
            "scaled_gradient_a": float(fields[7]),
            "scaled_gradient_b": float(fields[8]),
        }

    rows = []
    aggregate = 0.0
    for profile in profiles.values():
        profile_id = profile["profile_id"]
        result = results[profile_id]
        multiplicity = len(profile["families"])
        log_bound = result["log2_probability_upper_bound"]
        aggregate += multiplicity * (2.0 ** log_bound if log_bound < 0 else 1.0)
        rows.append(
            {
                "profile_id": profile_id,
                "packet_value_counts": {
                    str(value): count
                    for value, count in enumerate(profile["counts"], 1)
                    if count
                },
                "outer_families": profile["families"],
                "outer_family_count": multiplicity,
                **result,
            }
        )

    payload = {
        "schema": "riffle-parallelacc-g4-goal02-support33-moment-v1",
        "candidate": "Riffle ParallelAcc g=4",
        "evidence_label": "DIAGNOSTIC_NUMERICAL_UPPER_BOUND",
        "source_sha256": digest(SOURCE),
        "cpp_source_sha256": digest(CPP_SOURCE),
        "executable_sha256": digest(BINARY),
        "probability_space": (
            "For each fixed authenticated outer word, the 33 labeled nonzero packets "
            "and all zero packets undergo the uniform global packet permutation."
        ),
        "bound": {
            "packet_count": 524352,
            "binary_distance_threshold_inclusive": 188766,
            "outer_word_count": len(family_profiles),
            "unique_packet_multisets": len(rows),
            "aggregate_log2_upper_bound": math.log2(aggregate),
            "target_log2": -40.0,
            "passes_target_diagnostically": aggregate <= 2.0 ** -40,
        },
        "profiles": rows,
        "method": (
            "Two-parameter Chernoff and coefficient bound. A mixed-radix dynamic "
            "program sums the tilted prefix-XOR weight over every labeled ordering."
        ),
        "scope_limitation": (
            "Long-double evaluation and optimization are diagnostic. The receipt "
            "covers only the 26 authenticated support-33 outer words."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"outer_words={len(family_profiles)}")
    print(f"unique_profiles={len(rows)}")
    print(f"aggregate_log2_upper_bound={payload['bound']['aggregate_log2_upper_bound']:.12f}")
    print(f"passes_target_diagnostically={payload['bound']['passes_target_diagnostically']}")
    print(f"output={OUTPUT}")
    print("status=DIAGNOSTIC_NUMERICAL_UPPER_BOUND")


if __name__ == "__main__":
    main()
