#!/usr/bin/env python3
"""Bind the certified wrapped-EC parent to the exact k=2^20 wrapper."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys


WORKSTREAM = Path(__file__).resolve().parent
DEFAULT_EXPANDER_ROOT = Path(
    "C:/Users/peter/.codex/worktrees/30ff/permute_conv/expander_codes"
)
PARENT = WORKSTREAM / "expander_ec_d10_m15_k20_parent_candidate.json"
OUTPUT = WORKSTREAM / "expander_ec_d10_m15_k20_wrapper_outward.json"

TARGET_K = 1 << 20
TARGET_N = 1 << 21
TARGET_D = 228_590


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expander-root", type=Path, default=DEFAULT_EXPANDER_ROOT)
    parser.add_argument("--parent", type=Path, default=PARENT)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    scripts = args.expander_root.resolve() / "scripts"
    if not scripts.is_dir():
        parser.error(f"missing expander-code scripts directory: {scripts}")
    sys.path.insert(0, str(scripts))

    from binary_biregular_ec_certificate import verify_certificate
    from check_binary_biregular_ec_certificate import validate_certificate_structure

    parent = json.loads(args.parent.read_text(encoding="utf-8"))
    validate_certificate_structure(parent)
    expected = {
        "k": 1_048_585,
        "n": 2_097_170,
        "cutoff": 228_607,
        "left_degree": 10,
        "right_degree": 5,
        "memory": 15,
        "target_bits": 40,
    }
    if parent["parameters"] != expected:
        raise AssertionError("unexpected parent parameters")
    if len(parent["exact_weights"]) != 32:
        raise AssertionError("parent must use exact supports through 32")

    result = verify_certificate(parent)
    if not result.success or not bool(result.security_bits.lower() > 50):
        raise AssertionError("parent certificate does not exceed 50 bits")

    shortened_inputs = expected["k"] - TARGET_K
    punctured_outputs = expected["n"] - TARGET_N
    required_parent_distance = TARGET_D + punctured_outputs
    proved_parent_distance = expected["cutoff"] + 1
    if shortened_inputs != 9 or punctured_outputs != 18:
        raise AssertionError("unexpected wrapper dimensions")
    if proved_parent_distance != required_parent_distance:
        raise AssertionError("parent distance does not exactly cover puncturing loss")

    dependency_names = [
        "binary_biregular_ec_certificate.py",
        "check_binary_biregular_ec_certificate.py",
        "binary_biregular_diagnostic.py",
        "ea_certificate.py",
        "regular_ec_certificate.py",
        "regular_ec_diagnostic.py",
        "expander_bounds.py",
    ]
    dependencies = [
        {
            "path": str((scripts / name).resolve()),
            "sha256": sha256(scripts / name),
        }
        for name in dependency_names
    ]
    payload = {
        "schema": "expander-ec-d10-m15-k20-wrapper-outward-v1",
        "status": "OUTWARD_FINITE_DISTANCE_CERTIFICATE",
        "claim": {
            "message_bits": TARGET_K,
            "output_bits": TARGET_N,
            "minimum_distance_at_least": TARGET_D,
            "relative_distance_lower": TARGET_D / TARGET_N,
            "failure_probability_upper": str(result.total_bound.upper()),
            "failure_margin_bits_lower": str(result.security_bits.lower()),
            "comparison_to_2^-40": True,
        },
        "parent": {
            "parameters": expected,
            "minimum_distance_at_least_on_good_event": proved_parent_distance,
            "certificate": str(args.parent.resolve()),
            "certificate_sha256": sha256(args.parent),
            "outward_contributions": {
                "exact_supports_1_through_32": str(result.exact_bound.upper()),
                "outer_blocks": str(result.intermediate_bound.upper()),
                "central_blocks": str(result.dense_bound.upper()),
                "full_support": str(result.full_support_bound.upper()),
            },
        },
        "wrapper": {
            "shortened_input_coordinates": shortened_inputs,
            "punctured_output_coordinates": punctured_outputs,
            "distance_loss_upper": punctured_outputs,
            "dimension_argument": "positive wrapped distance implies that puncturing cannot create a kernel word",
            "randomness": "the wrapper is deterministic and adds no failure probability",
        },
        "ensemble": {
            "expander": "binary two-sided regular with left degree 10, right degree 5, and ten independently sampled regional permutations",
            "inner": "independently sampled wrapped binary convolution of memory 15",
            "reuse": "all expander permutations and convolution coefficients are sampled once and shared by every message",
            "linearity": "both the expander matrix and wrapped recurrence are binary linear maps",
            "encoding_complexity": "ordinary and transposed maps use O(k) word operations for fixed degrees and memory",
        },
        "verifier_dependencies": dependencies,
        "limitations": [
            "This is an Expand--Convolute code, not the structured SPIN construction.",
            "The theorem proves minimum distance and does not provide a decoding theorem.",
            "Linear ordinary and transposed encoding time is an operation-count statement; no implementation benchmark has been run.",
            "The verifier currently depends on the cited expander-code worktree paths and hashes.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
