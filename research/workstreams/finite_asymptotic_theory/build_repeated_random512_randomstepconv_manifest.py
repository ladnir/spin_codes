#!/usr/bin/env python3
"""Build the hash-bound repeated-random512 RandomStepConv manifest."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


WORKSTREAM = Path(__file__).resolve().parent
OUTPUT = WORKSTREAM / "FINITE_K20_REPEATED_RANDOM512_RANDOMSTEP_CONV_MANIFEST.json"
FILES = [
    "FINITE_K20_REPEATED_RANDOM512_RANDOMSTEP_CONV_CERTIFICATE.md",
    "evaluate_repeated_random512_randomstepconv_g1.py",
    "repeated_random512_randomstepconv_g1_s30_allq_d11.json",
    "certify_repeated_random512_randomstepconv_g1_outward.py",
    "build_repeated_random512_randomstepconv_manifest.py",
    "repeated_random512_randomstepconv_g1_s30_allq_outward_d11.json",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    receipt_path = WORKSTREAM / FILES[-1]
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    claim = receipt["claim"]
    event = receipt["constituent_event"]
    parameters = receipt["parameters"]
    assertions = {
        "outward_status": receipt.get("status") == "OUTWARD_DISTANCE_CERTIFICATE",
        "message_dimension": parameters.get("message_bits") == 1 << 20,
        "literal_eleven_percent": (
            parameters.get("distance_cutoff") == 233_165
            and claim.get("relative_distance_lower", 0.0) >= 0.11
        ),
        "all_occupations": claim.get("occupation_count") == 4096,
        "conditional_below_2^-108": (
            claim.get("conditional_bad_probability_log2_upper", 0.0) < -108.0
        ),
        "bounded_setup_below_2^-40": (
            claim.get("bounded_setup_attempts") == 28
            and claim.get("unconditional_abort_or_bad_log2_upper", 0.0) < -40.0
            and claim.get("unconditional_comparison_to_2^-40") is True
        ),
        "single_sample_success_above_0.63": (
            event.get("success_probability_lower", 0.0) > 0.63
        ),
        "one_repeated_constituent": (
            "repeat only that one code"
            in receipt["probability_space"]["outer_selection"]
        ),
    }
    if not all(assertions.values()):
        raise AssertionError(f"manifest assertion failed: {assertions}")

    files = []
    for name in FILES:
        path = WORKSTREAM / name
        files.append(
            {
                "path": name,
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    payload = {
        "schema": "finite-k20-repeated-random512-randomstepconv-manifest-v1",
        "status": "HASH_BOUND_OUTWARD_CERTIFICATE",
        "assertions": assertions,
        "claim": claim,
        "constituent_event": event,
        "files": files,
        "reproduction": [
            "python workstreams/finite_asymptotic_theory/evaluate_repeated_random512_randomstepconv_g1.py",
            "python workstreams/finite_asymptotic_theory/certify_repeated_random512_randomstepconv_g1_outward.py",
            "python workstreams/finite_asymptotic_theory/build_repeated_random512_randomstepconv_manifest.py",
        ],
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": payload["status"], "claim": claim}, indent=2))
    print(f"output={OUTPUT}")


if __name__ == "__main__":
    main()
