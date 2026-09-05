#!/usr/bin/env python3
"""Build the artifact manifest for the finite wrapped-EC certificate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


WORKSTREAM = Path(__file__).resolve().parent
OUTPUT = WORKSTREAM / "FINITE_K20_EXPANDER_EC_D10_M15_MANIFEST.json"
RECEIPT = WORKSTREAM / "expander_ec_d10_m15_k20_wrapper_outward.json"
ARTIFACT_NAMES = [
    "FINITE_K20_EXPANDER_EC_D10_M15_CERTIFICATE.md",
    "EXPANDER_CODE_ALTERNATIVE_AUDIT.md",
    "expander_ec_d10_m15_k20_parent_candidate.json",
    "expander_ec_d10_m15_k20_wrapper_outward.json",
    "generate_expander_ec_k20_candidate.py",
    "certify_expander_ec_k20_wrapper.py",
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    receipt = json.loads(RECEIPT.read_text(encoding="utf-8"))
    if receipt.get("status") != "OUTWARD_FINITE_DISTANCE_CERTIFICATE":
        raise AssertionError("wrapper receipt is not accepting")
    claim = receipt["claim"]
    if not claim.get("comparison_to_2^-40"):
        raise AssertionError("wrapper receipt does not close 40 bits")
    payload = {
        "schema": "finite-k20-expander-ec-d10-m15-manifest-v1",
        "status": "OUTWARD_FINITE_DISTANCE_CERTIFICATE",
        "claim": claim,
        "construction": {
            "family": "binary two-sided regular Expand--Convolute",
            "parent": "[2097170,1048585] ensemble with degrees 10/5 and wrapped-convolution memory 15",
            "wrapper": "shorten nine input coordinates and puncture eighteen output coordinates",
        },
        "artifacts": [
            {
                "path": f"workstreams/finite_asymptotic_theory/{name}",
                "sha256": sha256(WORKSTREAM / name),
                "size_bytes": (WORKSTREAM / name).stat().st_size,
            }
            for name in ARTIFACT_NAMES
        ],
        "external_verifier_dependencies": receipt["verifier_dependencies"],
        "reproduction": [
            "python workstreams/finite_asymptotic_theory/generate_expander_ec_k20_candidate.py",
            "python workstreams/finite_asymptotic_theory/certify_expander_ec_k20_wrapper.py",
            "python workstreams/finite_asymptotic_theory/build_expander_ec_k20_manifest.py",
        ],
        "limitations": receipt["limitations"],
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(f"output={OUTPUT}")


if __name__ == "__main__":
    main()
