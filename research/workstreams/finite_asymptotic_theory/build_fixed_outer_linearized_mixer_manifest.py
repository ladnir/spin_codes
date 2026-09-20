#!/usr/bin/env python3
"""Build the manifest for the fixed-outer linearized-mixer certificate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import flint


WORKSTREAM = Path(__file__).resolve().parent
REPO_ROOT = WORKSTREAM.parents[1]
OUTPUT = WORKSTREAM / "FINITE_K20_EBCH128X4_LINEARIZED_MIXER_RANDOMSTEP_CONV_MANIFEST.json"
FILES = [
    WORKSTREAM / "FINITE_K20_EBCH128X4_LINEARIZED_MIXER_RANDOMSTEP_CONV_CERTIFICATE.md",
    WORKSTREAM / "certify_fixed_outer_linearized_mixer_combined.py",
    WORKSTREAM / "certify_single_random_constituent_sparse_outward.py",
    WORKSTREAM / "certify_single_random_constituent_dense_outward.py",
    WORKSTREAM / "single_random_constituent_B512_sparse_outward_s22.json",
    WORKSTREAM / "single_random_constituent_B512_shared_two_band_dense_cover_outward_s22.json",
    WORKSTREAM / "fixed_outer_linearized_mixer_B512_combined_outward_s22.json",
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    receipt = json.loads(FILES[-1].read_text(encoding="utf-8"))
    payload = {
        "schema": "fixed-outer-linearized-mixer-randomstepconv-manifest-v1",
        "status": receipt["status"],
        "claim": receipt["claim"],
        "parameters": receipt["parameters"],
        "artifacts": [
            {
                "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
                "sha256": sha256(path),
                "size_bytes": path.stat().st_size,
            }
            for path in FILES
        ],
        "reproduction": [
            "python workstreams/finite_asymptotic_theory/certify_single_random_constituent_sparse_outward.py",
            "python workstreams/finite_asymptotic_theory/certify_single_random_constituent_dense_outward.py",
            "python workstreams/finite_asymptotic_theory/certify_fixed_outer_linearized_mixer_combined.py",
            "python workstreams/finite_asymptotic_theory/build_fixed_outer_linearized_mixer_manifest.py",
        ],
        "runtime": {
            "python_flint_version": flint.__version__,
            "arithmetic": "exact integers and rationals plus 192-bit Arb intervals",
        },
        "sampler": {
            "mixer_random_bits": 1024,
            "acceptance_test": None,
            "distribution": "two independent uniform elements of GF(2^512)",
        },
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"output={OUTPUT}")


if __name__ == "__main__":
    main()
