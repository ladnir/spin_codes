#!/usr/bin/env python3
"""Build the manifest for the one-shot random-constituent certificate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import flint


WORKSTREAM = Path(__file__).resolve().parent
REPO_ROOT = WORKSTREAM.parents[1]
OUTPUT = WORKSTREAM / "FINITE_K20_ONE_SHOT_RANDOM512_RANDOMSTEP_CONV_MANIFEST.json"
FILES = [
    WORKSTREAM / "FINITE_K20_ONE_SHOT_RANDOM512_RANDOMSTEP_CONV_CERTIFICATE.md",
    WORKSTREAM / "evaluate_single_random_constituent_highprob_renyi.py",
    WORKSTREAM / "evaluate_single_random_constituent_highprob_q1.py",
    WORKSTREAM / "evaluate_single_random_constituent_highprob_q2.py",
    WORKSTREAM / "evaluate_single_random_constituent_shared_two_band.py",
    WORKSTREAM / "evaluate_single_random_constituent_shared_two_band_dense.py",
    WORKSTREAM / "diagnose_single_random_constituent_shared_two_band_cover.py",
    WORKSTREAM / "certify_single_random_constituent_dense_outward.py",
    WORKSTREAM / "certify_single_random_constituent_sparse_outward.py",
    WORKSTREAM / "certify_single_random_constituent_combined.py",
    WORKSTREAM / "single_random_constituent_B512_highprob_q1_s22.json",
    WORKSTREAM / "single_random_constituent_B512_highprob_q2_s22.json",
    WORKSTREAM / "single_random_constituent_B512_shared_two_band_q3_16_s22.json",
    WORKSTREAM / "single_random_constituent_B512_shared_two_band_q17_32_s22.json",
    WORKSTREAM / "single_random_constituent_B512_shared_two_band_q33_64_s22.json",
    WORKSTREAM / "single_random_constituent_B512_shared_two_band_q65_128_s22.json",
    WORKSTREAM / "single_random_constituent_B512_shared_two_band_q129_159_s22.json",
    WORKSTREAM / "single_random_constituent_B512_shared_two_band_dense_cover_s22.json",
    WORKSTREAM / "single_random_constituent_B512_sparse_outward_s22.json",
    WORKSTREAM / "single_random_constituent_B512_shared_two_band_dense_cover_outward_s22.json",
    WORKSTREAM / "single_random_constituent_B512_combined_outward_s22.json",
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    receipt = json.loads(FILES[-1].read_text(encoding="utf-8"))
    payload = {
        "schema": "single-random-constituent-randomstepconv-manifest-v1",
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
            "python workstreams/finite_asymptotic_theory/certify_single_random_constituent_dense_outward.py",
            "python workstreams/finite_asymptotic_theory/certify_single_random_constituent_sparse_outward.py",
            "python workstreams/finite_asymptotic_theory/certify_single_random_constituent_combined.py",
            "python workstreams/finite_asymptotic_theory/build_single_random_constituent_manifest.py",
        ],
        "runtime": {
            "python_flint_version": flint.__version__,
            "arithmetic": "192-bit Arb plus directed positive binary64 recurrences",
        },
        "sampler": {
            "outer_random_bits": 256 * 512,
            "acceptance_test": None,
            "distribution": "uniform binary 256-by-512 matrix",
        },
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"output={OUTPUT}")


if __name__ == "__main__":
    main()
