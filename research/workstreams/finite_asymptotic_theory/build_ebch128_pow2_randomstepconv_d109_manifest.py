#!/usr/bin/env python3
"""Build the hash manifest for the power-of-two 10.9% certificate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import flint


WORKSTREAM = Path(__file__).resolve().parent
REPO_ROOT = WORKSTREAM.parents[1]
OUTPUT = WORKSTREAM / "FINITE_K20_EBCH128_POW2_RANDOMSTEP_CONV_D109_MANIFEST.json"
RECEIPT = WORKSTREAM / "ebch128_randomstepconv_g1_s20_pow2_outward_d109.json"
FILES = [
    REPO_ROOT / "scripts" / "ebch128x64_forward_paar_circuit.json",
    REPO_ROOT / "scripts" / "EBCH128_64.wd",
    WORKSTREAM / "FINITE_K20_EBCH128_OUTER_DEFINITION.md",
    WORKSTREAM / "FINITE_K20_EBCH128_POW2_RANDOMSTEP_CONV_D109_CERTIFICATE.md",
    WORKSTREAM / "evaluate_ebch128_randomstepconv_g1.py",
    WORKSTREAM / "evaluate_ebch128_randomstepconv_parity_pivot.py",
    WORKSTREAM / "evaluate_ebch128_randomstepconv_q1_exact.py",
    WORKSTREAM / "evaluate_ebch128_randomstepconv_dense_complement.py",
    WORKSTREAM / "certify_ebch128_randomstepconv_parity_pivot_outward.py",
    WORKSTREAM / "ebch128_randomstepconv_g1_s20_pow2_parity_pivot_d109_diagnostic.json",
    WORKSTREAM / "ebch128_randomstepconv_s20_pow2_dense_complement_r64_d109_u074.json",
    RECEIPT,
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    receipt = json.loads(RECEIPT.read_text(encoding="utf-8"))
    payload = {
        "schema": "ebch128-pow2-randomstepconv-d109-manifest-v1",
        "status": "OUTWARD_DISTANCE_CERTIFICATE",
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
            (
                "python workstreams/finite_asymptotic_theory/"
                "evaluate_ebch128_randomstepconv_parity_pivot.py "
                "--outer-rows 16384 --memory-bits 20 "
                "--distance-numerator 109 --distance-denominator 1000 "
                "--grid-min -10 --grid-max 1 --grid-step 0.25 "
                "--output workstreams/finite_asymptotic_theory/"
                "ebch128_randomstepconv_g1_s20_pow2_parity_pivot_d109_diagnostic.json"
            ),
            (
                "python workstreams/finite_asymptotic_theory/"
                "evaluate_ebch128_randomstepconv_dense_complement.py "
                "--memory-bits 20 --distance-numerator 109 "
                "--distance-denominator 1000 --maximum-inactive-rows 64 "
                "--maximum-region-zeros 4096 --log-surprisal 0.74 "
                "--log-r-min 3 --log-r-max 7 --log-r-step 0.0625 "
                "--output workstreams/finite_asymptotic_theory/"
                "ebch128_randomstepconv_s20_pow2_dense_complement_r64_d109_u074.json"
            ),
            (
                "python workstreams/finite_asymptotic_theory/"
                "certify_ebch128_randomstepconv_parity_pivot_outward.py "
                "--diagnostic workstreams/finite_asymptotic_theory/"
                "ebch128_randomstepconv_g1_s20_pow2_parity_pivot_d109_diagnostic.json "
                "--complement-diagnostic workstreams/finite_asymptotic_theory/"
                "ebch128_randomstepconv_s20_pow2_dense_complement_r64_d109_u074.json "
                "--outer-rows 16384 --memory-bits 20 "
                "--distance-numerator 109 --distance-denominator 1000 "
                "--certified-margin-bits 5 --output "
                "workstreams/finite_asymptotic_theory/"
                "ebch128_randomstepconv_g1_s20_pow2_outward_d109.json"
            ),
        ],
        "runtime": {
            "python_flint_version": flint.__version__,
            "arithmetic": "Arb plus directed positive binary64 recurrences",
        },
        "imported_premise": (
            "scripts/EBCH128_64.wd is the exact spectrum of the concrete "
            "generator fixed in FINITE_K20_EBCH128_OUTER_DEFINITION.md"
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"output={OUTPUT}")


if __name__ == "__main__":
    main()
