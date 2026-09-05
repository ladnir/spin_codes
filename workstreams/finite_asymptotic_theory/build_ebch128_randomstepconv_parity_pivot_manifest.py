#!/usr/bin/env python3
"""Build the hash manifest for the repeated-EBCH128 certificate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import flint


WORKSTREAM = Path(__file__).resolve().parent
REPO_ROOT = WORKSTREAM.parents[1]
OUTPUT = WORKSTREAM / "FINITE_K20_REPEATED_EBCH128_RANDOMSTEP_CONV_MANIFEST.json"
FILES = [
    REPO_ROOT / "scripts" / "ebch128x64_forward_paar_circuit.json",
    REPO_ROOT / "scripts" / "EBCH128_64.wd",
    WORKSTREAM / "evaluate_ebch128_randomstepconv_g1.py",
    WORKSTREAM / "evaluate_ebch128_randomstepconv_q1_exact.py",
    WORKSTREAM / "evaluate_ebch128_randomstepconv_parity_pivot.py",
    WORKSTREAM / "certify_ebch128_randomstepconv_parity_pivot_outward.py",
    WORKSTREAM / "ebch128_randomstepconv_g1_s30_parity_pivot_d11_diagnostic.json",
    WORKSTREAM / "ebch128_randomstepconv_g1_s30_parity_pivot_outward_d11.json",
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    receipt = json.loads(FILES[-1].read_text(encoding="utf-8"))
    payload = {
        "schema": "ebch128-randomstepconv-parity-pivot-manifest-v1",
        "status": "OUTWARD_DISTANCE_CERTIFICATE",
        "claim": receipt["claim"],
        "parameters": {
            **receipt["parameters"],
            "outer_cyclic_generator_polynomial_hex": "0xf4845518b9582a1f",
            "outer_extension_parity_coordinate": 127,
            "outer_generator_rows_sha256_little_endian_128": (
                "972cfc1c6de12e4ddc0c67680fd8e3cddedbd84d0a409c7061e09a86b818ca56"
            ),
        },
        "spectrum_generator_binding": {
            "status": "EXPLICIT_IMPORTED_PREMISE",
            "statement": (
                "The concrete generator row space has the weight enumerator "
                "in scripts/EBCH128_64.wd."
            ),
            "local_derivation_present": False,
        },
        "artifacts": [
            {
                "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
                "sha256": sha256(path),
                "size_bytes": path.stat().st_size,
            }
            for path in FILES
        ],
        "runtime": {
            "python_flint_version": flint.__version__,
            "arithmetic": "Arb plus directed positive binary64 recurrences",
        },
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"output={OUTPUT}")


if __name__ == "__main__":
    main()
