#!/usr/bin/env python3
"""Build the hash manifest for the partial finite k=2^20 certificate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


WORKSTREAM = Path(__file__).resolve().parent
REPO_ROOT = WORKSTREAM.parents[1]
OUTPUT = WORKSTREAM / "FINITE_K20_PARTIAL_CERTIFICATE_MANIFEST.json"

FILES = [
    WORKSTREAM / "FINITE_K20_INDEPENDENT_SETUP.md",
    WORKSTREAM / "FINITE_K20_COMPARISON.md",
    WORKSTREAM / "FINITE_K20_PROOF_OBLIGATIONS.md",
    WORKSTREAM / "FINITE_K20_DENSE_TRANSFER_TARGET.md",
    WORKSTREAM / "FINITE_K20_CONDITIONING_WINDOW.md",
    WORKSTREAM / "FINITE_K20_TWO_TRACK_DENSE_PLAN.md",
    WORKSTREAM / "certify_golay_ba_rm2sub_finite_one_active.py",
    WORKSTREAM / "golay_ba3_rm2sub_finite_B240_q1_outward_w23_217_k20_d11.json",
    WORKSTREAM / "certify_golay_ba_rm2sub_finite_q2_64.py",
    WORKSTREAM / "golay_ba3_rm2sub_finite_B240_q2_64_outward_w23_217_k20_d11.json",
    WORKSTREAM / "diagnose_finite_k20_conditioning_window.py",
    WORKSTREAM / "golay_ba3_rm2sub_finite_B240_conditioning_window_diagnostic.json",
    WORKSTREAM / "diagnose_finite_k20_band_pure.py",
    WORKSTREAM / "golay_ba3_rm2sub_finite_B240_band_pure_w23_217_d109_qL_diagnostic.json",
    WORKSTREAM / "diagnose_finite_k20_band_compositions_qL.py",
    WORKSTREAM / "golay_ba3_rm2sub_finite_B240_band_compositions_qL_d109_diagnostic.json",
    WORKSTREAM / "diagnose_finite_k20_qL_column_holder.py",
    WORKSTREAM / "golay_ba3_rm2sub_finite_B240_qL_column_holder_d11.json",
    WORKSTREAM / "diagnose_finite_k20_band_simplex_cover_qL.py",
    WORKSTREAM / "golay_ba3_rm2sub_finite_B240_k20_d11.json",
    WORKSTREAM / "golay_ba3_rm2sub_finite_B240_bulk_q2_64_k20_d11.json",
    REPO_ROOT
    / (
        "constructions/riffle_bchperm_transpose_bitshuffle_splitstate_"
        "preaddmul_rm2sub_t128_s20/receipts/min_state/"
        "s19_rm2sub_b_kernel_spectrum.json"
    ),
    REPO_ROOT
    / (
        "constructions/riffle_bchperm_transpose_bitshuffle_splitstate_"
        "preaddmul_rm2sub_t128_s20/receipts/min_state/"
        "s19_rm2sub_a_spectrum_audit.json"
    ),
]


def relative(path: Path) -> str:
    return path.relative_to(REPO_ROOT).as_posix()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    payload = {
        "schema": "spin-finite-k20-partial-certificate-manifest-v1",
        "status": "PARTIAL_CERTIFICATE_OPEN_Q65_8832",
        "construction": {
            "outer": "independent conditioned Golay--BA-3 rows",
            "outer_bits": 240,
            "outer_rows": 8832,
            "message_bits_after_shortening": 1 << 20,
            "parent_output_bits": 2_119_680,
            "bad_weight_at_most": 233_164,
            "permitted_outer_weights": [23, 217],
        },
        "proved_components": [
            {
                "occupations": [1, 1],
                "expected_bad_word_count_upper_hex": "0x1.6f3c66666f370p-47",
                "integer_margin_bits": 46,
            },
            {
                "occupations": [2, 64],
                "expected_bad_word_count_upper_scaled": {
                    "mantissa_binary64_hex": "0x1.4dc4b8b530941p-1",
                    "binary_exponent": -84,
                },
                "integer_margin_bits": 84,
            },
        ],
        "not_proved": [
            "occupations 65 through 8832",
            "efficient exact G_240^23 conditioning",
            "RM2Sub implementation equivalence and complete source audit",
            "the final 2^-40 setup-failure theorem",
            "ordinary encoding and conditioned-setup performance",
        ],
        "external_performance_receipt": {
            "sha256": "4b612e86c267a914b6cb5f6b231579f2ceb8dcfd7aa5981ee0b5807ce4fa01ee",
            "scope": "transposed online comparison only",
        },
        "files": [
            {
                "path": relative(path),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in FILES
        ],
        "reproduce": [
            "$env:SPIN_LOWER_WEIGHT='23'; $env:SPIN_UPPER_WEIGHT='217'; $env:SPIN_Q1_OUTPUT='golay_ba3_rm2sub_finite_B240_q1_outward_w23_217_k20_d11.json'; python workstreams/finite_asymptotic_theory/certify_golay_ba_rm2sub_finite_one_active.py",
            "$env:SPIN_LOWER_WEIGHT='23'; $env:SPIN_UPPER_WEIGHT='217'; $env:SPIN_Q2_64_OUTPUT='golay_ba3_rm2sub_finite_B240_q2_64_outward_w23_217_k20_d11.json'; python workstreams/finite_asymptotic_theory/certify_golay_ba_rm2sub_finite_q2_64.py",
            "python workstreams/finite_asymptotic_theory/build_finite_k20_partial_manifest.py",
        ],
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"manifest={OUTPUT}")


if __name__ == "__main__":
    main()
