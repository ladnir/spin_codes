#!/usr/bin/env python3
"""Build the hash-bound manifest for the 56-layer BCH250 certificate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path


WORKSTREAM = Path(__file__).resolve().parent
REPO_ROOT = WORKSTREAM.parents[1]
FILES = {
    "certificate": "FINITE_K20_BCH250_124_ROWLOCAL_FANOUT56_CERTIFICATE.md",
    "subcode_envelope": "shortened_bch250_124_subcode_constant_weight_envelope.json",
    "parent_packing_envelope": "shortened_bch250_125_constant_weight_envelope.json",
    "parameter_audit": "shortened_bch256_parameter_audit.json",
    "spectrum_outward": "bch250_124_parityfanout31x33_l56_cutoff12_spectrum_outward.json",
    "q1_outward": "bch250_124_parityfanout31x33_l56_q1_outward_d11.json",
    "q2_31_outward": "bch250_124_parityfanout31x33_l56_q2_31_outward_d11.json",
    "dense_outward": "bch250_124_parityfanout31x33_l56_cutoff12_dense_q32_8576_outward_d11.json",
    "q1_witnesses": "bch250_124_parityfanout31x33_l56_shell_sensitive_q1_d11_diagnostic.json",
    "q2_31_witnesses": "bch250_124_parityfanout31x33_l256_packing_expected_q2_31_d11_diagnostic.json",
    "dense_witnesses": "bch250_124_parityfanout31x33_l56_cutoff12_finite_interval_cover_q32_8576_d11_diagnostic.json",
    "spectrum_verifier": "certify_bch250_124_rowlocal_fanout_spectrum_arb.py",
    "q1_verifier": "certify_bch250_124_shell_sensitive_q1_outward.py",
    "q2_31_verifier": "certify_bch250_124_shell_sensitive_q2_31_outward.py",
    "dense_verifier": "certify_bch250_124_rowlocal_fanout_dense_intervals_outward.py",
    "dense_witness_builder": "diagnose_bch250_rowlocal_fanout_finite_interval_cover.py",
    "finite_cover_evaluator": "evaluate_bch250_rowlocal_fanout_finite_cover.py",
    "subcode_builder": "build_shortened_bch250_subcode_envelope.py",
    "manifest_builder": "build_bch250_124_rowlocal_fanout56_manifest.py",
}
EXTERNAL_FILES = {
    "rm2sub_selection": (
        "constructions/riffle_bchperm_transpose_bitshuffle_splitstate_"
        "preaddmul_rm2sub_t128_s20/receipts/min_state/"
        "s19_rm2sub_selection.json"
    ),
    "rm2sub_activation": (
        "constructions/riffle_bchperm_transpose_bitshuffle_splitstate_"
        "preaddmul_rm2sub_t128_s20/receipts/min_state/"
        "s19_rm2sub_b_kernel_spectrum.json"
    ),
    "rm2sub_live_spectrum": (
        "constructions/riffle_bchperm_transpose_bitshuffle_splitstate_"
        "preaddmul_rm2sub_t128_s20/receipts/min_state/"
        "s19_rm2sub_a_spectrum_audit.json"
    ),
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    paths = {name: WORKSTREAM / relative for name, relative in FILES.items()}
    external_paths = {
        name: REPO_ROOT / relative for name, relative in EXTERNAL_FILES.items()
    }
    for path in [*paths.values(), *external_paths.values()]:
        if not path.is_file():
            raise FileNotFoundError(path)

    spectrum = json.loads(paths["spectrum_outward"].read_text(encoding="utf-8"))
    q1 = json.loads(paths["q1_outward"].read_text(encoding="utf-8"))
    q2 = json.loads(paths["q2_31_outward"].read_text(encoding="utf-8"))
    dense = json.loads(paths["dense_outward"].read_text(encoding="utf-8"))
    if spectrum["status"] != "OUTWARD_SPECTRUM_CERTIFICATE":
        raise ValueError("spectrum receipt is not certified")
    if not spectrum["claim"]["tail_comparison_accepted"]:
        raise ValueError("tail-event receipt failed")
    if not q1["claim"]["comparison_to_2^-53"]:
        raise ValueError("Q=1 receipt misses 53 bits")
    if not q2["claim"]["comparison_to_2^-52"]:
        raise ValueError("Q=2..31 receipt misses 52 bits")
    if not dense["claim"]["comparison_to_2^-41"]:
        raise ValueError("dense receipt misses 41 bits")

    certified_sum = (
        Fraction(1, 1 << 45)
        + Fraction(1, 1 << 53)
        + Fraction(1, 1 << 52)
        + Fraction(1, 1 << 41)
    )
    if not certified_sum < Fraction(1, 1 << 40):
        raise ArithmeticError("coarse certified union does not reach 40 bits")

    displayed_margins = {
        "tail_event": -float(
            spectrum["claim"]["tail_expected_count_log2_upper_display"]
        ),
        "q1": float(q1["claim"]["union_margin_bits_display"]),
        "q2_through_q31_union": float(
            q2["claim"]["union_margin_bits_display"]
        ),
        "q32_through_q8576_union_on_central_event": float(
            dense["claim"]["union_margin_bits"]
        ),
    }
    displayed_sum = sum(math.exp2(-value) for value in displayed_margins.values())
    displayed_combined = -math.log2(displayed_sum)

    payload = {
        "schema": "bch250-124-rowlocal-fanout56-certificate-manifest-v1",
        "status": "HASH_BOUND_OUTWARD_PROOF_MODEL_CERTIFICATE",
        "parameters": {
            "requested_message_bits": 1 << 20,
            "parent_message_bits": 1_063_424,
            "shortened_input_coordinates": 14_848,
            "outer_bits": 250,
            "outer_dimension": 124,
            "outer_rows": 8_576,
            "output_bits": 2_144_000,
            "distance": 235_840,
            "relative_distance": 0.11,
            "fanout_layers_per_row": 56,
            "tail_weight_interval_low": [1, 12],
            "tail_weight_interval_high": [238, 250],
            "central_density_excess_bits": "389/1250",
        },
        "displayed_margins_bits": {
            **displayed_margins,
            "combined_union": displayed_combined,
        },
        "certified_failure_probability": {
            "tail_event_strictly_below": "2^-45",
            "q1_strictly_below": "2^-53",
            "q2_through_q31_strictly_below": "2^-52",
            "dense_on_central_event_strictly_below": "2^-41",
            "combined_strictly_below": "2^-40",
            "certified_integer_margin_bits": 40,
        },
        "fanout_operation_proxy": {
            "scalar_xors_per_row": 3_528,
            "scalar_xors_all_rows": 30_256_128,
            "reduction_from_256_layers_fraction": "25/32",
            "speedup_factor_against_256_layer_proxy": 256 / 56,
        },
        "files": [
            {"role": name, "path": FILES[name], "sha256": sha256(path)}
            for name, path in paths.items()
        ]
        + [
            {
                "role": name,
                "path": EXTERNAL_FILES[name],
                "sha256": sha256(path),
            }
            for name, path in external_paths.items()
        ],
        "blocking_obligations": [
            "implementation-equivalence audit",
            "optimized implementation and performance measurement",
        ],
    }
    output = WORKSTREAM / "FINITE_K20_BCH250_124_ROWLOCAL_FANOUT56_MANIFEST.json"
    output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"displayed_combined_margin_bits={displayed_combined:.12f}")
    print(f"output={output}")


if __name__ == "__main__":
    main()
