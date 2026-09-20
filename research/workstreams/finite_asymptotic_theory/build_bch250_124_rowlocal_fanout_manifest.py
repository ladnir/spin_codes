#!/usr/bin/env python3
"""Build the hash-bound certificate manifest for the BCH250-124 route."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path


WORKSTREAM = Path(__file__).resolve().parent
REPO_ROOT = WORKSTREAM.parents[1]
FILES = {
    "certificate": "FINITE_K20_BCH250_124_ROWLOCAL_FANOUT_CERTIFICATE.md",
    "status_note": "FINITE_K20_BCH250_124_ROWLOCAL_FANOUT_STATUS.md",
    "subcode_envelope": "shortened_bch250_124_subcode_constant_weight_envelope.json",
    "fanout_envelope": "bch250_124_parityfanout31x33_l256_packing_expected_envelope.json",
    "q1": "bch250_124_parityfanout31x33_l256_packing_expected_q1_d11_diagnostic.json",
    "q2_31": "bch250_124_parityfanout31x33_l256_packing_expected_q2_31_d11_diagnostic.json",
    "q32_8576": "bch250_124_parityfanout31x33_l256_finite_interval_cover_q32_8576_d11_diagnostic.json",
    "fanout_analyzer": "analyze_one_sampled_bch250_parityfanout.py",
    "dense_cover_builder": "diagnose_bch250_rowlocal_fanout_finite_interval_cover.py",
    "spectrum_outward": "bch250_124_parityfanout31x33_l256_spectrum_outward.json",
    "small_outward": "bch250_124_parityfanout31x33_l256_q1_31_outward_d11.json",
    "dense_outward": "bch250_124_parityfanout31x33_l256_dense_q32_8576_outward_d11.json",
    "spectrum_verifier": "certify_bch250_124_rowlocal_fanout_spectrum_arb.py",
    "small_verifier": "certify_bch250_124_rowlocal_fanout_q1_31_outward.py",
    "dense_verifier": "certify_bch250_124_rowlocal_fanout_dense_intervals_outward.py",
    "subcode_builder": "build_shortened_bch250_subcode_envelope.py",
    "parent_packing_envelope": "shortened_bch250_125_constant_weight_envelope.json",
    "parameter_audit": "shortened_bch256_parameter_audit.json",
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
    for path in paths.values():
        if not path.is_file():
            raise FileNotFoundError(path)
    external_paths = {
        name: REPO_ROOT / relative for name, relative in EXTERNAL_FILES.items()
    }
    for path in external_paths.values():
        if not path.is_file():
            raise FileNotFoundError(path)

    q1 = json.loads(paths["q1"].read_text(encoding="utf-8"))
    q2 = json.loads(paths["q2_31"].read_text(encoding="utf-8"))
    dense = json.loads(paths["q32_8576"].read_text(encoding="utf-8"))
    margins = [
        float(q1["aggregate_margin_bits"]),
        float(q2["partial_margin_bits"]),
        float(dense["union_margin_bits"]),
    ]
    maximum_log = max(-margin for margin in margins)
    combined_log = maximum_log + math.log2(
        sum(math.exp2(-margin - maximum_log) for margin in margins)
    )
    combined_margin = -combined_log
    if dense["rejected_singletons"]:
        raise ValueError("the dense interval cover has rejected singletons")
    if not dense["target_union_margin_met"]:
        raise ValueError("the dense interval cover misses its target")
    spectrum_outward = json.loads(
        paths["spectrum_outward"].read_text(encoding="utf-8")
    )
    small_outward = json.loads(
        paths["small_outward"].read_text(encoding="utf-8")
    )
    dense_outward = json.loads(
        paths["dense_outward"].read_text(encoding="utf-8")
    )
    if spectrum_outward["status"] != "OUTWARD_SPECTRUM_CERTIFICATE":
        raise ValueError("the spectrum receipt is not certified")
    if not small_outward["claim"]["comparison_to_2^-48"]:
        raise ValueError("the small-occupation receipt misses 48 bits")
    if not dense_outward["claim"]["comparison_to_2^-134"]:
        raise ValueError("the dense receipt misses 134 bits")

    payload = {
        "schema": "bch250-124-rowlocal-fanout-certificate-manifest-v1",
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
            "fanout_layers_per_row": 256,
        },
        "occupation_margins_bits": {
            "q1": margins[0],
            "q2_through_q31_union": margins[1],
            "q32_through_q8576_union": margins[2],
            "combined_union": combined_margin,
        },
        "comparison_to_2^-40_diagnostic": combined_margin > 40.0,
        "certified_failure_probability": {
            "small_occupations_strictly_below": "2^-48",
            "dense_occupations_strictly_below": "2^-134",
            "combined_strictly_below": "2^-47",
            "certified_integer_margin_bits": 47,
            "comparison_to_2^-40": True,
        },
        "files": [
            {
                "role": name,
                "path": FILES[name],
                "sha256": sha256(path),
            }
            for name, path in paths.items()
        ] + [
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
    output = WORKSTREAM / "FINITE_K20_BCH250_124_ROWLOCAL_FANOUT_MANIFEST.json"
    output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"combined_margin_bits={combined_margin:.12f}")
    print(f"output={output}")


if __name__ == "__main__":
    main()
