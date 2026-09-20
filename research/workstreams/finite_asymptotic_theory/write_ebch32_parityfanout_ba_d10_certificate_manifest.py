#!/usr/bin/env python3
"""Write the SHA-256 manifest for the finite 10% certificate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
WORKSTREAM = Path(__file__).resolve().parent
OUTPUT = WORKSTREAM / "FINITE_K20_D10_CERTIFICATE_MANIFEST.json"
FILES = (
    WORKSTREAM / "FINITE_K20_D10_CERTIFICATE.md",
    WORKSTREAM / "FINITE_K20_EBCH32_PARITYFANOUT_SETUP.md",
    WORKSTREAM / "certify_ebch32_parityfanout_ba_setup.py",
    WORKSTREAM / "ebch32_parityfanout31x33_ba3_B256_setup_outward.json",
    WORKSTREAM / "ebch32_16_delta8_spectrum.csv",
    WORKSTREAM / "certify_ebch32_parityfanout_ba_rm2sub_q1.py",
    WORKSTREAM / "ebch32_parityfanout31x33_ba3_rm2sub_B256_q1_outward_d11.json",
    WORKSTREAM / "certify_ebch32_parityfanout_ba_rm2sub_q2_64.py",
    WORKSTREAM / "ebch32_parityfanout31x33_ba3_rm2sub_B256_q2_64_outward_d11.json",
    WORKSTREAM / "diagnose_ebch32_parityfanout_ba_one_band_interval_cover_d10.py",
    WORKSTREAM / "ebch32_parityfanout31x33_ba3_B256_one_band_interval_cover_d10.json",
    WORKSTREAM / "certify_ebch32_parityfanout_ba_one_band_interval_cover_d10_outward.py",
    WORKSTREAM / "ebch32_parityfanout31x33_ba3_B256_one_band_interval_cover_d10_outward.json",
    WORKSTREAM / "certify_ebch32_parityfanout_dense_cells_outward.py",
    WORKSTREAM / "certify_ebch32_parityfanout_ba_combined_d10.py",
    WORKSTREAM / "ebch32_parityfanout31x33_ba3_B256_combined_outward_d10.json",
    WORKSTREAM / "certify_ebch32_parityfanout_ba_bounded_setup_d10.py",
    WORKSTREAM
    / "ebch32_parityfanout31x33_ba3_B256_bounded_setup_combined_outward_d10.json",
    ROOT
    / "constructions"
    / "riffle_bchperm_transpose_bitshuffle_splitstate_preaddmul_rm2sub_t128_s20"
    / "receipts"
    / "min_state"
    / "s19_rm2sub_selection.json",
    Path(__file__).resolve(),
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    missing = [str(path) for path in FILES if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing manifest inputs: " + ", ".join(missing))
    payload = {
        "schema": "finite-k20-ebch32-parityfanout-ba-d10-certificate-manifest-v2",
        "status": "COMPLETE_BOUNDED_SETUP_AND_DISTANCE_CERTIFICATE",
        "claim": {
            "message_bits": 2**20,
            "output_bits": 2**21,
            "bad_weight_at_most": 209_715,
            "conditional_bad_distance_probability_less_than_2^-51": True,
            "setup_abort_or_bad_distance_probability_less_than_2^-45": True,
            "setup_abort_or_bad_distance_probability_less_than_2^-40": True,
        },
        "files": [
            {
                "path": path.resolve().relative_to(ROOT.resolve()).as_posix(),
                "sha256": sha256(path),
            }
            for path in FILES
        ],
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"files={len(FILES)}")
    print(f"output={OUTPUT}")


if __name__ == "__main__":
    main()
