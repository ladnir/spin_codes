#!/usr/bin/env python3
"""Build the artifact manifest for the conditional Block Expand-5 result."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


WORKSTREAM = Path(__file__).resolve().parent
OUTPUT = WORKSTREAM / "FINITE_K20_BLOCK_EXPAND5_F2_CONDITIONAL_MANIFEST.json"
FILES = (
    "BLOCK_EXPAND_CONSTITUENT_ROUTE.md",
    "evaluate_block_expand_accumulate_spectrum.py",
    "block_expand_accumulate_512_256_d14_t0_5_spectrum.json",
    "evaluate_block_expand_cap_budget.py",
    "block_expand_cap_budget_diagnostic.json",
    "verify_block_expand_pair_region_small.py",
    "block_expand_pair_region_small_exact.json",
    "analyze_block_expand_pair_moments_small.py",
    "block_expand_pair_moments_K4_B8_exact.json",
    "block_expand_pair_moments_K6_B12_exact.json",
    "block_expand_pair_moments_K8_B16_exact.json",
    "certify_block_expand5_F2_conditional_transfer.py",
    "block_expand5_F2_conditional_transfer_outward.json",
    "certify_block_expand5_mean_cap_outward.py",
    "block_expand5_mean_cap_outward.json",
    "single_random_constituent_B512_sparse_outward_s22.json",
    "single_random_constituent_B512_shared_two_band_dense_cover_outward_s22.json",
)


def main() -> None:
    artifacts = []
    for name in FILES:
        path = WORKSTREAM / name
        if not path.is_file():
            raise FileNotFoundError(path)
        artifacts.append(
            {
                "file": name,
                "bytes": path.stat().st_size,
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
        )
    payload = {
        "schema": "finite-k20-block-expand5-F2-conditional-manifest-v1",
        "status": "CONDITIONAL_THEOREM_ARTIFACT_MANIFEST",
        "claim": {
            "message_bits": 1 << 20,
            "output_bits": 1 << 21,
            "minimum_distance_at_least": 228590,
            "relative_distance_at_least": 228590 / (1 << 21),
            "failure_margin_bits_lower_conditional": 42.6254,
        },
        "unproved_hypothesis": (
            "For the sampled degree-14 Block Expand-5 [512,256] constituent, "
            "Var(A_w) <= 2 E[A_w] for every positive shell w."
        ),
        "artifacts": artifacts,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(f"artifacts={len(artifacts)}")
    print(f"output={OUTPUT}")


if __name__ == "__main__":
    main()
