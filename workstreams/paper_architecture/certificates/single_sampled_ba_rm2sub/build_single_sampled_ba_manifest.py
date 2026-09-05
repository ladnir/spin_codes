#!/usr/bin/env python3
"""Build the hash manifest for the one-sampled BA/RM2Sub certificate."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


LOCAL_FILES = [
    "workstreams/finite_asymptotic_theory/build_single_sampled_ba_manifest.py",
    "workstreams/finite_asymptotic_theory/SINGLE_SAMPLED_BA_RM2SUB_D11.md",
    "workstreams/finite_asymptotic_theory/STRUCTURED_SPIN_THEOREM_TARGET.md",
    "workstreams/finite_asymptotic_theory/ASYMPTOTIC_SCALING.md",
    "workstreams/finite_asymptotic_theory/ARBITRARY_LENGTH_WRAPPER.md",
    "workstreams/finite_asymptotic_theory/certify_golay_ba_concave_majorant.py",
    "workstreams/finite_asymptotic_theory/golay_ba3_concave_majorant.json",
    "workstreams/finite_asymptotic_theory/certify_golay_ba_rm2sub_joint_interval.py",
    "workstreams/finite_asymptotic_theory/golay_ba3_rm2sub_joint_interval_d11.json",
    "workstreams/finite_asymptotic_theory/certify_golay_ba_rm2sub_sparse.py",
    "workstreams/finite_asymptotic_theory/golay_ba3_rm2sub_sparse_d11.json",
    "workstreams/finite_asymptotic_theory/WEIGHT_COUPLED_FIXED_OCCUPATION.md",
    "workstreams/finite_asymptotic_theory/certify_golay_ba_rm2sub_weight_coupled_fixed.py",
    "workstreams/finite_asymptotic_theory/golay_ba3_rm2sub_weight_coupled_fixed_d11.json",
    "workstreams/finite_asymptotic_theory/analyze_golay_ba_rm2sub_joint.py",
    "workstreams/finite_asymptotic_theory/analyze_rm2sub_dense_occupation.py",
    "workstreams/finite_asymptotic_theory/certify_rm2sub_dense_small.py",
    "workstreams/finite_asymptotic_theory/rm2sub_dense_small_exact_d11.json",
    "workstreams/finite_asymptotic_theory/certify_rm2sub_uniform_fixed_occupation.py",
    "workstreams/finite_asymptotic_theory/rm2sub_uniform_fixed_occupation_d11.json",
    "workstreams/finite_asymptotic_theory/RM2SUB_ONE_ACTIVE_CONTINUUM.md",
    "workstreams/finite_asymptotic_theory/RM2SUB_TWO_ACTIVE_CONTINUUM.md",
    "workstreams/finite_asymptotic_theory/RM2SUB_DENSE_OCCUPATION.md",
]

COMPANION_FILES = [
    "workstreams/linear_time_audit/ASYMPTOTIC_DISTANCE_CERTIFICATE.md",
    "workstreams/linear_time_audit/LINEAR_OUTER_CERTIFICATE.md",
    "workstreams/linear_time_audit/verify_outer_interval.py",
    "workstreams/linear_time_audit/verify_outer_sparse_constants.py",
    "workstreams/linear_time_audit/OUTER_INTERVAL_RECEIPT.json",
]

FROZEN_FILES = [
    "constructions/riffle_parityfanout31x33_bchperm_transpose_bitshuffle_"
    "splitstate_preaddmul_rm2sub_t128_s19/MAIN_CODE_FREEZE.md",
    "constructions/riffle_parityfanout31x33_bchperm_transpose_bitshuffle_"
    "splitstate_preaddmul_rm2sub_t128_s19/frozen_source/SOURCE_MANIFEST.json",
    "constructions/riffle_bchperm_transpose_bitshuffle_splitstate_preaddmul_"
    "rm2sub_t128_s20/receipts/min_state/s19_rm2sub_selection.json",
]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def records(root: Path, names: list[str]) -> list[dict[str, str]]:
    result = []
    for name in names:
        path = root / name
        if not path.is_file():
            raise FileNotFoundError(path)
        result.append({"path": name, "sha256": digest(path)})
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument(
        "--companion-root",
        type=Path,
        default=Path("C:/Users/peter/.codex/worktrees/3ef9/permute_conv"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "workstreams/finite_asymptotic_theory/"
            "SINGLE_SAMPLED_BA_RM2SUB_CERTIFICATE_MANIFEST.json"
        ),
    )
    args = parser.parse_args()

    payload = {
        "schema": "single-sampled-ba-rm2sub-certificate-manifest-v1",
        "claim": (
            "rate 1/2, relative distance 0.11, "
            "B=(39/4) log2(N)+O(1), "
            "and O(N) ordinary and transposed work"
        ),
        "local_files": records(args.repo_root, LOCAL_FILES),
        "companion_linear_time_files": records(
            args.companion_root, COMPANION_FILES
        ),
        "frozen_dependencies": records(args.repo_root, FROZEN_FILES),
        "companion_note": (
            "The companion paths are repository-relative.  This worktree "
            "reads them from --companion-root until integration merges both "
            "owned workstreams."
        ),
    }
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"status": "written", "output": str(args.output), "files": sum(len(payload[key]) for key in ("local_files", "companion_linear_time_files", "frozen_dependencies"))}, sort_keys=True))


if __name__ == "__main__":
    main()
