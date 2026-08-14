#!/usr/bin/env python3
"""Verify the frozen g=4 save-point artifacts and report invariants.

This is a fast integrity gate. It does not replace the complete outward replay
documented in PROOF_STATUS.md.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "G4_GLOBAL_LANE_CERTIFICATE_MANIFEST.json"
SUPPORTED_SCHEMAS = {
    "permute-conv.riffle-g4-savepoint.v1",
    "permute-conv.riffle-g4-global-lane-savepoint.v2",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require_equal(label: str, actual: Any, expected: Any) -> None:
    if actual != expected:
        raise SystemExit(f"FAIL {label}: expected {expected!r}, got {actual!r}")
    print(f"PASS {label}")


def resolve_artifact(path_text: str) -> Path:
    path = Path(path_text)
    return path if path.is_absolute() else ROOT / path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    args = parser.parse_args()

    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    schema = manifest.get("schema")
    if schema not in SUPPORTED_SCHEMAS:
        raise SystemExit(f"FAIL manifest schema: unsupported {schema!r}")
    print(f"PASS manifest schema {schema}")

    for name, metadata in manifest["artifacts"].items():
        path = resolve_artifact(metadata["path"])
        if not path.is_file():
            raise SystemExit(f"FAIL artifact {name}: missing {path}")
        require_equal(f"artifact {name} SHA-256", sha256_file(path), metadata["sha256"])

    report_path = resolve_artifact(manifest["artifacts"]["final_report"]["path"])
    report = json.loads(report_path.read_text(encoding="utf-8"))
    expected = manifest["final_report_assertions"]

    direct_fields = (
        "group_bits",
        "passed",
        "status",
        "certified_margin_bits",
        "bsp_replaced_root_cells",
        "bsp_leaves",
        "used_witnesses",
        "used_component_vertex_evaluations",
        "used_bsp_leaf_vertex_inequalities",
        "inequality_sha256",
        "witness_reports_sha256",
    )
    for field in direct_fields:
        require_equal(f"final report {field}", report[field], expected[field])

    require_equal(
        "final report union_log2_upper",
        report["union_log2_interval"][1],
        expected["union_log2_upper"],
    )
    require_equal(
        "final report feasible_strata",
        report["exact_geometry"]["feasible_strata"],
        expected["feasible_strata"],
    )
    require_equal(
        "final report covered_cells",
        report["exact_geometry"]["covered_cells"],
        expected["covered_cells"],
    )
    construction_expected = manifest.get("construction_assertions")
    if construction_expected is not None:
        construction = report.get("construction", {})
        require_equal(
            "construction rule",
            construction.get("rule"),
            construction_expected["rule"],
        )
        require_equal(
            "construction manifest SHA-256",
            construction.get("manifest_sha256"),
            construction_expected["manifest_sha256"],
        )
        require_equal(
            "construction lane-event independence",
            construction.get("lane_events_independent"),
            construction_expected["lane_events_independent"],
        )
        require_equal(
            "report source construction-manifest SHA-256",
            report["source_digests"]["construction_manifest"],
            construction_expected["manifest_sha256"],
        )
    print("G4 SAVEPOINT INTEGRITY: PASS")


if __name__ == "__main__":
    main()
