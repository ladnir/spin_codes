#!/usr/bin/env python3
"""Verify global-lane bindings before an explicit g=4 full rerun.

The default action is a lightweight, read-only preflight.  ``--invoke`` is an
explicit future-run switch; it executes the manifest's existing full-certifier
entrypoint only after every construction, input, and proof-code digest passes.
This script does not modify the shared certifier.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

import certify_global_lane_puncture_layout as layout


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "G4_GLOBAL_LANE_RERUN_MANIFEST.json"
EXPECTED_SCHEMA = "permute-conv.riffle-g4-global-lane-rerun.v1"
OLD_FINAL_REPORT = "out/g4_end_to_end_bsp_outward_certificate.json"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def resolve(path_text: str) -> Path:
    path = Path(path_text)
    return path if path.is_absolute() else ROOT / path


def require_equal(label: str, actual: Any, expected: Any) -> None:
    if actual != expected:
        raise ValueError(f"{label}: expected {expected!r}, observed {actual!r}")


def validate_binding(label: str, binding: dict[str, Any]) -> Path:
    if not isinstance(binding, dict):
        raise ValueError(f"{label} binding is not an object")
    path_text = binding.get("path")
    expected = binding.get("sha256")
    if not isinstance(path_text, str) or not isinstance(expected, str):
        raise ValueError(f"{label} binding lacks path/SHA-256")
    path = resolve(path_text)
    if not path.is_file():
        raise ValueError(f"{label} binding is missing: {path}")
    require_equal(f"{label} SHA-256", sha256_file(path), expected)
    return path


def validate_layout_certificate() -> dict[str, Any]:
    layout.certify_perfect_matchings()
    layout.certify_one_row_regrouping()
    layout.certify_exact_marginals()
    layout.certify_graph_hole_distribution()
    sample = layout.sample_global_lane_punctures(__import__("random").Random(0))
    layout.validate_sample(sample)
    if any(hole.data_lane != sample.retained_lane for hole in sample.holes):
        raise AssertionError("global-lane sample used more than one puncture lane")
    if len({hole.band_zero_tile for hole in sample.holes}) != layout.HOLES:
        raise AssertionError("global-lane sample repeated a band-zero tile")
    return {
        "perfect_matchings": layout.LANES,
        "blocks_per_matching": layout.TILES,
        "per_block_selection_marginal": "1/128",
        "per_band0_coordinate_marginal": "1/5376",
        "lane_events_independent": False,
        "deterministic_sampler_smoke_lane": sample.retained_lane,
    }


def preflight(manifest_path: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    require_equal("manifest schema", manifest.get("schema"), EXPECTED_SCHEMA)
    require_equal(
        "manifest status",
        manifest.get("status"),
        "STATIC_PREFLIGHT_READY_FULL_OUTWARD_RERUN_NOT_EXECUTED",
    )

    supersedes = manifest.get("supersedes")
    old_manifest = validate_binding("superseded savepoint", supersedes)
    old = json.loads(old_manifest.read_text(encoding="utf-8"))
    require_equal("old savepoint schema", old.get("schema"), "permute-conv.riffle-g4-savepoint.v1")
    if supersedes.get("disposition") is None:
        raise ValueError("superseded savepoint lacks its non-transfer disposition")

    construction = manifest.get("construction_rule")
    if not isinstance(construction, dict):
        raise ValueError("manifest lacks construction_rule")
    require_equal("construction name", construction.get("name"), "global-lane-puncture-v1")
    require_equal("per-block marginal", construction.get("per_block_selection_marginal"), "1/128")
    require_equal("per-coordinate marginal", construction.get("per_band0_coordinate_marginal"), "1/5376")
    require_equal("lane-event independence", construction.get("lane_events_independent"), False)
    require_equal(
        "one-row regrouping",
        construction.get("one_conditioned_row_regroups_into_complete_codewords"),
        True,
    )

    bindings = {"superseded_savepoint": supersedes}
    for section_name in (
        "construction_bindings",
        "theorem_data_bindings",
        "candidate_proof_inputs",
        "proof_code_bindings",
    ):
        section = manifest.get(section_name)
        if not isinstance(section, dict) or not section:
            raise ValueError(f"manifest lacks {section_name}")
        for name, binding in section.items():
            validate_binding(f"{section_name}.{name}", binding)
            bindings[f"{section_name}.{name}"] = binding

    certificate = manifest["construction_bindings"]["exact_layout_certificate"]
    require_equal(
        "layout certificate import path",
        resolve(certificate["path"]).resolve(),
        Path(layout.__file__).resolve(),
    )
    require_equal(
        "layout certificate required status",
        certificate.get("required_status"),
        "EXACT_GLOBAL_LANE_PUNCTURE_LAYOUT_CERTIFICATE",
    )
    layout_report = validate_layout_certificate()

    theorem = manifest["theorem_data_bindings"]
    required_theorem_bindings = {
        "global_lane_total_weight_lemma",
        "ebch_weight_spectrum",
        "ebch_weight_spectrum_source",
        "graph_weight_spectrum",
        "systematic_split_slices",
        "band01_split_spectrum",
        "band12_split_spectrum",
        "selected_global_lane_total_weight_witness",
    }
    require_equal("theorem-data binding names", set(theorem), required_theorem_bindings)
    selected_binding = theorem["selected_global_lane_total_weight_witness"]
    require_equal(
        "selected global-lane witness reference",
        selected_binding.get("reference"),
        "g4_support_lower30_failed8_retuned_local.json:0",
    )
    selected_source = json.loads(
        resolve(selected_binding["path"]).read_text(encoding="utf-8")
    )
    selected_rows = (
        selected_source
        if isinstance(selected_source, list)
        else selected_source.get("rows", [selected_source])
    )
    if not selected_rows:
        raise ValueError("selected global-lane witness source is empty")
    selected_row = selected_rows[0]
    require_equal("selected witness name", selected_row.get("name"), selected_binding.get("name"))
    require_equal(
        "selected witness outer type",
        selected_row.get("outer_type"),
        "exact_graph_puncture_total_weight",
    )

    rerun = manifest.get("rerun")
    if not isinstance(rerun, dict):
        raise ValueError("manifest lacks rerun command")
    entrypoint = resolve(str(rerun.get("entrypoint", "")))
    require_equal(
        "rerun entrypoint",
        entrypoint.resolve(),
        resolve(manifest["proof_code_bindings"]["full_certifier"]["path"]).resolve(),
    )
    arguments = rerun.get("arguments")
    if not isinstance(arguments, list) or not all(isinstance(value, str) for value in arguments):
        raise ValueError("rerun arguments must be a string array")
    output = str(rerun.get("output", ""))
    if output == OLD_FINAL_REPORT or "--output" not in arguments:
        raise ValueError("global-lane rerun would overwrite or omit its new output")
    output_index = arguments.index("--output") + 1
    if output_index >= len(arguments) or arguments[output_index] != output:
        raise ValueError("rerun output field and command arguments disagree")
    old_final = old.get("artifacts", {}).get("final_report", {}).get("path")
    if output == old_final:
        raise ValueError("global-lane output aliases the old final report")

    # Every candidate proof input must appear verbatim in the command.
    for binding in manifest["candidate_proof_inputs"].values():
        if binding["path"] == OLD_FINAL_REPORT:
            raise ValueError("old final report was included as a candidate proof input")
        if binding["path"] not in arguments:
            raise ValueError(f"rerun command omits candidate input {binding['path']}")

    report = {
        "schema": "permute-conv.riffle-g4-global-lane-preflight.v1",
        "status": "PASS_STATIC_GLOBAL_LANE_RERUN_PREFLIGHT",
        "manifest_path": str(manifest_path.resolve()),
        "manifest_sha256": sha256_file(manifest_path),
        "validated_bindings": len(bindings),
        "layout_certificate": layout_report,
        "theorem_data_bindings": sorted(theorem),
        "command": [sys.executable, str(entrypoint), *arguments],
        "output": str(resolve(output)),
        "scope_limit": "No full g=4 certifier was invoked by preflight.",
    }
    return manifest, report


def invoke(manifest_path: Path, manifest: dict[str, Any], report: dict[str, Any]) -> None:
    environment = os.environ.copy()
    environment["PERMUTE_CONV_CONSTRUCTION_MANIFEST"] = str(manifest_path.resolve())
    environment["PERMUTE_CONV_CONSTRUCTION_MANIFEST_SHA256"] = report["manifest_sha256"]
    command = report["command"]
    subprocess.run(command, cwd=ROOT, env=environment, check=True)
    output = Path(report["output"])
    if not output.is_file():
        raise RuntimeError("full certifier returned without producing the declared output")
    result = json.loads(output.read_text(encoding="utf-8"))
    if result.get("group_bits") != 4:
        raise RuntimeError("full certifier output has the wrong packet width")
    print(f"full_output={output}")
    print(f"full_output_sha256={sha256_file(output)}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument(
        "--invoke",
        action="store_true",
        help="after preflight, launch the full certifier (never used by static tests)",
    )
    args = parser.parse_args()
    manifest, report = preflight(args.manifest)
    print("g4 global-lane rerun preflight: PASS")
    print(f"manifest_sha256={report['manifest_sha256']}")
    print(f"validated_bindings={report['validated_bindings']}")
    print("lane_events_independent=NO_BY_CONSTRUCTION")
    print("command=" + subprocess.list2cmdline(report["command"]))
    if args.invoke:
        invoke(args.manifest, manifest, report)
    else:
        print("full_certifier_invoked=NO_STATIC_PREFLIGHT_ONLY")


if __name__ == "__main__":
    main()
