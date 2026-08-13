#!/usr/bin/env python3
"""Tune a resumable supplementary g=8 witness atlas at extracted targets.

The tool consumes the canonical tuning-target catalogue, tunes one paired
inner and conditioned-row witness at each selected exact profile, and writes
one content-bound checkpoint per target.  The final supplementary atlas is a
binary64 discovery artifact.  It is not an outward certificate or a coverage
claim.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from types import SimpleNamespace
from typing import Any

for _variable in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "BLIS_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_variable, "1")

from packet_group_native import load_point_caps, load_shared_drive_apply
from packet_group_profile_bound import split_cap_table
from probe_packet_group_conditioned_row_outer import load_split_spectrum
from survey_packet_group_g8_full_support_witness_reuse import optimize_seed_witness


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "G8_SUPPORT_MANIFEST.json"
DEFAULT_BASE_ATLAS = ROOT / "out" / "g8_support_seed_atlas.json"
DEFAULT_SPECTRUM01 = ROOT / "out" / "ebch85_band01_split_spectrum.csv"
DEFAULT_PUNCTURED01 = ROOT / "out" / "ebch84_punctured_band01_split_spectrum.csv"
DEFAULT_SPECTRUM12 = ROOT / "out" / "ebch86_band12_split_spectrum.csv"
EXPECTED_MANIFEST_SHA256 = (
    "cc446d0c6a0c986a8712ef78c4aa7c9d6073687665b7e0f315164339a5b81616"
)
EXPECTED_BASE_ATLAS_SHA256 = (
    "c48644fa62ad74619a4db30b07950b14aaf1419d040760f7d659f43288fca45e"
)
TARGET_SCHEMA = "permute-conv.packet-group-g8-tuning-target-catalogue.v1"
SUPPLEMENTARY_SCHEMA = "permute-conv.packet-group-g8-supplementary-witness-atlas.v1"
CHECKPOINT_SCHEMA = "permute-conv.packet-group-g8-supplementary-witness-checkpoint.v1"
MANIFEST_SCHEMA = "packet-group-g8-support-manifest-v1"
BASE_ATLAS_SCHEMA = "permute-conv.packet-group-g8-support-seed-atlas.v1"
CLASSES = 9
ATOMS = 262144

_WORKER_SPLIT_CAPS = None
_WORKER_SPECTRA = None
_WORKER_ARGS = None


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_path(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def write_canonical_atomic(path: Path, value: Any) -> str:
    payload = canonical_bytes(value) + b"\n"
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp-{os.getpid()}")
    temporary.write_bytes(payload)
    os.replace(temporary, path)
    return sha256_bytes(payload)


def parse_float_list(text: str) -> tuple[float, ...]:
    try:
        values = tuple(float(item) for item in text.split(",") if item)
    except ValueError as error:
        raise argparse.ArgumentTypeError("expected comma-separated floats") from error
    if not values or any(not math.isfinite(value) for value in values):
        raise argparse.ArgumentTypeError("float list must be finite and nonempty")
    return values


def validate_profile(value: Any) -> list[int]:
    if not isinstance(value, list) or len(value) != CLASSES:
        raise ValueError("target profile must contain nine integers")
    if any(isinstance(item, bool) for item in value):
        raise ValueError("target profile contains a Boolean")
    profile = [int(item) for item in value]
    if any(count < 0 for count in profile) or sum(profile) != ATOMS:
        raise ValueError("target profile is outside the g=8 profile domain")
    return profile


def validate_inputs(manifest_path: Path, base_atlas_path: Path, targets_path: Path):
    manifest_digest = sha256_path(manifest_path)
    base_digest = sha256_path(base_atlas_path)
    targets_digest = sha256_path(targets_path)
    if manifest_digest != EXPECTED_MANIFEST_SHA256:
        raise ValueError(
            f"manifest digest changed: {manifest_digest}; expected {EXPECTED_MANIFEST_SHA256}"
        )
    if base_digest != EXPECTED_BASE_ATLAS_SHA256:
        raise ValueError(
            f"base atlas digest changed: {base_digest}; expected {EXPECTED_BASE_ATLAS_SHA256}"
        )
    manifest = json.loads(manifest_path.read_bytes())
    base = json.loads(base_atlas_path.read_bytes())
    targets = json.loads(targets_path.read_bytes())
    if manifest.get("schema") != MANIFEST_SCHEMA:
        raise ValueError("support manifest schema changed")
    if base.get("schema") != BASE_ATLAS_SCHEMA:
        raise ValueError("base witness atlas schema changed")
    if targets.get("schema") != TARGET_SCHEMA:
        raise ValueError("input is not an extracted tuning-target catalogue")
    if targets.get("status") != "DIAGNOSTIC_TUNING_TARGETS_ONLY":
        raise ValueError("tuning-target catalogue has an unexpected status")
    if targets.get("manifest_sha256") != manifest_digest:
        raise ValueError("tuning targets bind another support manifest")
    if targets.get("atlas_sha256") != base_digest:
        raise ValueError("tuning targets bind another base atlas")
    sources = manifest.get("witness_sources", [])
    if len(sources) != 1 or sources[0].get("sha256") != base_digest:
        raise ValueError("support manifest does not bind the supplied base atlas")
    if targets_path.read_bytes() != canonical_bytes(targets) + b"\n":
        raise ValueError("tuning-target catalogue is not canonical sorted-compact LF JSON")
    return (
        manifest,
        base,
        targets,
        sources[0],
        manifest_digest,
        base_digest,
        targets_digest,
    )


def select_targets(catalogue: dict[str, Any], maximum: int) -> list[dict[str, Any]]:
    rows = catalogue.get("targets")
    if not isinstance(rows, list) or not rows:
        raise ValueError("tuning-target catalogue is empty")
    selected = rows[:maximum]
    seen_ids: set[str] = set()
    seen_profiles: set[tuple[int, ...]] = set()
    prior_rank = 0
    result = []
    for row in selected:
        rank = int(row.get("rank"))
        identifier = str(row.get("target_id", ""))
        profile = validate_profile(row.get("profile"))
        profile_digest = sha256_bytes(canonical_bytes(profile))
        if rank <= prior_rank or not identifier or identifier in seen_ids:
            raise ValueError("target ranks or identifiers are not unique and increasing")
        if row.get("profile_sha256") != profile_digest:
            raise ValueError(f"target profile digest changed at rank {rank}")
        key = tuple(profile)
        if key in seen_profiles:
            raise ValueError("target catalogue contains a duplicate exact profile")
        prior_rank = rank
        seen_ids.add(identifier)
        seen_profiles.add(key)
        result.append(row)
    return result


def configuration_record(args: argparse.Namespace, targets_digest: str) -> dict[str, Any]:
    sources = {
        "spectrum01": {"path": str(args.spectrum01), "sha256": sha256_path(args.spectrum01)},
        "punctured01": {
            "path": str(args.punctured01),
            "sha256": sha256_path(args.punctured01),
        },
        "spectrum12": {"path": str(args.spectrum12), "sha256": sha256_path(args.spectrum12)},
        "inner_evaluator": {
            "path": str(ROOT / "scripts" / "packet_group_profile_bound.py"),
            "sha256": sha256_path(ROOT / "scripts" / "packet_group_profile_bound.py"),
        },
        "conditioned_row_evaluator": {
            "path": str(ROOT / "scripts" / "probe_packet_group_conditioned_row_outer.py"),
            "sha256": sha256_path(
                ROOT / "scripts" / "probe_packet_group_conditioned_row_outer.py"
            ),
        },
    }
    record = {
        "targets_sha256": targets_digest,
        "inner_poles": list(args.inner_poles),
        "inner_exponents": list(args.inner_exponents),
        "inner_screen_iterations": args.inner_screen_iterations,
        "inner_final_iterations": args.inner_final_iterations,
        "outer_mode": args.outer_mode,
        "outer_max_iterations": args.outer_max_iterations,
        "outer_max_evaluations": args.outer_max_evaluations,
        "sources": sources,
    }
    record["configuration_digest"] = sha256_bytes(canonical_bytes(record))
    return record


def worker_configuration(record: dict[str, Any]) -> dict[str, Any]:
    return {
        "inner_poles": record["inner_poles"],
        "inner_exponents": record["inner_exponents"],
        "inner_screen_iterations": record["inner_screen_iterations"],
        "inner_final_iterations": record["inner_final_iterations"],
        "outer_mode": record["outer_mode"],
        "outer_max_iterations": record["outer_max_iterations"],
        "outer_max_evaluations": record["outer_max_evaluations"],
    }


def initialize_worker(
    spectrum01: str, punctured01: str, spectrum12: str, configuration: dict[str, Any]
) -> None:
    global _WORKER_SPLIT_CAPS, _WORKER_SPECTRA, _WORKER_ARGS
    native_apply = load_shared_drive_apply()
    native_caps = load_point_caps()
    if native_apply is None or native_caps is None:
        raise RuntimeError("supplementary tuner requires native inner apply and point-cap kernels")
    _WORKER_SPECTRA = (
        load_split_spectrum(Path(spectrum01)),
        load_split_spectrum(Path(punctured01), count_field="pair_count", divisor=42),
        load_split_spectrum(Path(spectrum12)),
    )
    _WORKER_SPLIT_CAPS = split_cap_table()
    _WORKER_ARGS = SimpleNamespace(**configuration)


def tune_worker(task: dict[str, Any]) -> dict[str, Any]:
    if _WORKER_SPLIT_CAPS is None or _WORKER_SPECTRA is None or _WORKER_ARGS is None:
        raise RuntimeError("supplementary tuning worker was not initialized")
    descriptor = {
        "name": f"supplementary_{task['target_id']}",
        "family": "highdim_unresolved_worst_vertex",
        "profile": task["profile"],
    }
    return optimize_seed_witness(
        descriptor, _WORKER_SPLIT_CAPS, _WORKER_SPECTRA, _WORKER_ARGS
    )


def target_binding(target: dict[str, Any], targets_digest: str) -> dict[str, Any]:
    return {
        "catalogue_sha256": targets_digest,
        "target_id": target["target_id"],
        "target_rank": int(target["rank"]),
        "target_row_sha256": sha256_bytes(canonical_bytes(target)),
        "profile_sha256": target["profile_sha256"],
        "provenance": target.get("provenance", []),
    }


def build_row(
    ordinal: int,
    target: dict[str, Any],
    targets_digest: str,
    tuned: dict[str, Any],
) -> dict[str, Any]:
    profile = validate_profile(target["profile"])
    support = [index for index, count in enumerate(profile) if count]
    if tuned.get("seed_profile") != profile:
        raise ValueError("tuned checkpoint profile does not match its target")
    affine = tuned.get("affine")
    if (
        not isinstance(affine, dict)
        or not math.isfinite(float(affine.get("constant_log2")))
        or len(affine.get("charge_log2", [])) != CLASSES
        or any(not math.isfinite(float(value)) for value in affine["charge_log2"])
    ):
        raise ValueError("tuned result has a malformed affine witness")
    return {
        "ordinal": ordinal,
        "name": f"supplementary_{target['target_id']}",
        "support": support,
        "support_mask": sum(1 << index for index in support),
        "dimension": len(support) - 1,
        "profile": profile,
        "physical_weight": sum(index * count for index, count in enumerate(profile)),
        "local_combined_log2": float(tuned["seed_combined_log2"]),
        "uniform_profile_target_margin_bits": float(
            tuned["seed_uniform_target_margin_bits"]
        ),
        "source_target": target_binding(target, targets_digest),
        "inner": tuned["inner"],
        "outer": tuned["outer"],
        "affine": affine,
        "elapsed_seconds": float(tuned["elapsed_seconds"]),
    }


def checkpoint_path(directory: Path, target: dict[str, Any]) -> Path:
    return directory / f"rank_{int(target['rank']):06d}_{target['profile_sha256'][:16]}.json"


def checkpoint_record(
    target: dict[str, Any],
    targets_digest: str,
    manifest_digest: str,
    base_source_id: str,
    base_digest: str,
    configuration_digest: str,
    tuned: dict[str, Any],
) -> dict[str, Any]:
    return {
        "schema": CHECKPOINT_SCHEMA,
        "status": "COMPLETE_DIAGNOSTIC_TUNING",
        "manifest_sha256": manifest_digest,
        "base_atlas_source_id": base_source_id,
        "base_atlas_sha256": base_digest,
        "configuration_digest": configuration_digest,
        "source_target": target_binding(target, targets_digest),
        "tuned": tuned,
    }


def load_checkpoint(
    path: Path,
    target: dict[str, Any],
    targets_digest: str,
    manifest_digest: str,
    base_source_id: str,
    base_digest: str,
    configuration_digest: str,
) -> dict[str, Any] | None:
    if not path.is_file():
        return None
    checkpoint = json.loads(path.read_bytes())
    expected = {
        "schema": CHECKPOINT_SCHEMA,
        "status": "COMPLETE_DIAGNOSTIC_TUNING",
        "manifest_sha256": manifest_digest,
        "base_atlas_source_id": base_source_id,
        "base_atlas_sha256": base_digest,
        "configuration_digest": configuration_digest,
    }
    if any(checkpoint.get(key) != value for key, value in expected.items()):
        raise ValueError(f"checkpoint binding mismatch: {path}")
    if checkpoint.get("source_target") != target_binding(target, targets_digest):
        raise ValueError(f"checkpoint target binding mismatch: {path}")
    if path.read_bytes() != canonical_bytes(checkpoint) + b"\n":
        raise ValueError(f"checkpoint is not canonical LF JSON: {path}")
    return checkpoint["tuned"]


def assemble_atlas(
    targets: list[dict[str, Any]],
    tuned_by_id: dict[str, dict[str, Any]],
    targets_path: Path,
    targets_digest: str,
    manifest_digest: str,
    base_source_id: str,
    base_digest: str,
    configuration: dict[str, Any],
    workers: int,
) -> dict[str, Any]:
    rows = [
        build_row(index, target, targets_digest, tuned_by_id[target["target_id"]])
        for index, target in enumerate(targets)
    ]
    row_bindings = [
        {
            "row": index,
            "target_id": row["source_target"]["target_id"],
            "row_sha256": sha256_bytes(canonical_bytes(row)),
        }
        for index, row in enumerate(rows)
    ]
    return {
        "schema": SUPPLEMENTARY_SCHEMA,
        "status": "DIAGNOSTIC_BINARY64_SUPPLEMENTARY_WITNESS_ATLAS",
        "scope": (
            "paired local tuning at a finite extracted target prefix; no outward replay "
            "or profile-domain coverage claim"
        ),
        "manifest_sha256": manifest_digest,
        "base_atlas_binding": {
            "source_id": base_source_id,
            "sha256": base_digest,
            "schema": BASE_ATLAS_SCHEMA,
        },
        "target_catalogue": {
            "path": str(targets_path),
            "sha256": targets_digest,
            "schema": TARGET_SCHEMA,
        },
        "configuration": configuration,
        "workers": workers,
        "target_count": len(targets),
        "row_bindings": row_bindings,
        "rows": rows,
    }


def synthetic_self_test() -> None:
    profile = [262136] + [1] * 8
    target = {
        "rank": 1,
        "target_id": "target-000001-synthetic",
        "profile": profile,
        "profile_sha256": sha256_bytes(canonical_bytes(profile)),
        "provenance": [{"node_id": "r1"}],
    }
    catalogue = {"targets": [target]}
    if select_targets(catalogue, 1) != [target]:
        raise AssertionError("synthetic target selection failed")
    tuned = {
        "seed_profile": profile,
        "seed_combined_log2": -1000.0,
        "seed_uniform_target_margin_bits": 800.0,
        "inner": {"pole": 0.2},
        "outer": {"pair_cauchy_theta": 0.5},
        "affine": {"constant_log2": 1.0, "charge_log2": [0.0] * CLASSES},
        "elapsed_seconds": 1.0,
    }
    checkpoint = checkpoint_record(
        target,
        "a" * 64,
        "synthetic-base",
        "b" * 64,
        "c" * 64,
        "d" * 64,
        tuned,
    )
    if checkpoint["source_target"]["target_row_sha256"] != sha256_bytes(
        canonical_bytes(target)
    ):
        raise AssertionError("synthetic target binding failed")
    atlas = assemble_atlas(
        [target],
        {target["target_id"]: tuned},
        Path("synthetic_targets.json"),
        "a" * 64,
        "b" * 64,
        "synthetic-base",
        "c" * 64,
        {"configuration_digest": "d" * 64},
        2,
    )
    if (
        atlas["target_count"] != 1
        or atlas["rows"][0]["profile"] != profile
        or atlas["row_bindings"][0]["row_sha256"]
        != sha256_bytes(canonical_bytes(atlas["rows"][0]))
    ):
        raise AssertionError("synthetic supplementary atlas assembly failed")
    print("synthetic_targets=1")
    print("synthetic_checkpoint_binding=PASS")
    print("synthetic_canonical_affine_row=PASS")
    print("status=SYNTHETIC_SUPPLEMENTARY_TUNER_SMOKE_PASS")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--targets", type=Path)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--base-atlas", type=Path, default=DEFAULT_BASE_ATLAS)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--checkpoint-dir", type=Path)
    parser.add_argument("--max-targets", type=int, default=8)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--inner-poles", type=parse_float_list, default=(0.2, 0.5, 0.8))
    parser.add_argument(
        "--inner-exponents", type=parse_float_list, default=(0.25, 0.5, 0.75)
    )
    parser.add_argument("--inner-screen-iterations", type=int, default=1)
    parser.add_argument("--inner-final-iterations", type=int, default=3)
    parser.add_argument(
        "--outer-mode", choices=("symmetric", "asymmetric"), default="asymmetric"
    )
    parser.add_argument("--outer-max-iterations", type=int, default=60)
    parser.add_argument("--outer-max-evaluations", type=int, default=800)
    parser.add_argument("--spectrum01", type=Path, default=DEFAULT_SPECTRUM01)
    parser.add_argument("--punctured01", type=Path, default=DEFAULT_PUNCTURED01)
    parser.add_argument("--spectrum12", type=Path, default=DEFAULT_SPECTRUM12)
    parser.add_argument("--synthetic-self-test", action="store_true")
    args = parser.parse_args()
    if args.synthetic_self_test:
        synthetic_self_test()
        return
    if args.targets is None or args.output is None or args.checkpoint_dir is None:
        parser.error("--targets, --output, and --checkpoint-dir are required")
    if (
        args.max_targets <= 0
        or args.workers <= 0
        or args.inner_screen_iterations <= 0
        or args.inner_final_iterations <= 0
        or args.outer_max_iterations <= 0
        or args.outer_max_evaluations <= 0
        or any(not 0.0 < value < 1.0 for value in args.inner_poles)
    ):
        parser.error("target, worker, iteration, and pole budgets must be positive")
    try:
        (
            _manifest,
            _base,
            catalogue,
            base_source,
            manifest_digest,
            base_digest,
            targets_digest,
        ) = validate_inputs(args.manifest, args.base_atlas, args.targets)
        targets = select_targets(catalogue, args.max_targets)
        configuration = configuration_record(args, targets_digest)
        configuration_digest = configuration["configuration_digest"]
        tuned_by_id: dict[str, dict[str, Any]] = {}
        missing = []
        for target in targets:
            path = checkpoint_path(args.checkpoint_dir, target)
            tuned = load_checkpoint(
                path,
                target,
                targets_digest,
                manifest_digest,
                base_source["source_id"],
                base_digest,
                configuration_digest,
            )
            if tuned is None:
                missing.append(target)
            else:
                tuned_by_id[target["target_id"]] = tuned
                print(f"checkpoint=resume target={target['target_id']}", flush=True)

        initializer_args = (
            str(args.spectrum01),
            str(args.punctured01),
            str(args.spectrum12),
            worker_configuration(configuration),
        )
        if missing and args.workers == 1:
            initialize_worker(*initializer_args)
            for index, target in enumerate(missing, 1):
                tuned = tune_worker(target)
                checkpoint = checkpoint_record(
                    target,
                    targets_digest,
                    manifest_digest,
                    base_source["source_id"],
                    base_digest,
                    configuration_digest,
                    tuned,
                )
                write_canonical_atomic(checkpoint_path(args.checkpoint_dir, target), checkpoint)
                tuned_by_id[target["target_id"]] = tuned
                print(
                    f"tuned={index}/{len(missing)} target={target['target_id']} "
                    f"combined={tuned['seed_combined_log2']:.9f}",
                    flush=True,
                )
        elif missing:
            with ProcessPoolExecutor(
                max_workers=args.workers,
                initializer=initialize_worker,
                initargs=initializer_args,
            ) as executor:
                futures = {executor.submit(tune_worker, target): target for target in missing}
                completed = 0
                for future in as_completed(futures):
                    target = futures[future]
                    tuned = future.result()
                    checkpoint = checkpoint_record(
                        target,
                        targets_digest,
                        manifest_digest,
                        base_source["source_id"],
                        base_digest,
                        configuration_digest,
                        tuned,
                    )
                    write_canonical_atomic(
                        checkpoint_path(args.checkpoint_dir, target), checkpoint
                    )
                    tuned_by_id[target["target_id"]] = tuned
                    completed += 1
                    print(
                        f"tuned={completed}/{len(missing)} target={target['target_id']} "
                        f"combined={tuned['seed_combined_log2']:.9f}",
                        flush=True,
                    )
        atlas = assemble_atlas(
            targets,
            tuned_by_id,
            args.targets,
            targets_digest,
            manifest_digest,
            base_source["source_id"],
            base_digest,
            configuration,
            args.workers,
        )
        digest = write_canonical_atomic(args.output, atlas)
    except (OSError, ValueError, RuntimeError, json.JSONDecodeError) as error:
        raise SystemExit(f"g=8 supplementary atlas tuner: {error}") from error
    print(f"targets={len(targets)} resumed={len(targets) - len(missing)} tuned={len(missing)}")
    print(f"output={args.output}")
    print(f"output_sha256={digest}")
    print("status=DIAGNOSTIC_BINARY64_SUPPLEMENTARY_WITNESS_ATLAS")


if __name__ == "__main__":
    main()
