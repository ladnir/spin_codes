#!/usr/bin/env python3
"""Build the finite-SPIN landscape database from explicitly catalogued results."""

from __future__ import annotations

import argparse
import csv
import decimal
import hashlib
import json
import math
import os
import pathlib
import re
import sqlite3
import tempfile
from typing import Any, Iterable


HERE = pathlib.Path(__file__).resolve().parent
SMALL_K = HERE.parent / "small_k_replay"
SCHEMA = HERE / "schema.sql"
CATALOG = HERE / "catalog.json"
DEFAULT_DB = HERE / "spin_landscape.sqlite3"
SCHEMA_VERSION = "4"


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def optional_int(value: Any) -> int | None:
    if value is None or value == "":
        return None
    return int(value)


def optional_float(value: Any) -> float | None:
    if value is None or value == "":
        return None
    return float(value)


def infer_family(label: str) -> str:
    lowered = label.lower()
    if lowered.startswith("random"):
        return "random"
    if "rm(" in lowered:
        return "rm"
    if "bch" in lowered:
        return "bch"
    raise ValueError(f"cannot infer outer family from {label!r}")


def parse_outer_parameters(label: str, row: dict[str, Any]) -> tuple[int, int, int | None]:
    block = optional_int(row.get("block_bits"))
    dimension = optional_int(row.get("dimension"))
    match = re.search(r"\[(\d+)\s*,\s*(\d+)(?:\s*,\s*(\d+))?\]", label)
    if block is None and match:
        block = int(match.group(1))
    if dimension is None and match:
        dimension = int(match.group(2))
    distance = int(match.group(3)) if match and match.group(3) else None
    if block is None or dimension is None:
        raise ValueError(f"missing outer parameters for {label!r}")
    return block, dimension, distance


def outer_model_kind(label: str, row: dict[str, Any]) -> tuple[str, str]:
    family = infer_family(label)
    model = str(row.get("outer_model", "")).lower()
    if model == 'random-simultaneous-spectrum-caps':
        if not row.get('setup_event_id') or row.get('setup_failure_bits') in (None, ''):
            raise ValueError('conditional random bounds require an explicit shared setup event')
        return 'random_setup_spectrum_caps', 'one reused random constituent; conditional on simultaneous shell caps'
    if family == "random" or "ensemble" in model:
        return "random_ensemble_expectation", "one sampled constituent reused; Q1 ensemble expectation"
    return "fixed_exact_spectrum", "one fixed constituent reused in every outer row"


def slug(text: str) -> str:
    value = re.sub(r"[^a-z0-9]+", "-", text.lower()).strip("-")
    return value or "item"


def register_source(
    db: sqlite3.Connection,
    path: pathlib.Path,
    ingestor: str,
    source_status: str,
    notes: str = "",
) -> int:
    relative = path.relative_to(HERE.parent).as_posix()
    db.execute(
        "INSERT INTO sources(path, sha256, ingestor, source_status, notes) VALUES(?,?,?,?,?)",
        (relative, sha256(path), ingestor, source_status, notes),
    )
    return int(db.execute("SELECT last_insert_rowid()").fetchone()[0])


def register_outer(db: sqlite3.Connection, label: str, row: dict[str, Any]) -> str:
    family = infer_family(label)
    block, dimension, distance = parse_outer_parameters(label, row)
    model_kind, reuse_policy = outer_model_kind(label, row)
    outer_id = f"{family}-b{block}-k{dimension}-{slug(model_kind)}"
    db.execute(
        """
        INSERT INTO outer_models(
            outer_id,label,family,model_kind,block_bits,dimension,
            minimum_distance,reuse_policy
        ) VALUES(?,?,?,?,?,?,?,?)
        ON CONFLICT(outer_id) DO UPDATE SET
            minimum_distance=COALESCE(outer_models.minimum_distance,excluded.minimum_distance)
        """,
        (outer_id, label, family, model_kind, block, dimension, distance, reuse_policy),
    )
    return outer_id


def register_inner(
    db: sqlite3.Connection,
    row: dict[str, Any],
    map_tag: str,
) -> str:
    step = int(row["step_bits"])
    state = int(row["state_bits"])
    persistence = optional_float(row.get("persistence_exponent"))
    if persistence is None:
        persistence = state + math.log2(step)
    a_distance = optional_int(row.get("a_minimum_distance"))
    kernel_distance = optional_int(row.get("kernel_minimum_distance"))
    kernel_weight_four = optional_int(row.get("kernel_weight_four"))
    metric_tag = f"a{a_distance}-kd{kernel_distance}-k4{kernel_weight_four}"
    inner_id = f"rm2sub-t{step}-s{state}-{slug(map_tag)}-{slug(metric_tag)}"
    db.execute(
        """
        INSERT OR IGNORE INTO inner_configs(
            inner_id,family,step_bits,state_bits,persistence_exponent,map_tag,
            a_minimum_distance,kernel_minimum_distance,kernel_weight_four
        ) VALUES(?,?,?,?,?,?,?,?,?)
        """,
        (
            inner_id,
            "rm2sub",
            step,
            state,
            persistence,
            map_tag,
            a_distance,
            kernel_distance,
            kernel_weight_four,
        ),
    )
    return inner_id


def classify_result(label: str, outer_model: str) -> str:
    if infer_family(label) == "random" or "ensemble" in outer_model.lower():
        return "reference"
    return "diagnostic"


def insert_result(
    db: sqlite3.Connection,
    *,
    source_id: int,
    source_locator: str,
    study: str,
    outer_id: str,
    inner_id: str,
    row: dict[str, Any],
    occupation_min: int,
    occupation_max: int,
    result_class: str,
    arithmetic: str,
    coverage_kind: str,
    margin_bits: float,
    margin_bits_text: str | None = None,
    failure_upper_text: str | None = None,
    notes: str = "",
) -> None:
    message_bits = int(row["message_bits"])
    message_exponent = optional_int(row.get("message_exponent"))
    if message_exponent is None and message_bits > 0 and message_bits & (message_bits - 1) == 0:
        message_exponent = message_bits.bit_length() - 1
    result_id = hashlib.sha256(f"{source_id}:{source_locator}".encode()).hexdigest()[:24]
    db.execute(
        """
        INSERT INTO results(
            result_id,source_id,source_locator,study,outer_id,inner_id,
            message_exponent,message_bits,output_bits,outer_rows,epochs_per_region,
            bad_weight,distance_target,occupation_min,occupation_max,coverage_kind,
            result_class,arithmetic,margin_bits,margin_bits_text,failure_upper_text,
            dominant_weight,dominant_second_weight,dominant_log_surprisal,dominant_witness_at_grid_edge,comparison_eligible,
            setup_event_id,setup_failure_bits,witness_shift,notes
        ) VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
        """,
        (
            result_id,
            source_id,
            source_locator,
            study,
            outer_id,
            inner_id,
            message_exponent,
            message_bits,
            int(row["output_bits"]),
            int(row["outer_rows"]),
            optional_int(row.get("epochs_per_region")),
            int(row["bad_weight"]),
            0.10,
            occupation_min,
            occupation_max,
            coverage_kind,
            result_class,
            arithmetic,
            margin_bits,
            margin_bits_text,
            failure_upper_text,
            optional_int(row.get("dominant_weight")),
            optional_int(row.get("dominant_second_weight")),
            optional_float(row.get("dominant_log_surprisal")),
            optional_int(row.get("dominant_witness_at_grid_edge")),
            int(row.get("comparison_eligible", 0)),
            row.get('setup_event_id') or None,
            optional_int(row.get('setup_failure_bits')),
            optional_float(row.get('witness_shift')),
            notes,
        ),
    )


def import_csv(db: sqlite3.Connection, spec: dict[str, Any]) -> None:
    path = (HERE if spec.get("base") == "landscape_db" else SMALL_K) / spec["path"]
    source_id = register_source(
        db,
        path,
        "csv_rows_v1",
        "binary64_diagnostic",
        "Rows with inner=RandomStepConv are rejected by the catalogued filter.",
    )
    review = spec.get("transfer_review_status", "historical_pending_review")
    if review == "activation_aware":
        manifest_path = HERE / spec["manifest"]
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        for relative, digest in manifest["source_sha256"].items():
            if sha256(HERE.parent.parent.parent / relative) != digest:
                raise ValueError(f"changed screen dependency: {relative}")
        if manifest["csv_sha256"] != sha256(path):
            raise ValueError(f"changed screen CSV: {path}")
    db.execute("UPDATE sources SET transfer_review_status=? WHERE source_id=?", (review, source_id))
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    required = {"message_bits", "output_bits", "outer_rows", "bad_weight", "margin_bits"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"unexpected CSV schema: {path}")
    imported = 0
    for index, row in enumerate(rows, start=2):
        row["comparison_eligible"] = int(
            review == "activation_aware"
            and int(row.get("message_exponent") or 0) not in spec.get("superseded_exponents", [])
        )
        required_inner = spec.get("require_inner")
        if required_inner is not None and row.get("inner") != required_inner:
            continue
        if row.get("inner") not in (None, "", "RM2Sub"):
            raise ValueError(f"non-RM2Sub row escaped filtering at {path}:{index}")
        label = row.get("series") or row.get("constituent")
        if not label:
            raise ValueError(f"missing constituent label at {path}:{index}")
        occupation = int(row.get("occupation") or spec.get("occupation") or 1)
        occupation_min = int(row.get('occupation_min') or occupation)
        occupation_max = int(row.get('occupation_max') or occupation)
        coverage_kind = row.get('coverage_kind') or ('single_occupation' if occupation_min==occupation_max else 'occupation_union')
        if not 1 <= occupation_min <= occupation_max <= int(row['outer_rows']):
            raise ValueError('invalid occupation range')
        if (coverage_kind=='single_occupation') != (occupation_min==occupation_max):
            raise ValueError('occupation range and coverage kind disagree')
        if coverage_kind=='full_distance' and (occupation_min!=1 or occupation_max!=int(row['outer_rows'])):
            raise ValueError('full distance requires every occupation')
        outer_id = register_outer(db, label, row)
        inner_id = register_inner(db, row, row.get("map_tag") or spec["map_tag"])
        result_class = classify_result(label, row.get("outer_model", ""))
        insert_result(
            db,
            source_id=source_id,
            source_locator=f"csv-line-{index}",
            study=spec["study"],
            outer_id=outer_id,
            inner_id=inner_id,
            row=row,
            occupation_min=occupation_min,
            occupation_max=occupation_max,
            result_class=result_class,
            arithmetic="nearest_binary64_diagnostic",
            coverage_kind=coverage_kind,
            margin_bits=float(row["margin_bits"]),
            notes=row.get("notes", ""),
        )
        imported += 1
    if imported == 0:
        raise ValueError(f"catalog filter imported no rows from {path}")


def import_q_ladder(db: sqlite3.Connection, spec: dict[str, Any]) -> None:
    path = SMALL_K / spec["path"]
    source_id = register_source(db, path, "q_ladder_v1", "binary64_diagnostic")
    payload = json.loads(path.read_text(encoding="utf-8"))
    parameters = payload["parameters"]
    row = {
        "block_bits": parameters["outer_block_bits"],
        "dimension": parameters["outer_dimension"],
        "message_exponent": parameters["message_exponent"],
        "message_bits": parameters["message_bits"],
        "outer_rows": parameters["outer_rows"],
        "output_bits": parameters["output_bits"],
        "bad_weight": parameters["bad_weight"],
        "step_bits": parameters["step_bits"],
        "state_bits": parameters["state_bits"],
        "persistence_exponent": parameters["persistence_exponent"],
        "epochs_per_region": parameters["epochs_per_region"],
        "a_minimum_distance": parameters.get("a_minimum_distance"),
        "kernel_minimum_distance": parameters.get("kernel_minimum_distance"),
        "kernel_weight_four": parameters.get("kernel_weight_four"),
    }
    label = "RM(4,9) [512,256,32]"
    outer_id = register_outer(db, label, row)
    inner_id = register_inner(db, row, spec["map_tag"])
    occupations = payload.get("occupation_rows")
    if not isinstance(occupations, list) or not occupations:
        raise ValueError(f"missing occupation_rows in {path}")
    for item in occupations:
        q = int(item["occupation"])
        insert_result(
            db,
            source_id=source_id,
            source_locator=f"occupation-{q}",
            study=spec["study"],
            outer_id=outer_id,
            inner_id=inner_id,
            row=row,
            occupation_min=q,
            occupation_max=q,
            result_class="diagnostic",
            arithmetic="nearest_binary64_diagnostic",
            coverage_kind="single_occupation",
            margin_bits=float(item["margin_bits_diagnostic"]),
            notes="Positive recurrence over regular and all-one RM rows.",
        )


INTERVAL_RE = re.compile(r"^\[([^ ]+) \+/- ([^\]]+)\]$")


def interval_upper(text: str) -> decimal.Decimal:
    match = INTERVAL_RE.match(text)
    if not match:
        raise ValueError(f"unsupported interval {text!r}")
    return decimal.Decimal(match.group(1)) + decimal.Decimal(match.group(2))


def decimal_margin_bits(probability_upper: decimal.Decimal) -> float:
    if probability_upper <= 0:
        raise ValueError(f"probability upper bound must be positive: {probability_upper}")
    with decimal.localcontext() as context:
        context.prec = 100
        return float(-(probability_upper.ln() / decimal.Decimal(2).ln()))


def import_certificate(db: sqlite3.Connection, spec: dict[str, Any]) -> None:
    path = SMALL_K / spec["path"]
    source_id = register_source(db, path, "full_certificate_v1", "outward_certificate")
    db.execute(
        "UPDATE sources SET transfer_review_status='activation_under_reaudit', notes=? WHERE source_id=?",
        ("Historical outward receipt; activation-state invariant requires re-audit. See PARAMETER_LANDSCAPE_HANDOFF.md.", source_id),
    )
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("status") != "OUTWARD_FULL_DISTANCE_CERTIFICATE":
        raise ValueError(f"not a full outward certificate: {path}")
    theorem = payload["theorem_parameters"]
    claim = payload["claim"]
    row = {
        "block_bits": 512,
        "dimension": 256,
        "message_bits": theorem["message_bits"],
        "output_bits": theorem["output_bits"],
        "outer_rows": 256,
        "epochs_per_region": 4,
        "bad_weight": theorem["bad_weight_upper"],
        "step_bits": 64,
        "state_bits": 14,
        "persistence_exponent": 20,
        "a_minimum_distance": 24,
        "kernel_minimum_distance": 6,
        "kernel_weight_four": 0,
    }
    outer_id = register_outer(db, "RM(4,9) [512,256,32]", row)
    inner_id = register_inner(db, row, spec["map_tag"])
    insert_result(
        db,
        source_id=source_id,
        source_locator="full-distance-union",
        study=spec["study"],
        outer_id=outer_id,
        inner_id=inner_id,
        row=row,
        occupation_min=int(payload["coverage"]["occupation_interval"][0]),
        occupation_max=int(payload["coverage"]["occupation_interval"][1]),
        result_class="certified",
        arithmetic="arb_256_outward",
        coverage_kind="full_distance",
        margin_bits=float(claim["margin_bits_lower"]),
        margin_bits_text=str(claim["margin_bits_lower"]),
        failure_upper_text=claim["failure_probability_upper_interval"],
        notes="Database row indexes the certificate; the source receipt remains authoritative.",
    )
    for index, component in enumerate(payload["components"]):
        q_min, q_max = map(int, component["occupation_interval"])
        failure_text = component["probability_upper_interval"]
        upper = interval_upper(failure_text)
        insert_result(
            db,
            source_id=source_id,
            source_locator=f"component-{index}-{component['role']}",
            study=spec["study"] + "_components",
            outer_id=outer_id,
            inner_id=inner_id,
            row=row,
            occupation_min=q_min,
            occupation_max=q_max,
            result_class="certified",
            arithmetic="arb_256_outward",
            coverage_kind="single_occupation" if q_min == q_max else "occupation_union",
            margin_bits=decimal_margin_bits(upper),
            failure_upper_text=failure_text,
            notes=f"Component role: {component['role']}.",
        )


def build(output: pathlib.Path) -> None:
    catalog = json.loads(CATALOG.read_text(encoding="utf-8"))
    if catalog.get("schema") != "spin-landscape-import-catalog-v1":
        raise ValueError("unsupported import catalog schema")
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        prefix=output.name + ".", suffix=".tmp", dir=output.parent, delete=False
    ) as handle:
        temporary = pathlib.Path(handle.name)
    try:
        db = sqlite3.connect(temporary)
        try:
            db.executescript(SCHEMA.read_text(encoding="utf-8"))
            db.executemany(
                "INSERT INTO schema_metadata(key,value) VALUES(?,?)",
                [
                    ("schema_version", SCHEMA_VERSION),
                    ("scope", catalog["scope"]),
                    ("catalog_sha256", sha256(CATALOG)),
                ],
            )
            for spec in catalog["csv_sources"]:
                import_csv(db, spec)
            for spec in catalog["q_ladder_sources"]:
                import_q_ladder(db, spec)
            for spec in catalog["certificate_sources"]:
                import_certificate(db, spec)
            # Append the new sources after every historical source, preserving
            # source IDs and result IDs in the original 414-row snapshot.
            for spec in catalog.get("activation_sources", []):
                import_csv(db, spec)
            db.commit()
            violations = db.execute("PRAGMA foreign_key_check").fetchall()
            if violations:
                raise ValueError(f"foreign-key violations: {violations}")
            non_rm2sub = db.execute(
                "SELECT COUNT(*) FROM inner_configs WHERE family <> 'rm2sub'"
            ).fetchone()[0]
            if non_rm2sub:
                raise ValueError("database contains a non-RM2Sub inner")
        finally:
            db.close()
        os.replace(temporary, output)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=pathlib.Path, default=DEFAULT_DB)
    return parser.parse_args()


if __name__ == "__main__":
    arguments = parse_args()
    build(arguments.output.resolve())
    print(arguments.output.resolve())
