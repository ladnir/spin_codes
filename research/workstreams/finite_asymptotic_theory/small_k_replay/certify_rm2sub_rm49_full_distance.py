#!/usr/bin/env python3
"""Assemble the authenticated Q=1,...,256 outward distance certificate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from flint import arb, ctx

import certify_rm2sub_rm49_q1_outward as q1


HERE = Path(__file__).resolve().parent
MANIFEST = HERE / "RM2SUB_RM49_FULL_CERTIFICATE_MANIFEST.json"
OUTPUT = HERE / "rm2sub_rm49_t64_s14_full_distance_outward.json"
CERTIFICATE_ROLES = (
    "q1_certificate",
    "q2_certificate",
    "q3_q4_certificate",
    "q05_q29_certificate",
    "q30_certificate",
    "q031_q256_certificate",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_manifest() -> tuple[dict[str, object], dict[str, Path]]:
    payload = json.loads(MANIFEST.read_text(encoding="utf-8"))
    paths: dict[str, Path] = {}
    for row in payload["files"]:
        path = (HERE / row["path"]).resolve()
        if sha256(path) != row["sha256"]:
            raise AssertionError(f"hash mismatch for {row['role']}")
        paths[str(row["role"])] = path
    return payload, paths


def covered_occupations(role: str, payload: dict[str, object]) -> list[int]:
    if role == "q1_certificate":
        return [1]
    if role == "q2_certificate":
        return [2]
    if role == "q3_q4_certificate":
        return [int(row["occupation"]) for row in payload["occupation_rows"]]
    if role == "q05_q29_certificate":
        return [int(row["occupation"]) for row in payload["occupation_rows"]]
    if role == "q30_certificate":
        return [int(payload["coverage"]["occupation"])]
    if role == "q031_q256_certificate":
        return [
            int(row["occupation"])
            for chunk in payload["chunk_rows"]
            for row in chunk["occupation_rows"]
        ]
    raise AssertionError(f"unknown certificate role {role}")


def probability_upper(payload: dict[str, object]) -> arb:
    claim = payload["claim"]
    if "failure_probability_upper_interval" in claim:
        return arb(str(claim["failure_probability_upper_interval"]))
    if "log2_failure_probability_upper" in claim:
        exponent = arb(float(claim["log2_failure_probability_upper"]))
        return arb(2) ** exponent
    raise AssertionError("component certificate has no probability upper bound")


def main() -> None:
    ctx.prec = 256
    manifest, paths = load_manifest()
    components = []
    coverage: list[int] = []
    total = arb(0)
    for role in CERTIFICATE_ROLES:
        payload = json.loads(paths[role].read_text(encoding="utf-8"))
        if not str(payload["status"]).startswith("OUTWARD_"):
            raise AssertionError(f"{role} is not an outward certificate")
        occupations = covered_occupations(role, payload)
        coverage.extend(occupations)
        probability = probability_upper(payload)
        total += probability
        components.append(
            {
                "role": role,
                "occupation_interval": [min(occupations), max(occupations)],
                "occupation_count": len(occupations),
                "probability_upper_interval": str(probability),
            }
        )

    if coverage != list(range(1, q1.OUTER_ROWS + 1)):
        raise AssertionError("component certificates do not cover Q=1,...,256 exactly")
    log2_total = q1.log2_interval(total)
    margin_lower = -q1.upper_float(log2_total)
    if not total < arb(1) / (1 << q1.TARGET_MARGIN_BITS):
        raise AssertionError("full occupation union does not clear 40 bits")

    payload = {
        "schema": "rm2sub-rm49-full-distance-outward-certificate-v1",
        "status": "OUTWARD_FULL_DISTANCE_CERTIFICATE",
        "theorem_parameters": {
            "message_bits": 1 << 16,
            "output_bits": 1 << 17,
            "outer": "one fixed RM(4,9) [512,256,32] constituent repeated in 256 rows",
            "inner": "fixed audited RM2Sub A/B pair with t=64 and s=14",
            "routing": "independent uniform row-coordinate permutations and independent uniform permutations in 512 transposed regions",
            "inner_randomness": "one independent uniform nonzero field scalar per RM2Sub epoch",
            "bad_weight_upper": q1.BAD_WEIGHT,
            "certified_minimum_distance_lower": q1.BAD_WEIGHT + 1,
        },
        "claim": {
            "event": "there exists a nonzero message whose encoded word has Hamming weight at most 13107",
            "failure_probability_upper_interval": str(total),
            "log2_failure_probability_interval": str(log2_total),
            "margin_bits_lower": margin_lower,
            "clears_40_bits": True,
        },
        "coverage": {
            "occupation_interval": [1, q1.OUTER_ROWS],
            "occupation_count": len(coverage),
            "complete_without_overlap": True,
        },
        "components": components,
        "arithmetic": {
            "library": "python-flint Arb",
            "precision_bits": ctx.prec,
            "union_rule": "sum the six authenticated positive component upper bounds",
        },
        "inputs": {
            "manifest": MANIFEST.name,
            "manifest_sha256": sha256(MANIFEST),
            "checker_sha256": sha256(Path(__file__).resolve()),
            "files": manifest["files"],
        },
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"status,{payload['status']}")
    print(f"minimum_distance_lower,{q1.BAD_WEIGHT + 1}")
    print(f"margin_bits_lower,{margin_lower:.12f}")
    print(f"wrote,{OUTPUT}")


if __name__ == "__main__":
    main()
