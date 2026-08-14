#!/usr/bin/env python3
"""Repair the one-conditioned-row outer for the frozen sloped layout.

The block ``(u,0)`` belongs to tile ``u`` in all three bands because every
band slope times lane zero is zero.  These 256 diagonal blocks form one
matching that covers every tile once.  No other lane has this property.

Condition the diagonal block in each unpunctured tile and apply the canonical
one-row pair-Cauchy bound.  This gives 256 copies of the normal-tile factor.
Do not move the conditioned row to a random punctured block.  Instead, treat
each of the 128 graph replacements pointwise.  Replacing one bit changes its
packet class from ``old`` to ``new`` with ``|new-old|<=1``.  For base-two
packet charges ``q_j=log2(t_j)``, add

    128 * max_{|new-old|<=1} (q_new-q_old)

to the unpunctured affine constant.  This is the audited ordered-adjacent
orientation.  Binary64 results are discovery diagnostics until an independent
outward evaluator certifies frozen parameters.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import tempfile
import time
from fractions import Fraction
from pathlib import Path
from typing import Any, Sequence

for _variable in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
):
    os.environ[_variable] = "1"

import numpy as np
from scipy.special import logsumexp

from packet_group_outer_profile import K, atom_count
from probe_packet_group_conditioned_row_outer import (
    evaluate_conditioned_outer,
    load_split_spectrum,
)


ROOT = Path(__file__).resolve().parents[1]
SPECTRUM01 = ROOT / "out/ebch85_band01_split_spectrum.csv"
PUNCTURED01 = ROOT / "out/ebch84_punctured_band01_split_spectrum.csv"
SPECTRUM12 = ROOT / "out/ebch86_band12_split_spectrum.csv"
G4_LEADER = ROOT / "out/g4_conditioned_row_residual_best140_asymmetric.json"
G8_TOP8 = ROOT / "out/g8_outer_family_price_probe_top8.json"
G8_V2_OUTER = ROOT / "out/g8_factorized_outer_components_v1.json"
TILE_CERTIFICATE = ROOT / "scripts/certify_three_band_tile_map.py"
CANONICAL_EVALUATOR = ROOT / "scripts/probe_packet_group_conditioned_row_outer.py"
DEFAULT_OUTPUT = ROOT / "out/g8_sloped_conditioned_row_repair_diagnostic.json"

PINNED_SHA256 = {
    SPECTRUM01: "10315518bf5aa72199a48bb201e9f11d159c61bd9a83aa6f0fc1948887cf1f06",
    PUNCTURED01: "78fd4730ce818fabf923c237f88115e7c9ad4ecc6588f347c8e653ff93839ce3",
    SPECTRUM12: "d7bf4816ab602a03bfdb526cabfa437ab412e5894e6ea64987ae736752f448e1",
    G4_LEADER: "20a967d9762411ae3d93b3c4e4815befcfab33f7bed20aac9e068084af9d6b13",
    G8_TOP8: "797e097f39fdfea4b30bd62a7c8c7efb11450beb0446776b5d3ccc07da85c6e5",
    G8_V2_OUTER: "ed47ab3fa6aec959dbcc6f1a2e44cdd574fd1b0dcdfafd299f5a08ff03494b1d",
    TILE_CERTIFICATE: "45b69f70d4a5cc459a881b3321004f942f8ac83cd7f1f45a3616f87e1bd1d2e3",
    CANONICAL_EVALUATOR: "c6ac971e857aa2f043ad0c61175cae151eab3d54c14442d100627760cd2d3774",
}
SLOPES = (0, 9, 20)
TILES = 256
GRAPH_REPLACEMENTS = 128
G4_LEADER_ROW = 132


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_sources() -> dict[str, dict[str, str]]:
    bindings = {}
    for path, expected in PINNED_SHA256.items():
        observed = sha256_file(path)
        if observed != expected:
            raise ValueError(
                f"sloped conditioned-row frozen input changed: {path}; "
                f"expected {expected}, observed {observed}"
            )
        bindings[path.name] = {"path": str(path), "sha256": observed}
    return bindings


def atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as handle:
            json.dump(payload, handle, sort_keys=True, separators=(",", ":"))
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, path)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def load_spectra() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    return (
        load_split_spectrum(SPECTRUM01),
        load_split_spectrum(PUNCTURED01, count_field="pair_count", divisor=42),
        load_split_spectrum(SPECTRUM12),
    )


def diagonal_matching_audit() -> dict:
    common_lanes = [
        lane
        for lane in range(256)
        if all((slope * lane) % 256 == 0 for slope in SLOPES)
    ]
    if common_lanes != [0] or math.gcd(9, 256) != 1:
        raise AssertionError("the frozen sloped layout has an unexpected common lane")
    covered = []
    for tile in range(256):
        images = tuple((tile + slope * 0) % 256 for slope in SLOPES)
        if images != (tile, tile, tile):
            raise AssertionError("diagonal lane does not align the three band tiles")
        covered.append(images[0])
    if sorted(covered) != list(range(256)):
        raise AssertionError("diagonal matching does not cover every tile exactly once")
    return {
        "slopes": list(SLOPES),
        "common_lane_classes_mod_256": common_lanes,
        "matching_blocks": 256,
        "tiles_covered_once_per_band": True,
        "proof": "9 is invertible modulo 256, so 9*l=0 implies l=0",
    }


def ordered_adjacent_tax(charge_log2: Sequence[float]) -> tuple[float, dict]:
    q = np.asarray(charge_log2, dtype=np.float64)
    if q.ndim != 1 or q.size < 2 or not np.all(np.isfinite(q)):
        raise ValueError("ordered-adjacent tax requires finite packet charges")
    candidates = [
        (float(q[new] - q[old]), old, new)
        for old in range(q.size)
        for new in range(q.size)
        if abs(new - old) <= 1
    ]
    gain, old, new = max(candidates, key=lambda row: (row[0], -row[1], -row[2]))
    if gain < 0.0:
        raise AssertionError("same-class replacement must make the maximum nonnegative")
    return GRAPH_REPLACEMENTS * gain, {
        "per_replacement_max_gain_bits": gain,
        "worst_old_packet_class": old,
        "worst_new_packet_class": new,
        "replacement_count": GRAPH_REPLACEMENTS,
        "orientation": "q_new-q_old",
        "constraint": "abs(new-old)<=1",
    }


def validate_profile(profile: Sequence[float], group_bits: int) -> np.ndarray:
    vector = np.asarray(profile, dtype=np.float64)
    if vector.shape != (group_bits + 1,) or not np.all(np.isfinite(vector)):
        raise ValueError("sloped repair received a malformed packet profile")
    if np.min(vector) < 0.0 or float(np.sum(vector)) != float(atom_count(group_bits)):
        raise ValueError("sloped repair profile has the wrong mass")
    return vector


def evaluate_repair(
    profile: Sequence[float],
    log_variables: Sequence[float],
    band_coefficients: Sequence[float],
    theta: float,
    group_bits: int,
    spectra: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> dict[str, Any]:
    vector = validate_profile(profile, group_bits)
    log_t = np.asarray(log_variables, dtype=np.float64)
    coefficients = tuple(float(value) for value in band_coefficients)
    if log_t.shape != (group_bits + 1,) or log_t[0] != 0.0 or not np.all(np.isfinite(log_t)):
        raise ValueError("sloped repair received malformed packet log variables")
    if len(coefficients) != 3:
        raise ValueError("sloped repair requires three band coefficients")
    spectrum01, punctured01, spectrum12 = spectra
    old_exact_graph_log2, old_details = evaluate_conditioned_outer(
        vector.tolist(),
        log_t,
        coefficients[1],
        theta,
        spectrum01,
        punctured01,
        spectrum12,
        conditioned_rows=1,
        band_coefficients=coefficients,
        group_bits=group_bits,
    )
    charge = log_t / math.log(2.0)
    unpunctured_constant = TILES * float(old_details["normal_tile_log"]) / math.log(2.0)
    replacement_tax, tax_details = ordered_adjacent_tax(charge)
    repaired_constant = unpunctured_constant + replacement_tax
    affine_value = repaired_constant - float(vector @ charge)
    direct_value = (
        TILES * float(old_details["normal_tile_log"])
        - float(vector @ log_t)
    ) / math.log(2.0) + replacement_tax
    identity_error = abs(affine_value - direct_value)
    if identity_error > 2e-9:
        raise AssertionError("sloped repaired scalar-affine identity failed")
    return {
        "group_bits": group_bits,
        "profile": vector.tolist(),
        "log_variables": log_t.tolist(),
        "band_coefficients": list(coefficients),
        "pair_cauchy_theta": theta,
        "old_exact_graph_log2": old_exact_graph_log2,
        "old_exact_graph_status": "LAYOUT_INVALID_REFERENCE_VALUE_ONLY",
        "unpunctured_diagonal_matching_log2": (
            unpunctured_constant - float(vector @ charge)
        ),
        "ordered_adjacent_replacement_tax_log2": replacement_tax,
        "repaired_log2": affine_value,
        "repaired_minus_old_bits": affine_value - old_exact_graph_log2,
        "affine": {
            "constant_log2": repaired_constant,
            "charge_log2": charge.tolist(),
            "scalar_from_affine_log2": affine_value,
            "scalar_from_direct_log2": direct_value,
            "absolute_identity_error_bits": identity_error,
            "orientation": "constant - profile dot charge",
        },
        "replacement_tax": tax_details,
        "normal_tile_log_natural": float(old_details["normal_tile_log"]),
        "proof_scope": (
            "condition the unique diagonal l=0 matching in all 256 unpunctured "
            "tiles; bound all 128 puncture/graph changes pointwise"
        ),
    }


def load_targets() -> list[dict[str, Any]]:
    g4_artifact = json.loads(G4_LEADER.read_bytes())
    g4 = g4_artifact["rows"][G4_LEADER_ROW]
    if g4.get("outer_type") != "conditioned_row_exact_graph_asymmetric":
        raise ValueError("canonical g4 leader is not the expected conditioned-row row")
    g4_details = g4["outer_details"]

    top8 = json.loads(G8_TOP8.read_bytes())
    row = top8["rows"][0]
    if row.get("node_id") != "h2:073" or row["conditioned_row_v2"].get("row") != 20:
        raise ValueError("frozen g8 h2:073 conditioned-row reference changed")
    v2 = json.loads(G8_V2_OUTER.read_bytes())
    v2_row = v2["rows"][20]
    if v2_row.get("origin_reference") != row["conditioned_row_v2"].get("reference"):
        raise ValueError("g8 h2:073 v2 conditioned-row origin no longer resolves")
    parameters = v2_row["parameters"]
    g8_profile = [float(Fraction(value)) for value in row["barycenter_exact"]]
    if sum(g8_profile) != atom_count(8):
        raise ValueError("g8 h2:073 barycenter mass changed")
    return [
        {
            "name": "canonical_g4_leader_row132",
            "group_bits": 4,
            "profile": g4["profile"],
            "log_variables": g4_details["log_variables"],
            "band_coefficients": g4_details["band_coefficients"],
            "theta": g4_details["pair_cauchy_theta"],
            "stored_old_exact_graph_log2": g4["outer_log2"],
        },
        {
            "name": "g8_h2_073_conditioned_row_v2_row20",
            "group_bits": 8,
            "profile": g8_profile,
            "profile_exact": row["barycenter_exact"],
            "log_variables": parameters["log_variables"],
            "band_coefficients": parameters["band_coefficients"],
            "theta": parameters["pair_cauchy_theta"],
            "stored_old_exact_graph_log2": row["conditioned_row_v2"]["price_log2"],
            "origin_reference": v2_row["origin_reference"],
        },
    ]


def spectrum_mass_audit(spectra: tuple[np.ndarray, np.ndarray, np.ndarray]) -> dict:
    labels = ("band01", "punctured_band01_average", "band12")
    result = {}
    for label, table in zip(labels, spectra):
        log_mass = float(logsumexp(table[:, 2]))
        error = abs(log_mass - 64 * math.log(2.0))
        if error > 2e-12:
            raise AssertionError(f"{label} split spectrum mass is not 2^64")
        result[label] = {"log_mass_natural": log_mass, "error_from_64ln2": error}
    return result


def static_self_test() -> None:
    validate_sources()
    matching = diagonal_matching_audit()
    spectra = load_spectra()
    spectrum_mass_audit(spectra)
    targets = load_targets()
    for group_bits in (4, 8):
        zero = evaluate_repair(
            [atom_count(group_bits)] + [0] * group_bits,
            [0.0] * (group_bits + 1),
            (0.4, 0.6, 0.6),
            0.5,
            group_bits,
            spectra,
        )
        if not math.isclose(zero["repaired_log2"], K, rel_tol=0.0, abs_tol=1e-6):
            raise AssertionError(f"g={group_bits} repaired mass identity failed")
        if zero["ordered_adjacent_replacement_tax_log2"] != 0.0:
            raise AssertionError("zero-tilt replacement tax is not zero")
    tax, details = ordered_adjacent_tax([0.0, 1.0, 2.0, 3.0, 4.0])
    if tax != 128.0 or details["worst_new_packet_class"] != details[
        "worst_old_packet_class"
    ] + 1:
        raise AssertionError("ordered-adjacent new-minus-old orientation failed")
    synthetic = evaluate_repair(
        [atom_count(4), 0, 0, 0, 0],
        [0.0, -0.2, 0.1, -0.3, 0.2],
        (0.43, 0.57, 0.61),
        0.37,
        4,
        spectra,
    )
    if synthetic["affine"]["absolute_identity_error_bits"] > 2e-9:
        raise AssertionError("synthetic repaired affine identity failed")
    print("sloped one-conditioned-row repair static self-test: PASS")
    print(f"diagonal_common_lanes={matching['common_lane_classes_mod_256']}")
    print("spectrum_masses=2^64,2^64,2^64")
    print("mass_identity_g4_g8=2^1048576")
    print("ordered_adjacent_orientation=q_new-q_old")
    print(f"frozen_targets={len(targets)}")


def build_report() -> dict[str, Any]:
    started = time.perf_counter()
    bindings = validate_sources()
    matching = diagonal_matching_audit()
    spectra = load_spectra()
    masses = spectrum_mass_audit(spectra)
    rows = []
    for target in load_targets():
        result = evaluate_repair(
            target["profile"],
            target["log_variables"],
            target["band_coefficients"],
            target["theta"],
            target["group_bits"],
            spectra,
        )
        stored = float(target["stored_old_exact_graph_log2"])
        replay_error = abs(result["old_exact_graph_log2"] - stored)
        if replay_error > 2e-8:
            raise ValueError(
                f"old exact-graph replay changed for {target['name']}: {replay_error}"
            )
        rows.append(
            {
                **target,
                **result,
                "stored_old_exact_graph_log2": stored,
                "old_exact_graph_replay_error_bits": replay_error,
            }
        )
    return {
        "schema": "permute-conv.packet-group-sloped-conditioned-row-repair.v1",
        "status": "DIAGNOSTIC_BINARY64_SOUND_SLOPED_LAYOUT_REPAIR_NO_OUTWARD_CLAIM",
        "configuration": {
            "serialized": True,
            "targets": 2,
            "conditioned_lane": 0,
            "unpunctured_tiles": 256,
            "ordered_adjacent_graph_replacements": 128,
            "graph_replacement_rule": (
                "128*max_{abs(new-old)<=1}(q_new-q_old)"
            ),
        },
        "diagonal_matching_audit": matching,
        "spectrum_mass_audit": masses,
        "source_bindings": bindings,
        "rows": rows,
        "summary": {
            "maximum_repair_penalty_vs_old_bits": max(
                row["repaired_minus_old_bits"] for row in rows
            ),
            "minimum_repair_penalty_vs_old_bits": min(
                row["repaired_minus_old_bits"] for row in rows
            ),
            "maximum_affine_identity_error_bits": max(
                row["affine"]["absolute_identity_error_bits"] for row in rows
            ),
        },
        "runtime_seconds": time.perf_counter() - started,
        "scope_limit": (
            "Sound layout repair with binary64 frozen-parameter evaluation; "
            "shared outward certifier remains unchanged"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.self_test:
        static_self_test()
        return
    report = build_report()
    atomic_json(args.output, report)
    print("sloped conditioned-row repair diagnostic: COMPLETE")
    for row in report["rows"]:
        print(
            f"target={row['name']} old={row['old_exact_graph_log2']:.12f} "
            f"unpunctured={row['unpunctured_diagonal_matching_log2']:.12f} "
            f"tax={row['ordered_adjacent_replacement_tax_log2']:.12f} "
            f"repaired={row['repaired_log2']:.12f} "
            f"delta={row['repaired_minus_old_bits']:.12f}"
        )
    print(f"runtime_seconds={report['runtime_seconds']:.6f}")
    print(f"output={args.output}")
    print(f"sha256={sha256_file(args.output)}")


if __name__ == "__main__":
    main()
