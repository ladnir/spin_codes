#!/usr/bin/env python3
"""Compare three frozen/canonical g=8 outer families at leaf dual barycenters."""

from __future__ import annotations

import argparse
import json
import math
import os
import time
from fractions import Fraction
from pathlib import Path
from typing import Any

for _variable in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
):
    os.environ[_variable] = "1"

import numpy as np

import extract_packet_group_g8_component_pricing_target as pricing
import produce_packet_group_g8_highdim_adaptive_shards as coordinate
import probe_packet8_outer_profile_enumerator as total_spectrum
import probe_packet8_three_band_linear_bl2 as three_band
import run_packet_group_g8_factorized_geometry_persistence_wave as geometry_wave


ROOT = Path(__file__).resolve().parents[1]
GEOMETRY = ROOT / "out" / "g8_factorized_geometry_persistence_wave.json"
GEOMETRY_SHA256 = "636103bd034315c9c4381f3319d3d908f98be93e607a8efe3d1fd4fcff6bfa4e"
V2_MANIFEST = ROOT / "G8_FACTORIZED_MANIFEST_V2.json"
V2_INNER = ROOT / "out" / "g8_factorized_inner_components_v1.json"
V2_OUTER = ROOT / "out" / "g8_factorized_outer_components_v1.json"
DEFAULT_OUTPUT = ROOT / "out" / "g8_outer_family_price_probe_top8.json"
SCHEMA = "permute-conv.packet-group-g8-outer-family-price-probe.v1"
SAVING_GATE = 50000.0
TOP_LEAVES = 8


def active_top_distinct(state: dict[str, Any]) -> list[dict[str, Any]]:
    ranked = sorted(
        (node for node in state["nodes"].values() if node["state"] == "ACTIVE_LEAF"),
        key=lambda node: (-float(node["contribution_log2"]), node["node_id"]),
    )
    selected, seen = [], set()
    for node in ranked:
        cell = int(node["h2_cell"])
        if cell in seen:
            continue
        seen.add(cell)
        selected.append(node)
        if len(selected) == TOP_LEAVES:
            return selected
    raise RuntimeError("could not select eight distinct h2 leaves")


def current_dual(node: dict[str, Any], bank: geometry_wave.FactorizedBank):
    matrix = np.asarray(node["exact_vertices"], dtype=np.float64)
    normalization = coordinate.diagnostic_normalization(matrix, 262144)
    inner = bank.inner_constants[None, :] - matrix @ bank.inner_charges.T
    outer = bank.outer_constants[None, :] - matrix @ bank.outer_charges.T
    lp = pricing.solve_component_lp(inner, outer, normalization)
    return pricing.rationalize_distribution(
        -np.asarray(lp.ineqlin.marginals, dtype=np.float64),
        geometry_wave.factorized.one_shot.DUAL_DENOMINATOR,
    )


def exact_barycenter(
    node: dict[str, Any], dual: tuple[tuple[int, Fraction], ...]
) -> tuple[Fraction, ...]:
    if sum((weight for _index, weight in dual), Fraction()) != 1:
        raise ValueError("stored leaf dual mass is not one")
    profiles = tuple(tuple(map(int, row)) for row in node["exact_vertices"])
    return pricing.exact_barycenter(profiles, dual)


def reference_key(reference: dict[str, Any]) -> tuple[str, int, str]:
    return (
        str(reference["source_sha256"]), int(reference["row"]),
        str(reference["row_sha256"]),
    )


def v2_conditioned_pool(args: argparse.Namespace, factorized_state: dict[str, Any]):
    base, _bindings = geometry_wave.component.load_component_bank(
        args.manifest, args.atlas, args.old_supplementary, args.new_supplementary
    )
    full = geometry_wave.build_bank(base, factorized_state["candidates"])
    v2 = json.loads(args.v2_outer.read_bytes())
    lookup = {
        reference_key(reference): index
        for index, reference in enumerate(full.outer_references)
    }
    indices = []
    for row in v2["rows"]:
        key = reference_key(row["origin_reference"])
        if key not in lookup:
            raise ValueError("v2 outer origin does not resolve in full conditioned bank")
        indices.append(lookup[key])
    if len(indices) != 34 or len(set(indices)) != 34:
        raise ValueError("unexpected v2 conditioned outer bank cardinality")
    pool = (
        full.outer_constants[np.asarray(indices, dtype=np.int64)],
        full.outer_charges[np.asarray(indices, dtype=np.int64)],
        tuple(full.outer_references[index] for index in indices),
    )
    return pool, full


def conditioned_price(profile: np.ndarray, pool) -> dict[str, Any]:
    constants, charges, references = pool
    values = constants - charges @ profile
    index = int(np.argmin(values))
    return {"price_log2": float(values[index]), "row": index, "reference": references[index]}


def spectrum_price(profile: np.ndarray) -> dict[str, Any]:
    value, result, parts = total_spectrum.optimize_outer(profile.tolist())
    _score, log_alpha, log_enumerator, envelope = parts
    log_t = np.concatenate(([0.0], result.x[:8]))
    active = np.flatnonzero(envelope >= np.max(envelope) - 1e-7)
    return {
        "price_log2": float(value / math.log(2.0)),
        "parameters": {
            "log_variables": log_t.tolist(), "log_beta": float(result.x[8]),
            "log_alpha": float(log_alpha), "log_local_enumerator": float(log_enumerator),
            "active_envelope_weights": active.tolist(),
        },
        "optimizer": {"success": bool(result.success), "iterations": int(result.nit), "evaluations": int(result.nfev), "message": str(result.message)},
        "family": "unpunctured-exact-ebch-total-spectrum-envelope-v1",
    }


def three_band_price(profile: np.ndarray) -> dict[str, Any]:
    value, result = three_band.optimize_outer(
        profile.tolist(), optimize_band_coefficients=True, log_fugacity_bound=40.0
    )
    band1 = float(result.x[8])
    return {
        "price_log2": float(value / math.log(2.0)),
        "parameters": {
            "log_variables": np.concatenate(([0.0], result.x[:8])).tolist(),
            "band_coefficients": [1.0 - band1, band1, band1],
        },
        "optimizer": {"success": bool(result.success), "iterations": int(result.nit), "evaluations": int(result.nfev), "message": str(result.message)},
        "family": "unpunctured-three-band-linear-bl-direct-sum-v1",
    }


def build(args: argparse.Namespace) -> dict[str, Any]:
    started = time.perf_counter()
    bindings = {
        "geometry_sha256": coordinate.sha256_path(args.geometry),
        "v2_manifest_sha256": coordinate.sha256_path(args.v2_manifest),
        "v2_inner_components_sha256": coordinate.sha256_path(args.v2_inner),
        "v2_outer_components_sha256": coordinate.sha256_path(args.v2_outer),
        "factorized_checkpoint_sha256": coordinate.sha256_path(geometry_wave.FACTORIZED),
    }
    if bindings["geometry_sha256"] != GEOMETRY_SHA256:
        raise ValueError("frozen 197-leaf geometry changed")
    state = json.loads(args.geometry.read_bytes())
    factorized_state = json.loads(geometry_wave.FACTORIZED.read_bytes())
    pool, full_bank = v2_conditioned_pool(args, factorized_state)
    rows = []
    for node in active_top_distinct(state):
        dual = current_dual(node, full_bank)
        barycenter = exact_barycenter(node, dual)
        profile = np.asarray([float(value) for value in barycenter], dtype=np.float64)
        if abs(float(np.sum(profile)) - 262144.0) > 1e-8:
            raise ValueError("dual barycenter mass changed")
        conditioned = conditioned_price(profile, pool)
        spectrum = spectrum_price(profile)
        band = three_band_price(profile)
        for candidate in (spectrum, band):
            candidate["saving_vs_conditioned_bits"] = conditioned["price_log2"] - candidate["price_log2"]
            candidate["saves_at_least_50000"] = candidate["saving_vs_conditioned_bits"] >= SAVING_GATE
        rows.append({
            "node_id": node["node_id"], "h2_cell": int(node["h2_cell"]),
            "contribution_log2": float(node["contribution_log2"]),
            "vertex_dual": [
                {"vertex": index, "weight": pricing.fraction_text(weight)}
                for index, weight in dual
            ],
            "barycenter_exact": [pricing.fraction_text(value) for value in barycenter],
            "conditioned_row_v2": conditioned, "total_spectrum": spectrum,
            "three_band_linear_bl": band,
        })
    best = max(
        (
            (row[family]["saving_vs_conditioned_bits"], row["node_id"], family)
            for row in rows for family in ("total_spectrum", "three_band_linear_bl")
        ),
        key=lambda item: (item[0], item[1], item[2]),
    )
    return {
        "schema": SCHEMA, "status": "DIAGNOSTIC_BINARY64_OUTER_FAMILY_PRICE_COMPARISON",
        "source_bindings": bindings,
        "configuration": {
            "top_distinct_h2_leaves": TOP_LEAVES, "serialized": True,
            "conditioned_rows": len(pool[2]), "saving_gate_bits": SAVING_GATE,
            "total_spectrum_wrapper": "probe_packet8_outer_profile_enumerator.optimize_outer",
            "three_band_wrapper": "probe_packet8_three_band_linear_bl2.optimize_outer(optimize_band_coefficients=True)",
        },
        "rows": rows,
        "summary": {
            "any_family_saves_at_least_50000": best[0] >= SAVING_GATE,
            "best_saving_bits": best[0], "best_node_id": best[1], "best_family": best[2],
            "total_spectrum_pass_count": sum(row["total_spectrum"]["saves_at_least_50000"] for row in rows),
            "three_band_pass_count": sum(row["three_band_linear_bl"]["saves_at_least_50000"] for row in rows),
        },
        "runtime_seconds": time.perf_counter() - started,
        "scope_limit": "Binary64 pricing diagnostic at fixed dual barycenters; no outward or completeness claim.",
    }


def self_test(args: argparse.Namespace) -> None:
    state = json.loads(args.geometry.read_bytes())
    cells = [int(node["h2_cell"]) for node in active_top_distinct(state)]
    if cells != [73, 39, 67, 79, 22, 95, 46, 60]:
        raise SystemExit(f"frozen top8 distinct-h2 ranking changed: {cells}")
    factorized_state = json.loads(geometry_wave.FACTORIZED.read_bytes())
    pool, full_bank = v2_conditioned_pool(args, factorized_state)
    for node in active_top_distinct(state):
        barycenter = exact_barycenter(node, current_dual(node, full_bank))
        if sum(barycenter, Fraction()) != 262144:
            raise SystemExit("exact dual barycenter mass failed")
    if len(pool[2]) != 34:
        raise SystemExit("v2 conditioned pool static check failed")
    print("frozen_top8_distinct_h2_073_039_067_079_022_095_046_060=PASS")
    print("exact_dual_barycenter_mass=PASS")
    print("v2_conditioned_outer_rows_34=PASS")
    print("status=PASS_G8_OUTER_FAMILY_PRICE_PROBE_SELF_TEST")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geometry", type=Path, default=GEOMETRY)
    parser.add_argument("--v2-manifest", type=Path, default=V2_MANIFEST)
    parser.add_argument("--v2-inner", type=Path, default=V2_INNER)
    parser.add_argument("--v2-outer", type=Path, default=V2_OUTER)
    parser.add_argument("--manifest", type=Path, default=geometry_wave.component.h2.DEFAULT_MANIFEST)
    parser.add_argument("--atlas", type=Path, default=geometry_wave.component.h2.DEFAULT_ATLAS)
    parser.add_argument("--old-supplementary", type=Path, default=geometry_wave.component.h2.DEFAULT_SUPPLEMENTARY)
    parser.add_argument("--new-supplementary", type=Path, default=geometry_wave.component.NEW_ATLAS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test(args)
        return
    report = build(args)
    digest = geometry_wave.engine.atomic_write(args.output, report)
    print(f"any_saves_50000={report['summary']['any_family_saves_at_least_50000']}")
    print(f"best_saving_bits={report['summary']['best_saving_bits']:.12f}")
    print(f"best_node={report['summary']['best_node_id']}")
    print(f"best_family={report['summary']['best_family']}")
    print(f"runtime_seconds={report['runtime_seconds']:.3f}")
    print(f"output={args.output}")
    print(f"sha256={digest}")


if __name__ == "__main__":
    main()
