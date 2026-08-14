#!/usr/bin/env python3
"""Replay fixed g=8 leaves with independent inner and outer mixtures."""

from __future__ import annotations

import argparse
import json
import math
import os
import time
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any

for _variable in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
):
    os.environ[_variable] = "1"

import numpy as np
from scipy.optimize import linprog

import evaluate_packet_group_g8_cumulative_boxes as h2
import plan_packet_group_g8_lowdim_root_cover as minimax
import produce_packet_group_g8_highdim_adaptive_shards as coordinate
import produce_packet_group_g8_highdim_dominance_shards as dominance
import replay_packet_group_g8_adaptive_checkpoint_with_atlas as replay


ROOT = Path(__file__).resolve().parents[1]
NEW_ATLAS = ROOT / "out" / "g8_adaptive_cumulative_supplementary_diverse16.json"
NEW_ATLAS_SHA256 = "530f042f0f4be3fd64ad166c68380302dbb5abb745ffb3745a4fd4ce88b1236c"
DEFAULT_OUTPUT = ROOT / "out" / "g8_full_support_adaptive_independent_component_replay.json"
SCHEMA = "permute-conv.packet-group-g8-independent-component-replay.v1"
INNER_EXPONENT = 9 * (1 << 21) // 100


@dataclass(frozen=True)
class ComponentBank:
    inner_constants: np.ndarray
    inner_charges: np.ndarray
    outer_constants: np.ndarray
    outer_charges: np.ndarray
    references: tuple[dict[str, Any], ...]
    sources: tuple[dict[str, Any], ...]


def row_components(row: dict[str, Any]) -> tuple[float, np.ndarray, float, np.ndarray]:
    inner = row.get("inner", {})
    outer = row.get("outer", {})
    fugacities = np.asarray(inner.get("fugacities"), dtype=np.float64)
    log_variables = np.asarray(outer.get("log_variables"), dtype=np.float64)
    if (
        fugacities.shape != (9,)
        or log_variables.shape != (9,)
        or np.any(fugacities <= 0.0)
        or not np.all(np.isfinite(fugacities))
        or not np.all(np.isfinite(log_variables))
    ):
        raise ValueError("witness row has malformed inner or outer component data")
    inner_charge = np.log2(fugacities)
    outer_charge = log_variables / math.log(2.0)
    if "constant_log2" in inner:
        inner_constant = float(inner["constant_log2"])
    else:
        inner_constant = float(inner["inner_mgf_log2"]) - INNER_EXPONENT * math.log2(
            float(inner["pole"])
        )
    if "constant_log2" in outer:
        outer_constant = float(outer["constant_log2"])
    else:
        profile = np.asarray(row["profile"], dtype=np.float64)
        outer_constant = float(outer["outer_log2"]) + float(profile @ outer_charge)
    combined = row["affine"]
    combined_constant = float(combined["constant_log2"])
    combined_charge = np.asarray(combined["charge_log2"], dtype=np.float64)
    if (
        not math.isfinite(inner_constant)
        or not math.isfinite(outer_constant)
        or not math.isclose(
            inner_constant + outer_constant,
            combined_constant,
            rel_tol=0.0,
            abs_tol=5e-8,
        )
        or not np.allclose(
            inner_charge + outer_charge,
            combined_charge,
            rtol=0.0,
            atol=5e-12,
        )
    ):
        raise ValueError("inner/outer decomposition does not reproduce combined affine row")
    return inner_constant, inner_charge, outer_constant, outer_charge


def load_component_bank(
    manifest: Path, base: Path, old: Path, new: Path
) -> tuple[ComponentBank, dict[str, Any]]:
    if coordinate.sha256_path(new) != NEW_ATLAS_SHA256:
        raise ValueError("new16 atlas digest changed")
    (
        _manifest,
        base_artifact,
        base_source,
        base_bindings,
        manifest_digest,
        base_digest,
    ) = coordinate.load_frozen_inputs(manifest, base)
    old_loaded = dominance.load_supplementary_atlas(old, manifest_digest, base_digest)
    new_loaded = dominance.load_supplementary_atlas(new, manifest_digest, base_digest)
    artifacts = [base_artifact, json.loads(old.read_bytes()), json.loads(new.read_bytes())]
    references = []
    sources = []
    base_refs = coordinate.build_atlas_bank(base_artifact, base_source, base_bindings).references
    references.extend(base_refs)
    sources.append(
        {
            "source_id": base_source["source_id"],
            "path": str(base),
            "sha256": base_digest,
            "row_count": len(base_artifact["rows"]),
            "manifest_bound": True,
        }
    )
    for loaded, path in ((old_loaded, old), (new_loaded, new)):
        references.extend(loaded[2])
        source = dict(loaded[3])
        source["path"] = str(path)
        source["manifest_bound"] = False
        sources.append(source)
    rows = [row for artifact in artifacts for row in artifact["rows"]]
    if len(rows) != 558 or len(references) != 558:
        raise ValueError("component bank does not contain base510+old32+new16")
    components = [row_components(row) for row in rows]
    return ComponentBank(
        inner_constants=np.asarray([row[0] for row in components]),
        inner_charges=np.asarray([row[1] for row in components]),
        outer_constants=np.asarray([row[2] for row in components]),
        outer_charges=np.asarray([row[3] for row in components]),
        references=tuple(references),
        sources=tuple(sources),
    ), {
        "manifest_sha256": manifest_digest,
        "base_atlas_sha256": base_digest,
        "old_supplementary_sha256": coordinate.sha256_path(old),
        "new_supplementary_sha256": coordinate.sha256_path(new),
    }


def rationalize(raw: np.ndarray, denominator: int) -> tuple[tuple[int, Fraction], ...]:
    indices = tuple(range(len(raw)))
    return minimax.rationalize_weights(raw, indices, denominator)


def independent_minimax(
    inner_values: np.ndarray,
    outer_values: np.ndarray,
    normalization: np.ndarray,
    denominator: int,
) -> dict[str, Any]:
    vertex_count, inner_count = inner_values.shape
    if outer_values.shape[0] != vertex_count or normalization.shape != (vertex_count,):
        raise ValueError("component matrix dimensions do not agree")
    outer_count = outer_values.shape[1]
    objective = np.concatenate((np.zeros(inner_count + outer_count), [1.0]))
    inequalities = np.hstack(
        (inner_values, outer_values, -np.ones((vertex_count, 1)))
    )
    equalities = np.zeros((2, inner_count + outer_count + 1))
    equalities[0, :inner_count] = 1.0
    equalities[1, inner_count : inner_count + outer_count] = 1.0
    result = linprog(
        objective,
        A_ub=inequalities,
        b_ub=normalization,
        A_eq=equalities,
        b_eq=np.ones(2),
        bounds=[(0.0, None)] * (inner_count + outer_count) + [(None, None)],
        method="highs-ds",
    )
    if not result.success:
        raise RuntimeError(f"independent component LP failed: {result.message}")
    inner_rational = rationalize(result.x[:inner_count], denominator)
    outer_rational = rationalize(
        result.x[inner_count : inner_count + outer_count], denominator
    )
    vertex_values = []
    for vertex in range(vertex_count):
        inner_value = math.fsum(
            float(weight) * float(inner_values[vertex, index])
            for index, weight in inner_rational
        )
        outer_value = math.fsum(
            float(weight) * float(outer_values[vertex, index])
            for index, weight in outer_rational
        )
        vertex_values.append(inner_value + outer_value - float(normalization[vertex]))
    binary64 = (
        inner_values @ result.x[:inner_count]
        + outer_values @ result.x[inner_count : inner_count + outer_count]
        - normalization
    )
    return {
        "inner": inner_rational,
        "outer": outer_rational,
        "score": max(vertex_values),
        "vertex_values": tuple(vertex_values),
        "binary64_lp_score": float(np.max(binary64)),
        "rationalization_loss_bits": max(vertex_values) - float(np.max(binary64)),
    }


def component_selector(
    rows: tuple[tuple[int, Fraction], ...], references: tuple[dict[str, Any], ...]
) -> dict[str, Any]:
    return {
        "components": [
            {"weight": f"{weight.numerator}/{weight.denominator}", "witness": references[index]}
            for index, weight in rows
        ]
    }


def replay_leaf(
    node: dict[str, Any], bank: ComponentBank, denominator: int
) -> dict[str, Any]:
    profiles = replay.validate_leaf_payload(node)
    matrix = np.asarray(profiles, dtype=np.float64)
    normalization = coordinate.diagnostic_normalization(matrix, 262144)
    inner_values = bank.inner_constants[None, :] - matrix @ bank.inner_charges.T
    outer_values = bank.outer_constants[None, :] - matrix @ bank.outer_charges.T
    mixed = independent_minimax(inner_values, outer_values, normalization, denominator)
    score = float(mixed["score"])
    count = int(node["exact_count"])
    contribution = math.log2(count) + score
    worst = int(np.argmax(np.asarray(mixed["vertex_values"])))
    return {
        "node_id": node["node_id"],
        "h2_cell": int(node["h2_cell"]),
        "depth": int(node["depth"]),
        "exact_count": node["exact_count"],
        "exact_vertex_count": len(profiles),
        "exact_vertex_sha256": node["exact_vertex_sha256"],
        "old_candidate_upper_log2": float(node["candidate_upper_log2"]),
        "candidate_upper_log2": score,
        "upper_improvement_bits": float(node["candidate_upper_log2"]) - score,
        "old_contribution_log2": float(node["contribution_log2"]),
        "contribution_log2": contribution,
        "contribution_improvement_bits": float(node["contribution_log2"]) - contribution,
        "worst_vertex": list(profiles[worst]),
        "inner_selector": component_selector(mixed["inner"], bank.references),
        "outer_selector": component_selector(mixed["outer"], bank.references),
        "inner_support_size": len(mixed["inner"]),
        "outer_support_size": len(mixed["outer"]),
        "diagnostic": {
            "binary64_lp_upper_log2": mixed["binary64_lp_score"],
            "rationalization_loss_bits": mixed["rationalization_loss_bits"],
        },
    }


def build(args: argparse.Namespace) -> dict[str, Any]:
    started = time.perf_counter()
    checkpoint, checkpoint_digest = replay.load_checkpoint(args.checkpoint)
    bank, bindings = load_component_bank(
        args.manifest, args.atlas, args.old_supplementary, args.new_supplementary
    )
    checkpoint_sources = {
        source["sha256"] for source in checkpoint["source_bindings"]["witness_sources"]
    }
    if checkpoint_sources != {bindings["base_atlas_sha256"], bindings["old_supplementary_sha256"]}:
        raise ValueError("checkpoint source bindings differ from base+old32")
    leaves = replay.active_leaves(checkpoint)
    if sum(int(node["exact_count"]) for node in leaves) != math.comb(262143, 8):
        raise ValueError("active leaf counts do not recover the exact census")
    rows = [replay_leaf(node, bank, args.mixture_denominator) for node in leaves]
    old_aggregate = replay.log2_sum([float(node["contribution_log2"]) for node in leaves])
    aggregate = replay.log2_sum([row["contribution_log2"] for row in rows])
    worst = max(rows, key=lambda row: (row["contribution_log2"], row["node_id"]))
    inner_sizes = [row["inner_support_size"] for row in rows]
    outer_sizes = [row["outer_support_size"] for row in rows]
    return {
        "schema": SCHEMA,
        "status": "DIAGNOSTIC_INDEPENDENT_COMPONENT_MIXTURE_REQUIRES_OUTWARD_REPLAY",
        "checkpoint_binding": {
            "path": str(args.checkpoint),
            "sha256": checkpoint_digest,
            "schema": checkpoint["schema"],
            "completed_waves": int(checkpoint["completed_waves"]),
        },
        **bindings,
        "witness_sources": list(bank.sources),
        "component_model": {
            "inner_mixture_mass": "1/1",
            "outer_mixture_mass": "1/1",
            "shared_profile_normalization_subtracted_once": True,
            "binary64_lp": "min max_v sum_i alpha_i I_i(v)+sum_j beta_j O_j(v)-N(v)",
            "rational_denominator": args.mixture_denominator,
        },
        "leaves": rows,
        "summary": {
            "active_leaf_count": len(rows),
            "eligible_inner_count": len(bank.references),
            "eligible_outer_count": len(bank.references),
            "old_aggregate_log2_union_diagnostic": old_aggregate,
            "aggregate_log2_union_diagnostic": aggregate,
            "aggregate_improvement_bits": old_aggregate - aggregate,
            "worst_leaf_id": worst["node_id"],
            "worst_contribution_log2": worst["contribution_log2"],
            "worst_improvement_bits": (
                float(checkpoint["global_diagnostic"]["worst_contribution_log2"])
                - worst["contribution_log2"]
            ),
            "inner_support_size_minimum": min(inner_sizes),
            "inner_support_size_maximum": max(inner_sizes),
            "inner_support_size_mean": math.fsum(inner_sizes) / len(inner_sizes),
            "outer_support_size_minimum": min(outer_sizes),
            "outer_support_size_maximum": max(outer_sizes),
            "outer_support_size_mean": math.fsum(outer_sizes) / len(outer_sizes),
            "runtime_seconds": time.perf_counter() - started,
        },
        "scope_limit": (
            "Fixed stored vertices and counts; no geometry or split changes. Binary64 "
            "component-mixture discovery requires independent outward replay."
        ),
    }


def explicit_cross_product(
    inner: np.ndarray, outer: np.ndarray, normalization: np.ndarray, denominator: int
) -> dict[str, Any]:
    combined = np.column_stack(
        [inner[:, i] + outer[:, j] - normalization for i in range(inner.shape[1]) for j in range(outer.shape[1])]
    )
    return minimax.optimize_minimax_mixture(
        combined,
        0.0,
        seed_per_vertex=1,
        reduced_cost_tolerance=1e-10,
        weight_tolerance=1e-12,
        denominator=denominator,
    )


def run_self_test() -> None:
    inner = np.asarray([[0.0, 2.0], [2.0, 0.0], [1.0, 1.0]])
    outer = np.asarray([[3.0, 0.0], [0.0, 3.0], [1.0, 1.0]])
    normalization = np.asarray([0.0, 0.0, 0.0])
    independent = independent_minimax(inner, outer, normalization, 1 << 20)
    explicit = explicit_cross_product(inner, outer, normalization, 1 << 20)
    if (
        abs(independent["binary64_lp_score"] - explicit["binary64_lp_score"]) > 1e-10
        or abs(independent["score"] - explicit["score"]) > 1e-5
        or sum((weight for _index, weight in independent["inner"]), Fraction()) != 1
        or sum((weight for _index, weight in independent["outer"]), Fraction()) != 1
    ):
        raise SystemExit("independent component LP cross-product self-test failed")
    print(f"independent_lp_score={independent['binary64_lp_score']:.12f}")
    print(f"explicit_cross_product_score={explicit['binary64_lp_score']:.12f}")
    print("separate_rational_mass_constraints=PASS")
    print("status=PASS_G8_INDEPENDENT_COMPONENT_MIXTURE_SELF_TEST")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=h2.DEFAULT_MANIFEST)
    parser.add_argument("--atlas", type=Path, default=h2.DEFAULT_ATLAS)
    parser.add_argument("--old-supplementary", type=Path, default=h2.DEFAULT_SUPPLEMENTARY)
    parser.add_argument("--new-supplementary", type=Path, default=NEW_ATLAS)
    parser.add_argument("--checkpoint", type=Path, default=replay.DEFAULT_CHECKPOINT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--mixture-denominator", type=int, default=1 << 30)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        run_self_test()
        return
    if args.mixture_denominator < 1:
        parser.error("mixture denominator must be positive")
    report = build(args)
    digest = coordinate.write_canonical(args.output, report)
    summary = report["summary"]
    print(f"active_leaves={summary['active_leaf_count']}")
    print(f"aggregate_log2={summary['aggregate_log2_union_diagnostic']:.12f}")
    print(f"aggregate_improvement_bits={summary['aggregate_improvement_bits']:.12f}")
    print(f"worst_improvement_bits={summary['worst_improvement_bits']:.12f}")
    print(
        f"inner_support={summary['inner_support_size_minimum']}:"
        f"{summary['inner_support_size_maximum']} "
        f"outer_support={summary['outer_support_size_minimum']}:"
        f"{summary['outer_support_size_maximum']}"
    )
    print(f"runtime_seconds={summary['runtime_seconds']:.3f}")
    print(f"output={args.output}")
    print(f"sha256={digest}")
    print("status=DIAGNOSTIC_INDEPENDENT_COMPONENT_MIXTURE_REQUIRES_OUTWARD_REPLAY")


if __name__ == "__main__":
    main()
