#!/usr/bin/env python3
"""Plan direct root certificates for g=8 support dimensions zero through five.

The output is a binary64 discovery plan, not a proof artifact.  For every
selected exact support, the script reconstructs the exact rational vertices
of the shifted simplex with the physical-weight cut.  It then selects the
frozen atlas witness with the smallest diagnostic maximum over those vertices.
An independent outward evaluator must replay that selector before a shard can
claim ``CERTIFIED_LEAF`` or ``COMPLETE``.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import linprog
from scipy.special import gammaln


GROUP_BITS = 8
CLASSES = GROUP_BITS + 1
M = 262144
MINIMUM_PROFILE_WEIGHT = 21
FEASIBLE_PROFILE_COUNT = 553169839211945865258921061892182603726
UNIFORM_TARGET_LOG2 = -40.0 - math.log2(FEASIBLE_PROFILE_COUNT)
EXPECTED_SUPPORT_COUNTS = (8, 36, 84, 126, 126, 84, 36, 9, 1)
EXPECTED_PROFILE_COUNTS = (
    8,
    9437096,
    2886184992415,
    378293710105607109,
    24791478291771024604728,
    866508443723541949201514342,
    16224627887200131391413327495204,
    151895545730963601519675551057510391,
    553017927440720481221689864611271442433,
)
ATLAS_SCHEMA = "permute-conv.packet-group-g8-support-seed-atlas.v1"
ATLAS_SHA256 = "c48644fa62ad74619a4db30b07950b14aaf1419d040760f7d659f43288fca45e"
MANIFEST_SCHEMA = "packet-group-g8-support-manifest-v1"
MANIFEST_SHA256 = "cc446d0c6a0c986a8712ef78c4aa7c9d6073687665b7e0f315164339a5b81616"
SOURCE_ID = "atlas-3934dae"
ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ATLAS = ROOT / "out" / "g8_support_seed_atlas.json"
DEFAULT_MANIFEST = ROOT / "G8_SUPPORT_MANIFEST.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def support_mask(support: tuple[int, ...]) -> int:
    return sum(1 << index for index in support)


def excluded_below_weight(support: tuple[int, ...]) -> int:
    """Count exact-support profiles below physical weight 21."""

    if 0 not in support:
        return 0
    positive = tuple(index for index in support if index)
    threshold = MINIMUM_PROFILE_WEIGHT - sum(positive)
    if threshold <= 0:
        return 0
    counts = [0] * threshold
    counts[0] = 1
    for weight in positive:
        for total in range(weight, threshold):
            counts[total] += counts[total - weight]
    return sum(counts)


def exact_root_count(support: tuple[int, ...]) -> int:
    if support == (0,):
        return 0
    return math.comb(M - 1, len(support) - 1) - excluded_below_weight(support)


def root_vertices(support: tuple[int, ...]) -> tuple[tuple[Fraction, ...], ...]:
    """Return exact vertices of the continuous root proof polytope.

    Active coordinates have lower bound one, inactive coordinates are zero,
    all coordinates sum to M, and physical weight is at least 21.
    """

    size = len(support)
    residual = M - size
    base = [Fraction(1 if index in support else 0) for index in range(CLASSES)]

    def pure(index: int) -> tuple[Fraction, ...]:
        point = base.copy()
        point[index] += residual
        return tuple(point)

    if 0 not in support:
        return tuple(pure(index) for index in support)
    positive = tuple(index for index in support if index)
    if not positive:
        return ()
    missing_weight = MINIMUM_PROFILE_WEIGHT - sum(positive)
    if missing_weight <= 0:
        return tuple(pure(index) for index in support)
    vertices = [pure(index) for index in positive]
    for index in positive:
        positive_mass = Fraction(missing_weight, index)
        if not 0 <= positive_mass <= residual:
            raise RuntimeError(f"physical cut misses edge for support {support}")
        point = base.copy()
        point[index] += positive_mass
        point[0] += residual - positive_mass
        vertices.append(tuple(point))
    return tuple(vertices)


def normalization_log2(vertices: tuple[tuple[Fraction, ...], ...]) -> np.ndarray:
    rows = np.asarray([[float(value) for value in vertex] for vertex in vertices])
    classes = np.asarray([math.comb(GROUP_BITS, index) for index in range(CLASSES)])
    return (
        gammaln(M + 1)
        - np.sum(gammaln(rows + 1), axis=1)
        + rows @ np.log(classes)
    ) / math.log(2.0)


def load_atlas(path: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    actual = sha256(path)
    if actual != ATLAS_SHA256:
        raise ValueError(f"frozen atlas digest mismatch: {actual}")
    artifact = json.loads(path.read_text(encoding="utf-8"))
    if artifact.get("schema") != ATLAS_SCHEMA:
        raise ValueError("frozen atlas schema mismatch")
    rows = artifact.get("rows")
    if not isinstance(rows, list) or len(rows) != 510:
        raise ValueError("frozen atlas must contain 510 rows")
    ordinals = [int(row["ordinal"]) for row in rows]
    if ordinals != list(range(510)):
        raise ValueError("frozen atlas ordinals are not canonical")
    for row in rows:
        if len(row.get("affine", {}).get("charge_log2", [])) != CLASSES:
            raise ValueError(f"atlas row {row['ordinal']} has malformed charge")
        if any(float(value) <= 0.0 for value in row["inner"]["fugacities"]):
            raise ValueError(f"atlas row {row['ordinal']} has a zero fugacity")
    return artifact, rows


def load_manifest(path: Path) -> tuple[dict[str, Any], dict[int, str], dict[int, dict]]:
    actual = sha256(path)
    if actual != MANIFEST_SHA256:
        raise ValueError(f"frozen manifest digest mismatch: {actual}")
    manifest = json.loads(path.read_text(encoding="utf-8"))
    if manifest.get("schema") != MANIFEST_SCHEMA:
        raise ValueError("frozen manifest schema mismatch")
    sources = [
        source
        for source in manifest.get("witness_sources", [])
        if source.get("source_id") == SOURCE_ID
    ]
    if len(sources) != 1 or sources[0].get("sha256") != ATLAS_SHA256:
        raise ValueError("manifest does not bind the frozen atlas")
    row_hashes = {
        int(row["row"]): str(row["row_sha256"]) for row in sources[0]["rows"]
    }
    if len(row_hashes) != 510:
        raise ValueError("manifest does not bind all atlas rows")
    census = {
        int(row["mask"], 16): row for row in manifest["census"]["feasible_masks"]
    }
    if len(census) != 510:
        raise ValueError("manifest does not contain all feasible support masks")
    return manifest, row_hashes, census


def candidate_rows(
    atlas_rows: list[dict[str, Any]], support: tuple[int, ...], scope: str
) -> list[dict[str, Any]]:
    if scope == "local":
        selected = [row for row in atlas_rows if tuple(row["support"]) == support]
    elif scope == "same-dimension":
        selected = [row for row in atlas_rows if len(row["support"]) == len(support)]
    else:
        selected = atlas_rows
    if not selected:
        raise RuntimeError(f"candidate scope {scope!r} is empty for {support}")
    return selected


def witness_value_matrix(
    vertices: tuple[tuple[Fraction, ...], ...], candidates: list[dict[str, Any]]
) -> np.ndarray:
    profiles = np.asarray(
        [[float(value) for value in vertex] for vertex in vertices], dtype=np.float64
    )
    normalizations = normalization_log2(vertices)
    constants = np.asarray(
        [float(row["affine"]["constant_log2"]) for row in candidates]
    )
    charges = np.asarray(
        [row["affine"]["charge_log2"] for row in candidates], dtype=np.float64
    )
    values = constants[:, None] - charges @ profiles.T - normalizations[None, :]
    return values.T


def best_singleton(values: np.ndarray) -> dict[str, Any]:
    maxima = np.max(values, axis=0)
    index = int(np.argmin(maxima))
    return {
        "index": index,
        "score": float(maxima[index]),
        "vertex_values": values[:, index].tolist(),
    }


def seed_columns(values: np.ndarray, per_vertex: int) -> np.ndarray:
    finite = np.flatnonzero(np.all(np.isfinite(values), axis=0))
    if not len(finite):
        raise RuntimeError("root has no finite eligible witness")
    selected: set[int] = set()
    for vertex in range(values.shape[0]):
        order = finite[np.lexsort((finite, values[vertex, finite]))]
        selected.update(map(int, order[:per_vertex]))
    return np.asarray(sorted(selected), dtype=np.int64)


def rationalize_weights(
    raw_weights: np.ndarray, indices: tuple[int, ...], denominator: int
) -> tuple[tuple[int, Fraction], ...]:
    if denominator < 1 or len(raw_weights) != len(indices):
        raise ValueError("invalid mixture rationalization input")
    clipped = np.maximum(np.asarray(raw_weights, dtype=np.float64), 0.0)
    mass = float(np.sum(clipped))
    if not math.isfinite(mass) or mass <= 0.0:
        raise RuntimeError("mixture has nonpositive binary64 mass")
    normalized = clipped / mass
    scaled = normalized * denominator
    units = np.floor(scaled).astype(object)
    missing = denominator - sum(map(int, units))
    order = sorted(
        range(len(indices)),
        key=lambda local: (-(scaled[local] - int(units[local])), indices[local]),
    )
    for local in order[:missing]:
        units[local] = int(units[local]) + 1
    result = tuple(
        (index, Fraction(int(unit), denominator))
        for index, unit in zip(indices, units)
        if int(unit) > 0
    )
    if sum((weight for _index, weight in result), Fraction(0)) != 1:
        raise RuntimeError("rationalized mixture weights do not sum to one")
    return result


def optimize_minimax_mixture(
    values: np.ndarray,
    target: float,
    *,
    seed_per_vertex: int,
    reduced_cost_tolerance: float,
    weight_tolerance: float,
    denominator: int,
) -> dict[str, Any]:
    """Solve the finite vertex game by safe full-scan column generation."""

    eligible = np.flatnonzero(np.all(np.isfinite(values), axis=0))
    if not len(eligible):
        raise RuntimeError("root has no finite eligible witness")
    selected = seed_columns(values, seed_per_vertex)

    def solve(columns: np.ndarray):
        shifted = values[:, columns] - target
        count = len(columns)
        result = linprog(
            np.append(np.zeros(count), 1.0),
            A_ub=np.hstack((shifted, -np.ones((values.shape[0], 1)))),
            b_ub=np.zeros(values.shape[0]),
            A_eq=np.asarray([[1.0] * count + [0.0]]),
            b_eq=np.asarray([1.0]),
            bounds=[(0.0, None)] * count + [(None, None)],
            method="highs-ds",
        )
        if not result.success:
            raise RuntimeError(f"lowdim minimax mixture LP failed: {result.message}")
        return result

    rounds = 0
    full_scans = 0
    result = solve(selected)
    minimum_reduced_cost = None
    maximum_pricing_gap = None
    shifted_all = values[:, eligible] - target
    while True:
        dual = -np.asarray(result.ineqlin.marginals, dtype=np.float64)
        if (
            len(dual) != values.shape[0]
            or np.min(dual) < -1e-10
            or abs(float(np.sum(dual)) - 1.0) > 1e-7
        ):
            raise RuntimeError("minimax mixture LP returned an invalid vertex dual")
        restricted_prices = dual @ (values[:, selected] - target)
        restricted_gap = float(result.fun) - float(np.min(restricted_prices))
        if abs(restricted_gap) > 1e-7 * max(1.0, abs(float(result.fun))):
            raise RuntimeError("minimax mixture LP failed restricted primal-dual KKT")
        prices = dual @ shifted_all
        reduced = prices - float(result.fun)
        full_scans += 1
        minimum_reduced_cost = float(np.min(reduced))
        maximum_pricing_gap = float(result.fun) - float(np.min(prices))
        selected_set = set(map(int, selected))
        violating = [
            (float(reduced[local]), int(index))
            for local, index in enumerate(eligible)
            if int(index) not in selected_set
            and float(reduced[local]) < -reduced_cost_tolerance
        ]
        if not violating:
            break
        violating.sort(key=lambda item: (item[0], item[1]))
        additions = [index for _cost, index in violating[: values.shape[0]]]
        selected = np.asarray(sorted(selected_set.union(additions)), dtype=np.int64)
        result = solve(selected)
        rounds += 1

    active_local = np.flatnonzero(result.x[: len(selected)] > weight_tolerance)
    if not len(active_local):
        raise RuntimeError("minimax mixture LP returned empty support")
    active_indices = tuple(int(selected[local]) for local in active_local)
    rational = rationalize_weights(
        result.x[active_local], active_indices, denominator
    )
    mixed = tuple(
        math.fsum(float(weight) * float(values[vertex, index]) for index, weight in rational)
        for vertex in range(values.shape[0])
    )
    binary_mixed = tuple(
        math.fsum(
            float(result.x[local]) * float(values[vertex, selected[local]])
            for local in range(len(selected))
        )
        for vertex in range(values.shape[0])
    )
    return {
        "indices": tuple(index for index, _weight in rational),
        "weights": tuple(weight for _index, weight in rational),
        "score": max(mixed),
        "vertex_values": mixed,
        "binary64_lp_score": max(binary_mixed),
        "rationalization_loss_bits": max(mixed) - max(binary_mixed),
        "candidate_count": int(len(selected)),
        "column_generation_rounds": rounds,
        "full_column_scans": full_scans,
        "minimum_reduced_cost": minimum_reduced_cost,
        "maximum_pricing_gap": maximum_pricing_gap,
        "reduced_cost_tolerance": reduced_cost_tolerance,
        "sparse_basic_solution": len(rational) <= values.shape[0],
        "rational_denominator": denominator,
    }


def witness_reference(row: dict[str, Any], row_hashes: dict[int, str]) -> dict[str, Any]:
    ordinal = int(row["ordinal"])
    return {
        "source_id": SOURCE_ID,
        "source_sha256": ATLAS_SHA256,
        "row": ordinal,
        "row_sha256": row_hashes[ordinal],
    }


def mixture_selector(
    result: dict[str, Any], candidates: list[dict[str, Any]], row_hashes: dict[int, str]
) -> dict[str, Any]:
    if len(result["indices"]) == 1:
        return {
            "kind": "witness",
            "witness": witness_reference(candidates[result["indices"][0]], row_hashes),
        }
    return {
        "kind": "mixture",
        "components": [
            {
                "weight": fraction_text(weight),
                "witness": witness_reference(candidates[index], row_hashes),
            }
            for index, weight in zip(result["indices"], result["weights"])
        ],
    }


def run_self_test() -> None:
    # Neither singleton passes target zero: both maxima equal one.  Their
    # one-half mixture equals minus one at both vertices and does pass.
    values = np.asarray([[1.0, -3.0], [-3.0, 1.0]], dtype=np.float64)
    singleton = best_singleton(values)
    mixture = optimize_minimax_mixture(
        values,
        0.0,
        seed_per_vertex=1,
        reduced_cost_tolerance=1e-10,
        weight_tolerance=1e-12,
        denominator=1024,
    )
    if not (
        singleton["score"] > 0.0
        and mixture["score"] < 0.0
        and mixture["weights"] == (Fraction(1, 2), Fraction(1, 2))
        and sum(mixture["weights"], Fraction(0)) == 1
        and mixture["minimum_reduced_cost"] >= -1e-10
        and mixture["maximum_pricing_gap"] <= 1e-10
    ):
        raise SystemExit("lowdim mixture self-test failed")
    selector = mixture_selector(
        mixture,
        [{"ordinal": 0}, {"ordinal": 1}],
        {0: "0" * 64, 1: "1" * 64},
    )
    if (
        selector.get("kind") != "mixture"
        or [row.get("weight") for row in selector.get("components", [])]
        != ["1/2", "1/2"]
    ):
        raise SystemExit("lowdim mixture selector self-test failed")
    # The balanced third column is not a per-vertex winner.  A safe pricing
    # scan must still add it and improve the two-winner mixture from 50 to 49.
    hidden = np.asarray([[0.0, 100.0, 49.0], [100.0, 0.0, 49.0]])
    priced = optimize_minimax_mixture(
        hidden,
        -100.0,
        seed_per_vertex=1,
        reduced_cost_tolerance=1e-10,
        weight_tolerance=1e-12,
        denominator=1024,
    )
    if not (
        abs(priced["score"] - 49.0) < 1e-10
        and priced["indices"] == (2,)
        and priced["column_generation_rounds"] >= 1
        and priced["minimum_reduced_cost"] >= -1e-10
        and priced["maximum_pricing_gap"] <= 1e-10
    ):
        raise SystemExit("lowdim mixture pricing self-test failed")
    print("synthetic_singleton_score=1.0")
    print("synthetic_rational_mixture=1/2,1/2")
    print("synthetic_mixture_score=-1.0")
    print("synthetic_hidden_column_score=49.0")
    print("status=PASS_LOWDIM_MINIMAX_MIXTURE_SELF_TEST")


def validate_census() -> None:
    support_counts = [0] * CLASSES
    profile_counts = [0] * CLASSES
    for size in range(1, CLASSES + 1):
        for support in itertools.combinations(range(CLASSES), size):
            count = exact_root_count(support)
            if count:
                support_counts[size - 1] += 1
                profile_counts[size - 1] += count
    if tuple(support_counts) != EXPECTED_SUPPORT_COUNTS:
        raise RuntimeError(f"support census changed: {support_counts}")
    if tuple(profile_counts) != EXPECTED_PROFILE_COUNTS:
        raise RuntimeError(f"profile census changed: {profile_counts}")
    if sum(profile_counts) != FEASIBLE_PROFILE_COUNT:
        raise RuntimeError("profile census total changed")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", type=Path, default=DEFAULT_ATLAS)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--min-dimension", type=int, default=0)
    parser.add_argument("--max-dimension", type=int, default=2)
    parser.add_argument(
        "--candidate-scope",
        choices=("local", "same-dimension", "all"),
        default="all",
    )
    parser.add_argument(
        "--selector-mode",
        choices=("singleton", "minimax-mixture"),
        default="minimax-mixture",
    )
    parser.add_argument("--mixture-seed-per-vertex", type=int, default=1)
    parser.add_argument("--mixture-denominator", type=int, default=1 << 30)
    parser.add_argument("--reduced-cost-tolerance", type=float, default=1e-9)
    parser.add_argument("--weight-tolerance", type=float, default=1e-12)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.self_test:
        run_self_test()
        return
    if not 0 <= args.min_dimension <= args.max_dimension <= 5:
        parser.error("require 0 <= min-dimension <= max-dimension <= 5")
    if args.output is None:
        parser.error("--output is required unless --self-test is selected")
    if (
        args.mixture_seed_per_vertex < 1
        or args.mixture_denominator < 1
        or args.reduced_cost_tolerance < 0.0
        or args.weight_tolerance < 0.0
    ):
        parser.error("invalid mixture discovery limits")

    validate_census()
    manifest, row_hashes, manifest_census = load_manifest(args.manifest)
    atlas, atlas_rows = load_atlas(args.atlas)
    plans = []
    for size in range(args.min_dimension + 1, args.max_dimension + 2):
        for support in itertools.combinations(range(CLASSES), size):
            count = exact_root_count(support)
            if not count:
                continue
            mask = support_mask(support)
            census_row = manifest_census.get(mask)
            if (
                census_row is None
                or census_row["active_classes"] != list(support)
                or int(census_row["dimension"]) != size - 1
                or int(census_row["count"]) != count
            ):
                raise RuntimeError(f"manifest census mismatch for support {support}")
            vertices = root_vertices(support)
            candidates = candidate_rows(atlas_rows, support, args.candidate_scope)
            values = witness_value_matrix(vertices, candidates)
            singleton = best_singleton(values)
            candidate = candidates[singleton["index"]]
            singleton_selector = {
                "kind": "witness",
                "witness": witness_reference(candidate, row_hashes),
            }
            mixture = optimize_minimax_mixture(
                values,
                UNIFORM_TARGET_LOG2,
                seed_per_vertex=args.mixture_seed_per_vertex,
                reduced_cost_tolerance=args.reduced_cost_tolerance,
                weight_tolerance=args.weight_tolerance,
                denominator=args.mixture_denominator,
            )
            fixed_mixture_selector = mixture_selector(mixture, candidates, row_hashes)
            if args.selector_mode == "minimax-mixture" and mixture["score"] <= singleton["score"]:
                selector = fixed_mixture_selector
                upper = float(mixture["score"])
                vertex_values = list(map(float, mixture["vertex_values"]))
                selected_kind = fixed_mixture_selector["kind"]
            else:
                selector = singleton_selector
                upper = float(singleton["score"])
                vertex_values = list(map(float, singleton["vertex_values"]))
                selected_kind = "witness"
            contribution = math.log2(count) + upper
            count_record = {
                "kind": "exact",
                "method": "root-census-v1",
                "value": str(count),
            }
            count_derivation = {
                "method": "exact-support-positive-compositions-weight-cut-v1",
                "parameters": {
                    "active_classes": list(support),
                    "active_lower_bounds": ["1"] * size,
                    "active_upper_bounds": [str(M)] * size,
                    "total_mass": str(M),
                    "minimum_profile_weight": str(MINIMUM_PROFILE_WEIGHT),
                    "weight_cut_exclusions": str(excluded_below_weight(support)),
                },
            }
            plans.append(
                {
                    "manifest_sha256": MANIFEST_SHA256,
                    "support_mask": f"0x{mask:03x}",
                    "active_classes": list(support),
                    "dimension": size - 1,
                    "root_count": str(count),
                    "root_count_record": count_record,
                    "diagnostic_root_count_derivation": count_derivation,
                    "root_vertices": [
                        [fraction_text(value) for value in vertex]
                        for vertex in vertices
                    ],
                    "candidate_selector": selector,
                    "selected_selector_kind": selected_kind,
                    "diagnostic_best_singleton": {
                        "selector": singleton_selector,
                        "root_upper_log2": float(singleton["score"]),
                        "vertex_values_log2": singleton["vertex_values"],
                        "uniform_target_pass": singleton["score"] <= UNIFORM_TARGET_LOG2,
                    },
                    "diagnostic_minimax_mixture": {
                        "selector": fixed_mixture_selector,
                        "root_upper_log2": float(mixture["score"]),
                        "vertex_values_log2": list(map(float, mixture["vertex_values"])),
                        "binary64_lp_score": float(mixture["binary64_lp_score"]),
                        "rationalization_loss_bits": float(
                            mixture["rationalization_loss_bits"]
                        ),
                        "candidate_count": int(mixture["candidate_count"]),
                        "components": len(mixture["indices"]),
                        "column_generation_rounds": int(
                            mixture["column_generation_rounds"]
                        ),
                        "full_column_scans": int(mixture["full_column_scans"]),
                        "minimum_reduced_cost": mixture["minimum_reduced_cost"],
                        "maximum_pricing_gap": mixture["maximum_pricing_gap"],
                        "reduced_cost_tolerance": mixture["reduced_cost_tolerance"],
                        "sparse_basic_solution": mixture["sparse_basic_solution"],
                        "rational_denominator": mixture["rational_denominator"],
                        "uniform_target_pass": mixture["score"] <= UNIFORM_TARGET_LOG2,
                    },
                    "diagnostic_vertex_values_log2": vertex_values,
                    "diagnostic_root_upper_log2": upper,
                    "diagnostic_contribution_log2": contribution,
                    "diagnostic_uniform_target_pass": upper <= UNIFORM_TARGET_LOG2,
                    "proof_state": "REQUIRES_OUTWARD_REPLAY",
                    "proposed_shard_after_successful_outward_replay": {
                        "schema": "packet-group-g8-support-shard-v1",
                        "manifest_sha256": MANIFEST_SHA256,
                        "support_mask": f"0x{mask:03x}",
                        "active_classes": list(support),
                        "root_count": str(count),
                        "root_node": "r",
                        "nodes": [
                            {
                                "node_id": "r",
                                "state": "CERTIFIED_LEAF",
                                "selector": selector,
                                "count": count_record,
                            }
                        ],
                        "aggregation_requested": "collapsed",
                        "state": "COMPLETE",
                    },
                }
            )

    worst = sorted(plans, key=lambda row: row["diagnostic_root_upper_log2"], reverse=True)
    contribution_order = sorted(
        plans, key=lambda row: row["diagnostic_contribution_log2"], reverse=True
    )
    report = {
        "schema": "permute-conv.packet-group-g8-lowdim-root-cover-plan.v1",
        "status": "DIAGNOSTIC_BINARY64_REQUIRES_OUTWARD_REPLAY",
        "scope": (
            "direct one-leaf root candidates for selected exact supports; no proof, "
            "completeness, or certified coverage claim"
        ),
        "frozen_atlas_commit": "3934dae73836e9051638db2e51b8b66a07055b3c",
        "manifest_commit": "78774ee32efd81b176399ba506f823c10bb99267",
        "manifest": {
            "path": str(args.manifest),
            "sha256": MANIFEST_SHA256,
            "schema": MANIFEST_SCHEMA,
            "run_id": manifest.get("run_id"),
            "canonical_json": manifest.get("arithmetic", {}).get("canonical_json"),
        },
        "atlas": {
            "source_id": SOURCE_ID,
            "path": str(args.atlas),
            "sha256": ATLAS_SHA256,
            "schema": ATLAS_SCHEMA,
            "configuration_digest": atlas.get("configuration_digest"),
        },
        "parameters": {
            "group_bits": GROUP_BITS,
            "M": str(M),
            "minimum_profile_weight": str(MINIMUM_PROFILE_WEIGHT),
            "uniform_profile_target_log2": UNIFORM_TARGET_LOG2,
        },
        "configuration": {
            "min_dimension": args.min_dimension,
            "max_dimension": args.max_dimension,
            "candidate_scope": args.candidate_scope,
            "selector_mode": args.selector_mode,
            "mixture_seed_per_vertex": args.mixture_seed_per_vertex,
            "mixture_denominator": args.mixture_denominator,
            "reduced_cost_tolerance": args.reduced_cost_tolerance,
            "weight_tolerance": args.weight_tolerance,
        },
        "method": {
            "geometry": "shifted-simplex-exact-vertices-with-weight-cut-v1",
            "selection": (
                "full-scan column-generation minimax over exact root vertices, "
                "followed by deterministic exact rational weight rounding"
                if args.selector_mode == "minimax-mixture"
                else "minimum singleton binary64 maximum over exact root vertices"
            ),
            "proof_condition": (
                "independent outward evaluation of the fixed selector at every "
                "reconstructed root vertex"
            ),
            "counting": "manifest-bound root-census-v1",
            "aggregation": "collapsed exact root count times root vertex maximum",
        },
        "supports_planned": len(plans),
        "diagnostic_singleton_uniform_target_passes": sum(
            row["diagnostic_best_singleton"]["uniform_target_pass"] for row in plans
        ),
        "diagnostic_mixture_uniform_target_passes": sum(
            row["diagnostic_minimax_mixture"]["uniform_target_pass"] for row in plans
        ),
        "diagnostic_uniform_target_passes": sum(
            row["diagnostic_uniform_target_pass"] for row in plans
        ),
        "selected_fixed_mixtures": sum(
            row["selected_selector_kind"] == "mixture" for row in plans
        ),
        "worst_root_candidates": worst[:10],
        "largest_contribution_candidates": contribution_order[:10],
        "plans": plans,
        "proof_blockers": [
            "binary64 constants, charges, log-gamma values, and maxima are not outward bounds",
            "a final shard may claim COMPLETE only after independent outward replay",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(
        f"supports={len(plans)} "
        f"singleton_passes={report['diagnostic_singleton_uniform_target_passes']} "
        f"mixture_passes={report['diagnostic_mixture_uniform_target_passes']} "
        f"selected_passes={report['diagnostic_uniform_target_passes']} "
        f"worst={worst[0]['diagnostic_root_upper_log2']:.9f} "
        f"largest_contribution={contribution_order[0]['diagnostic_contribution_log2']:.9f}"
    )
    print(f"output={args.output}")
    print("status=DIAGNOSTIC_BINARY64_REQUIRES_OUTWARD_REPLAY")


if __name__ == "__main__":
    main()
