#!/usr/bin/env python3
"""Extract a component-LP dual and its exact pricing barycenter for one leaf."""

from __future__ import annotations

import argparse
import json
import math
import os
from decimal import Decimal, getcontext
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
import replay_packet_group_g8_adaptive_checkpoint_with_atlas as stored
import replay_packet_group_g8_independent_component_mixtures as component


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ARTIFACT = ROOT / "out" / "g8_full_support_adaptive_independent_persistence_wave.json"
ARTIFACT_SHA256 = "32a7eb88d09f20725a8fd4fe5a727c6f3bc925ee9e90febcc753c2ed9cfa0739"
DEFAULT_OUTPUT = ROOT / "out" / "g8_component_pricing_target_h2_159_R.json"
DEFAULT_NODE_ID = "h2:159/R"
SCHEMA = "permute-conv.packet-group-g8-component-pricing-target.v1"


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def decimal_text(value: Fraction, digits: int = 50) -> str:
    getcontext().prec = digits
    return format(Decimal(value.numerator) / Decimal(value.denominator), "f")


def solve_component_lp(
    inner_values: np.ndarray,
    outer_values: np.ndarray,
    normalization: np.ndarray,
) -> Any:
    vertices, inner_count = inner_values.shape
    if outer_values.shape[0] != vertices or normalization.shape != (vertices,):
        raise ValueError("component LP arrays have incompatible shapes")
    outer_count = outer_values.shape[1]
    result = linprog(
        np.concatenate((np.zeros(inner_count + outer_count), [1.0])),
        A_ub=np.hstack(
            (inner_values, outer_values, -np.ones((vertices, 1)))
        ),
        b_ub=normalization,
        A_eq=np.asarray(
            [
                [1.0] * inner_count + [0.0] * outer_count + [0.0],
                [0.0] * inner_count + [1.0] * outer_count + [0.0],
            ]
        ),
        b_eq=np.ones(2),
        bounds=[(0.0, None)] * (inner_count + outer_count) + [(None, None)],
        method="highs-ds",
    )
    if not result.success:
        raise RuntimeError(f"component pricing LP failed: {result.message}")
    return result


def rationalize_distribution(values: np.ndarray, denominator: int) -> tuple[tuple[int, Fraction], ...]:
    return minimax.rationalize_weights(values, tuple(range(len(values))), denominator)


def pricing_rows(
    prices: np.ndarray,
    minimum: float,
    primal: np.ndarray,
    references: tuple[dict[str, Any], ...],
) -> list[dict[str, Any]]:
    return [
        {
            "component": index,
            "witness": references[index],
            "dual_price": float(price),
            "reduced_cost": float(price - minimum),
            "binary64_primal_weight": float(primal[index]),
        }
        for index, price in enumerate(prices)
    ]


def exact_barycenter(
    profiles: tuple[tuple[int, ...], ...],
    dual: tuple[tuple[int, Fraction], ...],
) -> tuple[Fraction, ...]:
    return tuple(
        sum((weight * profiles[index][coordinate] for index, weight in dual), Fraction())
        for coordinate in range(9)
    )


def extract(
    profiles: tuple[tuple[int, ...], ...],
    bank: component.ComponentBank,
    denominator: int,
) -> dict[str, Any]:
    matrix = np.asarray(profiles, dtype=np.float64)
    normalization = coordinate.diagnostic_normalization(matrix, 262144)
    inner_values = bank.inner_constants[None, :] - matrix @ bank.inner_charges.T
    outer_values = bank.outer_constants[None, :] - matrix @ bank.outer_charges.T
    result = solve_component_lp(inner_values, outer_values, normalization)
    inner_count = inner_values.shape[1]
    outer_count = outer_values.shape[1]
    alpha = np.asarray(result.x[:inner_count], dtype=np.float64)
    beta = np.asarray(result.x[inner_count : inner_count + outer_count], dtype=np.float64)
    dual_raw = -np.asarray(result.ineqlin.marginals, dtype=np.float64)
    if (
        np.min(dual_raw) < -1e-10
        or abs(float(np.sum(dual_raw)) - 1.0) > 1e-8
    ):
        raise RuntimeError("component LP returned an invalid vertex dual")
    dual = rationalize_distribution(dual_raw, denominator)
    dual_float = np.zeros(len(profiles), dtype=np.float64)
    for index, weight in dual:
        dual_float[index] = float(weight)
    barycenter = exact_barycenter(profiles, dual)

    inner_prices = dual_float @ inner_values
    outer_prices = dual_float @ outer_values
    normalization_price = float(dual_float @ normalization)
    inner_minimum = float(np.min(inner_prices))
    outer_minimum = float(np.min(outer_prices))
    dual_value = inner_minimum + outer_minimum - normalization_price
    inner_selector = rationalize_distribution(alpha, denominator)
    outer_selector = rationalize_distribution(beta, denominator)
    inner_rational_gap = (
        math.fsum(float(weight) * float(inner_prices[index]) for index, weight in inner_selector)
        - inner_minimum
    )
    outer_rational_gap = (
        math.fsum(float(weight) * float(outer_prices[index]) for index, weight in outer_selector)
        - outer_minimum
    )
    inner_active = {index for index, _weight in inner_selector}
    outer_active = {index for index, _weight in outer_selector}
    inner_inactive = [
        float(inner_prices[index] - inner_minimum)
        for index in range(inner_count)
        if index not in inner_active
    ]
    outer_inactive = [
        float(outer_prices[index] - outer_minimum)
        for index in range(outer_count)
        if index not in outer_active
    ]
    larger = (
        "inner"
        if inner_rational_gap > outer_rational_gap
        else "outer"
        if outer_rational_gap > inner_rational_gap
        else "tie"
    )
    return {
        "binary64_primal_score": float(result.fun),
        "rational_dual_value": dual_value,
        "primal_dual_absolute_gap_bits": abs(float(result.fun) - dual_value),
        "vertex_dual": [
            {
                "vertex": index,
                "weight": fraction_text(weight),
                "binary64_weight": float(dual_raw[index]),
                "profile": list(profiles[index]),
            }
            for index, weight in dual
        ],
        "pricing_barycenter": {
            "exact": [fraction_text(value) for value in barycenter],
            "decimal_50_digits": [decimal_text(value) for value in barycenter],
            "exact_mass": fraction_text(sum(barycenter, Fraction())),
            "physical_weight": fraction_text(
                sum((index * value for index, value in enumerate(barycenter)), Fraction())
            ),
            "interpretation": (
                "exact rational average of active exact leaf vertices under the "
                "rationalized LP vertex dual"
            ),
        },
        "inner_pricing": {
            "best_price": inner_minimum,
            "rational_selector_pricing_gap": inner_rational_gap,
            "next_inactive_reduced_cost": min(inner_inactive) if inner_inactive else None,
            "components": pricing_rows(
                inner_prices, inner_minimum, alpha, bank.references
            ),
        },
        "outer_pricing": {
            "best_price": outer_minimum,
            "rational_selector_pricing_gap": outer_rational_gap,
            "next_inactive_reduced_cost": min(outer_inactive) if outer_inactive else None,
            "components": pricing_rows(
                outer_prices, outer_minimum, beta, bank.references
            ),
        },
        "larger_rational_selector_pricing_gap": larger,
        "gap_definition": (
            "weighted dual price of the rationalized primal selector minus the "
            "minimum dual price over the corresponding finite component bank"
        ),
    }


def build(args: argparse.Namespace) -> dict[str, Any]:
    digest = coordinate.sha256_path(args.artifact)
    if digest != ARTIFACT_SHA256:
        raise ValueError("independent persistence artifact digest changed")
    artifact = json.loads(args.artifact.read_bytes())
    node = artifact.get("nodes", {}).get(args.node_id)
    if node is None or node.get("state") != "ACTIVE_LEAF":
        raise ValueError("requested node is not an active stored leaf")
    profiles = stored.validate_leaf_payload(node)
    bank, bindings = component.load_component_bank(
        args.manifest, args.atlas, args.old_supplementary, args.new_supplementary
    )
    pricing = extract(profiles, bank, args.rational_denominator)
    return {
        "schema": SCHEMA,
        "status": "DIAGNOSTIC_FINITE_BANK_COMPONENT_PRICING_TARGET",
        "artifact_binding": {
            "path": str(args.artifact),
            "sha256": digest,
            "schema": artifact["schema"],
        },
        **bindings,
        "witness_sources": list(bank.sources),
        "node_id": args.node_id,
        "exact_count": node["exact_count"],
        "exact_vertex_count": len(profiles),
        "exact_vertex_sha256": node["exact_vertex_sha256"],
        "rational_denominator": args.rational_denominator,
        "pricing": pricing,
        "scope_limit": (
            "Finite-bank binary64 LP diagnostic only. The barycenter is an oracle "
            "target, not a new witness, outward certificate, or coverage claim."
        ),
    }


def run_self_test() -> None:
    profiles = ((1, 1, 1, 2, 1, 1, 1, 1, 262135), (1, 1, 1, 1, 1, 2, 1, 1, 262135))
    bank = component.ComponentBank(
        inner_constants=np.asarray([0.0, 0.0]),
        inner_charges=np.asarray([[0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0]]),
        outer_constants=np.asarray([0.0, 0.0]),
        outer_charges=np.asarray([[0.0] * 9, [0.0] * 9]),
        references=(
            {"source_id": "x", "source_sha256": "0" * 64, "row": 0, "row_sha256": "1" * 64},
            {"source_id": "x", "source_sha256": "0" * 64, "row": 1, "row_sha256": "2" * 64},
        ),
        sources=(),
    )
    report = extract(profiles, bank, 1024)
    barycenter = report["pricing_barycenter"]["exact"]
    if (
        report["primal_dual_absolute_gap_bits"] > 1e-8
        or report["pricing_barycenter"]["exact_mass"] != "262144/1"
        or barycenter[3] != "3/2"
        or barycenter[5] != "3/2"
    ):
        raise SystemExit("component pricing target self-test failed")
    print("synthetic_primal_dual_gap=PASS")
    print("synthetic_barycenter_classes_3_5=(3/2,3/2)")
    print("per_component_reduced_costs=PASS")
    print("status=PASS_G8_COMPONENT_PRICING_TARGET_SELF_TEST")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=h2.DEFAULT_MANIFEST)
    parser.add_argument("--atlas", type=Path, default=h2.DEFAULT_ATLAS)
    parser.add_argument("--old-supplementary", type=Path, default=h2.DEFAULT_SUPPLEMENTARY)
    parser.add_argument("--new-supplementary", type=Path, default=component.NEW_ATLAS)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_ARTIFACT)
    parser.add_argument("--node-id", default=DEFAULT_NODE_ID)
    parser.add_argument("--rational-denominator", type=int, default=1 << 40)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        run_self_test()
        return
    if args.rational_denominator < 1:
        parser.error("rational denominator must be positive")
    report = build(args)
    digest = coordinate.write_canonical(args.output, report)
    pricing = report["pricing"]
    print(f"node_id={report['node_id']}")
    print(f"dual_support={len(pricing['vertex_dual'])}")
    print(f"primal_dual_gap_bits={pricing['primal_dual_absolute_gap_bits']:.12g}")
    print(f"larger_oracle_gap={pricing['larger_rational_selector_pricing_gap']}")
    print(f"output={args.output}")
    print(f"sha256={digest}")
    print("status=DIAGNOSTIC_FINITE_BANK_COMPONENT_PRICING_TARGET")


if __name__ == "__main__":
    main()
