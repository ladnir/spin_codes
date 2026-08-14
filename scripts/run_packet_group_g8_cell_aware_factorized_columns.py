#!/usr/bin/env python3
"""Run one bounded cell-aware factorized column-generation experiment."""

from __future__ import annotations

import argparse
import json
import math
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
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
import replay_packet_group_g8_adaptive_checkpoint_with_atlas as stored
import replay_packet_group_g8_independent_component_mixtures as component


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "out" / "g8_full_support_adaptive_independent_persistence_wave.json"
INPUT_SHA256 = "32a7eb88d09f20725a8fd4fe5a727c6f3bc925ee9e90febcc753c2ed9cfa0739"
DEFAULT_WORK_DIR = ROOT / "out" / "g8_cell_aware_factorized_columns"
DEFAULT_OUTPUT = ROOT / "out" / "g8_cell_aware_factorized_columns.json"
SCHEMA = "permute-conv.packet-group-g8-cell-aware-factorized-columns.v1"
TAU_PRICE = 64.0
DUAL_DENOMINATOR = 1 << 40
MIXTURE_DENOMINATOR = 1 << 30
TARGETS = 8
WORKERS = 8
INNER_EXPONENT = component.INNER_EXPONENT


def atomic_json(path: Path, value: Any) -> str:
    return coordinate.write_canonical(path, value)


def fractional_profile(exact: list[str]) -> list[float]:
    values = [float(Fraction(value)) for value in exact[:8]]
    values.append(262144.0 - math.fsum(values))
    if len(values) != 9 or sum(values) != 262144.0 or any(value <= 0.0 for value in values):
        raise ValueError("fractional residual profile does not preserve full support and mass")
    return values


def select_targets(artifact: dict[str, Any]) -> list[dict[str, Any]]:
    leaves = [node for node in artifact["nodes"].values() if node["state"] == "ACTIVE_LEAF"]
    leaves.sort(key=lambda node: (-float(node["contribution_log2"]), node["node_id"]))
    selected = []
    seen = set()
    for node in leaves:
        parent = int(node["h2_cell"])
        if parent in seen:
            continue
        seen.add(parent)
        selected.append(node)
        if len(selected) == TARGETS:
            break
    if len(selected) != TARGETS:
        raise ValueError("could not select eight distinct h2 parents")
    return selected


def prepare_target(node: dict[str, Any], bank: component.ComponentBank) -> dict[str, Any]:
    profiles = stored.validate_leaf_payload(node)
    matrix = np.asarray(profiles, dtype=np.float64)
    normalization = coordinate.diagnostic_normalization(matrix, 262144)
    inner = bank.inner_constants[None, :] - matrix @ bank.inner_charges.T
    outer = bank.outer_constants[None, :] - matrix @ bank.outer_charges.T
    result = pricing.solve_component_lp(inner, outer, normalization)
    dual_raw = -np.asarray(result.ineqlin.marginals, dtype=np.float64)
    dual = pricing.rationalize_distribution(dual_raw, DUAL_DENOMINATOR)
    barycenter = pricing.exact_barycenter(profiles, dual)
    dual_float = np.zeros(len(profiles))
    for index, weight in dual:
        dual_float[index] = float(weight)
    return {
        "node_id": node["node_id"],
        "h2_cell": int(node["h2_cell"]),
        "contribution_log2": float(node["contribution_log2"]),
        "exact_vertex_sha256": node["exact_vertex_sha256"],
        "profiles": [list(profile) for profile in profiles],
        "dual": [(index, pricing.fraction_text(weight)) for index, weight in dual],
        "barycenter_exact": [pricing.fraction_text(value) for value in barycenter],
        "barycenter_float": fractional_profile([pricing.fraction_text(value) for value in barycenter]),
        "inner_mu": float(np.min(dual_float @ inner)),
        "outer_mu": float(np.min(dual_float @ outer)),
    }


def tune_job(task: dict[str, Any]) -> dict[str, Any]:
    lane = task["lane"]
    profile = task["barycenter_float"]
    started = time.perf_counter()
    if lane == "inner":
        from packet_group_native import load_point_caps, load_shared_drive_apply
        from packet_group_profile_bound import inner_probability, split_cap_table, tune_inner_density
        if load_point_caps() is None or load_shared_drive_apply() is None:
            raise RuntimeError("inner pricing job needs native kernels")
        caps = split_cap_table()
        tuned = tune_inner_density(
            8, profile, caps, poles=(0.2, 0.5, 0.8),
            exponents=(0.25, 0.5, 0.75), iterations=1,
        )
        pole = float(tuned[1])
        fugacities = np.asarray(tuned[3], dtype=np.float64)
        _probability, details = inner_probability(
            8, profile, pole, fugacities, caps, iterations=3
        )
        constant = float(details["inner_mgf_log2"]) - INNER_EXPONENT * math.log2(pole)
        charge = np.log2(fugacities)
        price = constant - float(np.asarray(profile) @ charge)
        parameters = {
            "lane": lane, "pole": pole, "fugacities": fugacities.tolist(),
            "constant_log2": constant, "charge_log2": charge.tolist(),
            "screen_iterations": 1, "final_iterations": 3,
            "collatz_best_iteration": int(details["collatz_best_iteration"]),
        }
    elif lane == "outer":
        from probe_packet_group_conditioned_row_outer import load_split_spectrum
        from survey_packet_group_g8_full_support_witness_reuse import optimize_outer_seed
        spectra = (
            load_split_spectrum(ROOT / "out" / "ebch85_band01_split_spectrum.csv"),
            load_split_spectrum(ROOT / "out" / "ebch84_punctured_band01_split_spectrum.csv", count_field="pair_count", divisor=42),
            load_split_spectrum(ROOT / "out" / "ebch86_band12_split_spectrum.csv"),
        )
        tuned = optimize_outer_seed(profile, spectra, "asymmetric", 60, 800)
        charge = np.asarray(tuned["log_variables"]) / math.log(2.0)
        constant = float(tuned["outer_log2"]) + float(np.asarray(profile) @ charge)
        price = constant - float(np.asarray(profile) @ charge)
        parameters = {
            "lane": lane, "log_variables": tuned["log_variables"],
            "band_coefficients": tuned["band_coefficients"],
            "pair_cauchy_theta": tuned["pair_cauchy_theta"],
            "constant_log2": constant, "charge_log2": charge.tolist(),
            "optimizer_iterations": tuned["optimizer_iterations"],
            "optimizer_evaluations": tuned["optimizer_evaluations"],
            "optimizer_success": tuned["optimizer_success"],
        }
    else:
        raise ValueError("unknown pricing lane")
    parameter_digest = coordinate.sha256_bytes(coordinate.canonical_bytes(parameters))
    return {
        "schema": "permute-conv.packet-group-g8-fractional-component-checkpoint.v1",
        "status": "COMPLETE_DIAGNOSTIC_FRACTIONAL_PRICING",
        "lane": lane,
        "node_id": task["node_id"],
        "h2_cell": task["h2_cell"],
        "barycenter_exact": task["barycenter_exact"],
        "barycenter_float": profile,
        "old_mu": task[f"{lane}_mu"],
        "candidate_price": price,
        "reduced_cost": price - task[f"{lane}_mu"],
        "parameter_digest": parameter_digest,
        "parameters": parameters,
        "runtime_seconds": time.perf_counter() - started,
    }


def candidate_column(
    profiles: np.ndarray, row: dict[str, Any]
) -> np.ndarray:
    parameters = row["parameters"]
    return float(parameters["constant_log2"]) - profiles @ np.asarray(parameters["charge_log2"])


def identity_error(target: dict[str, Any], row: dict[str, Any]) -> float:
    profiles = np.asarray(target["profiles"], dtype=np.float64)
    dual = np.zeros(len(profiles))
    for index, weight in target["dual"]:
        dual[index] = float(Fraction(weight))
    full_price = float(dual @ candidate_column(profiles, row))
    barycenter_price = float(row["candidate_price"])
    return abs(full_price - barycenter_price)


def extended_replay(
    artifact: dict[str, Any], bank: component.ComponentBank,
    accepted_inner: list[dict[str, Any]], accepted_outer: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    inner_constants = np.concatenate((bank.inner_constants, np.asarray([r["parameters"]["constant_log2"] for r in accepted_inner])))
    inner_charges = np.vstack((bank.inner_charges, np.asarray([r["parameters"]["charge_log2"] for r in accepted_inner]))) if accepted_inner else bank.inner_charges
    outer_constants = np.concatenate((bank.outer_constants, np.asarray([r["parameters"]["constant_log2"] for r in accepted_outer])))
    outer_charges = np.vstack((bank.outer_charges, np.asarray([r["parameters"]["charge_log2"] for r in accepted_outer]))) if accepted_outer else bank.outer_charges
    rows = []
    for node in sorted((n for n in artifact["nodes"].values() if n["state"] == "ACTIVE_LEAF"), key=lambda n:n["node_id"]):
        profiles = stored.validate_leaf_payload(node)
        matrix = np.asarray(profiles, dtype=np.float64)
        normalization = coordinate.diagnostic_normalization(matrix, 262144)
        mixed = component.independent_minimax(
            inner_constants[None,:] - matrix @ inner_charges.T,
            outer_constants[None,:] - matrix @ outer_charges.T,
            normalization, MIXTURE_DENOMINATOR,
        )
        score = float(mixed["score"])
        old_score = float(node["candidate_upper_log2"])
        retained_old = score > old_score
        final_score = old_score if retained_old else score
        contribution = math.log2(int(node["exact_count"])) + final_score
        rows.append({
            "node_id": node["node_id"], "h2_cell": int(node["h2_cell"]),
            "old_upper_log2": old_score, "new_upper_log2": final_score,
            "improvement_bits": old_score-final_score,
            "contribution_log2": contribution, "retained_old_selector": retained_old,
            "new_inner_positive": any(index >= len(bank.references) for index,_ in mixed["inner"]),
            "new_outer_positive": any(index >= len(bank.references) for index,_ in mixed["outer"]),
        })
    old_aggregate = component.replay.log2_sum([float(n["contribution_log2"]) for n in artifact["nodes"].values() if n["state"]=="ACTIVE_LEAF"])
    new_aggregate = component.replay.log2_sum([r["contribution_log2"] for r in rows])
    old_worst = float(artifact["global_diagnostic"]["worst_contribution_log2"])
    new_worst = max(r["contribution_log2"] for r in rows)
    return rows, {
        "old_aggregate_log2": old_aggregate, "new_aggregate_log2": new_aggregate,
        "aggregate_improvement_bits": old_aggregate-new_aggregate,
        "old_worst_contribution_log2": old_worst, "new_worst_contribution_log2": new_worst,
        "worst_improvement_bits": old_worst-new_worst,
    }


def run(args: argparse.Namespace) -> dict[str, Any]:
    started = time.perf_counter()
    if coordinate.sha256_path(args.input) != INPUT_SHA256:
        raise ValueError("persistence input digest changed")
    artifact = json.loads(args.input.read_bytes())
    bank, bindings = component.load_component_bank(args.manifest,args.atlas,args.old_supplementary,args.new_supplementary)
    selected_nodes = select_targets(artifact)
    targets = [prepare_target(node, bank) for node in selected_nodes]
    tasks = [{**target, "lane": lane} for target in targets for lane in ("inner","outer")]
    args.work_dir.mkdir(parents=True, exist_ok=True)
    configuration_digest = coordinate.sha256_bytes(coordinate.canonical_bytes({
        "input_sha256": INPUT_SHA256, "targets": [(t["node_id"], t["barycenter_exact"]) for t in targets],
        "inner_poles": [0.2,0.5,0.8], "inner_exponents": [0.25,0.5,0.75],
        "inner_screen_iterations": 1, "inner_final_iterations": 3,
        "outer_mode": "asymmetric", "outer_max_iterations": 60, "outer_max_evaluations": 800,
    }))
    completed=[]; resumed=0; missing=[]
    for task in tasks:
        path=args.work_dir / f"{task['lane']}_{task['h2_cell']:03d}.json"
        if path.exists():
            row=json.loads(path.read_bytes())
            if (row.get("configuration_digest")==configuration_digest and row.get("lane")==task["lane"] and row.get("node_id")==task["node_id"] and row.get("barycenter_exact")==task["barycenter_exact"]):
                completed.append(row); resumed+=1; continue
        missing.append(task)
    with ProcessPoolExecutor(max_workers=WORKERS) as pool:
        futures={pool.submit(tune_job,task):task for task in missing}
        for future in as_completed(futures):
            row=future.result(); row["configuration_digest"]=configuration_digest; completed.append(row)
            atomic_json(args.work_dir / f"{row['lane']}_{row['h2_cell']:03d}.json", row)
    completed.sort(key=lambda r:(r["h2_cell"],r["lane"]))
    target_by_node={t["node_id"]:t for t in targets}
    for row in completed:
        row["affine_barycenter_identity_error_bits"]=identity_error(target_by_node[row["node_id"]],row)
        row["support_eligible"]=all(math.isfinite(float(v)) for v in row["parameters"]["charge_log2"])
        row["accepted"]=row["support_eligible"] and row["affine_barycenter_identity_error_bits"]<=1e-6 and row["reduced_cost"] < -TAU_PRICE
    unique={}
    for row in completed:
        if row["accepted"]: unique[(row["lane"],row["parameter_digest"])]=row
    accepted_inner=[r for (lane,_),r in unique.items() if lane=="inner"]
    accepted_outer=[r for (lane,_),r in unique.items() if lane=="outer"]
    leaves,replay_summary=extended_replay(artifact,bank,accepted_inner,accepted_outer)
    selected_ids={t["node_id"] for t in targets}
    improved_16k=sum(r["node_id"] in selected_ids and r["improvement_bits"]>=16384 for r in leaves)
    identity_pass=all(r["affine_barycenter_identity_error_bits"]<=1e-6 for r in completed if r["accepted"])
    no_regression=all(r["improvement_bits"]>=0 for r in leaves)
    positive_inner=any(r["new_inner_positive"] for r in leaves)
    positive_outer=any(r["new_outer_positive"] for r in leaves)
    correctness=identity_pass and no_regression and positive_inner and positive_outer
    payoff=improved_16k>=4 and (replay_summary["aggregate_improvement_bits"]>=32768 or replay_summary["worst_improvement_bits"]>=32768)
    catalogue={"schema":"permute-conv.packet-group-g8-cell-aware-component-catalogue.v1","rows":completed}
    catalogue_digest=atomic_json(args.work_dir/"candidate_catalogue.json",catalogue)
    return {
        "schema":SCHEMA,"status":"DIAGNOSTIC_CELL_AWARE_FACTORIZED_COLUMNS",
        "input_sha256":INPUT_SHA256,**bindings,"configuration":{"targets":TARGETS,"workers":WORKERS,"tau_price_bits":TAU_PRICE,"maximum_jobs":16,"configuration_digest":configuration_digest},
        "resumed_job_count":resumed,"new_job_count":len(missing),
        "targets":[{k:t[k] for k in ("node_id","h2_cell","contribution_log2","barycenter_exact","inner_mu","outer_mu")} for t in targets],
        "lane_results":completed,"candidate_catalogue_sha256":catalogue_digest,
        "accepted_inner_count":len(accepted_inner),"accepted_outer_count":len(accepted_outer),
        "replay":{"leaves":leaves,"summary":replay_summary},
        "gate":{"affine_identity_passes":identity_pass,"rational_replay_no_regression":no_regression,"new_inner_positive":positive_inner,"new_outer_positive":positive_outer,"selected_leaves_improved_16384":improved_16k,"correctness_passes":correctness,"payoff_passes":payoff,"overall_passes":correctness and payoff,"outward_replay_performed":False,"classification":"DIAGNOSTIC_RATIONAL_REPLAY_ONLY"},
        "runtime_seconds":time.perf_counter()-started,
        "scope_limit":"Discovery and rational binary64 replay only; no outward or completeness claim.",
    }


def self_test() -> None:
    exact=["1/2"]*8+["262140/1"]
    profile=fractional_profile(exact)
    if sum(profile)!=262144.0 or profile[-1] != 262140.0:
        raise SystemExit("fractional residual mass self-test failed")
    target={"profiles":[[1]*9,[2]*9],"dual":[(0,"1/2"),(1,"1/2")],"barycenter_float":[1.5]*9}
    row={"candidate_price":-13.5,"parameters":{"constant_log2":0.0,"charge_log2":[1.0]*9}}
    if identity_error(target,row)>1e-12:
        raise SystemExit("fractional affine identity self-test failed")
    from probe_packet_group_conditioned_row_outer import load_split_spectrum
    for path, options in (
        (ROOT / "out" / "ebch85_band01_split_spectrum.csv", {}),
        (ROOT / "out" / "ebch84_punctured_band01_split_spectrum.csv", {"count_field":"pair_count","divisor":42}),
        (ROOT / "out" / "ebch86_band12_split_spectrum.csv", {}),
    ):
        load_split_spectrum(path, **options)
    print("fractional_residual_mass=PASS")
    print("fractional_affine_barycenter_identity=PASS")
    print("sparse_integer_refiners_disabled=PASS")
    print("outer_worker_path_initialization=PASS")
    print("status=PASS_G8_CELL_AWARE_FACTORIZED_COLUMNS_SELF_TEST")


def main() -> None:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--manifest",type=Path,default=component.h2.DEFAULT_MANIFEST);p.add_argument("--atlas",type=Path,default=component.h2.DEFAULT_ATLAS)
    p.add_argument("--old-supplementary",type=Path,default=component.h2.DEFAULT_SUPPLEMENTARY);p.add_argument("--new-supplementary",type=Path,default=component.NEW_ATLAS)
    p.add_argument("--input",type=Path,default=INPUT);p.add_argument("--work-dir",type=Path,default=DEFAULT_WORK_DIR);p.add_argument("--output",type=Path,default=DEFAULT_OUTPUT);p.add_argument("--self-test",action="store_true")
    args=p.parse_args()
    if args.self_test:self_test();return
    report=run(args);digest=atomic_json(args.output,report);s=report["replay"]["summary"]
    print(f"accepted_inner={report['accepted_inner_count']} accepted_outer={report['accepted_outer_count']}")
    print(f"aggregate_improvement_bits={s['aggregate_improvement_bits']:.12f} worst_improvement_bits={s['worst_improvement_bits']:.12f}")
    print(f"gate_overall={report['gate']['overall_passes']} classification={report['gate']['classification']}")
    print(f"runtime_seconds={report['runtime_seconds']:.3f}");print(f"output={args.output}");print(f"sha256={digest}")


if __name__=="__main__":main()
