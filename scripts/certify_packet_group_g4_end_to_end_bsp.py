#!/usr/bin/env python3
"""Certify the complete ``g=4`` support mesh with exact BSP replacements.

The lower-support ledger supplies supports 1 through 30.  The full-support
ledger supplies support 31.  A batch BSP artifact may replace selected root
cells of support 31.  The verifier audits the combined exact geometry, replays
every BSP partition over exact rationals, outward-hardens every used witness,
and performs one cell-local union bound over all positive supports.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from decimal import Decimal
from fractions import Fraction
from pathlib import Path
from typing import Any

import certify_packet_group_g4_anchor_mesh as g4
import certify_packet_group_g4_cell_bsp as bsp
import certify_packet_group_triangle_ledger as g2
from outward_log2 import Interval, log2_int, self_check


BATCH_SCHEMA = "packet-group-g4-cell-dominance-bsp-batch-v1"
STATUS = "OUTWARD_CERTIFIED_G4_END_TO_END_BSP_2^-40"


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compact_hash(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"certificate input {path} is not a JSON object")
    return value


def combine_support_ledgers(
    lower: dict[str, Any], full: dict[str, Any]
) -> dict[str, Any]:
    for label, ledger in (("lower", lower), ("full", full)):
        if ledger.get("schema") != g4.SUPPORT_SCHEMA:
            raise ValueError(f"{label} ledger has an unsupported schema")
        if int(ledger.get("group_bits", -1)) != 4:
            raise ValueError(f"{label} ledger has the wrong group size")
    lower_by_mask = {
        int(row.get("support_mask", -1)): row for row in lower.get("support_strata", [])
    }
    full_by_mask = {
        int(row.get("support_mask", -1)): row for row in full.get("support_strata", [])
    }
    if set(lower_by_mask) != set(range(1, 32)) or set(full_by_mask) != set(range(1, 32)):
        raise ValueError("both ledgers must carry placeholder records for supports 1..31")
    combined = dict(full)
    combined["complete_support_stratified_cover"] = True
    combined["support_strata"] = [
        lower_by_mask[mask] if mask < 31 else full_by_mask[mask]
        for mask in range(1, 32)
    ]
    return combined


def resolved_witness_paths(
    inputs: tuple[tuple[Path, dict[str, Any]], ...]
) -> list[Path]:
    by_basename: dict[str, tuple[str, Path]] = {}
    for artifact_path, artifact in inputs:
        for path in g4.source_paths(artifact_path, artifact):
            digest = file_sha256(path)
            previous = by_basename.get(path.name)
            if previous is not None and previous[0] != digest:
                raise ValueError(
                    f"witness basename {path.name!r} resolves to different contents"
                )
            by_basename[path.name] = (digest, path)
    return [by_basename[name][1] for name in sorted(by_basename)]


def maximum_interval(left: Interval | None, right: Interval) -> Interval:
    if left is None:
        return right
    return Interval(max(left.lo, right.lo), max(left.hi, right.hi))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lower-ledger", type=Path, required=True)
    parser.add_argument("--full-ledger", type=Path, required=True)
    parser.add_argument("--bsp-batch", type=Path, action="append", required=True)
    parser.add_argument("--iterations", type=int, default=40)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--checkpoint-dir", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if not 1 <= args.workers <= 8:
        parser.error("--workers must lie in 1..8")
    if args.iterations < 1:
        parser.error("--iterations must be positive")
    self_check()

    lower = load_json(args.lower_ledger)
    full = load_json(args.full_ledger)
    batches = [(path, load_json(path)) for path in args.bsp_batch]
    for _path, batch in batches:
        if batch.get("schema") != BATCH_SCHEMA or int(batch.get("group_bits", -1)) != 4:
            raise ValueError("unsupported g4 BSP batch schema")
        source_ledger = batch.get("source_ledger")
        if not isinstance(source_ledger, dict):
            raise ValueError("BSP batch lacks its source-ledger record")
        source_digest = str(source_ledger.get("sha256", ""))
        if source_digest not in {
            file_sha256(args.lower_ledger), file_sha256(args.full_ledger)
        }:
            raise ValueError("BSP batch was not built from a supplied support ledger")

    combined = combine_support_ledgers(lower, full)
    support_mode = g4.support_stratum_contexts(combined)
    if support_mode is None:
        raise ValueError("combined ledger did not enter support-stratified mode")
    contexts, geometry = support_mode
    by_identifier = {context["identifier"]: context for context in contexts}
    if len(by_identifier) != len(contexts):
        raise ValueError("combined root-cell identifiers are not unique")

    replacements: dict[str, dict[str, Any]] = {}
    replacement_leaves: dict[str, tuple[list[dict[str, Any]], dict[str, Any]]] = {}
    for _batch_path, batch in batches:
        raw_cells = batch.get("cells")
        if not isinstance(raw_cells, list) or not raw_cells:
            raise ValueError("BSP batch has no replacement cells")
        batch_masks = batch.get("support_mask")
        if isinstance(batch_masks, list):
            allowed_masks = {int(value) for value in batch_masks}
        else:
            allowed_masks = {int(batch_masks)}
        before = len(replacements)
        for cell in raw_cells:
            if not isinstance(cell, dict):
                raise ValueError("BSP cell records must be objects")
            support_mask = int(cell.get("support_mask", next(iter(allowed_masks))))
            if support_mask not in allowed_masks:
                raise ValueError("BSP cell support mask is outside its batch declaration")
            raw_identifier = str(cell.get("cell_id", ""))
            identifier = f"s{support_mask:02x}:{raw_identifier}"
            if not raw_identifier or identifier in replacements:
                raise ValueError("BSP cell identifiers are empty or duplicated")
            context = by_identifier.get(identifier)
            if context is None:
                raise ValueError(f"BSP replacement {identifier!r} is not a root cell")
            if not bool(cell.get("complete_exact_partition_claimed", False)):
                raise ValueError(f"BSP replacement {identifier!r} is not claimed complete")
            row = context["row"]
            anchors = context["anchors"]
            indices = g4.stratum_cell_indices(row, len(anchors), context["dimension"])
            if list(indices) != list(map(int, cell.get("root_anchor_indices", []))):
                raise ValueError(f"BSP replacement {identifier!r} has wrong root indices")
            root_vertices = tuple(anchors[index] for index in indices)
            bsp.configure_support(support_mask)
            if root_vertices != bsp.parse_serialized_vertices(cell.get("root_anchor_profiles")):
                raise ValueError(f"BSP replacement {identifier!r} has wrong root profiles")
            root_constraints = bsp.simplex_constraints(root_vertices)
            if bsp.enumerate_vertices(root_constraints) != tuple(sorted(root_vertices)):
                raise ValueError(f"BSP replacement {identifier!r} root reconstruction failed")
            local_count, _omitted, _ranges = g4.simplex_integer_profile_count_upper(
                indices, anchors
            )
            if local_count != int(row.get("integer_profile_count_upper", -1)):
                raise ValueError(f"source root {identifier!r} has the wrong count")
            if local_count != int(cell.get("root_integer_profile_count_upper", -1)):
                raise ValueError(f"BSP replacement {identifier!r} has the wrong count")
            leaves, vertices = bsp.replay_tree(cell, root_constraints)
            replacements[identifier] = cell
            replacement_leaves[identifier] = (leaves, vertices)
        if int(batch.get("diagnostic", {}).get("processed_cells", -1)) != len(replacements) - before:
            raise ValueError("BSP batch processed-cell count mismatch")

    mixtures: dict[str, tuple[tuple[str, Fraction], ...]] = {}
    used: set[str] = set()
    for context in contexts:
        identifier = context["identifier"]
        if identifier in replacements:
            leaves, _vertices = replacement_leaves[identifier]
            used.update(str(leaf["owner_witness"]) for leaf in leaves)
            continue
        mixture = g4.parse_mixture(context["row"])
        mixtures[identifier] = mixture
        used.update(name for name, _weight in mixture)

    paths = resolved_witness_paths(
        ((args.lower_ledger, lower), (args.full_ledger, full), *batches)
    )
    witnesses = g2.load_witnesses(paths)
    missing = sorted(name for name in used if name != "full_bijection" and name not in witnesses)
    if missing:
        raise ValueError("missing witness references: " + ", ".join(missing))

    g4.configure_generic_hardener()
    selected = sorted(name for name in used if name != "full_bijection")
    hardened = g2.harden_selected_witnesses(
        selected,
        witnesses,
        args.iterations,
        args.workers,
        args.checkpoint_dir,
        group_bits=4,
    )
    if "full_bijection" in used:
        hardened["full_bijection"] = g4.harden_full_bijection()

    for context in contexts:
        identifier = context["identifier"]
        if identifier in replacements:
            leaves, _vertices = replacement_leaves[identifier]
            for leaf in leaves:
                owner = str(leaf["owner_witness"])
                g4.require_support_eligible(
                    identifier, context["support"], ((owner, Fraction(1)),), hardened
                )
        else:
            g4.require_support_eligible(
                identifier, context["support"], mixtures[identifier], hardened
            )

    inequality_digest = hashlib.sha256()
    terms: list[Interval] = []
    top_terms: list[tuple[Decimal, dict[str, Any]]] = []
    component_evaluations = 0
    root_cells = 0
    bsp_leaves = 0
    bsp_vertex_inequalities = 0
    maximum: Interval | None = None
    worst: tuple[str, str] | None = None

    for context in contexts:
        identifier = context["identifier"]
        row = context["row"]
        anchors = context["anchors"]
        indices = g4.stratum_cell_indices(row, len(anchors), context["dimension"])
        local_count, omitted, ranges = g4.simplex_integer_profile_count_upper(indices, anchors)
        if local_count != int(row.get("integer_profile_count_upper", local_count)):
            raise ValueError(f"cell {identifier!r} integer-profile count mismatch")
        cell_maximum: Interval | None = None
        owner_description = "fixed_mixture"
        local_evaluations = 0
        if identifier in replacements:
            owner_description = "exact_bsp_partition"
            leaves, vertices_by_leaf = replacement_leaves[identifier]
            bsp_leaves += len(leaves)
            for leaf in sorted(leaves, key=lambda value: str(value["node_id"])):
                leaf_id = str(leaf["node_id"])
                owner = str(leaf["owner_witness"])
                mixture = ((owner, Fraction(1)),)
                for vertex_index, profile in enumerate(vertices_by_leaf[leaf_id]):
                    value, components = g4.evaluate_profile(profile, mixture, hardened)
                    component_evaluations += len(components)
                    local_evaluations += 1
                    bsp_vertex_inequalities += 1
                    inequality_digest.update(
                        json.dumps(
                            [identifier, leaf_id, vertex_index, owner,
                             bsp.serialized_vertices((profile,))[0],
                             str(value.lo), str(value.hi)],
                            separators=(",", ":"),
                        ).encode()
                    )
                    cell_maximum = maximum_interval(cell_maximum, value)
                    if maximum is None or value.hi > maximum.hi:
                        maximum = value
                        worst = (identifier, f"{leaf_id}:{vertex_index}")
        else:
            root_cells += 1
            mixture = mixtures[identifier]
            for vertex in indices:
                value, components = g4.evaluate_profile(anchors[vertex], mixture, hardened)
                component_evaluations += len(components)
                local_evaluations += 1
                inequality_digest.update(
                    json.dumps(
                        [identifier, vertex,
                         [[name, str(weight), str(part.lo), str(part.hi)]
                          for name, weight, part in components],
                         str(value.lo), str(value.hi)],
                        separators=(",", ":"),
                    ).encode()
                )
                cell_maximum = maximum_interval(cell_maximum, value)
                if maximum is None or value.hi > maximum.hi:
                    maximum = value
                    worst = (identifier, str(vertex))
        if cell_maximum is None:
            raise ValueError(f"cell {identifier!r} has no verified vertices")
        if local_count <= 0:
            continue
        term = cell_maximum + log2_int(local_count)
        terms.append(term)
        summary = {
            "cell_id": identifier,
            "owner_kind": owner_description,
            "integer_profile_count_upper": local_count,
            "omitted_coordinate": omitted,
            "coordinate_ranges": [list(pair) for pair in ranges],
            "evaluated_vertices": local_evaluations,
            "maximum_branch_log2_interval": [str(cell_maximum.lo), str(cell_maximum.hi)],
            "union_term_log2_interval": [str(term.lo), str(term.hi)],
        }
        top_terms.append((term.hi, summary))

    if maximum is None or worst is None or not terms:
        raise ValueError("complete union evaluation produced no terms")
    # The accelerated inner hardener deliberately carries a very coarse
    # finite lower endpoint.  The proof comparison uses only the upper
    # endpoint, so evaluate that endpoint from exact point intervals.  The
    # maximum of the valid term lower endpoints remains a valid lower bound
    # on the logarithm of their sum.
    union_upper = g2._log2_sum_exp(Interval.exact(term.hi) for term in terms)
    union = Interval(max(term.lo for term in terms), union_upper.hi)
    margin = Decimal(-40) - union.hi
    passed = margin >= 0
    top_terms.sort(key=lambda item: item[0], reverse=True)
    source_digests = {
        "lower_ledger": file_sha256(args.lower_ledger),
        "full_ledger": file_sha256(args.full_ledger),
        "bsp_batches": {path.name: file_sha256(path) for path, _batch in batches},
        "witness_sources": {path.name: file_sha256(path) for path in paths},
    }
    witness_reports_digest = compact_hash(
        [hardened[name]["report"] for name in sorted(hardened)]
    )
    report = {
        "status": STATUS if passed else "OUTWARD_G4_END_TO_END_BSP_DID_NOT_CLOSE",
        "passed": passed,
        "group_bits": 4,
        "scope": "all feasible positive-support profiles for the frozen g=4 construction",
        "source_digests": source_digests,
        "exact_geometry": geometry,
        "profile_owner_rule": g4.SUPPORT_OWNER_RULE,
        "bsp_partition_rule": (
            "each selected support-mesh root is replaced by its exact rational BSP; "
            "the root count upper multiplies the maximum bound over all BSP leaves"
        ),
        "root_cells_without_bsp": root_cells,
        "bsp_replaced_root_cells": len(replacements),
        "bsp_leaves": bsp_leaves,
        "used_witnesses": len(used),
        "used_component_vertex_evaluations": component_evaluations,
        "used_bsp_leaf_vertex_inequalities": bsp_vertex_inequalities,
        "inequality_sha256": inequality_digest.hexdigest(),
        "witness_reports_sha256": witness_reports_digest,
        "witness_kind_counts": dict(
            sorted(Counter(str(hardened[name]["report"].get("kind", "combined"))
                           for name in hardened).items())
        ),
        "worst_branch_cell": worst[0],
        "worst_branch_vertex": worst[1],
        "maximum_branch_log2_interval": [str(maximum.lo), str(maximum.hi)],
        "cell_local_union_terms": len(terms),
        "union_accounting": (
            "outward log2 sum over exact support-mesh root cells; BSP roots use "
            "their exact leaf maximum and the original verifier-recomputed root count upper"
        ),
        "union_log2_interval": [str(union.lo), str(union.hi)],
        "certified_margin_bits": str(margin),
        "top_cell_union_contributions": [row for _value, row in top_terms[:20]],
    }
    rendered = json.dumps(report, indent=2, sort_keys=True)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    if not passed:
        raise SystemExit("g4 end-to-end BSP verifier: outward 40-bit union did not close")


if __name__ == "__main__":
    main()
