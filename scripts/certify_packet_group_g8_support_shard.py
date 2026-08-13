#!/usr/bin/env python3
"""Independently replay exact-support g=8 slab/BSP shard certificates.

Producer floating-point scores, vertices, counts, and completion claims are
diagnostic only.  This verifier reconstructs the manifest census, exact
rational geometry, content-addressed witness rows, outward witness bounds,
and union aggregation.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import subprocess
import sys
from collections import Counter
from decimal import Decimal
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable

import certify_packet_group_triangle_ledger as outward
from certify_packet_group_profile import outward_normalization
from outward_log2 import Interval, LN2, ln_factorial, log2_int, self_check
from packet_group_drive_stratified import profile_classes


MANIFEST_SCHEMA = "packet-group-g8-support-manifest-v1"
SHARD_SCHEMA = "packet-group-g8-support-shard-v1"
REPLAY_SCHEMA = "packet-group-g8-support-replay-v1"
CANONICAL_JSON = "sorted-compact-json-v1"
SPLIT_RULE = "primitive-affine-gap-v1"
CLASSES = 9
GROUP_BITS = 8
EXPECTED_MASKS = set(range(2, 1 << CLASSES))
DEFAULT_MANIFEST_SHA256 = (
    "dafff5c978d51515355960293832223a8a4736ca4b2c0405496c89d395ac3890"
)
FROZEN_COMMIT = "3934dae73836e9051638db2e51b8b66a07055b3c"
FROZEN_TREE = "b7a7d3cbdc0b0361c5f55931e6e09afeddf762e8"

Profile = tuple[Fraction, ...]
Constraint = tuple[tuple[Fraction, ...], Fraction]


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")


def canonical_sha256(value: Any) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def parse_int(value: Any, label: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{label} must be an integer")
    try:
        parsed = int(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} must be an integer") from error
    if str(value).strip() not in (str(parsed), f"+{parsed}"):
        raise ValueError(f"{label} must use a canonical integer encoding")
    return parsed


def parse_mask(value: Any, label: str) -> int:
    if not isinstance(value, str) or len(value) != 5 or not value.startswith("0x"):
        raise ValueError(f"{label} must be a three-digit hexadecimal mask")
    try:
        parsed = int(value, 16)
    except ValueError as error:
        raise ValueError(f"{label} is malformed") from error
    if value != f"0x{parsed:03x}":
        raise ValueError(f"{label} is not canonical")
    return parsed


def parse_fraction(value: Any, label: str) -> Fraction:
    if not isinstance(value, str):
        raise ValueError(f"{label} must be a rational string")
    try:
        result = Fraction(value)
    except (ValueError, ZeroDivisionError) as error:
        raise ValueError(f"{label} is not a rational") from error
    canonical = f"{result.numerator}/{result.denominator}"
    if value != canonical:
        raise ValueError(f"{label} is not reduced and canonical")
    return result


def mask_text(mask: int) -> str:
    return f"0x{mask:03x}"


def support_of(mask: int) -> tuple[int, ...]:
    return tuple(index for index in range(CLASSES) if mask & (1 << index))


def excluded_support_profiles(support: tuple[int, ...], cutoff: int) -> int:
    """Count exact-support profiles below the physical-weight cutoff."""

    if 0 not in support:
        return 0
    positive = tuple(index for index in support if index)
    coefficients = [0] * cutoff
    coefficients[0] = 1
    for weight in positive:
        updated = [0] * cutoff
        for total, count in enumerate(coefficients):
            if not count:
                continue
            for multiplicity in range(1, (cutoff - 1 - total) // weight + 1):
                updated[total + multiplicity * weight] += count
        coefficients = updated
    return sum(coefficients)


def exact_support_count(
    support: tuple[int, ...], mass: int, cutoff: int
) -> tuple[int, int]:
    if not support:
        return 0, 0
    unconstrained = math.comb(mass - 1, len(support) - 1)
    excluded = excluded_support_profiles(support, cutoff)
    return unconstrained - excluded, excluded


def resolve_source(
    source: dict[str, Any], manifest_path: Path, overrides: dict[str, Path]
) -> Path:
    identifier = str(source["source_id"])
    if identifier in overrides:
        path = overrides[identifier]
    else:
        # The path is only a locator.  The digest is the authority.
        locator = Path(str(source.get("path", "")))
        path = locator if locator.is_absolute() else manifest_path.parent / locator
    if not path.is_file():
        raise ValueError(f"cannot resolve source_id {identifier!r}")
    actual = file_sha256(path)
    if actual != str(source["sha256"]):
        raise ValueError(f"source_id {identifier!r} digest mismatch")
    return path


def validate_manifest(
    manifest: dict[str, Any], manifest_path: Path, overrides: dict[str, Path]
) -> dict[str, Any]:
    if manifest.get("schema") != MANIFEST_SCHEMA:
        raise ValueError("unsupported g8 support manifest schema")
    parameters = manifest.get("parameters", {})
    expected = {
        "group_bits": GROUP_BITS,
        "N": 1 << 21,
        "K": 1 << 20,
        "M": 1 << 18,
        "minimum_profile_weight": 21,
    }
    parsed = {
        key: parse_int(parameters.get(key), f"parameters.{key}")
        for key in expected
    }
    if parsed != expected:
        raise ValueError(f"manifest parameter mismatch: {parsed!r}")
    if parse_fraction(parameters.get("probability_target_log2"), "probability target") != -40:
        raise ValueError("manifest probability target is not -40")
    arithmetic = manifest.get("arithmetic", {})
    if arithmetic.get("canonical_json") != CANONICAL_JSON:
        raise ValueError("manifest canonical JSON convention changed")
    if arithmetic.get("integer_split_rule") != SPLIT_RULE:
        raise ValueError("manifest integer split rule changed")
    if arithmetic.get("inner_point_cap_policy") != "pure-python-reference-v1":
        raise ValueError("manifest must require the pure-Python point-cap reference")
    closure = manifest.get("code_closure", {})
    if closure != {
        "kind": "git-tree-v1",
        "repository_commit": FROZEN_COMMIT,
        "repository_tree": FROZEN_TREE,
        "scope": "all imported discovery and outward proof code",
    }:
        raise ValueError("manifest code closure does not match the frozen proof tree")
    verify_git_code_closure(manifest_path)

    census = manifest.get("census", {})
    rows = census.get("feasible_masks")
    if not isinstance(rows, list) or len(rows) != len(EXPECTED_MASKS):
        raise ValueError("manifest must contain all 510 feasible masks")
    by_mask: dict[int, dict[str, Any]] = {}
    total = 0
    for row in rows:
        mask = parse_mask(row.get("mask"), "census mask")
        if mask in by_mask:
            raise ValueError(f"duplicate census mask {mask_text(mask)}")
        support = support_of(mask)
        if list(support) != row.get("active_classes"):
            raise ValueError(f"active classes mismatch for {mask_text(mask)}")
        if parse_int(row.get("dimension"), "census dimension") != len(support) - 1:
            raise ValueError(f"dimension mismatch for {mask_text(mask)}")
        count, excluded = exact_support_count(support, expected["M"], 21)
        if parse_int(row.get("count"), "census count") != count:
            raise ValueError(f"support count mismatch for {mask_text(mask)}")
        if parse_int(row.get("weight_cut_exclusions"), "weight exclusions") != excluded:
            raise ValueError(f"weight exclusions mismatch for {mask_text(mask)}")
        by_mask[mask] = row
        total += count
    if set(by_mask) != EXPECTED_MASKS:
        raise ValueError("manifest feasible-mask ownership is incomplete")
    infeasible = {
        parse_mask(value, "infeasible mask") for value in census.get("infeasible_masks", [])
    }
    if infeasible != {0, 1}:
        raise ValueError("manifest must reject exactly the empty mask and {0}")
    if parse_int(census.get("feasible_mask_count"), "feasible mask count") != 510:
        raise ValueError("manifest feasible mask count changed")
    if parse_int(census.get("feasible_profile_count"), "feasible profile count") != total:
        raise ValueError("manifest support counts do not sum to the stated census")

    source_paths: dict[str, Path] = {}
    source_records: dict[str, dict[str, Any]] = {}
    all_sources = manifest.get("witness_sources", []) + manifest.get("fixed_sources", [])
    for source in all_sources:
        identifier = str(source.get("source_id"))
        if not identifier or identifier in source_paths:
            raise ValueError(f"duplicate or empty source_id {identifier!r}")
        source_paths[identifier] = resolve_source(source, manifest_path, overrides)
        source_records[identifier] = source
    evaluator = str(arithmetic.get("outward_evaluator_source_id", ""))
    if evaluator not in source_paths:
        raise ValueError("manifest outward evaluator source is absent")
    required_runtime_sources = {
        "spectrum-01": Path(outward.__file__).resolve().parent.parent
        / "out"
        / "ebch85_band01_split_spectrum.csv",
        "spectrum-12": Path(outward.__file__).resolve().parent.parent
        / "out"
        / "ebch86_band12_split_spectrum.csv",
        "graph24-spectrum": Path(outward.__file__).resolve().with_name(
            "ebch128_graph24_spectrum.csv"
        ),
        "ebch-weight-distribution": Path(outward.__file__).resolve().with_name(
            "EBCH128_64.wd"
        ),
        "systematic-split-slices": Path(outward.__file__).resolve().with_name(
            "ebch128_systematic_split_slices.csv"
        ),
    }
    for identifier, runtime_path in required_runtime_sources.items():
        if identifier not in source_records:
            raise ValueError(f"manifest omits transitive proof input {identifier!r}")
        if file_sha256(runtime_path) != source_records[identifier]["sha256"]:
            raise ValueError(f"runtime proof input {identifier!r} digest mismatch")

    witness_rows: dict[tuple[str, int], tuple[dict[str, Any], str]] = {}
    for source in manifest.get("witness_sources", []):
        identifier = str(source["source_id"])
        artifact = json.loads(source_paths[identifier].read_text(encoding="utf-8"))
        if artifact.get("schema") != source.get("schema"):
            raise ValueError(f"witness source schema mismatch for {identifier!r}")
        rows_value = artifact.get("rows")
        declarations = source.get("rows")
        if not isinstance(rows_value, list) or not isinstance(declarations, list):
            raise ValueError(f"witness source {identifier!r} has no row catalogue")
        if len(rows_value) != len(declarations):
            raise ValueError(f"witness row count mismatch for {identifier!r}")
        for declaration in declarations:
            index = parse_int(declaration.get("row"), "witness row index")
            if not 0 <= index < len(rows_value):
                raise ValueError(f"witness row index outside {identifier!r}")
            row = rows_value[index]
            digest = canonical_sha256(row)
            if digest != declaration.get("row_sha256"):
                raise ValueError(f"witness row digest mismatch at {identifier}:{index}")
            row_mask = parse_mask(declaration.get("support_mask"), "witness support mask")
            if int(row.get("support_mask", -1)) != row_mask:
                raise ValueError(f"witness support mask mismatch at {identifier}:{index}")
            witness_rows[(identifier, index)] = (row, digest)
    return {
        "manifest_path": manifest_path,
        "parameters": parsed,
        "census": by_mask,
        "profile_count": total,
        "source_paths": source_paths,
        "source_records": source_records,
        "witness_rows": witness_rows,
        "hardened_witnesses": {},
    }


def verify_git_code_closure(manifest_path: Path) -> None:
    """Check every currently imported repository module against the frozen tree."""

    repository = manifest_path.resolve().parent
    git_directory = repository / ".git"
    if not git_directory.exists():
        # A verifier bundle without Git relies on the exact manifest fields
        # and individually hashed proof entry points/data files.
        return
    try:
        tree = subprocess.check_output(
            ["git", "-C", str(repository), "rev-parse", f"{FROZEN_COMMIT}^{{tree}}"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (OSError, subprocess.CalledProcessError) as error:
        raise ValueError("cannot resolve the frozen Git code closure") from error
    if tree != FROZEN_TREE:
        raise ValueError("frozen Git commit does not resolve to the manifest tree")
    scripts_root = (repository / "scripts").resolve()
    verifier_path = Path(__file__).resolve()
    imported_paths = set()
    for module in tuple(sys.modules.values()):
        raw_path = getattr(module, "__file__", None)
        if not raw_path:
            continue
        path = Path(raw_path).resolve()
        if path == verifier_path or scripts_root not in path.parents or path.suffix != ".py":
            continue
        imported_paths.add(path)
    for path in sorted(imported_paths):
        relative = path.relative_to(repository).as_posix()
        tracked = subprocess.run(
            ["git", "-C", str(repository), "cat-file", "-e", f"{FROZEN_COMMIT}:{relative}"],
            check=False,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        if tracked.returncode != 0:
            raise ValueError(f"imported proof module is absent from frozen tree: {relative}")
        result = subprocess.run(
            ["git", "-C", str(repository), "diff", "--quiet", FROZEN_COMMIT, "--", relative],
            check=False,
        )
        if result.returncode != 0:
            raise ValueError(f"imported proof module differs from frozen tree: {relative}")


def solve_linear(
    matrix: Iterable[Iterable[Fraction]], right: Iterable[Fraction]
) -> tuple[Fraction, ...] | None:
    rows = [list(row) + [value] for row, value in zip(matrix, right)]
    n = len(rows)
    if n == 0 or any(len(row) != n + 1 for row in rows):
        raise ValueError("malformed exact linear system")
    for column in range(n):
        pivot = next((index for index in range(column, n) if rows[index][column]), None)
        if pivot is None:
            return None
        rows[column], rows[pivot] = rows[pivot], rows[column]
        scale = rows[column][column]
        rows[column] = [value / scale for value in rows[column]]
        for index in range(n):
            if index == column or not rows[index][column]:
                continue
            factor = rows[index][column]
            rows[index] = [
                left - factor * pivot_value
                for left, pivot_value in zip(rows[index], rows[column])
            ]
    return tuple(row[-1] for row in rows)


def root_constraints(support: tuple[int, ...], mass: int, cutoff: int) -> tuple[Constraint, ...]:
    if len(support) < 2:
        return ()
    chart = support[:-1]
    dependent = support[-1]
    dimension = len(chart)
    constraints: list[Constraint] = []
    for coordinate in range(dimension):
        left = [Fraction(0)] * dimension
        left[coordinate] = -1
        constraints.append((tuple(left), Fraction(-1)))
    constraints.append((tuple(Fraction(1) for _ in chart), Fraction(mass - 1)))
    # sum_j j*a_j >= cutoff after eliminating the dependent coordinate.
    constraints.append(
        (
            tuple(Fraction(dependent - index) for index in chart),
            Fraction(dependent * mass - cutoff),
        )
    )
    return tuple(constraints)


def split_constraint(
    coefficients: tuple[int, ...], threshold: int, support: tuple[int, ...], left: bool
) -> Constraint:
    chart = support[:-1]
    dependent = support[-1]
    reduced = tuple(Fraction(coefficients[index] - coefficients[dependent]) for index in chart)
    if left:
        return reduced, Fraction(threshold - coefficients[dependent] * (1 << 18))
    return tuple(-value for value in reduced), Fraction(
        coefficients[dependent] * (1 << 18) - threshold - 1
    )


def full_profile(support: tuple[int, ...], chart_value: tuple[Fraction, ...]) -> Profile:
    profile = [Fraction(0)] * CLASSES
    for index, value in zip(support[:-1], chart_value):
        profile[index] = value
    profile[support[-1]] = Fraction(1 << 18) - sum(chart_value)
    return tuple(profile)


def enumerate_vertices(
    support: tuple[int, ...], constraints: tuple[Constraint, ...], maximum_systems: int
) -> tuple[Profile, ...]:
    dimension = len(support) - 1
    if dimension == 0:
        profile = [Fraction(0)] * CLASSES
        profile[support[0]] = 1 << 18
        candidate = tuple(profile)
        if sum(index * value for index, value in enumerate(candidate)) < 21:
            return ()
        return (candidate,)
    systems = math.comb(len(constraints), dimension)
    if systems > maximum_systems:
        raise ValueError(
            f"vertex reconstruction needs {systems} systems, above limit {maximum_systems}"
        )
    vertices: set[Profile] = set()
    for active in itertools.combinations(range(len(constraints)), dimension):
        solution = solve_linear(
            [constraints[index][0] for index in active],
            [constraints[index][1] for index in active],
        )
        if solution is None:
            continue
        if all(
            sum(coefficient * value for coefficient, value in zip(left, solution)) <= right
            for left, right in constraints
        ):
            vertices.add(full_profile(support, solution))
    return tuple(sorted(vertices))


def parse_witness_reference(
    value: Any,
    context: dict[str, Any],
    expected_manifest_sources: dict[str, dict[str, Any]],
) -> tuple[str, int, dict[str, Any]]:
    if not isinstance(value, dict):
        raise ValueError("witness reference must be an object")
    identifier = str(value.get("source_id", ""))
    index = parse_int(value.get("row"), "witness row")
    key = (identifier, index)
    if key not in context["witness_rows"]:
        raise ValueError(f"unmanifested witness reference {identifier}:{index}")
    source = expected_manifest_sources[identifier]
    if value.get("source_sha256") != source.get("sha256"):
        raise ValueError(f"witness source digest mismatch for {identifier}:{index}")
    row, digest = context["witness_rows"][key]
    if value.get("row_sha256") != digest:
        raise ValueError(f"witness row digest mismatch for {identifier}:{index}")
    return identifier, index, row


def parse_selector(
    value: Any, context: dict[str, Any]
) -> tuple[tuple[str, Fraction], ...]:
    if not isinstance(value, dict):
        raise ValueError("certified leaf selector must be an object")
    kind = value.get("kind")
    sources = context["source_records"]
    if kind == "witness":
        identifier, index, _row = parse_witness_reference(value.get("witness"), context, sources)
        return ((f"{identifier}:{index}", Fraction(1)),)
    if kind != "mixture" or not isinstance(value.get("components"), list):
        raise ValueError("selector must be one witness or one fixed mixture")
    combined: dict[str, Fraction] = {}
    for component in value["components"]:
        if not isinstance(component, dict):
            raise ValueError("mixture component must be an object")
        weight = parse_fraction(component.get("weight"), "mixture weight")
        if weight < 0:
            raise ValueError("mixture weight must be nonnegative")
        identifier, index, _row = parse_witness_reference(
            component.get("witness"), context, sources
        )
        name = f"{identifier}:{index}"
        combined[name] = combined.get(name, Fraction(0)) + weight
    if sum(combined.values(), Fraction(0)) != 1:
        raise ValueError("fixed mixture weights must sum exactly to one")
    return tuple(sorted((name, weight) for name, weight in combined.items() if weight))


def validate_split(value: Any) -> tuple[tuple[int, ...], int]:
    if not isinstance(value, dict) or not isinstance(value.get("coefficients"), list):
        raise ValueError("split must contain an integer coefficient vector")
    coefficients = tuple(
        parse_int(item, "split coefficient") for item in value["coefficients"]
    )
    if len(coefficients) != CLASSES or not any(coefficients):
        raise ValueError("split coefficient vector must have nine entries and be nonzero")
    if math.gcd(*map(abs, coefficients)) != 1:
        raise ValueError("split coefficient vector must be primitive")
    if next(item for item in coefficients if item) < 0:
        raise ValueError("split coefficient vector has noncanonical sign")
    return coefficients, parse_int(value.get("threshold"), "split threshold")


def replay_tree(
    shard: dict[str, Any], support: tuple[int, ...], context: dict[str, Any], maximum_systems: int
) -> tuple[
    list[dict[str, Any]],
    dict[str, tuple[Profile, ...]],
    dict[str, tuple[Profile, ...]],
    Counter,
]:
    raw_nodes = shard.get("nodes")
    if not isinstance(raw_nodes, list) or not raw_nodes:
        raise ValueError("shard node list is empty")
    nodes: dict[str, dict[str, Any]] = {}
    for node in raw_nodes:
        if not isinstance(node, dict) or not isinstance(node.get("node_id"), str):
            raise ValueError("every node needs a string node_id")
        identifier = node["node_id"]
        if identifier in nodes:
            raise ValueError(f"duplicate node_id {identifier!r}")
        nodes[identifier] = node
    root = shard.get("root_node")
    if root not in nodes:
        raise ValueError("root_node is absent")
    visited: set[str] = set()
    leaves: list[dict[str, Any]] = []
    vertices_by_leaf: dict[str, tuple[Profile, ...]] = {}
    vertices_by_node: dict[str, tuple[Profile, ...]] = {}
    counts: Counter = Counter()

    def visit(identifier: str, constraints: tuple[Constraint, ...]) -> None:
        if identifier in visited:
            raise ValueError(f"node {identifier!r} has multiple parents or a cycle")
        visited.add(identifier)
        node = nodes.get(identifier)
        if node is None:
            raise ValueError(f"missing node {identifier!r}")
        state = node.get("state")
        if state not in ("EMPTY", "UNRESOLVED", "SPLIT", "CERTIFIED_LEAF"):
            raise ValueError(f"node {identifier!r} has invalid state {state!r}")
        counts[state.lower()] += 1
        vertices = enumerate_vertices(support, constraints, maximum_systems)
        vertices_by_node[identifier] = vertices
        if state == "EMPTY":
            if node.get("emptiness_method") != "exact-rational-polytope-v1":
                raise ValueError(f"empty node {identifier!r} uses an unsupported proof")
            if vertices:
                raise ValueError(f"empty node {identifier!r} has a nonempty real polytope")
            return
        if not vertices:
            raise ValueError(f"nonempty node {identifier!r} has an empty real polytope")
        if state == "UNRESOLVED":
            return
        if state == "CERTIFIED_LEAF":
            node["_mixture"] = parse_selector(node.get("selector"), context)
            leaves.append(node)
            vertices_by_leaf[identifier] = vertices
            return
        coefficients, threshold = validate_split(node.get("split"))
        left = node.get("left")
        right = node.get("right")
        if not isinstance(left, str) or not isinstance(right, str) or left == right:
            raise ValueError(f"split node {identifier!r} has invalid children")
        visit(left, constraints + (split_constraint(coefficients, threshold, support, True),))
        visit(right, constraints + (split_constraint(coefficients, threshold, support, False),))

    visit(str(root), root_constraints(support, 1 << 18, 21))
    if visited != set(nodes):
        raise ValueError(f"unreachable nodes: {sorted(set(nodes) - visited)!r}")
    return leaves, vertices_by_leaf, vertices_by_node, counts


def iter_support_profiles(support: tuple[int, ...], mass: int, cutoff: int):
    """Enumerate a small exact-support root; callers enforce a size limit."""

    profile = [0] * CLASSES

    def visit(position: int, remaining: int):
        if position == len(support) - 1:
            profile[support[position]] = remaining
            if remaining >= 1 and sum(index * count for index, count in enumerate(profile)) >= cutoff:
                yield tuple(profile)
            profile[support[position]] = 0
            return
        index = support[position]
        minimum_tail = len(support) - position - 1
        for count in range(1, remaining - minimum_tail + 1):
            profile[index] = count
            yield from visit(position + 1, remaining - count)
        profile[index] = 0

    yield from visit(0, mass)


def owned_path(profile: tuple[int, ...], shard: dict[str, Any]) -> tuple[str, ...]:
    nodes = {str(node["node_id"]): node for node in shard["nodes"]}
    identifier = str(shard["root_node"])
    path = []
    while True:
        path.append(identifier)
        node = nodes[identifier]
        state = node["state"]
        if state == "CERTIFIED_LEAF":
            return tuple(path)
        if state == "EMPTY":
            raise ValueError(f"integer profile {profile!r} is owned by an EMPTY node")
        if state == "UNRESOLVED":
            raise ValueError("unresolved node owns an integer profile")
        coefficients, threshold = validate_split(node["split"])
        value = sum(left * right for left, right in zip(coefficients, profile))
        identifier = str(node["left"] if value <= threshold else node["right"])


def bounded_leaf_counts(
    support: tuple[int, ...], root_count: int, shard: dict[str, Any], limit: int
) -> dict[str, int] | None:
    if root_count > limit:
        return None
    counts: Counter = Counter()
    seen = 0
    for profile in iter_support_profiles(support, 1 << 18, 21):
        for identifier in owned_path(profile, shard):
            counts[identifier] += 1
        seen += 1
    if seen != root_count:
        raise ValueError(f"bounded enumeration count mismatch: {seen} != {root_count}")
    return dict(counts)


def box_count(lower: tuple[int, ...], upper: tuple[int, ...], support: tuple[int, ...]) -> int:
    if len(lower) != CLASSES or len(upper) != CLASSES:
        raise ValueError("coordinate box needs nine lower and upper bounds")
    for index in range(CLASSES):
        expected_positive = index in support
        if lower[index] < (1 if expected_positive else 0) or upper[index] < lower[index]:
            raise ValueError("coordinate box bounds violate exact support")
        if not expected_positive and (lower[index] or upper[index]):
            raise ValueError("coordinate box activates an absent class")
    active_lower = [lower[index] for index in support]
    widths = [upper[index] - lower[index] for index in support]
    remaining = (1 << 18) - sum(active_lower)
    total = 0
    for subset in range(1 << len(support)):
        removed = sum(
            widths[index] + 1 for index in range(len(support)) if subset & (1 << index)
        )
        top = remaining - removed + len(support) - 1
        term = 0 if top < len(support) - 1 else math.comb(top, len(support) - 1)
        total += -term if subset.bit_count() & 1 else term
    # The excluded region has physical weight at most 20 and is tiny.
    excluded = 0
    if 0 in support:
        positive = tuple(index for index in support if index)
        profile = [0] * CLASSES

        def enumerate_low(position: int, physical: int):
            nonlocal excluded
            if position == len(positive):
                profile[0] = (1 << 18) - sum(profile)
                if all(lower[index] <= profile[index] <= upper[index] for index in range(CLASSES)):
                    excluded += 1
                profile[0] = 0
                return
            index = positive[position]
            for count in range(1, (20 - physical) // index + 1):
                profile[index] = count
                enumerate_low(position + 1, physical + index * count)
            profile[index] = 0

        enumerate_low(0, 0)
    return total - excluded


def parse_bound_vector(value: Any, label: str) -> tuple[int, ...]:
    if not isinstance(value, list) or len(value) != CLASSES:
        raise ValueError(f"{label} must contain nine integers")
    return tuple(parse_int(item, label) for item in value)


def verify_upper_count_record(
    record: dict[str, Any],
    vertices: tuple[Profile, ...],
    support: tuple[int, ...],
    root_count: int,
) -> int:
    method = record.get("method")
    claimed = parse_int(record.get("value"), "count value")
    if claimed < 0:
        raise ValueError("count upper bound must be nonnegative")
    if method == "root-census-upper-v1":
        expected = root_count
    elif method == "coordinate-box-ie-v1":
        lower = parse_bound_vector(record.get("lower"), "coordinate-box lower")
        upper = parse_bound_vector(record.get("upper"), "coordinate-box upper")
        for vertex in vertices:
            if any(
                vertex[index] < lower[index] or vertex[index] > upper[index]
                for index in range(CLASSES)
            ):
                raise ValueError("coordinate box does not contain the complete leaf polytope")
        expected = box_count(lower, upper, support)
    elif method == "verified-min-v1":
        bounds = record.get("bounds")
        if not isinstance(bounds, list) or len(bounds) < 2:
            raise ValueError("verified-min-v1 needs at least two parameterized bounds")
        values = []
        for bound in bounds:
            if not isinstance(bound, dict) or bound.get("kind") != "upper":
                raise ValueError("verified-min-v1 components must be upper count records")
            if bound.get("method") == "verified-min-v1":
                raise ValueError("nested verified-min-v1 count records are not canonical")
            values.append(verify_upper_count_record(bound, vertices, support, root_count))
        expected = min(values)
    else:
        raise ValueError(f"unsupported upper-count method {method!r}")
    if claimed != expected:
        raise ValueError(f"count upper bound mismatch: {claimed} != {expected}")
    return expected


def verify_node_count_records(
    nodes: list[dict[str, Any]],
    vertices_by_node: dict[str, tuple[Profile, ...]],
    exact_node_counts: dict[str, int] | None,
    support: tuple[int, ...],
    root_count: int,
    root_node: str,
) -> None:
    by_id = {str(node["node_id"]): node for node in nodes}
    for node in nodes:
        identifier = str(node["node_id"])
        record = node.get("count")
        if record is None:
            continue
        if not isinstance(record, dict):
            raise ValueError(f"leaf {identifier!r} count record must be an object")
        kind = record.get("kind")
        if kind == "exact":
            method = record.get("method")
            claimed = parse_int(record.get("value"), "exact count value")
            if method == "root-census-v1":
                if identifier != root_node:
                    raise ValueError("root-census-v1 is valid only on the root node")
                expected = root_count
            elif method == "enumeration-v1":
                if exact_node_counts is None:
                    raise ValueError("enumeration-v1 exceeds the verifier enumeration limit")
                expected = exact_node_counts.get(identifier, 0)
            else:
                raise ValueError(f"unsupported exact-count method {method!r}")
            if claimed != expected:
                raise ValueError(f"exact leaf count mismatch: {claimed} != {expected}")
        elif kind == "upper":
            upper = verify_upper_count_record(
                record, vertices_by_node[identifier], support, root_count
            )
            if exact_node_counts is not None and upper < exact_node_counts.get(identifier, 0):
                raise ValueError("claimed count upper bound is below exact enumeration")
        else:
            raise ValueError(f"leaf {identifier!r} has invalid count kind {kind!r}")
    for node in nodes:
        if node.get("state") != "SPLIT":
            continue
        parent_record = node.get("count")
        left_record = by_id[str(node["left"])].get("count")
        right_record = by_id[str(node["right"])].get("count")
        records = (parent_record, left_record, right_record)
        if all(isinstance(record, dict) and record.get("kind") == "exact" for record in records):
            parent, left, right = (
                parse_int(record["value"], "exact conservation count") for record in records
            )
            if left + right != parent:
                raise ValueError(f"exact child counts do not conserve split node {node['node_id']!r}")


def scale_interval(value: Interval, factor: Fraction) -> Interval:
    return value * (
        Interval.exact(factor.numerator) / Interval.exact(factor.denominator)
    )


def fractional_normalization(profile: Profile) -> Interval:
    natural = ln_factorial(1 << 18)
    for count in profile:
        if count.denominator == 1:
            gamma = ln_factorial(count.numerator)
        else:
            floor = count.numerator // count.denominator
            fraction = count - floor
            chord = scale_interval(ln_factorial(floor), 1 - fraction) + scale_interval(
                ln_factorial(floor + 1), fraction
            )
            gamma = Interval(Decimal(-1), chord.hi)
        natural = natural - gamma
    result = natural / LN2
    for count, classes in zip(profile, profile_classes(GROUP_BITS)):
        if count:
            result = result + scale_interval(log2_int(classes), count)
    return result


def atlas_row_to_hardening(row: dict[str, Any]) -> dict[str, Any]:
    inner = row.get("inner", {})
    outer_row = row.get("outer", {})
    return {
        "fugacities": inner.get("fugacities"),
        "pole": inner.get("pole"),
        "outer_type": "conditioned_row_exact_graph_asymmetric",
        "outer_details": {
            "log_variables": outer_row.get("log_variables"),
            "band_coefficients": outer_row.get("band_coefficients"),
            "pair_cauchy_theta": outer_row.get("pair_cauchy_theta"),
        },
    }


def harden_used_witnesses(
    names: set[str], context: dict[str, Any], iterations: int
) -> dict[str, dict[str, Any]]:
    # The frozen manifest does not content-address a native library.  Force
    # the independent pure-Python outward recurrence; native diagnostics may
    # still propose a positive Collatz vector, which cannot affect soundness.
    outward.load_native_point_caps = lambda: None
    cached = context["hardened_witnesses"]
    missing = sorted(names - set(cached))
    if not missing:
        return {name: cached[name] for name in names}
    outward._initialize_hardening_worker(GROUP_BITS)
    split_caps = outward._WORKER_SPLIT_CAPS
    if split_caps is None:
        raise RuntimeError("outward hardening was not initialized")
    for name in missing:
        identifier, raw_index = name.rsplit(":", 1)
        index = int(raw_index)
        row, _digest = context["witness_rows"][(identifier, index)]
        source_path = context["source_paths"][identifier]
        cached[name] = outward.harden_witness(
            name,
            atlas_row_to_hardening(row),
            source_path,
            split_caps,
            iterations,
        )
    return {name: cached[name] for name in names}


def evaluate_profile(
    profile: Profile,
    mixture: tuple[tuple[str, Fraction], ...],
    hardened: dict[str, dict[str, Any]],
) -> Interval:
    total = Interval.exact(0)
    normalization = fractional_normalization(profile)
    for name, weight in mixture:
        witness = hardened[name]
        value = witness["constant"]
        for count, charge in zip(profile, witness["charges"]):
            if not count:
                continue
            if charge is None:
                raise ValueError(f"witness {name!r} has zero fugacity on an active class")
            value = value - scale_interval(charge, count)
        if witness.get("subtract_normalization", True):
            value = value - normalization
        total = total + scale_interval(value, weight)
    return total


def maximum_interval(left: Interval | None, right: Interval) -> Interval:
    if left is None:
        return right
    return Interval(max(left.lo, right.lo), max(left.hi, right.hi))


def verify_shard(
    shard_path: Path,
    manifest_sha256: str,
    context: dict[str, Any],
    iterations: int,
    maximum_systems: int,
    enumeration_limit: int,
    structure_only: bool,
) -> dict[str, Any]:
    verifier_sha256 = file_sha256(Path(__file__).resolve())
    shard = json.loads(shard_path.read_text(encoding="utf-8"))
    if shard.get("schema") != SHARD_SCHEMA:
        raise ValueError(f"unsupported shard schema in {shard_path}")
    if shard.get("manifest_sha256") != manifest_sha256:
        raise ValueError(f"manifest digest mismatch in {shard_path}")
    mask = parse_mask(shard.get("support_mask"), "shard support mask")
    if mask not in context["census"]:
        raise ValueError(f"shard owns infeasible mask {mask_text(mask)}")
    support = support_of(mask)
    if shard.get("active_classes") != list(support):
        raise ValueError(f"shard active classes mismatch for {mask_text(mask)}")
    root_count = int(context["census"][mask]["count"])
    if parse_int(shard.get("root_count"), "shard root count") != root_count:
        raise ValueError(f"shard root count mismatch for {mask_text(mask)}")
    leaves, vertices_by_leaf, vertices_by_node, node_counts = replay_tree(
        shard, support, context, maximum_systems
    )
    if node_counts["unresolved"]:
        raise ValueError(f"shard {mask_text(mask)} contains unresolved nodes")
    if shard.get("state") != "COMPLETE":
        raise ValueError(f"shard {mask_text(mask)} does not claim COMPLETE")
    exact_node_counts = bounded_leaf_counts(support, root_count, shard, enumeration_limit)
    leaf_ids = {str(leaf["node_id"]) for leaf in leaves}
    if exact_node_counts is not None and sum(
        exact_node_counts.get(identifier, 0) for identifier in leaf_ids
    ) != root_count:
        raise ValueError("exact leaf counts do not conserve the root")
    verify_node_count_records(
        shard["nodes"],
        vertices_by_node,
        exact_node_counts,
        support,
        root_count,
        str(shard["root_node"]),
    )

    used = {
        name for leaf in leaves for name, _weight in leaf["_mixture"]
    }
    if structure_only:
        return {
            "schema": REPLAY_SCHEMA,
            "manifest_sha256": manifest_sha256,
            "verifier_source_sha256": verifier_sha256,
            "shard_sha256": file_sha256(shard_path),
            "support_mask": mask_text(mask),
            "status": "STRUCTURE_ONLY_COMPLETE",
            "root_count": str(root_count),
            "node_counts": {
                state: node_counts[state]
                for state in ("empty", "certified_leaf", "split", "unresolved")
            },
            "used_witnesses": len(used),
        }

    hardened = harden_used_witnesses(used, context, iterations)
    # Capture Python modules imported lazily by the hardening path as well.
    verify_git_code_closure(context["manifest_path"])
    for name in used:
        for index in support:
            if hardened[name]["charges"][index] is None:
                raise ValueError(f"witness {name!r} is ineligible for {mask_text(mask)}")
    leaf_maxima: dict[str, Interval] = {}
    inequality_digest = hashlib.sha256()
    for leaf in sorted(leaves, key=lambda row: str(row["node_id"])):
        identifier = str(leaf["node_id"])
        maximum = None
        for profile in vertices_by_leaf[identifier]:
            value = evaluate_profile(profile, leaf["_mixture"], hardened)
            maximum = maximum_interval(maximum, value)
            inequality_digest.update(
                canonical_bytes(
                    [
                        identifier,
                        [str(item) for item in profile],
                        [[name, str(weight)] for name, weight in leaf["_mixture"]],
                        str(value.lo),
                        str(value.hi),
                    ]
                )
            )
        if maximum is None:
            raise ValueError(f"certified leaf {identifier!r} has no vertices")
        leaf_maxima[identifier] = maximum

    global_maximum = None
    for value in leaf_maxima.values():
        global_maximum = maximum_interval(global_maximum, value)
    assert global_maximum is not None
    collapsed = global_maximum + log2_int(root_count)
    expanded = None
    if exact_node_counts is not None:
        terms = [
            leaf_maxima[identifier] + log2_int(count)
            for identifier, count in exact_node_counts.items()
            if identifier in leaf_ids
            if count
        ]
        if terms:
            upper = outward._log2_sum_exp(Interval.exact(term.hi) for term in terms)
            expanded = Interval(max(term.lo for term in terms), upper.hi)
    requested = shard.get("aggregation_requested")
    if requested not in ("collapsed", "expanded", "minimum"):
        raise ValueError("shard aggregation_requested is invalid")
    if requested == "expanded" and expanded is None:
        raise ValueError("expanded aggregation requires independently exact leaf counts")
    if requested == "collapsed" or expanded is None:
        aggregation_used, contribution = "collapsed", collapsed
    elif requested == "expanded":
        aggregation_used, contribution = "expanded", expanded
    elif expanded.hi < collapsed.hi:
        aggregation_used, contribution = "expanded", expanded
    else:
        aggregation_used, contribution = "collapsed", collapsed
    return {
        "schema": REPLAY_SCHEMA,
        "manifest_sha256": manifest_sha256,
        "verifier_source_sha256": verifier_sha256,
        "shard_sha256": file_sha256(shard_path),
        "support_mask": mask_text(mask),
        "status": "VERIFIED_COMPLETE",
        "root_count": str(root_count),
        "node_counts": {
            state: node_counts[state]
            for state in ("empty", "certified_leaf", "split", "unresolved")
        },
        "used_witnesses": len(used),
        "vertex_inequalities": sum(len(vertices_by_leaf[str(row["node_id"])]) for row in leaves),
        "inequality_sha256": inequality_digest.hexdigest(),
        "aggregation_used": aggregation_used,
        "collapsed_log2_interval": [str(collapsed.lo), str(collapsed.hi)],
        "expanded_log2_interval": (
            None if expanded is None else [str(expanded.lo), str(expanded.hi)]
        ),
        "contribution_log2_interval": [str(contribution.lo), str(contribution.hi)],
        "_contribution": contribution,
    }


def structural_self_test() -> None:
    assert exact_support_count((0, 1), 1 << 18, 21) == (262123, 20)
    assert exact_support_count((1,), 1 << 18, 21) == (1, 0)
    support = (0, 1)
    constraints = root_constraints(support, 1 << 18, 21)
    vertices = enumerate_vertices(support, constraints, 100)
    assert vertices == (
        (Fraction(1), Fraction(262143), *([Fraction(0)] * 7)),
        (Fraction(262123), Fraction(21), *([Fraction(0)] * 7)),
    )
    left = constraints + (split_constraint((1, 0, 0, 0, 0, 0, 0, 0, 0), 100, support, True),)
    right = constraints + (split_constraint((1, 0, 0, 0, 0, 0, 0, 0, 0), 100, support, False),)
    assert max(vertex[0] for vertex in enumerate_vertices(support, left, 100)) == 100
    assert min(vertex[0] for vertex in enumerate_vertices(support, right, 100)) == 101


def parse_overrides(values: list[str]) -> dict[str, Path]:
    result = {}
    for value in values:
        if "=" not in value:
            raise ValueError("--source must have the form source_id=PATH")
        identifier, path = value.split("=", 1)
        if not identifier or identifier in result:
            raise ValueError(f"duplicate or empty source override {identifier!r}")
        result[identifier] = Path(path)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=Path("G8_SUPPORT_MANIFEST.json"))
    parser.add_argument("--manifest-sha256", default=DEFAULT_MANIFEST_SHA256)
    parser.add_argument("--shard", type=Path, action="append", required=True)
    parser.add_argument("--source", action="append", default=[])
    parser.add_argument("--iterations", type=int, default=3)
    parser.add_argument("--maximum-vertex-systems", type=int, default=2_000_000)
    parser.add_argument("--enumeration-limit", type=int, default=1_000_000)
    parser.add_argument("--structure-only", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.iterations < 1 or args.maximum_vertex_systems < 1 or args.enumeration_limit < 0:
        parser.error("iteration and reconstruction limits are invalid")
    self_check()
    structural_self_test()
    actual_manifest_sha256 = file_sha256(args.manifest)
    if actual_manifest_sha256 != args.manifest_sha256:
        raise SystemExit(
            f"g8 shard verifier: manifest digest mismatch: {actual_manifest_sha256}"
        )
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    context = validate_manifest(manifest, args.manifest, parse_overrides(args.source))
    receipts = [
        verify_shard(
            path,
            actual_manifest_sha256,
            context,
            args.iterations,
            args.maximum_vertex_systems,
            args.enumeration_limit,
            args.structure_only,
        )
        for path in args.shard
    ]
    masks = [parse_mask(receipt["support_mask"], "receipt mask") for receipt in receipts]
    if len(masks) != len(set(masks)):
        raise ValueError("duplicate support shards were supplied")
    complete = set(masks) == EXPECTED_MASKS
    if complete and sum(int(receipt["root_count"]) for receipt in receipts) != context["profile_count"]:
        raise ValueError("integrated shard root counts do not sum to the manifest census")
    report: dict[str, Any] = {
        "schema": "packet-group-g8-support-integration-v1",
        "manifest_sha256": actual_manifest_sha256,
        "verifier_source_sha256": file_sha256(Path(__file__).resolve()),
        "verified_shards": len(receipts),
        "complete_support_census": complete,
        "receipts": receipts,
    }
    if complete and not args.structure_only:
        contributions = [receipt.pop("_contribution") for receipt in receipts]
        upper = outward._log2_sum_exp(
            Interval.exact(contribution.hi) for contribution in contributions
        )
        union = Interval(max(value.lo for value in contributions), upper.hi)
        report.update(
            {
                "status": "VERIFIED_COMPLETE" if union.hi <= -40 else "VERIFIED_COMPLETE_DID_NOT_CLOSE",
                "global_union_log2_interval": [str(union.lo), str(union.hi)],
                "passed": union.hi <= -40,
            }
        )
    else:
        for receipt in receipts:
            receipt.pop("_contribution", None)
        report.update(
            {
                "status": "PARTIAL_REPLAY",
                "passed": False,
                "scope": "not all 510 exact-support shards were supplied",
            }
        )
    rendered = json.dumps(report, indent=2, sort_keys=True)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
