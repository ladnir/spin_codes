#!/usr/bin/env python3
"""Independent gates for g=8 factorized inner/outer selectors.

This module does not consume producer affine constants or charges.  It
reconstructs each component from frozen parameters, evaluates two rational
marginals, and subtracts the packet-profile normalization exactly once.
It also replays the exact cumulative-prefix geometry used by the 189-leaf
full-support certificate shape.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math
from functools import lru_cache
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np

import certify_packet_group_triangle_ledger as outward
from certify_packet8_hard_face_drive_inner import (
    exact_float,
    outward_transport_image_decimal_upper,
)
from outward_log2 import Interval, log2_fraction


MANIFEST_SCHEMA = "packet-group-g8-factorized-manifest-v2"
CERTIFICATE_SCHEMA = "packet-group-g8-factorized-full-support-certificate-v1"
SELECTOR_KIND = "independent-sum-mixture-v1"
NORMALIZATION_RULE = "packet-profile-orbit-v1"
INNER_INTERFACE = "uniform-conditional-on-complete-outer-word-v1"
COMBINATION_RULE = "sum-marginals-subtract-normalization-once-v1"
GROUP_BITS = 8
CLASSES = 9
PROFILE_MASS = 1 << 18
FREE_MASS = PROFILE_MASS - CLASSES
H2_BINS = 4
EXPECTED_ROOTS = math.comb(H2_BINS + CLASSES - 2, CLASSES - 1)
EXPECTED_LEAVES = 189

Profile = tuple[Fraction, ...]
Reference = tuple[str, int]
Marginal = tuple[tuple[Reference, Fraction], ...]


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_sha256(value: Any) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def parse_int(value: Any, label: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{label} must be an integer")
    if isinstance(value, int):
        return value
    if isinstance(value, str) and value and value.lstrip("-").isdigit():
        parsed = int(value)
        if str(parsed) != value:
            raise ValueError(f"{label} is not canonical")
        return parsed
    raise ValueError(f"{label} must be an integer")


def parse_fraction(value: Any, label: str) -> Fraction:
    if not isinstance(value, str) or value.count("/") != 1:
        raise ValueError(f"{label} must be a reduced p/q string")
    numerator_text, denominator_text = value.split("/")
    numerator = parse_int(numerator_text, f"{label} numerator")
    denominator = parse_int(denominator_text, f"{label} denominator")
    if denominator <= 0:
        raise ValueError(f"{label} denominator must be positive")
    result = Fraction(numerator, denominator)
    if f"{result.numerator}/{result.denominator}" != value:
        raise ValueError(f"{label} must be reduced")
    return result


def balanced_bins(total: int = FREE_MASS, bins: int = H2_BINS) -> tuple[tuple[int, int], ...]:
    points = total + 1
    return tuple(
        (index * points // bins, (index + 1) * points // bins - 1)
        for index in range(bins)
    )


def close_bounds(
    lower: tuple[int, ...], upper: tuple[int, ...]
) -> tuple[tuple[int, ...], tuple[int, ...]] | None:
    if len(lower) != CLASSES - 1 or len(upper) != CLASSES - 1:
        raise ValueError("cumulative bound has the wrong dimension")
    low = list(lower)
    high = list(upper)
    for index in range(1, len(low)):
        low[index] = max(low[index], low[index - 1])
    for index in range(len(high) - 2, -1, -1):
        high[index] = min(high[index], high[index + 1])
    if any(left > right for left, right in zip(low, high)):
        return None
    return tuple(low), tuple(high)


def child_bounds(
    lower: tuple[int, ...], upper: tuple[int, ...], coordinate: int,
    threshold: int, left: bool,
) -> tuple[tuple[int, ...], tuple[int, ...]] | None:
    if not 0 <= coordinate < CLASSES - 1:
        raise ValueError("cumulative split coordinate is invalid")
    low = list(lower)
    high = list(upper)
    if left:
        high[coordinate] = min(high[coordinate], threshold)
    else:
        low[coordinate] = max(low[coordinate], threshold + 1)
    return close_bounds(tuple(low), tuple(high))


@lru_cache(maxsize=None)
def exact_chain_count(lower: tuple[int, ...], upper: tuple[int, ...]) -> int:
    """Count a bounded nondecreasing integer chain with exact integer DP."""

    closed = close_bounds(lower, upper)
    if closed is None:
        return 0
    lower, upper = closed
    previous_low = lower[0]
    previous_high = upper[0]
    previous = [1] * (previous_high - previous_low + 1)
    for low, high in zip(lower[1:], upper[1:]):
        prefixes = list(itertools.accumulate(previous))
        total = prefixes[-1]
        if low > previous_high:
            current = [total] * (high - low + 1)
        else:
            start = max(0, low - previous_low)
            current = prefixes[start : high - previous_low + 1]
            if high > previous_high:
                current.extend([total] * (high - previous_high))
        previous = current
        previous_low = low
        previous_high = high
    return sum(previous)


def compositions(total: int):
    if total == 0:
        yield ()
        return
    for first in range(1, total + 1):
        for rest in compositions(total - first):
            yield (first, *rest)


def chain_vertices(
    lower: tuple[int, ...], upper: tuple[int, ...]
) -> tuple[tuple[int, ...], ...]:
    """Reconstruct exact vertices of a bounded chain order polytope."""

    closed = close_bounds(lower, upper)
    if closed is None:
        return ()
    lower, upper = closed
    result: set[tuple[int, ...]] = set()
    for lengths in compositions(len(lower)):
        blocks = []
        start = 0
        for length in lengths:
            end = start + length
            anchors = sorted(
                {value for index in range(start, end) for value in (lower[index], upper[index])}
            )
            blocks.append((start, end, anchors))
            start = end
        for selected in itertools.product(*(block[2] for block in blocks)):
            chain = tuple(
                value
                for (start, end, _anchors), value in zip(blocks, selected)
                for _ in range(end - start)
            )
            if (
                all(low <= value <= high for value, low, high in zip(chain, lower, upper))
                and all(left <= right for left, right in zip(chain, chain[1:]))
            ):
                result.add(chain)
    return tuple(sorted(result))


def profile_from_chain(chain: tuple[int, ...]) -> tuple[int, ...]:
    shifted = []
    previous = 0
    for value in chain:
        shifted.append(value - previous)
        previous = value
    shifted.append(FREE_MASS - previous)
    profile = tuple(value + 1 for value in shifted)
    if len(profile) != CLASSES or sum(profile) != PROFILE_MASS or min(profile) < 1:
        raise ValueError("cumulative vertex left the full-support profile domain")
    return profile


def parse_reference(
    value: Any, expected_role: str, context: dict[str, Any]
) -> Reference:
    if not isinstance(value, dict):
        raise ValueError("component reference must be an object")
    identifier = str(value.get("source_id", ""))
    index = parse_int(value.get("row"), "component row")
    key = (identifier, index)
    if key not in context["component_rows"]:
        raise ValueError(f"unmanifested component reference {identifier}:{index}")
    source = context["component_sources"][identifier]
    if source["component_role"] != expected_role:
        raise ValueError("component reference crosses marginal roles")
    if value.get("source_sha256") != source["sha256"]:
        raise ValueError("component source digest mismatch")
    _row, row_digest = context["component_rows"][key]
    if value.get("row_sha256") != row_digest:
        raise ValueError("component row digest mismatch")
    return key


def parse_marginal(value: Any, role: str, context: dict[str, Any]) -> Marginal:
    if not isinstance(value, list) or not value:
        raise ValueError(f"{role} marginal must be a nonempty array")
    result = []
    seen: set[Reference] = set()
    total = Fraction(0)
    for component in value:
        if not isinstance(component, dict):
            raise ValueError("marginal component must be an object")
        weight = parse_fraction(component.get("weight"), "marginal weight")
        if weight <= 0:
            raise ValueError("serialized marginal weights must be positive")
        reference = parse_reference(component.get("component"), role, context)
        if reference in seen:
            raise ValueError("marginal references must be unique")
        seen.add(reference)
        total += weight
        result.append((reference, weight))
    if total != 1:
        raise ValueError(f"{role} marginal weights must sum exactly to one")
    return tuple(result)


def parse_selector(value: Any, context: dict[str, Any]) -> tuple[Marginal, Marginal]:
    if not isinstance(value, dict) or value.get("kind") != SELECTOR_KIND:
        raise ValueError(f"factorized selector kind must be {SELECTOR_KIND!r}")
    if value.get("normalization") != NORMALIZATION_RULE:
        raise ValueError("factorized selector normalization rule changed")
    inner = parse_marginal(value.get("inner_marginal"), "inner", context)
    outer = parse_marginal(value.get("outer_marginal"), "outer", context)
    return inner, outer


def validate_manifest(
    manifest: dict[str, Any], manifest_path: Path,
    overrides: dict[str, Path] | None = None,
) -> dict[str, Any]:
    if manifest.get("schema") != MANIFEST_SCHEMA:
        raise ValueError("unsupported factorized manifest schema")
    parameters = manifest.get("parameters", {})
    if (parse_int(parameters.get("group_bits"), "group_bits"), parse_int(parameters.get("M"), "M")) != (
        GROUP_BITS,
        PROFILE_MASS,
    ):
        raise ValueError("factorized manifest parameter mismatch")
    contract = manifest.get("recombination_contract", {})
    expected = {
        "selector_kind": SELECTOR_KIND,
        "inner_interface": INNER_INTERFACE,
        "combination_rule": COMBINATION_RULE,
        "normalization": NORMALIZATION_RULE,
    }
    if contract != expected:
        raise ValueError("factorized recombination contract mismatch")
    overrides = overrides or {}
    sources: dict[str, dict[str, Any]] = {}
    rows: dict[Reference, tuple[dict[str, Any], str]] = {}
    for source in manifest.get("component_sources", []):
        if not isinstance(source, dict):
            raise ValueError("component source declaration must be an object")
        identifier = str(source.get("source_id", ""))
        role = source.get("component_role")
        if not identifier or identifier in sources or role not in ("inner", "outer"):
            raise ValueError("component source identifier or role is invalid")
        path = overrides.get(identifier)
        if path is None:
            locator = Path(str(source.get("path", "")))
            path = locator if locator.is_absolute() else manifest_path.parent / locator
        if sha256_path(path) != source.get("sha256"):
            raise ValueError(f"component source digest mismatch for {identifier!r}")
        artifact = json.loads(path.read_text(encoding="utf-8"))
        if artifact.get("schema") != source.get("schema"):
            raise ValueError("component source schema mismatch")
        declarations = source.get("rows")
        if not isinstance(declarations, list) or len(declarations) != len(artifact.get("rows", [])):
            raise ValueError("component row catalogue is incomplete")
        sources[identifier] = {**source, "path": path}
        for declaration in declarations:
            index = parse_int(declaration.get("row"), "component row index")
            if not 0 <= index < len(artifact["rows"]):
                raise ValueError("component row index is outside its source")
            row = artifact["rows"][index]
            digest = canonical_sha256(row)
            if declaration.get("row_sha256") != digest:
                raise ValueError("component row digest mismatch")
            if row.get("component_role") != role:
                raise ValueError("component row role mismatch")
            rows[(identifier, index)] = (row, digest)
    if not sources or {source["component_role"] for source in sources.values()} != {"inner", "outer"}:
        raise ValueError("manifest must bind both component roles")
    return {"component_sources": sources, "component_rows": rows, "hardened": {}}


def _require_float_vector(value: Any, length: int, label: str) -> np.ndarray:
    result = np.asarray(value, dtype=np.float64)
    if result.shape != (length,) or not np.all(np.isfinite(result)):
        raise ValueError(f"{label} must contain {length} finite binary64 values")
    return result


def harden_inner_component(name: str, row: dict[str, Any]) -> dict[str, Any]:
    parameters = row.get("parameters", {})
    fugacities = _require_float_vector(parameters.get("fugacities"), CLASSES, "inner fugacities")
    if np.any(fugacities < 0):
        raise ValueError("inner fugacities must be nonnegative")
    pole = float(parameters.get("pole"))
    if not math.isfinite(pole) or not 0 < pole < 1:
        raise ValueError("inner pole must lie in (0,1)")
    values = _require_float_vector(parameters.get("collatz_vector"), 65, "Collatz vector")
    if np.any(values <= 0):
        raise ValueError("Collatz vector must be positive")

    outward.load_native_point_caps = lambda: None
    outward._initialize_hardening_worker(GROUP_BITS)
    split_caps = outward._WORKER_SPLIT_CAPS
    if split_caps is None:
        raise RuntimeError("outward split caps were not initialized")
    histograms = outward._block_histograms_outward(fugacities)
    caps = outward._point_caps_outward(fugacities)
    image, statistics = outward_transport_image_decimal_upper(
        values, histograms, caps, split_caps, pole
    )
    value_fractions = [exact_float(float(value)) for value in values]
    eigenvalue = max(entry / value for entry, value in zip(image, value_fractions))
    domination = max(Fraction(1) / value for value in value_fractions)
    mgf = (
        log2_fraction(domination)
        + log2_fraction(eigenvalue).times_int(outward.INNER_BLOCKS)
        + log2_fraction(value_fractions[0])
    )
    constant = mgf - log2_fraction(exact_float(pole)).times_int(outward.D)
    charges = tuple(
        None if value == 0 else log2_fraction(exact_float(float(value)))
        for value in fugacities
    )
    return {"name": name, "role": "inner", "constant": constant, "charges": charges, "report": statistics}


def harden_outer_component(name: str, row: dict[str, Any]) -> dict[str, Any]:
    parameters = row.get("parameters", {})
    log_variables = _require_float_vector(
        parameters.get("log_variables"), CLASSES, "outer log variables"
    )
    if log_variables[0] != 0:
        raise ValueError("outer class-zero log variable must be exactly zero")
    coefficients = _require_float_vector(
        parameters.get("band_coefficients"), 3, "outer BL coefficients"
    )
    theta = float(parameters.get("pair_cauchy_theta"))
    if not math.isfinite(theta):
        raise ValueError("outer Cauchy split must be finite")
    constant, charges, report = outward.conditioned_row_exact_graph_outer_outward(
        GROUP_BITS,
        log_variables.tolist(),
        coefficients.tolist(),
        theta,
    )
    return {"name": name, "role": "outer", "constant": constant, "charges": charges, "report": report}


def harden_reference(reference: Reference, context: dict[str, Any]) -> dict[str, Any]:
    cached = context["hardened"]
    if reference in cached:
        return cached[reference]
    row, _digest = context["component_rows"][reference]
    name = f"{reference[0]}:{reference[1]}"
    role = context["component_sources"][reference[0]]["component_role"]
    hardened = (
        harden_inner_component(name, row)
        if role == "inner"
        else harden_outer_component(name, row)
    )
    cached[reference] = hardened
    return hardened


def scale_interval(value: Interval, weight: Fraction) -> Interval:
    return value * (Interval.exact(weight.numerator) / Interval.exact(weight.denominator))


def evaluate_affine(profile: Profile, component: dict[str, Any]) -> Interval:
    value = component["constant"]
    for count, charge in zip(profile, component["charges"]):
        if not count:
            continue
        if charge is None:
            raise ValueError("component has a zero variable on an active class")
        value = value - scale_interval(charge, count)
    return value


def evaluate_selector(
    profile: Profile, selector: tuple[Marginal, Marginal], context: dict[str, Any]
) -> Interval:
    if len(profile) != CLASSES or sum(profile, Fraction(0)) != PROFILE_MASS:
        raise ValueError("profile has the wrong dimension or mass")
    total = Interval.exact(0)
    for marginal in selector:
        subtotal = Interval.exact(0)
        for reference, weight in marginal:
            subtotal = subtotal + scale_interval(
                evaluate_affine(profile, harden_reference(reference, context)), weight
            )
        total = total + subtotal
    # Import lazily so this module's synthetic parser tests do not initialize
    # the complete support verifier or any discovery code.
    from certify_packet_group_g8_support_shard import fractional_normalization

    return total - fractional_normalization(profile)


def _node_map(certificate: dict[str, Any]) -> dict[str, dict[str, Any]]:
    raw = certificate.get("nodes")
    values = list(raw.values()) if isinstance(raw, dict) else raw
    if not isinstance(values, list):
        raise ValueError("factorized certificate nodes must be an array or object")
    result = {}
    for node in values:
        if not isinstance(node, dict):
            raise ValueError("factorized node must be an object")
        identifier = str(node.get("node_id", ""))
        if not identifier or identifier in result:
            raise ValueError("duplicate or empty factorized node identifier")
        result[identifier] = node
    return result


def replay_cumulative_geometry(
    certificate: dict[str, Any], *, require_selectors: bool = True,
    expected_leaves: int = EXPECTED_LEAVES,
) -> dict[str, Any]:
    """Replay the full-support h=2 roots and every cumulative-prefix split."""

    if certificate.get("schema") != CERTIFICATE_SCHEMA:
        raise ValueError("unsupported factorized full-support certificate schema")
    nodes = _node_map(certificate)
    roots = sorted(
        (node for node in nodes.values() if node.get("parent_id") is None),
        key=lambda node: parse_int(node.get("h2_cell"), "h2 cell"),
    )
    owners = tuple(itertools.combinations_with_replacement(range(H2_BINS), CLASSES - 1))
    if len(roots) != EXPECTED_ROOTS or len(owners) != EXPECTED_ROOTS:
        raise ValueError("factorized geometry does not contain all 165 h2 roots")
    bins = balanced_bins()
    visited: set[str] = set()
    leaves: list[dict[str, Any]] = []

    def visit(identifier: str, expected: tuple[tuple[int, ...], tuple[int, ...]]) -> int:
        if identifier in visited or identifier not in nodes:
            raise ValueError("factorized tree has a cycle, duplicate parent, or missing child")
        visited.add(identifier)
        node = nodes[identifier]
        lower = tuple(parse_int(value, "lower bound") for value in node.get("lower", []))
        upper = tuple(parse_int(value, "upper bound") for value in node.get("upper", []))
        closed = close_bounds(lower, upper)
        if closed is None or closed != expected:
            raise ValueError("stored cumulative bounds do not match exact ownership")
        count = exact_chain_count(*closed)
        if count <= 0 or parse_int(node.get("exact_count"), "node exact count") != count:
            raise ValueError("factorized node exact count mismatch")
        state = node.get("state")
        if state == "CERTIFIED_LEAF":
            vertices = chain_vertices(*closed)
            if not vertices:
                raise ValueError("certified cumulative leaf has no exact vertices")
            if require_selectors and not isinstance(node.get("selector"), dict):
                raise ValueError("certified cumulative leaf lacks a factorized selector")
            leaves.append({"node": node, "vertices": tuple(profile_from_chain(row) for row in vertices)})
            return count
        if state != "SPLIT":
            raise ValueError("complete factorized tree contains a nonterminal unresolved state")
        split = node.get("split", {})
        coordinate = parse_int(split.get("coordinate"), "cumulative split coordinate")
        threshold = parse_int(split.get("threshold"), "cumulative split threshold")
        coefficients = split.get("profile_coefficients")
        expected_coefficients = [1 if index <= coordinate else 0 for index in range(CLASSES)]
        expected_profile_threshold = threshold + coordinate + 1
        if coefficients != expected_coefficients or parse_int(
            split.get("profile_threshold"), "profile split threshold"
        ) != expected_profile_threshold:
            raise ValueError("cumulative split and exact profile ownership disagree")
        left_bounds = child_bounds(*closed, coordinate, threshold, True)
        right_bounds = child_bounds(*closed, coordinate, threshold, False)
        if left_bounds is None or right_bounds is None:
            raise ValueError("cumulative split has an empty child")
        left = str(node.get("left", ""))
        right = str(node.get("right", ""))
        left_count = visit(left, left_bounds)
        right_count = visit(right, right_bounds)
        if left_count + right_count != count:
            raise ValueError("exact cumulative child counts do not conserve their parent")
        return count

    root_total = 0
    for ordinal, (root, owner) in enumerate(zip(roots, owners)):
        if parse_int(root.get("h2_cell"), "h2 cell") != ordinal or root.get("h2_bin_indices") != list(owner):
            raise ValueError("h2 root owner catalogue is noncanonical")
        expected = (
            tuple(bins[index][0] for index in owner),
            tuple(bins[index][1] for index in owner),
        )
        root_total += visit(str(root["node_id"]), expected)
    if visited != set(nodes):
        raise ValueError("factorized certificate contains unreachable nodes")
    if root_total != math.comb(PROFILE_MASS - 1, CLASSES - 1):
        raise ValueError("h2 root counts do not recover the full-support census")
    if len(leaves) != expected_leaves:
        raise ValueError(f"factorized tree has {len(leaves)} leaves, expected {expected_leaves}")
    return {"root_count": root_total, "root_count_count": len(roots), "leaf_count": len(leaves), "leaves": leaves}


def synthetic_self_test() -> None:
    inner_source = {
        "component_role": "inner",
        "sha256": "1" * 64,
    }
    outer_source = {
        "component_role": "outer",
        "sha256": "2" * 64,
    }
    inner_row = {"component_role": "inner", "parameters": {}}
    outer_row = {"component_role": "outer", "parameters": {}}
    context = {
        "component_sources": {"i": inner_source, "o": outer_source},
        "component_rows": {
            ("i", 0): (inner_row, canonical_sha256(inner_row)),
            ("o", 0): (outer_row, canonical_sha256(outer_row)),
        },
        "hardened": {
            ("i", 0): {
                "role": "inner",
                "constant": Interval.exact(7),
                "charges": tuple(Interval.exact(0) for _ in range(CLASSES)),
            },
            ("o", 0): {
                "role": "outer",
                "constant": Interval.exact(11),
                "charges": tuple(Interval.exact(0) for _ in range(CLASSES)),
            },
        },
    }
    selector_wire = {
        "kind": SELECTOR_KIND,
        "normalization": NORMALIZATION_RULE,
        "inner_marginal": [{
            "weight": "1/1",
            "component": {"source_id": "i", "source_sha256": "1" * 64, "row": 0,
                          "row_sha256": canonical_sha256(inner_row)},
        }],
        "outer_marginal": [{
            "weight": "1/1",
            "component": {"source_id": "o", "source_sha256": "2" * 64, "row": 0,
                          "row_sha256": canonical_sha256(outer_row)},
        }],
    }
    selector = parse_selector(selector_wire, context)
    profile = tuple([Fraction(1)] * 8 + [Fraction(PROFILE_MASS - 8)])
    from certify_packet_group_g8_support_shard import fractional_normalization

    normalization = fractional_normalization(profile)
    value = evaluate_selector(profile, selector, context)
    expected = Interval.exact(18) - normalization
    if value != expected:
        raise AssertionError("factorized selector did not subtract normalization exactly once")

    invalid = json.loads(json.dumps(selector_wire))
    invalid["outer_marginal"][0]["component"]["source_id"] = "i"
    invalid["outer_marginal"][0]["component"]["source_sha256"] = "1" * 64
    invalid["outer_marginal"][0]["component"]["row_sha256"] = canonical_sha256(inner_row)
    try:
        parse_selector(invalid, context)
    except ValueError:
        pass
    else:
        raise AssertionError("cross-role component unexpectedly passed")

    toy_lower = (0,) * (CLASSES - 1)
    toy_upper = (2,) * (CLASSES - 1)
    parent = exact_chain_count(toy_lower, toy_upper)
    left = child_bounds(toy_lower, toy_upper, 3, 0, True)
    right = child_bounds(toy_lower, toy_upper, 3, 0, False)
    if left is None or right is None or exact_chain_count(*left) + exact_chain_count(*right) != parent:
        raise AssertionError("exact cumulative split did not conserve its parent")


if __name__ == "__main__":
    synthetic_self_test()
    print("factorized_selector_normalization_once=PASS")
    print("factorized_selector_cross_role_rejection=PASS")
    print("cumulative_exact_split_conservation=PASS")
    print("status=PASS_G8_FACTORIZED_SELECTOR_SYNTHETIC")
