#!/usr/bin/env python3
"""Validate the schema and support coverage of a binary biregular EC certificate.

This checker uses only the Python standard library.  It does not import the
optimizer, the transfer-matrix implementation, SciPy, or Arb.  Its purpose is
to reject malformed parameters, markers, and support partitions before the
numerical verifier evaluates the proved bounds.
"""

from __future__ import annotations

import argparse
import json
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Any


DEGREE_FIVE_SCHEMA = "binary-biregular-ec-degree-five-v1"
ODD_DEGREE_SCHEMA = "binary-biregular-ec-odd-degree-v1"
SCHEMA = DEGREE_FIVE_SCHEMA


def _require_keys(value: dict[str, Any], expected: set[str], where: str) -> None:
    actual = set(value)
    if actual != expected:
        missing = sorted(expected - actual)
        extra = sorted(actual - expected)
        raise ValueError(f"{where} keys differ: missing={missing}, extra={extra}")


def _integer(value: Any, where: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{where} must be an integer")
    return value


def _marker(value: Any, where: str, *, at_most_one: bool) -> None:
    if not isinstance(value, str):
        raise ValueError(f"{where} must be a decimal string")
    try:
        parsed = Decimal(value)
    except InvalidOperation as error:
        raise ValueError(f"{where} is not a decimal number") from error
    if not parsed.is_finite() or parsed <= 0:
        raise ValueError(f"{where} must be finite and positive")
    if at_most_one and parsed > 1:
        raise ValueError(f"{where} must be at most one")


def validate_certificate_structure(certificate: dict[str, Any]) -> None:
    """Validate all non-numerical-proof obligations in ``certificate``."""
    if not isinstance(certificate, dict):
        raise ValueError("certificate root must be an object")
    _require_keys(
        certificate,
        {
            "schema",
            "parameters",
            "verification",
            "exact_weights",
            "outer_blocks",
            "central_blocks",
            "full_support",
        },
        "certificate",
    )
    if certificate["schema"] not in {DEGREE_FIVE_SCHEMA, ODD_DEGREE_SCHEMA}:
        raise ValueError(f"unsupported certificate schema: {certificate['schema']!r}")

    parameters = certificate["parameters"]
    if not isinstance(parameters, dict):
        raise ValueError("parameters must be an object")
    _require_keys(
        parameters,
        {
            "k",
            "n",
            "cutoff",
            "left_degree",
            "right_degree",
            "memory",
            "target_bits",
        },
        "parameters",
    )
    k = _integer(parameters["k"], "parameters.k")
    n = _integer(parameters["n"], "parameters.n")
    cutoff = _integer(parameters["cutoff"], "parameters.cutoff")
    left_degree = _integer(parameters["left_degree"], "parameters.left_degree")
    right_degree = _integer(parameters["right_degree"], "parameters.right_degree")
    memory = _integer(parameters["memory"], "parameters.memory")
    target_bits = _integer(parameters["target_bits"], "parameters.target_bits")
    if min(k, n, left_degree, right_degree, memory) <= 0:
        raise ValueError("lengths, degrees, and memory must be positive")
    if right_degree % 2 != 1 or right_degree < 3:
        raise ValueError("the local-limit proof requires odd right degree at least three")
    if certificate["schema"] == DEGREE_FIVE_SCHEMA and right_degree != 5:
        raise ValueError("the legacy degree-five schema requires right degree five")
    if k % right_degree:
        raise ValueError("right_degree must divide k")
    if n != left_degree * (k // right_degree):
        raise ValueError("the two expander degrees do not balance")
    if n != 2 * k:
        raise ValueError("this schema requires rate one half")
    if k % 2 != 1:
        raise ValueError("the complement partition requires odd k")
    if not 0 <= cutoff < n // 2:
        raise ValueError("cutoff must lie in [0,n/2)")
    if target_bits < 0:
        raise ValueError("target_bits must be nonnegative")

    verification = certificate["verification"]
    if not isinstance(verification, dict):
        raise ValueError("verification must be an object")
    _require_keys(verification, {"precision_bits"}, "verification")
    precision_bits = _integer(
        verification["precision_bits"], "verification.precision_bits"
    )
    if precision_bits < 64:
        raise ValueError("verification precision must be at least 64 bits")

    covered: list[tuple[int, int, str]] = []
    exact = certificate["exact_weights"]
    if not isinstance(exact, list) or not exact:
        raise ValueError("exact_weights must be a nonempty list")
    for index, item in enumerate(exact):
        if not isinstance(item, dict):
            raise ValueError(f"exact_weights[{index}] must be an object")
        _require_keys(item, {"r", "z"}, f"exact_weights[{index}]")
        r = _integer(item["r"], f"exact_weights[{index}].r")
        if r != index + 1:
            raise ValueError("exact_weights must list a contiguous prefix from one")
        _marker(item["z"], f"exact_weights[{index}].z", at_most_one=True)
        covered.append((r, r, f"exact weight {r}"))
    central_start = k // right_degree
    if len(exact) >= central_start:
        raise ValueError("exact_weights must end before the central range")

    outer = certificate["outer_blocks"]
    if not isinstance(outer, list) or not outer:
        raise ValueError("outer_blocks must be a nonempty list")
    previous_outer_hi = 0
    for index, item in enumerate(outer):
        if not isinstance(item, dict):
            raise ValueError(f"outer_blocks[{index}] must be an object")
        _require_keys(item, {"lo", "hi", "x", "z"}, f"outer_blocks[{index}]")
        lo = _integer(item["lo"], f"outer_blocks[{index}].lo")
        hi = _integer(item["hi"], f"outer_blocks[{index}].hi")
        if lo <= previous_outer_hi:
            raise ValueError("outer_blocks must be ordered and disjoint")
        previous_outer_hi = hi
        _marker(item["x"], f"outer_blocks[{index}].x", at_most_one=False)
        _marker(item["z"], f"outer_blocks[{index}].z", at_most_one=True)
        covered.append((lo, hi, f"outer block {index}"))

    central = certificate["central_blocks"]
    if not isinstance(central, list) or not central:
        raise ValueError("central_blocks must be a nonempty list")
    half = (k - 1) // 2
    expected_central = central_start
    for index, item in enumerate(central):
        if not isinstance(item, dict):
            raise ValueError(f"central_blocks[{index}] must be an object")
        _require_keys(
            item, {"lo", "hi", "include_complement"}, f"central_blocks[{index}]"
        )
        lo = _integer(item["lo"], f"central_blocks[{index}].lo")
        hi = _integer(item["hi"], f"central_blocks[{index}].hi")
        if item["include_complement"] is not True:
            raise ValueError("every central block must include its complement")
        if not central_start <= lo <= hi <= half:
            raise ValueError(
                f"central block {index} lies outside [k/d_R,(k-1)/2]"
            )
        if lo != expected_central:
            raise ValueError("central_blocks must be ordered and contiguous")
        expected_central = hi + 1
        covered.append((lo, hi, f"central block {index}"))
        covered.append((k - hi, k - lo, f"central complement {index}"))
    if expected_central != half + 1:
        raise ValueError("central_blocks do not reach (k-1)/2")

    full = certificate["full_support"]
    if not isinstance(full, dict):
        raise ValueError("full_support must be an object")
    _require_keys(full, {"r", "z"}, "full_support")
    full_r = _integer(full["r"], "full_support.r")
    _marker(full["z"], "full_support.z", at_most_one=True)
    if full_r != k:
        raise ValueError("full_support.r must equal k")
    covered.append((full_r, full_r, "full support"))

    expected = 1
    for lo, hi, label in sorted(covered):
        if lo != expected:
            relation = "overlaps" if lo < expected else "leaves a gap before"
            raise ValueError(f"{label} {relation} support {expected}")
        if hi < lo or hi > k:
            raise ValueError(f"{label} has invalid endpoints [{lo},{hi}]")
        expected = hi + 1
    if expected != k + 1:
        raise ValueError(f"support coverage stops at {expected - 1}, expected {k}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("certificate", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    certificate = json.loads(args.certificate.read_text(encoding="utf-8"))
    validate_certificate_structure(certificate)
    parameters = certificate["parameters"]
    print(
        "valid certificate: "
        f"supports 1..{parameters['k']}, "
        f"{len(certificate['outer_blocks'])} outer blocks, "
        f"{len(certificate['central_blocks'])} complemented central blocks"
    )


if __name__ == "__main__":
    main()
