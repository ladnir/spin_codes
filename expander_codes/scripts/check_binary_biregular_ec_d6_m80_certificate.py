#!/usr/bin/env python3
"""Check the structure of a degree-three binary biregular EC certificate.

This checker uses only the Python standard library.  It verifies parameters,
decimal marker domains, and the partition of all nonzero message supports.
It does not import the optimizer, transfer matrices, SciPy, or Arb.
"""

from __future__ import annotations

import argparse
import json
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Any


SCHEMA = "binary-biregular-ec-degree-three-v1"


def _require_keys(value: dict[str, Any], expected: set[str], where: str) -> None:
    actual = set(value)
    if actual != expected:
        raise ValueError(
            f"{where} keys differ: "
            f"missing={sorted(expected - actual)}, extra={sorted(actual - expected)}"
        )


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
    """Check every non-numerical obligation in the certificate."""
    if not isinstance(certificate, dict):
        raise ValueError("certificate root must be an object")
    _require_keys(
        certificate,
        {
            "schema",
            "parameters",
            "verification",
            "exact_blocks",
            "outer_blocks",
            "central_blocks",
            "full_support",
        },
        "certificate",
    )
    if certificate["schema"] != SCHEMA:
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
    if (left_degree, right_degree) != (6, 3):
        raise ValueError("this schema requires left/right degrees 6/3")
    if k % right_degree:
        raise ValueError("right_degree must divide k")
    if n != left_degree * (k // right_degree) or n != 2 * k:
        raise ValueError("parameters must define a rate-one-half biregular graph")
    if k % 2 != 1:
        raise ValueError("the complemented central partition requires odd k")
    if not 0 <= cutoff < n // 2:
        raise ValueError("cutoff must lie in [0,n/2)")
    if target_bits < 0:
        raise ValueError("target_bits must be nonnegative")

    verification = certificate["verification"]
    if not isinstance(verification, dict):
        raise ValueError("verification must be an object")
    _require_keys(verification, {"precision_bits"}, "verification")
    if _integer(verification["precision_bits"], "verification.precision_bits") < 64:
        raise ValueError("verification precision must be at least 64 bits")

    covered: list[tuple[int, int, str]] = []
    exact = certificate["exact_blocks"]
    if not isinstance(exact, list) or not exact:
        raise ValueError("exact_blocks must be a nonempty list")
    expected_exact = 1
    for index, item in enumerate(exact):
        if not isinstance(item, dict):
            raise ValueError(f"exact_blocks[{index}] must be an object")
        _require_keys(item, {"lo", "hi", "z"}, f"exact_blocks[{index}]")
        lo = _integer(item["lo"], f"exact_blocks[{index}].lo")
        hi = _integer(item["hi"], f"exact_blocks[{index}].hi")
        if lo != expected_exact or hi < lo:
            raise ValueError("exact_blocks must form a contiguous prefix from one")
        _marker(item["z"], f"exact_blocks[{index}].z", at_most_one=True)
        covered.append((lo, hi, f"exact block {index}"))
        expected_exact = hi + 1

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
        if lo <= previous_outer_hi or hi < lo:
            raise ValueError("outer_blocks must be ordered and disjoint")
        previous_outer_hi = hi
        _marker(item["x"], f"outer_blocks[{index}].x", at_most_one=False)
        _marker(item["z"], f"outer_blocks[{index}].z", at_most_one=True)
        covered.append((lo, hi, f"outer block {index}"))

    central = certificate["central_blocks"]
    if not isinstance(central, list) or not central:
        raise ValueError("central_blocks must be a nonempty list")
    half = (k - 1) // 2
    expected_central: int | None = None
    for index, item in enumerate(central):
        if not isinstance(item, dict):
            raise ValueError(f"central_blocks[{index}] must be an object")
        _require_keys(
            item, {"lo", "hi", "include_complement"},
            f"central_blocks[{index}]",
        )
        lo = _integer(item["lo"], f"central_blocks[{index}].lo")
        hi = _integer(item["hi"], f"central_blocks[{index}].hi")
        if item["include_complement"] is not True:
            raise ValueError("every central block must include its complement")
        if not 1 <= lo <= hi <= half:
            raise ValueError("central block lies outside the low half")
        if expected_central is not None and lo != expected_central:
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
    if _integer(full["r"], "full_support.r") != k:
        raise ValueError("full_support.r must equal k")
    _marker(full["z"], "full_support.z", at_most_one=True)
    covered.append((k, k, "full support"))

    expected = 1
    for lo, hi, label in sorted(covered):
        if lo != expected:
            relation = "overlaps" if lo < expected else "leaves a gap before"
            raise ValueError(f"{label} {relation} support {expected}")
        if hi > k:
            raise ValueError(f"{label} exceeds support {k}")
        expected = hi + 1
    if expected != k + 1:
        raise ValueError(f"support coverage stops at {expected - 1}, expected {k}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("certificate", type=Path)
    args = parser.parse_args()
    certificate = json.loads(args.certificate.read_text(encoding="utf-8"))
    validate_certificate_structure(certificate)
    print(
        "valid degree-three certificate structure: "
        f"supports 1..{certificate['parameters']['k']}, "
        f"{len(certificate['exact_blocks'])} exact blocks, "
        f"{len(certificate['outer_blocks'])} outer blocks, "
        f"{len(certificate['central_blocks'])} complemented central blocks"
    )


if __name__ == "__main__":
    main()
