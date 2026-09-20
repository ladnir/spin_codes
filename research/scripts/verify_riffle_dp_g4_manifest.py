#!/usr/bin/env python3
"""Verify the frozen Riffle DP g=4 construction manifest.

This verifier checks the exact derived parameters, the field coefficient
schedule, the local BCH generator rows, and the hashes of the implementation
sources. It does not turn the finite-seed implementation sampler into the
uniform distribution declared by the proof ensemble.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = REPOSITORY / "explorations" / "riffle_dp_g4_construction_manifest.json"
DEFAULT_LIBOTE_ROOT = Path(r"C:\Users\peter\.codex\worktrees\permute-conv-code\libOTe")
MASK64 = (1 << 64) - 1
FIELD_REDUCTION = 0x1B
FIELD_ORDER_FACTORS = (3, 5, 17, 257, 641, 65537, 6700417)


def fail(message: str) -> None:
    raise RuntimeError(message)


def multiply(left: int, right: int) -> int:
    result = 0
    for _ in range(64):
        if right & 1:
            result ^= left
        carry = left >> 63
        left = (left << 1) & MASK64
        if carry:
            left ^= FIELD_REDUCTION
        right >>= 1
    return result


def power(value: int, exponent: int) -> int:
    result = 1
    while exponent:
        if exponent & 1:
            result = multiply(result, value)
        value = multiply(value, value)
        exponent >>= 1
    return result


def generator_rows(generator: int) -> list[int]:
    rows = []
    for row in range(64):
        low = (generator << row) & MASK64
        high = (generator >> (64 - row) if row else 0) | (1 << 63)
        rows.append(low | (high << 64))
    return rows


def rows_digest(rows: list[int]) -> str:
    material = b"".join(row.to_bytes(16, "little") for row in rows)
    return hashlib.sha256(material).hexdigest()


def source_digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require_source(pattern: str, source: str, description: str) -> None:
    if re.search(pattern, source, re.MULTILINE) is None:
        fail(f"source assertion failed: {description}")


def verify(manifest: dict, libote_root: Path) -> None:
    claim = manifest["claim_under_study"]
    outer = manifest["outer_code"]
    local = manifest["local_code"]
    layout = manifest["layout"]
    inner = manifest["inner_map"]
    terminal = manifest["derived_terminal_event"]

    data_symbols = claim["message_dimension_bits"] // outer["symbol_bits"]
    if data_symbols != outer["data_symbol_count"]:
        fail("data-symbol count is inconsistent")
    outer_words = data_symbols + outer["parity_symbol_count"]
    if outer_words != layout["outer_word_count"]:
        fail("outer-word count is inconsistent")
    binary_length = outer_words * local["length"] - layout["punctured_coordinates"]
    if binary_length != claim["binary_length"]:
        fail("binary length is inconsistent")
    packet_count, packet_remainder = divmod(binary_length, layout["packet_width_bits"])
    if packet_remainder or packet_count != layout["packet_count"]:
        fail("packet count is inconsistent")
    node_count, node_remainder = divmod(binary_length, inner["node_width_bits"])
    if node_remainder or node_count != inner["node_count"]:
        fail("inner-node count is inconsistent")
    threshold = 9 * binary_length // 100
    if threshold != claim["failure_weight_inclusive"]:
        fail("failure threshold is not floor(0.09 N)")
    terminal_nodes = threshold // inner["node_width_bits"]
    if terminal_nodes != terminal["terminal_node_count"]:
        fail("terminal-node count is inconsistent")
    terminal_packets = terminal_nodes * inner["packets_per_node"]
    if terminal_packets != terminal["terminal_packet_count"]:
        fail("terminal-packet count is inconsistent")
    if local["length"] // layout["packet_width_bits"] != layout["packets_per_local_word"]:
        fail("packets per local word are inconsistent")

    order = (1 << 64) - 1
    if power(2, order) != 1:
        fail("x does not have multiplicative order dividing 2^64-1")
    for prime in FIELD_ORDER_FACTORS:
        if power(2, order // prime) == 1:
            fail(f"x is not primitive; failed factor {prime}")
    coefficients = []
    value = 1
    for _ in range(data_symbols):
        coefficients.append(value)
        value = multiply(value, 2)
    if 0 in coefficients or len(set(coefficients)) != data_symbols:
        fail("coefficient schedule is not distinct and nonzero")

    generator = int(local["cyclic_generator_polynomial_hex"], 16)
    if rows_digest(generator_rows(generator)) != local["generator_rows_sha256_little_endian_128"]:
        fail("local BCH generator-row digest is inconsistent")

    binding = manifest["source_binding"]
    for relative, expected in binding["files"].items():
        path = libote_root / relative
        if not path.is_file():
            fail(f"bound source file does not exist: {path}")
        actual = source_digest(path)
        if actual != expected:
            fail(f"source hash mismatch for {relative}: {actual} != {expected}")

    verification = manifest["verification"]
    test_source = REPOSITORY / verification["equivalence_test_source"]
    if source_digest(test_source) != verification["equivalence_test_source_sha256"]:
        fail("implementation-equivalence test source hash mismatch")

    dp_source = (libote_root / "libOTe/Tools/RiffleCode/RiffleDoubleParityG4.h").read_text()
    inner_source = (libote_root / "libOTe/Tools/RiffleCode/RifflePacketInnerChain.h").read_text()
    require_source(r"packetWidth\s*=\s*4", dp_source, "packet width four")
    require_source(r"fieldReduction\s*=\s*0x1b", dp_source, "GF(2^64) reduction")
    require_source(r"mInner\.initUniform\(codeSize", dp_source, "finite-seed setup adapter")
    require_source(r"mInner\.initOrder\(codeSize", dp_source, "explicit-order theorem interface")
    require_source(r"void\s+initOrder\(", inner_source, "explicit packet-order initializer")
    if "composeNodePacketOrder" in dp_source:
        fail("frozen candidate unexpectedly composes per-node packet orders")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--libote-root", type=Path, default=DEFAULT_LIBOTE_ROOT)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    verify(manifest, args.libote_root)
    print("candidate=Riffle DP g=4")
    print("revision=g0-v1")
    print("status=EXACT_CONSTRUCTION_MANIFEST_VERIFIED")
    print("implementation_sampler=FINITE_SEED_ADAPTER_NOT_PROOF_DISTRIBUTION")


if __name__ == "__main__":
    main()
