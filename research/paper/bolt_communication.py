"""Reproducible Bolt-max communication scenario, not a serialized Bolt proof.

The outer calculation follows Hypercat's results/bolt-proof-size/account.py.
Inner byte counts come from the same verified proxy runs as our time estimate.
See artifact/BOLT_PCS_ESTIMATE.md for scope, receipts, and omitted components.
No benchmark runs and no files are written by this module.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

MIB = 2**20
PROXY_CSV = "tools/standalone/opening-results/27.csv"
PROXY_CSV_SHA256 = "e421aa8167a9570a533cea846a5f4d5202985b39fe3ce68bc91bb8f8e3164723"
PROXY_BYTES = {"inner_message_proxy": 0x6E270,
               "inner_syndrome_proxy": 0x2B430}


def queries(distance, divisor, bits=100):
    return math.ceil(bits * math.log(2) / -math.log1p(-distance / divisor))


def sibling_count(indices, depth):
    """Digest count for deduplicated queries, including odd-tree dummy siblings."""
    active = set(indices)
    count = 0
    for _ in range(depth):
        count += sum((i ^ 1) not in active for i in active)
        active = {i >> 1 for i in active}
    return count


def expected_siblings(n, q):
    """Exact expectation (floating-point evaluation), uniform distinct queries."""
    if n < 1 or n & (n - 1) or not 1 <= q <= n:
        raise ValueError("Require a power-of-two leaf count and 1 <= q <= n")
    if n == 1:
        return 0.0
    occupied = []
    for h in range(1, n.bit_length() - 1):
        size = 1 << h
        if size > n - q:
            occupied.append(n / size)
        else:
            # P(empty) = binom(n-size,q) / binom(n,q).
            log_empty = math.fsum(math.log1p(-size / (n - i)) for i in range(q))
            occupied.append((n / size) * -math.expm1(log_empty))
    return math.fsum(occupied) + 2 - q


def outer_payload(log_symbols=27):
    k = 1 << (log_symbols - 7)
    qx, qy = queries(.013, 3), queries(.5, 2)
    hashes = expected_siblings(k, qx) + expected_siblings(k // 4, qy) + 2
    leaves = (qx + qy) * 128 * 4
    return {"message_queries": qx, "syndrome_queries": qy,
            "opened_column_bytes": leaves, "expected_merkle_bytes": 32 * hashes,
            "expected_outer_bytes": leaves + 32 * hashes}


def communication_estimate():
    outer = outer_payload()
    # Same pair of complete, unconstrained Ligerito proxies used for timing.
    # They include roots, query data, Merkle proofs, and their own transcripts.
    # Their use here is a model substitution, not a Bolt protocol implementation.
    components = {"outer_columns": outer["opened_column_bytes"],
                  "outer_merkle": outer["expected_merkle_bytes"],
                  **PROXY_BYTES,
                  "row_evaluations": 128 * 16,
                  "two_sumchecks": 2 * 20 * 3 * 16,
                  "scalar_allowance": 8 * 16}
    total = math.fsum(components.values())
    large_outer = outer_payload(30)["expected_outer_bytes"]
    # A secondary, paper-anchored check, NOT an independent full-proof receipt.
    # Table 1 reports 12.92 MB; preserve both possible byte-unit conventions.
    drop = large_outer - outer["expected_outer_bytes"]
    anchored = [(12.92 * unit - drop) / MIB for unit in (10**6, MIB)]
    return {"scope": "Amortized-limit communication proxy; not a bound or measured Bolt proof",
            "input_bytes": 2**29, "components_bytes": components,
            "query_model": "Uniform distinct queries in each of two independent pools",
            **outer, "estimated_bytes": total, "estimated_mib": total / MIB,
            "paper_anchored_mib_scenarios": anchored}


def validate_proxy_receipt(source):
    path = Path(source) / PROXY_CSV
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != PROXY_CSV_SHA256:
        raise ValueError(f"Changed proxy receipt: {path}")
    rows = list(csv.DictReader(raw.decode().splitlines()))
    for phase, expected in PROXY_BYTES.items():
        subset = [r for r in rows if r["phase"] == phase]
        if ([int(r["rep"]) for r in subset] != list(range(6))
                or any(int(r["log_n"]) != 27 for r in subset)
                or any(int(r["checksum"], 16) != expected for r in subset)):
            raise ValueError(f"Invalid size receipts for {phase}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, help="Optional pinned bolt-standalone checkout")
    args = parser.parse_args()
    if args.source:
        validate_proxy_receipt(args.source)
    print(json.dumps(communication_estimate(), indent=2))
