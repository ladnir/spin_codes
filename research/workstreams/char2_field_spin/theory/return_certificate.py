"""Exact two-scan return bound for a characteristic-two GRS[8,4] outer.

This adapts the companion's characteristic-independent return-rank bound.
The encoder uses input-weighted accumulators and any Q distinct nonzero
coefficients. Discovery uses floating proposals; verify() checks rational
endpoint bounds, exact sparse sums, and complete integer-domain coverage.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys


REPO = Path(__file__).resolve().parents[4]
COMPANION = REPO.parent / "gen-BAA" / "field-universality" / "expander-lift"
DEPENDENCIES = (
    "distance_ladder.py", "signed_scan_certificate.py", "signed_scan_design.py",
    "block_support_precode.py", "repeat_accumulate_precode.py",
    "sparse_mixing_optimization.py",
)


def load_bounds(companion: Path, n: int):
    sys.path.insert(0, str(companion.resolve()))
    import distance_ladder as bounds
    import signed_scan_certificate as intervals

    for module in (bounds, intervals):
        module.N, module.K, module.L = n, n // 2, n // 8
    return bounds


def verify(data: dict, companion: Path = COMPANION) -> dict:
    n, q, e = data["n"], data["alphabet_size"], data["forbidden_weight"]
    assert n > 0 and n % 8 == 0 and data["k"] == n // 2
    assert data["local_code"] == [8, 4, 5] and data["scans"] == 2
    assert q >= 2 and 128 < 4 * e < n and data["sparse_cutoff"] == 128
    assert data["scan"] == "independent input-weighted accumulator"
    assert data["routing"] == "two independent uniform symbol permutations"
    bounds = load_bounds(companion, n)
    t, sparse = bounds.weighted_sparse(e, q, F(data["sparse_tau"]), F(data["sparse_t"]))
    assert str(t) == data["sparse_t"] and str(sparse) == data["sparse_failure"]
    following = 129
    for record in data["intervals"]:
        assert record["lower"] == following
        actual = bounds.interval_certificate(
            record["lower"], record["upper"],
            *(F(record[name]) for name in ("x", "y", "v", "z", "tau", "t")), q, e,
        )
        assert actual == record and F(actual["interval_log2_upper"]) < -80
        following = record["upper"] + 1
    assert following == 4 * e + 1
    total = sparse + F(len(data["intervals"]), 2**80)
    failure_log = bounds.cached_logs(total)[1]
    assert total < F(1, 2**data["required_margin"])
    return {
        "n": n, "k": n // 2, "local_code": [8, 4, 5], "alphabet_size": q,
        "minimum_distance_at_least": e + 1,
        "relative_distance": str(F(e + 1, n)),
        "failure_log2_upper": str(failure_log),
        "margin_bits": -float(failure_log),
        "sparse_log2_upper": str(bounds.cached_logs(sparse)[1]),
        "bulk_intervals": len(data["intervals"]),
        "complete_distance_certificate": True,
        "field_condition": "a field containing at least eight elements and Q distinct nonzero elements",
        "scope": "All nonzero messages; independent uniform permutations and independent small coefficients; no odd-characteristic refinement.",
        "verifier_sha256": {
            name: hashlib.sha256((companion / name).read_bytes()).hexdigest()
            for name in DEPENDENCIES
        },
    }


def generate(n: int, q: int, e: int, tau: int, margin: int, companion: Path) -> dict:
    bounds = load_bounds(companion, n)
    t, sparse = bounds.weighted_sparse(e, q, F(tau))
    print("SPARSE", -float(bounds.cached_logs(sparse)[1]), flush=True)
    if sparse >= F(1, 2**margin):
        raise ValueError("Sparse bound exceeds the requested margin")
    intervals, lower = [], 129
    while lower <= 4 * e:
        upper = min(4 * e, max(lower, int(lower * 1.08)))
        record = bounds.propose(lower, upper, q, e)
        while F(record["interval_log2_upper"]) >= -80 and upper > lower:
            upper = (lower + upper) // 2
            record = bounds.propose(lower, upper, q, e)
        if F(record["interval_log2_upper"]) >= -80:
            raise ValueError(f"Unclosed bulk source weight {lower}")
        intervals.append(record)
        lower = upper + 1
        if len(intervals) % 20 == 0:
            print("BULK", len(intervals), "through", lower - 1, flush=True)
    return {
        "n": n, "k": n // 2, "local_code": [8, 4, 5], "scans": 2,
        "alphabet_size": q, "forbidden_weight": e,
        "scan": "independent input-weighted accumulator",
        "routing": "two independent uniform symbol permutations",
        "sparse_cutoff": 128, "sparse_tau": str(tau), "sparse_t": str(t),
        "sparse_failure": str(sparse), "intervals": intervals,
        "required_margin": margin,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--generate", action="store_true")
    mode.add_argument("--verify", type=Path)
    parser.add_argument("--companion", type=Path, default=COMPANION)
    parser.add_argument("--log-k", type=int, default=16)
    parser.add_argument("--alphabet", type=int, default=128)
    parser.add_argument("--distance", default="0.10")
    parser.add_argument("--tau", type=int, default=512)
    parser.add_argument("--margin", type=int, default=40)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.verify:
        data = json.loads(args.verify.read_text(encoding="utf-8"))
    else:
        if args.output is None:
            parser.error("--generate requires --output")
        n = 2 << args.log_k
        data = generate(n, args.alphabet, int(F(args.distance) * n), args.tau, args.margin, args.companion)
    result = verify(data, args.companion)
    if args.generate:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
        args.output.with_suffix(".verified.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
