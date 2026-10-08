"""Replay the companion's generic two-scan certificate in characteristic two.

Only the cardinality of the nonzero coefficient alphabet and the MDS
shortening dimensions enter this bound.  The companion's prime-specific
outer implementation and arithmetic-cost annotations do not transfer.
See ../SMALL_COEFFICIENTS.md for the construction and transfer argument.

Discovery may use floating optimization.  Acceptance re-evaluates every
leaf using the companion's exact rational/outward-logarithm verifier.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from functools import partial
import hashlib
import json
from pathlib import Path
import sys


REPO = Path(__file__).resolve().parents[4]
DEFAULT_COMPANION = REPO.parent / "gen-BAA" / "field-universality" / "expander-lift"
DEPENDENCIES = (
    "convolution_width_one_finite.py",
    "convolution_design_space.py",
    "repeat_accumulate_precode.py",
    "sparse_mixing_optimization.py",
)


def load_verifier(companion: Path, n: int):
    """Override all N-bound defaults, not just the module-level constant."""
    sys.path.insert(0, str(companion.resolve()))
    import convolution_width_one_finite as finite

    finite.N = n
    original_entropy = finite.entropy_terms
    finite.entropy_terms = lambda weight, n=n: original_entropy(weight, n)
    finite.sparse_bound = partial(finite.sparse_bound, n=n)
    finite.propose = partial(finite.propose, n=n)
    return finite


def transferred_result(result: dict, witness: Path, companion: Path) -> dict:
    return {
        "construction": "GRS outer; two independent uniform symbol permutations; two independently input-weighted accumulators",
        "fields": [8, 32, 64, 128],
        "field_size_condition": "at least block length and alphabet size plus one",
        "coefficient_alphabet": "any Q distinct nonzero field elements",
        "n": result["n"],
        "k": result["k"],
        "block_n": 2 * result["block_k"],
        "block_k": result["block_k"],
        "alphabet_size": result["q"],
        "minimum_distance_at_least": result["minimum_distance_at_least"],
        "relative_distance": str(Fraction(result["minimum_distance_at_least"], result["n"])),
        "failure_upper": result["failure_upper"],
        "failure_log2_upper": result["failure_log2_upper"],
        "sparse_intervals": result["sparse_intervals"],
        "dense_leaves": result["dense_leaves"],
        "complete_distance_certificate": result["complete_distance_certificate"],
        "witness": str(witness.resolve()),
        "witness_sha256": hashlib.sha256(witness.read_bytes()).hexdigest(),
        "verifier_sha256": {
            name: hashlib.sha256((companion / name).read_bytes()).hexdigest()
            for name in DEPENDENCIES
        },
        "scope": "Generic quotient-rank bound only; no odd-characteristic pair refinement; no benchmark or Newton arithmetic-cost claim.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--companion", type=Path, default=DEFAULT_COMPANION)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--verify", type=Path)
    mode.add_argument("--generate", action="store_true")
    parser.add_argument("--log-k", type=int, default=16)
    parser.add_argument("--block-k", type=int, default=4)
    parser.add_argument("--alphabet", type=int, default=128)
    parser.add_argument("--distance", default="0.10")
    parser.add_argument("--source-min", type=int, default=32)
    parser.add_argument("--max-nodes", type=int, default=1500)
    parser.add_argument("--sparse-bits", type=int, default=48)
    parser.add_argument("--output-dir", type=Path, default=REPO / "tmp" / "char2-field-spin")
    args = parser.parse_args()
    if args.verify:
        data = json.loads(args.verify.read_text(encoding="utf-8"))
        finite = load_verifier(args.companion, data["n"])
        result = finite.verify(data)
        witness = args.verify
    else:
        n = 2 << args.log_k
        if n % (2 * args.block_k):
            parser.error("message length must be divisible by the outer dimension")
        e = int(Fraction(args.distance) * n)
        args.output_dir.mkdir(parents=True, exist_ok=True)
        finite = load_verifier(args.companion, n)
        finite.ROOT = args.output_dir
        result = finite.generate(
            args.alphabet, args.block_k, e, args.source_min,
            args.max_nodes, args.sparse_bits,
        )
        witness = args.output_dir / f"convolution-width-one-direct-r{args.block_k}-q{args.alphabet}-e{e}.json"
    print(json.dumps(transferred_result(result, witness, args.companion), indent=2))


if __name__ == "__main__":
    main()
