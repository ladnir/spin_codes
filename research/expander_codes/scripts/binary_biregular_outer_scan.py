#!/usr/bin/env python3
"""Scan one exact outer-support block for binary two-sided regular EC.

The scan reuses one family of uniform-slice transfer matrices across every
support in the requested block.  It uses binary64 arithmetic and is therefore
a parameter diagnostic, not a certificate.  Run blocks serially because the
slice family can require substantial time and memory.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict

from binary_biregular_diagnostic import exact_biregular_fixed_marker_block


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--k", type=int, default=1_048_575)
    parser.add_argument("--left-degree", type=int, default=6)
    parser.add_argument("--right-degree", type=int, default=3)
    parser.add_argument("--cutoff", type=int, default=230_729)
    parser.add_argument("--memory", type=int, default=80)
    parser.add_argument("--support-start", type=int, required=True)
    parser.add_argument("--support-limit", type=int, required=True)
    parser.add_argument("--output-marker", type=float, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    result = exact_biregular_fixed_marker_block(
        code="ec",
        k=args.k,
        left_degree=args.left_degree,
        right_degree=args.right_degree,
        cutoff=args.cutoff,
        memory=args.memory,
        support_start=args.support_start,
        support_limit=args.support_limit,
        output_marker=args.output_marker,
    )
    print(json.dumps(asdict(result), indent=2))


if __name__ == "__main__":
    main()
