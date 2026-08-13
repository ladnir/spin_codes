#!/usr/bin/env python3
"""LP envelope cover for ordered packet-profile simplices.

On one ordered chamber, barycentric interpolation of the concave
``log2 Q(a)`` is a rigorous lower bound.  After this replacement, maximizing
the pointwise minimum of all affine witness parts is a linear program with ten
variables.  Unlike a fixed or mixed witness, this bound permits the winning
witness to change anywhere inside the simplex.

This first probe handles chambers whose first class is nonzero.  Class-zero
chambers are truncated by the minimum-weight constraint and require a small
triangulation handled separately.  Arithmetic and LP discovery are binary64.
"""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

from probe_packet8_profile_ordered_chambers import add_full_bijection, uniform_subset
from probe_packet8_profile_simplex_landscape import PROFILE_COUNT_LOG2, normalization_log2
from probe_packet8_shared_witness_io import load_shared_witness_arrays


def chamber_vertices(permutation: tuple[int, ...]) -> np.ndarray:
    subset = 0
    rows = []
    for packet_class in permutation:
        subset |= 1 << packet_class
        rows.append(uniform_subset(subset))
    return np.vstack(rows)


def envelope_upper(
    vertices: np.ndarray,
    constants: np.ndarray,
    charges: np.ndarray,
):
    normalizations = normalization_log2(vertices)
    # Variables are nine barycentric coordinates followed by the lower
    # affine-envelope value t.  Each witness imposes t <= C_i - <a,q_i>.
    inequalities = np.hstack((charges @ vertices.T, np.ones((len(constants), 1))))
    objective = np.append(normalizations, -1.0)
    result = linprog(
        objective,
        A_ub=inequalities,
        b_ub=constants,
        A_eq=np.asarray([[1.0] * len(vertices) + [0.0]]),
        b_eq=np.asarray([1.0]),
        bounds=[(0.0, None)] * len(vertices) + [(None, None)],
        method="highs",
    )
    if not result.success:
        raise RuntimeError(result.message)
    # linprog minimized chord(normalization)-t.
    return -float(result.fun), result.x[:-1]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--upgraded-cache", type=Path, required=True)
    parser.add_argument(
        "--extra-shared-report", type=Path, action="append", default=[]
    )
    parser.add_argument("--skip", type=int, default=0)
    parser.add_argument("--limit", type=int, default=100)
    parser.add_argument("--progress-every", type=int, default=25)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.skip < 0 or args.limit <= 0:
        raise SystemExit("simplex envelope cover: invalid skip/limit")

    names, constants, charges = load_shared_witness_arrays(
        args.upgraded_cache, args.extra_shared_report
    )
    full = add_full_bijection([])[0]
    names.append(full.name)
    constants = np.append(constants, full.constant_log2)
    charges = np.vstack((charges, full.linear_charge))
    target = -40.0 - PROFILE_COUNT_LOG2

    selected = itertools.islice(
        (row for row in itertools.permutations(range(9)) if row[0] != 0),
        args.skip,
        args.skip + args.limit,
    )
    rows = []
    for index, permutation in enumerate(selected, 1):
        score, barycentric = envelope_upper(
            chamber_vertices(permutation), constants, charges
        )
        rows.append(
            {
                "permutation": list(permutation),
                "score_log2": score,
                "margin_bits": target - score,
                "maximizer_barycentric": barycentric.tolist(),
            }
        )
        if args.progress_every and index % args.progress_every == 0:
            print(
                f"processed={index} covered={sum(row['score_log2'] <= target for row in rows)} "
                f"worst={max(row['score_log2'] for row in rows):.9f}",
                flush=True,
            )

    covered = sum(row["score_log2"] <= target for row in rows)
    worst = max(rows, key=lambda row: row["score_log2"])
    report = {
        "status": "DIAGNOSTIC_BINARY64_ORDERED_SIMPLEX_ENVELOPE_LP",
        "target_log2": target,
        "skip": args.skip,
        "tested": len(rows),
        "covered": covered,
        "worst": worst,
        "rows": rows,
    }
    print("packet-8 ordered-simplex affine-envelope cover")
    print(f"tested={len(rows)} covered={covered} uncovered={len(rows)-covered}")
    print(f"worst_score_log2={worst['score_log2']:.12f}")
    print(f"target_log2={target:.12f}")
    print("worst_order=" + ",".join(map(str, worst["permutation"])))
    print(
        "tested_simplex_envelope_cover="
        + ("PASS" if covered == len(rows) else "NO")
    )
    print("status=DIAGNOSTIC_BINARY64_ORDERED_SIMPLEX_ENVELOPE_LP")
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
