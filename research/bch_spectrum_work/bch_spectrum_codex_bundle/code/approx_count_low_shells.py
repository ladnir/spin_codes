#!/usr/bin/env python3
"""PAC approximate counts for the seven relevant BCH-256 weight shells.

The Boolean formula can represent P, Q, or one concrete five-dimensional
intermediate code Q <= C <= P.  Counts for P and Q can use their 2-design shell
symmetry; the exact quotient identity then recovers the spectrum of C.

ApproxMC returns a multiplicative (1+epsilon) approximation except with
probability delta.  The result is therefore probabilistic evidence with a
quantified guarantee, not a deterministic weight-enumerator certificate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from array import array
from pathlib import Path

import pyapproxmc
import pycryptosat
from pysat.card import CardEnc, EncType

from bch_quotient import (
    LENGTH,
    binary_poly_degree,
    generator_polynomial,
)


ROOT = Path(__file__).resolve().parents[1]
N = 256


def extended_row(punctured: int) -> int:
    """Append the parity coordinate at bit 255."""

    return punctured | ((punctured.bit_count() & 1) << LENGTH)


def generator_rows(code_name: str) -> list[int]:
    """Build a systematic basis for P, Q, or one intermediate C."""

    p_generator = generator_polynomial(37)
    q_generator = generator_polynomial(39)
    assert binary_poly_degree(p_generator) == 124
    assert binary_poly_degree(q_generator) == 132

    if code_name == "p":
        rows = [extended_row(p_generator << shift) for shift in range(131)]
    elif code_name == "q":
        rows = [extended_row(q_generator << shift) for shift in range(123)]
    elif code_name == "c":
        q_rows = [extended_row(q_generator << shift) for shift in range(123)]
        quotient_rows = [extended_row(p_generator << shift) for shift in range(5)]
        rows = q_rows + quotient_rows
    else:
        raise ValueError(code_name)
    rows = systematic_generator_rows(rows)
    assert all(row.bit_count() % 2 == 0 for row in rows)
    assert gf2_rank(rows) == len(rows)
    return rows


def gf2_rank(rows: list[int]) -> int:
    pivots: dict[int, int] = {}
    for original in rows:
        value = original
        while value:
            pivot = value.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = value
                break
            value ^= pivots[pivot]
    return len(pivots)


def systematic_generator_rows(rows: list[int]) -> list[int]:
    """Return an equivalent basis with 128 identity pivot columns."""

    reduced = list(rows)
    pivot_row = 0
    for coordinate in range(N):
        selected = next(
            (
                row
                for row in range(pivot_row, len(reduced))
                if (reduced[row] >> coordinate) & 1
            ),
            None,
        )
        if selected is None:
            continue
        reduced[pivot_row], reduced[selected] = reduced[selected], reduced[pivot_row]
        pivot = reduced[pivot_row]
        for row in range(len(reduced)):
            if row != pivot_row and ((reduced[row] >> coordinate) & 1):
                reduced[row] ^= pivot
        pivot_row += 1
        if pivot_row == len(reduced):
            break
    if pivot_row != len(reduced):
        raise ArithmeticError("generator rows are not independent")
    return reduced


def add_xor_definition(
    clauses: list[list[int]], inputs: list[int], output: int, next_var: int
) -> int:
    """Add a deterministic CNF definition output = XOR(inputs)."""

    if not inputs:
        clauses.append([-output])
        return next_var
    if len(inputs) == 1:
        clauses.extend(([-output, inputs[0]], [output, -inputs[0]]))
        return next_var

    accumulator = inputs[0]
    for position, right in enumerate(inputs[1:], start=1):
        last = position == len(inputs) - 1
        target = output if last else next_var
        if not last:
            next_var += 1
        # target <-> accumulator XOR right
        clauses.extend(
            (
                [accumulator, right, -target],
                [-accumulator, -right, -target],
                [accumulator, -right, target],
                [-accumulator, right, target],
            )
        )
        accumulator = target
    return next_var


def build_formula(
    rows: list[int], weight: int, fix_pair: bool
) -> tuple[list[list[int]], int, list[int]]:
    dimension = len(rows)
    message_variables = list(range(1, dimension + 1))
    output_variables = list(range(dimension + 1, dimension + N + 1))
    next_var = dimension + N + 1
    clauses: list[list[int]] = []

    for coordinate, output in enumerate(output_variables):
        inputs = [
            message_variables[row]
            for row in range(len(rows))
            if (rows[row] >> coordinate) & 1
        ]
        next_var = add_xor_definition(clauses, inputs, output, next_var)

    cardinality = CardEnc.equals(
        lits=output_variables,
        bound=weight,
        top_id=next_var - 1,
        encoding=EncType.seqcounter,
    )
    clauses.extend(cardinality.clauses)
    if fix_pair:
        clauses.extend(([output_variables[0]], [output_variables[1]]))
    top_var = max(next_var - 1, cardinality.nv)
    return clauses, top_var, message_variables


def flattened_clauses(clauses: list[list[int]]) -> array:
    flat = array("i")
    for clause in clauses:
        flat.extend(clause)
        flat.append(0)
    return flat


def count_shell(
    rows: list[int],
    weight: int,
    epsilon: float,
    delta: float,
    seed: int,
    fix_pair: bool,
) -> dict:
    started = time.perf_counter()
    clauses, top_var, projection = build_formula(rows, weight, fix_pair)
    build_seconds = time.perf_counter() - started

    counter = pyapproxmc.Counter(
        seed=seed,
        epsilon=epsilon,
        delta=delta,
    )
    counter.add_clauses(flattened_clauses(clauses))
    count_started = time.perf_counter()
    mantissa, exponent = counter.count(projection)
    count_seconds = time.perf_counter() - count_started
    estimate = int(mantissa) << int(exponent)
    factor = 1.0 + epsilon
    result = {
        "weight": weight,
        "estimate_mantissa": int(mantissa),
        "estimate_exponent": int(exponent),
        "estimate": str(estimate),
        "estimate_log2": math.log2(estimate) if estimate else -math.inf,
        "pac_lower": estimate / factor,
        "pac_upper": estimate * factor,
        "epsilon": epsilon,
        "delta": delta,
        "seed": seed,
        "variables": top_var,
        "clauses": len(clauses),
        "build_seconds": build_seconds,
        "count_seconds": count_seconds,
    }
    if fix_pair:
        numerator = math.comb(N, 2)
        denominator = math.comb(weight, 2)
        result.update(
            {
                "counted_object": "shell words containing coordinates 0 and 1",
                "two_design_recovery_numerator": numerator,
                "two_design_recovery_denominator": denominator,
                "recovered_full_shell_estimate": estimate * numerator / denominator,
                "recovered_full_shell_log2": (
                    math.log2(estimate * numerator / denominator)
                    if estimate
                    else -math.inf
                ),
            }
        )
    return result


def native_xor_sat_probe(rows: list[int], weight: int, fix_pair: bool) -> dict:
    """Find one shell word using native BCH parity equations and cardinality CNF."""

    pivot_coordinates = [(row & -row).bit_length() - 1 for row in rows]
    if len(set(pivot_coordinates)) != len(rows):
        raise ArithmeticError("systematic pivots are not unique")
    pivot_set = set(pivot_coordinates)

    solver = pycryptosat.Solver()
    for coordinate in range(N):
        if coordinate in pivot_set:
            continue
        xor_variables = [coordinate + 1]
        xor_variables.extend(
            pivot_coordinates[row] + 1
            for row in range(len(rows))
            if (rows[row] >> coordinate) & 1
        )
        solver.add_xor_clause(xor_variables, False)

    cardinality = CardEnc.equals(
        lits=list(range(1, N + 1)),
        bound=weight,
        top_id=N,
        encoding=EncType.seqcounter,
    )
    for clause in cardinality.clauses:
        solver.add_clause(clause)
    if fix_pair:
        solver.add_clause([1])
        solver.add_clause([2])

    started = time.perf_counter()
    satisfiable, model = solver.solve()
    seconds = time.perf_counter() - started
    observed_weight = None
    if satisfiable:
        observed_weight = sum(bool(model[index]) for index in range(1, N + 1))
        if observed_weight != weight:
            raise AssertionError((observed_weight, weight))
    return {
        "weight": weight,
        "satisfiable": bool(satisfiable),
        "observed_weight": observed_weight,
        "seconds": seconds,
        "fixed_coordinate_pair": fix_pair,
        "native_xor_constraints": N - len(rows),
        "cardinality_clauses": len(cardinality.clauses),
        "variables": cardinality.nv,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--weights",
        type=int,
        nargs="+",
        default=[38],
        choices=list(range(38, 52, 2)),
    )
    parser.add_argument("--epsilon", type=float, default=0.8)
    parser.add_argument("--delta", type=float, default=0.2)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--code", choices=("p", "q", "c"), default="p")
    parser.add_argument(
        "--no-fix-pair",
        action="store_true",
        help="disable the 2-design fixed-pair reduction (required for C counts)",
    )
    parser.add_argument(
        "--native-xor-sat-only",
        action="store_true",
        help="find one shell word with native XORs; do not run ApproxMC",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "generated" / "approxmc_low_shells.json",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.epsilon <= 0.0:
        raise SystemExit("epsilon must be positive")
    if not 0.0 < args.delta <= 1.0:
        raise SystemExit("delta must lie in (0,1]")

    rows = generator_rows(args.code)
    fix_pair = not args.no_fix_pair
    if args.code == "c" and fix_pair and not args.native_xor_sat_only:
        raise SystemExit("C need not be a 2-design; pass --no-fix-pair for C counts")
    digest = hashlib.sha256(
        b"".join(row.to_bytes(32, "little") for row in rows)
    ).hexdigest()
    if args.native_xor_sat_only:
        probes = []
        for weight in args.weights:
            print(f"native-XOR SAT probe at weight {weight}", flush=True)
            probe = native_xor_sat_probe(rows, weight, fix_pair)
            probes.append(probe)
            print(json.dumps(probe, indent=2), flush=True)
        return

    results = []
    for offset, weight in enumerate(args.weights):
        print(f"counting weight {weight}", flush=True)
        results.append(
            count_shell(
                rows,
                weight,
                epsilon=args.epsilon,
                delta=args.delta,
                seed=args.seed + offset,
                fix_pair=fix_pair,
            )
        )
        print(
            f"weight {weight}: approximately {results[-1]['estimate']} "
            f"({results[-1]['count_seconds']:.2f} s)",
            flush=True,
        )

    payload = {
        "classification": "ApproxMC PAC approximate model count",
        "code": {
            "name": args.code.upper(),
            "parameters": [N, len(rows)],
            "construction": {
                "p": "extended primitive BCH with designed distance 37",
                "q": "extended primitive BCH with designed distance 39",
                "c": (
                    "Q plus representatives p_generator*x^j for j=0..4"
                ),
            }[args.code],
            "generator_rows_sha256": digest,
            "rank": gf2_rank(rows),
            "fixed_coordinate_pair": fix_pair,
        },
        "guarantee": (
            "For each invocation, with probability at least 1-delta, the exact "
            "projected count lies between estimate/(1+epsilon) and "
            "estimate*(1+epsilon)."
        ),
        "results": results,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {args.output.resolve()}")


if __name__ == "__main__":
    main()
