#!/usr/bin/env python3
"""Construct one EBCH32--PF31x33--BA row and test an excluded weight tail.

The checker accepts a query only when the selected exact solver proves that
no codeword lies in the requested low- or high-weight tail.  An unknown
result is not an acceptance.
"""

from __future__ import annotations

import argparse
from collections import Counter
import json
import random
import time
from pathlib import Path

import z3


WORKSTREAM = Path(__file__).resolve().parent
OUTPUT = WORKSTREAM / "ebch32_parityfanout31x33_ba3_row_min_distance_check.json"
B = 256
K = 128
CONSTITUENT_SPECTRUM = {
    0: 1,
    8: 620,
    12: 13_888,
    16: 36_518,
    20: 13_888,
    24: 620,
    32: 1,
}


def gf_multiply(left: int, right: int) -> int:
    result = 0
    value = left
    multiplier = right
    while multiplier:
        if multiplier & 1:
            result ^= value
        multiplier >>= 1
        value <<= 1
        if value & 0x20:
            value ^= 0x25  # x^5 + x^2 + 1
    return result


def gf_power(exponent: int) -> int:
    result = 1
    alpha = 2
    power = exponent % 31
    while power:
        if power & 1:
            result = gf_multiply(result, alpha)
        power >>= 1
        if power:
            alpha = gf_multiply(alpha, alpha)
    return result


def ebch32_generator_rows() -> list[int]:
    exponents = (1, 2, 4, 8, 16, 3, 6, 12, 24, 17, 5, 10, 20, 9, 18)
    polynomial = [1]
    for exponent in exponents:
        root = gf_power(exponent)
        updated = [0] * (len(polynomial) + 1)
        for degree, coefficient in enumerate(polynomial):
            updated[degree] ^= gf_multiply(coefficient, root)
            updated[degree + 1] ^= coefficient
        polynomial = updated
    if any(coefficient not in (0, 1) for coefficient in polynomial):
        raise ArithmeticError("BCH generator polynomial did not descend to GF(2)")
    generator = sum(coefficient << degree for degree, coefficient in enumerate(polynomial))
    rows = []
    for shift in range(16):
        cyclic = generator << shift
        parity = cyclic.bit_count() & 1
        rows.append(cyclic | (parity << 31))
    return rows


def audit_constituent(rows: list[int]) -> None:
    spectrum: Counter[int] = Counter()
    word = 0
    previous_gray = 0
    spectrum[0] += 1
    for value in range(1, 1 << 16):
        gray = value ^ (value >> 1)
        changed = gray ^ previous_gray
        word ^= rows[changed.bit_length() - 1]
        spectrum[word.bit_count()] += 1
        previous_gray = gray
    if dict(sorted(spectrum.items())) != CONSTITUENT_SPECTRUM:
        raise ArithmeticError(f"constituent spectrum mismatch: {dict(spectrum)}")


def direct_sum_rows(local_rows: list[int]) -> list[int]:
    return [
        row << (32 * block)
        for block in range(8)
        for row in local_rows
    ]


def parity_fanout(rows: list[int], source: list[int], target: list[int]) -> list[int]:
    source_mask = sum(1 << index for index in source)
    target_mask = sum(1 << index for index in target)
    return [
        row ^ (target_mask if (row & source_mask).bit_count() & 1 else 0)
        for row in rows
    ]


def permute_and_accumulate(rows: list[int], permutation: list[int]) -> list[int]:
    transformed = []
    for row in rows:
        prefix = 0
        output = 0
        for destination, source in enumerate(permutation):
            prefix ^= (row >> source) & 1
            output |= prefix << destination
        transformed.append(output)
    return transformed


def sample_row(seed: int) -> tuple[list[int], dict[str, object]]:
    rng = random.Random(seed)
    local = ebch32_generator_rows()
    audit_constituent(local)
    rows = direct_sum_rows(local)
    source = rng.sample(range(B), 31)
    remaining = [index for index in range(B) if index not in set(source)]
    target = rng.sample(remaining, 33)
    rows = parity_fanout(rows, source, target)
    first = list(range(B))
    second = list(range(B))
    rng.shuffle(first)
    rng.shuffle(second)
    rows = permute_and_accumulate(rows, first)
    rows = permute_and_accumulate(rows, second)
    return rows, {
        "seed": seed,
        "source": source,
        "target": target,
        "first_interleaver": first,
        "second_interleaver": second,
    }


def xor_all(terms: list[z3.BoolRef]) -> z3.BoolRef:
    """Return a balanced XOR expression for any number of terms."""
    if not terms:
        return z3.BoolVal(False)
    layer = terms
    while len(layer) > 1:
        layer = [
            z3.Xor(layer[index], layer[index + 1])
            if index + 1 < len(layer)
            else layer[index]
            for index in range(0, len(layer), 2)
        ]
    return layer[0]


def check_z3(rows: list[int], threshold: int, timeout_ms: int, coset: bool):
    messages = [z3.Bool(f"m_{index}") for index in range(K)]
    outputs = [z3.Bool(f"y_{index}") for index in range(B)]
    solver = z3.Solver()
    if timeout_ms:
        solver.set(timeout=timeout_ms)
    solver.add(z3.Or(messages))
    for coordinate in range(B):
        terms = [
            messages[row]
            for row in range(K)
            if (rows[row] >> coordinate) & 1
        ]
        solver.add(outputs[coordinate] == xor_all(terms))
    tail_bits = [z3.Not(output) for output in outputs] if coset else outputs
    solver.add(z3.PbLe([(output, 1) for output in tail_bits], threshold))
    started = time.perf_counter()
    result = solver.check()
    elapsed = time.perf_counter() - started
    witness = None
    if result == z3.sat:
        model = solver.model()
        message = sum(
            (1 << index) if z3.is_true(model.evaluate(value, model_completion=True)) else 0
            for index, value in enumerate(messages)
        )
        codeword = 0
        for index, row in enumerate(rows):
            if (message >> index) & 1:
                codeword ^= row
        witness = {
            "message_hex": hex(message),
            "codeword_hex": hex(codeword),
            "weight": codeword.bit_count(),
        }
        tail_weight = B - codeword.bit_count() if coset else codeword.bit_count()
        witness["tail_weight"] = tail_weight
        if not tail_weight <= threshold or (not coset and codeword == 0):
            raise ArithmeticError("solver witness does not satisfy the query")
    reason_unknown = solver.reason_unknown() if result == z3.unknown else ""
    return str(result), elapsed, witness, reason_unknown


def check_gurobi(rows: list[int], threshold: int, timeout_seconds: float, coset: bool):
    """Solve the same exact 0--1 feasibility problem through integer parity lifts."""
    import gurobipy as gp

    model = gp.Model("ebch32_pf_ba_tail")
    model.Params.OutputFlag = 0
    model.Params.TimeLimit = timeout_seconds
    messages = model.addVars(K, vtype=gp.GRB.BINARY, name="m")
    outputs = model.addVars(B, vtype=gp.GRB.BINARY, name="y")
    halves = model.addVars(B, lb=0, ub=64, vtype=gp.GRB.INTEGER, name="half")
    model.addConstr(gp.quicksum(messages[index] for index in range(K)) >= 1)
    for coordinate in range(B):
        terms = gp.quicksum(
            messages[row]
            for row in range(K)
            if (rows[row] >> coordinate) & 1
        )
        model.addConstr(terms == outputs[coordinate] + 2 * halves[coordinate])
    tail_weight = (
        gp.quicksum(1 - outputs[index] for index in range(B))
        if coset
        else gp.quicksum(outputs[index] for index in range(B))
    )
    model.addConstr(tail_weight <= threshold)
    started = time.perf_counter()
    model.optimize()
    elapsed = time.perf_counter() - started
    if model.Status == gp.GRB.INFEASIBLE:
        return "unsat", elapsed, None, ""
    if model.SolCount:
        message = sum(
            (1 << index) if messages[index].X > 0.5 else 0
            for index in range(K)
        )
        codeword = 0
        for index, row in enumerate(rows):
            if (message >> index) & 1:
                codeword ^= row
        actual_tail = B - codeword.bit_count() if coset else codeword.bit_count()
        if actual_tail > threshold or (not coset and codeword == 0):
            raise ArithmeticError("solver witness does not satisfy the query")
        return (
            "sat",
            elapsed,
            {
                "message_hex": hex(message),
                "codeword_hex": hex(codeword),
                "weight": codeword.bit_count(),
                "tail_weight": actual_tail,
            },
            "",
        )
    return "unknown", elapsed, None, f"Gurobi status {model.Status}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--threshold", type=int, default=16)
    parser.add_argument("--timeout-seconds", type=float, default=60.0)
    parser.add_argument("--solver", choices=("z3", "gurobi"), default="z3")
    parser.add_argument(
        "--coset",
        action="store_true",
        help="exclude codewords within threshold of the all-ones word",
    )
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    rows, setup = sample_row(args.seed)
    if args.solver == "z3":
        result, elapsed, witness, unknown = check_z3(
            rows, args.threshold, round(1000 * args.timeout_seconds), args.coset
        )
        solver_version = z3.get_version_string()
    else:
        result, elapsed, witness, unknown = check_gurobi(
            rows, args.threshold, args.timeout_seconds, args.coset
        )
        import gurobipy as gp

        solver_version = ".".join(str(value) for value in gp.gurobi.version())
    payload = {
        "schema": "ebch32-parityfanout31x33-ba3-row-min-distance-check-v1",
        "status": (
            "ACCEPTED_EXACT_SOLVER_RESULT"
            if result == "unsat"
            else "REJECTED_LOW_WEIGHT_WITNESS"
            if result == "sat"
            else "INCONCLUSIVE"
        ),
        "claim": {
            "minimum_distance_greater_than": args.threshold if result == "unsat" else None,
            "solver_result": result,
            "elapsed_seconds": elapsed,
            "witness": witness,
            "reason_unknown": unknown,
        },
        "parameters": {
            "row_length": B,
            "row_dimension": K,
            "threshold": args.threshold,
            "tail": "all_ones_coset" if args.coset else "zero_codeword",
            "timeout_seconds": args.timeout_seconds,
            "constituent_spectrum_audit": "passed",
        },
        "setup": setup,
        "solver": {
            "name": args.solver,
            "version": solver_version,
            "acceptance_rule": "Accept only an unsat result.",
        },
        "scope": (
            "This is an implementation feasibility checker. An independently "
            "checkable UNSAT proof format is not yet emitted."
        ),
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
