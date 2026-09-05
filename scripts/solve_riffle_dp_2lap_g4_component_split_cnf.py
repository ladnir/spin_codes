#!/usr/bin/env python3
"""Exact CNF prototype for the component-split 24-node nearest-pair problem."""

from __future__ import annotations

import argparse
import hashlib
import json
import threading
import time
from pathlib import Path

from pysat.card import CardEnc, EncType
from pysat.solvers import Solver

from analyze_systematic_group_kernel import systematic_state_columns
from analyze_riffle_dp_g4_autonomous_map import accumulate, binary_rank, build_apply
from audit_riffle_dp_2lap_g4_component_endpoint_dimension30 import (
    observations,
    support_degrees,
    support_dimension,
)
from probe_riffle_dp_g4_component_mixing import transpose_columns
from solve_riffle_dp_2lap_g4_component_split_sat import (
    balanced_split,
    component_basis,
    generator_words,
)


ROOT = Path(__file__).resolve().parents[1]
LENGTH = 384
WINDOW = 24
REQUIRED_DISTANCE = 97
ENCODINGS = {
    "seqcounter": EncType.seqcounter,
    "cardnetwrk": EncType.cardnetwrk,
    "totalizer": EncType.totalizer,
    "kmtotalizer": EncType.kmtotalizer,
}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def xor_output_clauses(
    inputs: list[int], output: int, first_auxiliary: int
) -> tuple[list[list[int]], int]:
    literals = inputs + [output]
    if not literals:
        return [], first_auxiliary
    if len(literals) == 1:
        return [[-literals[0]]], first_auxiliary
    clauses = []
    accumulator = literals[0]
    next_variable = first_auxiliary
    for literal in literals[1:]:
        result = next_variable
        next_variable += 1
        # result = accumulator xor literal.
        clauses.extend(
            (
                [accumulator, literal, -result],
                [-accumulator, -literal, -result],
                [accumulator, -literal, result],
                [-accumulator, literal, result],
            )
        )
        accumulator = result
    clauses.append([-accumulator])
    return clauses, next_variable


def solve_side(
    *,
    side: str,
    words: tuple[int, ...],
    timeout_seconds: float,
    solver_name: str,
    encoding_name: str,
) -> dict:
    dimension = len(words)
    information_variables = list(range(1, dimension + 1))
    output_variables = list(range(dimension + 1, dimension + LENGTH + 1))
    cardinality_literals = (
        output_variables if side == "weight" else [-variable for variable in output_variables]
    )
    cardinality = CardEnc.atmost(
        lits=cardinality_literals,
        bound=REQUIRED_DISTANCE - 1,
        top_id=output_variables[-1],
        encoding=ENCODINGS[encoding_name],
    )
    clauses = list(cardinality.clauses)
    clauses.append(information_variables)
    next_variable = cardinality.nv + 1
    for coordinate, output in enumerate(output_variables):
        inputs = [
            index + 1
            for index, word in enumerate(words)
            if (word >> coordinate) & 1
        ]
        parity_clauses, next_variable = xor_output_clauses(
            inputs, output, next_variable
        )
        clauses.extend(parity_clauses)

    started = time.perf_counter()
    with Solver(name=solver_name, bootstrap_with=clauses) as solver:
        timer = threading.Timer(timeout_seconds, solver.interrupt)
        timer.start()
        try:
            satisfiable = solver.solve_limited(expect_interrupt=True)
        finally:
            timer.cancel()
        elapsed = time.perf_counter() - started
        stats = solver.accum_stats()
        model = solver.get_model() if satisfiable is True else None

    row = {
        "side": side,
        "timeout_seconds": timeout_seconds,
        "elapsed_seconds": elapsed,
        "solver": solver_name,
        "cardinality_encoding": encoding_name,
        "variable_count": next_variable - 1,
        "clause_count": len(clauses),
        "solver_statistics": stats,
        "result": (
            "COUNTEREXAMPLE" if satisfiable is True else "PASS" if satisfiable is False else "UNKNOWN"
        ),
    }
    if model is not None:
        positive = {literal for literal in model if literal > 0}
        information = sum(
            (1 << index)
            for index in range(dimension)
            if index + 1 in positive
        )
        row["information_hex"] = hex(information)
    return row


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--support", type=lambda value: int(value, 0), required=True)
    parser.add_argument("--packet-value", type=int, default=1)
    parser.add_argument("--timeout-seconds", type=float, default=60.0)
    parser.add_argument("--solver", default="cd19")
    parser.add_argument("--encoding", choices=tuple(ENCODINGS), default="seqcounter")
    parser.add_argument("--side", choices=("weight", "complement", "both"), default="both")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not 1 <= args.support <= 0x7F:
        raise SystemExit("component-split CNF: support must lie in 1..127")
    if not 1 <= args.packet_value <= 14:
        raise SystemExit("component-split CNF: packet value must lie in 1..14")
    if args.timeout_seconds <= 0:
        raise SystemExit("component-split CNF: timeout must be positive")

    apply_state = build_apply(systematic_state_columns())

    def step(state: int) -> int:
        return apply_state(accumulate(state))

    transpose_step = build_apply(
        transpose_columns(tuple(step(1 << bit) for bit in range(64)))
    )
    left_mask, right_mask = balanced_split(args.support)
    left_basis = component_basis(transpose_step, left_mask)
    right_basis = component_basis(transpose_step, right_mask)
    basis = left_basis + right_basis
    dimension = support_dimension(args.support)
    if len(basis) != dimension or binary_rank(basis) != dimension:
        raise RuntimeError("component-split CNF: basis mismatch")
    states = observations(step, args.packet_value, WINDOW)
    words = generator_words(states, basis)
    sides = ("weight", "complement") if args.side == "both" else (args.side,)
    rows = [
        solve_side(
            side=side,
            words=words,
            timeout_seconds=args.timeout_seconds,
            solver_name=args.solver,
            encoding_name=args.encoding,
        )
        for side in sides
    ]

    for row in rows:
        if row["result"] != "COUNTEREXAMPLE":
            continue
        information = int(row.pop("information_hex"), 16)
        character = 0
        for index, basis_character in enumerate(basis):
            if (information >> index) & 1:
                character ^= basis_character
        weight = sum((character & state).bit_count() & 1 for state in states)
        side_weight = weight if row["side"] == "weight" else LENGTH - weight
        if not character or side_weight >= REQUIRED_DISTANCE:
            raise RuntimeError("component-split CNF: model replay failed")
        row["character_hex"] = hex(character)
        row["weight"] = weight
        row["complement_weight"] = LENGTH - weight

    result = (
        "COUNTEREXAMPLE"
        if any(row["result"] == "COUNTEREXAMPLE" for row in rows)
        else "PASS"
        if len(rows) == 2 and all(row["result"] == "PASS" for row in rows)
        else rows[0]["result"]
        if len(rows) == 1
        else "UNKNOWN"
    )
    payload = {
        "schema": "riffle-dp-2lap-g4-component-split-cnf-v1",
        "candidate": "Riffle DP-2Lap g=4",
        "evidence_label": "EXACT_CNF_RESULT_OR_TIMEOUT_DIAGNOSTIC",
        "source_sha256": digest(Path(__file__).resolve()),
        "support_mask_hex": hex(args.support),
        "support_degrees": list(support_degrees(args.support)),
        "code_dimension": dimension,
        "packet_value": args.packet_value,
        "window_nodes": WINDOW,
        "code_length": LENGTH,
        "required_two_sided_distance": REQUIRED_DISTANCE,
        "left_support_mask_hex": hex(left_mask),
        "left_support_degrees": list(support_degrees(left_mask)),
        "left_dimension": len(left_basis),
        "right_support_mask_hex": hex(right_mask),
        "right_support_degrees": list(support_degrees(right_mask)),
        "right_dimension": len(right_basis),
        "rows": rows,
        "result": result,
        "scope_limitation": (
            "PASS or COUNTEREXAMPLE is an exact solver conclusion, and each "
            "counterexample is replayed from the reconstructed code. UNKNOWN "
            "reports only a timeout. This prototype does not retain a proof log."
        ),
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"support={hex(args.support)} split={hex(left_mask)}+{hex(right_mask)}")
    for row in rows:
        print(
            f"side={row['side']} result={row['result']} "
            f"elapsed={row['elapsed_seconds']:.6f} conflicts={row['solver_statistics'].get('conflicts')}"
        )
    print(f"output={args.output}")
    print(f"status={result}")


if __name__ == "__main__":
    main()
