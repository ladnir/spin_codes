#!/usr/bin/env python3
"""Independently solve Goal 05 with a pure-CNF observation-code encoding."""

from __future__ import annotations

import argparse
import hashlib
import json
import threading
import time
from pathlib import Path

from pysat.card import CardEnc, EncType
from pysat.solvers import Solver


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / "constructions" / "riffle_packetmul_2lap_g4"
GOAL04 = CANDIDATE / "receipts" / "goal04_four_node_primary.json"
OUTPUT = CANDIDATE / "receipts" / "goal05_distance_independent.json"
GENERATOR = 0xF4845518B9582A1F
MASK64 = (1 << 64) - 1
DIMENSION = 64
BLOCKS = 4
LENGTH = DIMENSION * BLOCKS
REJECTED_WEIGHT = 41


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def apply_columns(columns: tuple[int, ...], value: int) -> int:
    result = 0
    while value:
        bit = value & -value
        result ^= columns[bit.bit_length() - 1]
        value ^= bit
    return result


def inverse_columns(columns: tuple[int, ...]) -> tuple[int, ...]:
    rows = []
    for output in range(DIMENSION):
        left = sum(
            ((columns[input_bit] >> output) & 1) << input_bit
            for input_bit in range(DIMENSION)
        )
        rows.append(left | (1 << (DIMENSION + output)))
    for column in range(DIMENSION):
        pivot = next(row for row in range(column, DIMENSION) if (rows[row] >> column) & 1)
        rows[column], rows[pivot] = rows[pivot], rows[column]
        for row in range(DIMENSION):
            if row != column and ((rows[row] >> column) & 1):
                rows[row] ^= rows[column]
    inverse_rows = tuple(row >> DIMENSION for row in rows)
    return tuple(
        sum(
            ((inverse_rows[input_bit] >> output_bit) & 1) << input_bit
            for input_bit in range(DIMENSION)
        )
        for output_bit in range(DIMENSION)
    )


def systematic_state_columns() -> tuple[int, ...]:
    generator_rows = tuple((GENERATOR << row) | (1 << 127) for row in range(64))
    left = tuple(row & MASK64 for row in generator_rows)
    right = tuple((row >> 64) & MASK64 for row in generator_rows)
    inverse_left = inverse_columns(left)
    result = tuple(apply_columns(right, inverse_left[column]) for column in range(64))
    if any(apply_columns(left, inverse_left[column]) != 1 << column for column in range(64)):
        raise RuntimeError("Goal 05 pure CNF: systematic reconstruction failed")
    return result


def accumulate_by_bits(value: int) -> int:
    running = 0
    result = 0
    for bit in range(64):
        running ^= (value >> bit) & 1
        result |= running << bit
    return result


def observation_generators() -> tuple[int, ...]:
    columns = systematic_state_columns()

    def step(value: int) -> int:
        return accumulate_by_bits(apply_columns(columns, value))

    words = []
    for bit in range(64):
        state = 1 << bit
        word = 0
        for block in range(BLOCKS):
            word |= state << (64 * block)
            state = step(state)
        words.append(word)
    return tuple(words)


def add_xor_equivalence(
    clauses: list[list[int]],
    inputs: list[int],
    output: int,
    next_variable: int,
) -> int:
    if not inputs:
        clauses.append([-output])
        return next_variable
    if len(inputs) == 1:
        clauses.extend(([-inputs[0], output], [inputs[0], -output]))
        return next_variable

    def add_gate(left: int, right: int, result: int) -> None:
        clauses.extend(
            (
                [left, right, -result],
                [left, -right, result],
                [-left, right, result],
                [-left, -right, -result],
            )
        )

    previous = next_variable
    next_variable += 1
    add_gate(inputs[0], inputs[1], previous)
    for index in range(2, len(inputs)):
        result = output if index == len(inputs) - 1 else next_variable
        if result == next_variable:
            next_variable += 1
        add_gate(previous, inputs[index], result)
        previous = result
    if len(inputs) == 2:
        clauses.extend(([-previous, output], [previous, -output]))
    return next_variable


def replay_goal04(generators: tuple[int, ...]) -> dict:
    receipt = json.loads(GOAL04.read_text())
    witness = receipt["minimum_witness"]
    initial = int(witness["left2_hex"], 16)
    codeword = 0
    for bit, word in enumerate(generators):
        if (initial >> bit) & 1:
            codeword ^= word
    blocks = [
        (codeword >> (64 * block)) & MASK64 for block in range(BLOCKS)
    ]
    expected = [
        int(witness["left2_hex"], 16),
        int(witness["left1_hex"], 16),
        int(witness["anchor_hex"], 16),
        int(witness["right1_hex"], 16),
    ]
    if blocks != expected or codeword.bit_count() != 42:
        raise RuntimeError("Goal 05 pure CNF: weight-42 witness replay failed")
    return {
        "initial_output_hex": hex(initial),
        "block_hex": [hex(block) for block in blocks],
        "block_weights": [block.bit_count() for block in blocks],
        "total_weight": 42,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--timeout-seconds", type=float, default=900.0)
    parser.add_argument("--conflict-budget", type=int, default=100_000)
    args = parser.parse_args()
    if args.timeout_seconds <= 0:
        raise SystemExit("Goal 05 pure CNF: timeout must be positive")
    if args.conflict_budget <= 0:
        raise SystemExit("Goal 05 pure CNF: conflict budget must be positive")

    generators = observation_generators()
    witness = replay_goal04(generators)
    information_variables = list(range(1, DIMENSION + 1))
    output_variables = list(range(DIMENSION + 1, DIMENSION + LENGTH + 1))
    top_variable = output_variables[-1]
    global_cardinality = CardEnc.atmost(
        lits=output_variables,
        bound=REJECTED_WEIGHT,
        top_id=top_variable,
        encoding=EncType.seqcounter,
    )
    clauses = list(global_cardinality.clauses)
    top_variable = global_cardinality.nv
    for block in range(BLOCKS):
        literals = output_variables[64 * block : 64 * (block + 1)]
        lower = CardEnc.atleast(
            lits=literals,
            bound=9,
            top_id=top_variable,
            encoding=EncType.seqcounter,
        )
        clauses.extend(lower.clauses)
        top_variable = lower.nv
        upper = CardEnc.atmost(
            lits=literals,
            bound=14,
            top_id=top_variable,
            encoding=EncType.seqcounter,
        )
        clauses.extend(upper.clauses)
        top_variable = upper.nv
    clauses.append(information_variables)

    xor_clause_start = len(clauses)
    next_variable = top_variable + 1
    for coordinate, output_variable in enumerate(output_variables):
        inputs = [
            bit + 1
            for bit, word in enumerate(generators)
            if (word >> coordinate) & 1
        ]
        next_variable = add_xor_equivalence(
            clauses, inputs, output_variable, next_variable
        )
    top_variable = next_variable - 1

    started = time.perf_counter()
    with Solver(name="cadical195", bootstrap_with=clauses) as solver:
        solver.conf_budget(args.conflict_budget)
        timer = threading.Timer(args.timeout_seconds, solver.interrupt)
        timer.start()
        try:
            satisfiable = solver.solve_limited(expect_interrupt=True)
        finally:
            timer.cancel()
        elapsed = time.perf_counter() - started
        statistics = solver.accum_stats()
        model = solver.get_model() if satisfiable is True else None

    result = "SAT" if satisfiable is True else "UNSAT" if satisfiable is False else "UNKNOWN"
    payload = {
        "schema": "riffle-packetmul-2lap-g4-goal05-pure-cnf-v1",
        "candidate": "Riffle PacketMul-2Lap g=4",
        "evidence_label": "EXACT_SAT_RESULT_OR_TIMEOUT_DIAGNOSTIC",
        "source_sha256": digest(Path(__file__).resolve()),
        "goal04_primary_sha256": digest(GOAL04),
        "code": {
            "definition": "C4={(o,U(o),U^2(o),U^3(o)):o in F_2^64}",
            "length": LENGTH,
            "dimension": DIMENSION,
            "systematic": True,
        },
        "authenticated_upper_bound_witness": witness,
        "decision_problem": {
            "information_word_nonzero": True,
            "output_weight_upper_bound": REJECTED_WEIGHT,
            "proved_block_weight_range": [9, 14],
        },
        "encoding": {
            "solver": "CaDiCaL 1.9.5 through python-sat",
            "linear_constraints": "Tseitin XOR chains in pure CNF",
            "cardinality_encoding": "sequential counters",
            "cnf_variables": top_variable,
            "cnf_clauses": len(clauses),
            "xor_chain_clauses": len(clauses) - xor_clause_start,
        },
        "timeout_seconds": args.timeout_seconds,
        "conflict_budget": args.conflict_budget,
        "elapsed_seconds": elapsed,
        "solver_statistics": statistics,
        "result": result,
        "scope_limitation": (
            "UNSAT is an exact independent solver conclusion. UNKNOWN records only a "
            "timeout. This receipt does not contain a DRAT proof."
        ),
    }
    if model is not None:
        positive = {literal for literal in model if literal > 0}
        information = sum(
            1 << bit for bit in range(DIMENSION) if bit + 1 in positive
        )
        codeword = 0
        for bit, word in enumerate(generators):
            if (information >> bit) & 1:
                codeword ^= word
        weight = codeword.bit_count()
        if information == 0 or weight > REJECTED_WEIGHT:
            raise RuntimeError("Goal 05 pure CNF: SAT model replay failed")
        payload["counterexample"] = {
            "information_hex": hex(information),
            "block_hex": [
                hex((codeword >> (64 * block)) & MASK64) for block in range(BLOCKS)
            ],
            "block_weights": [
                ((codeword >> (64 * block)) & MASK64).bit_count()
                for block in range(BLOCKS)
            ],
            "total_weight": weight,
        }

    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"authenticated_weight_42_state={witness['initial_output_hex']}", flush=True)
    print(f"result={result}", flush=True)
    print(f"elapsed_seconds={elapsed:.6f}", flush=True)
    print(f"output={OUTPUT.relative_to(ROOT)}", flush=True)


if __name__ == "__main__":
    main()
