#!/usr/bin/env python3
"""Decide whether the eight-node autonomous response code has weight at most 92."""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

from pysat.card import CardEnc, EncType
from pysat.solvers import Solver


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / "constructions" / "riffle_packetmul_wrapmul_2lap_g4"
TAIL = CANDIDATE / "receipts" / "goal02_tail_list_primary.json"
OUTPUT = CANDIDATE / "receipts" / "goal02_c8_native_xor.json"
GENERATOR = 0xF4845518B9582A1F
MASK64 = (1 << 64) - 1
DIMENSION = 64
NODES = 8
LENGTH = DIMENSION * NODES
REJECTED_WEIGHT = 92


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def apply_columns(columns: tuple[int, ...], value: int) -> int:
    result = 0
    while value:
        low = value & -value
        result ^= columns[low.bit_length() - 1]
        value ^= low
    return result


def inverse_columns(columns: tuple[int, ...]) -> tuple[int, ...]:
    rows = []
    for output in range(64):
        left = sum(
            ((columns[input_bit] >> output) & 1) << input_bit
            for input_bit in range(64)
        )
        rows.append(left | (1 << (64 + output)))
    for column in range(64):
        pivot = next(row for row in range(column, 64) if (rows[row] >> column) & 1)
        rows[column], rows[pivot] = rows[pivot], rows[column]
        for row in range(64):
            if row != column and ((rows[row] >> column) & 1):
                rows[row] ^= rows[column]
    inverse_rows = tuple(row >> 64 for row in rows)
    return tuple(
        sum(
            ((inverse_rows[input_bit] >> output_bit) & 1) << input_bit
            for input_bit in range(64)
        )
        for output_bit in range(64)
    )


def systematic_right_columns() -> tuple[int, ...]:
    rows = tuple((GENERATOR << row) | (1 << 127) for row in range(64))
    left = tuple(row & MASK64 for row in rows)
    right = tuple((row >> 64) & MASK64 for row in rows)
    inverse = inverse_columns(left)
    if any(apply_columns(left, inverse[bit]) != 1 << bit for bit in range(64)):
        raise RuntimeError("Goal 02 C8: systematic reconstruction failed")
    return tuple(apply_columns(right, inverse[bit]) for bit in range(64))


def accumulate(value: int) -> int:
    value ^= value << 1
    value ^= value << 2
    value ^= value << 4
    value ^= value << 8
    value ^= value << 16
    value ^= value << 32
    return value & MASK64


def observation_generators() -> tuple[int, ...]:
    parity = systematic_right_columns()

    def step(output: int) -> int:
        return accumulate(apply_columns(parity, output))

    words = []
    for bit in range(DIMENSION):
        output = 1 << bit
        word = 0
        for node in range(NODES):
            word |= output << (64 * node)
            output = step(output)
        words.append(word)
    return tuple(words)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--timeout-seconds", type=float, default=300.0)
    args = parser.parse_args()
    if args.timeout_seconds <= 0:
        raise SystemExit("Goal 02 C8: timeout must be positive")

    generators = observation_generators()
    information_variables = list(range(1, DIMENSION + 1))
    output_variables = list(range(DIMENSION + 1, DIMENSION + LENGTH + 1))
    top_variable = output_variables[-1]
    cardinality = CardEnc.atmost(
        lits=output_variables,
        bound=REJECTED_WEIGHT,
        top_id=top_variable,
        encoding=EncType.totalizer,
    )
    clauses = list(cardinality.clauses)
    top_variable = cardinality.nv
    for block in range(2):
        literals = output_variables[256 * block : 256 * (block + 1)]
        lower = CardEnc.atleast(
            lits=literals,
            bound=36,
            top_id=top_variable,
            encoding=EncType.totalizer,
        )
        clauses.extend(lower.clauses)
        top_variable = lower.nv
        upper = CardEnc.atmost(
            lits=literals,
            bound=56,
            top_id=top_variable,
            encoding=EncType.totalizer,
        )
        clauses.extend(upper.clauses)
        top_variable = upper.nv
    clauses.append(information_variables)

    started = time.perf_counter()
    with Solver(name="cms", bootstrap_with=clauses) as solver:
        xor_widths = []
        for coordinate, output_variable in enumerate(output_variables):
            literals = [
                bit + 1
                for bit, word in enumerate(generators)
                if (word >> coordinate) & 1
            ]
            xor_widths.append(len(literals))
            solver.add_xor_clause(literals + [output_variable], value=False)
        solver.solver.time_budget(args.timeout_seconds)
        solver.conf_budget((1 << 31) - 1)
        satisfiable = solver.solve_limited()
        elapsed = time.perf_counter() - started
        try:
            statistics = solver.accum_stats()
        except NotImplementedError:
            statistics = {}
        model = solver.get_model() if satisfiable is True else None

    result = "SAT" if satisfiable is True else "UNSAT" if satisfiable is False else "UNKNOWN"
    payload = {
        "schema": "riffle-packetmul-wrapmul-2lap-g4-goal02-c8-native-xor-v1",
        "candidate": "Riffle PacketMul-WrapMul-2Lap g=4",
        "candidate_id": "riffle_packetmul_wrapmul_2lap_g4",
        "evidence_label": "EXACT_SAT_RESULT_OR_TIMEOUT_DIAGNOSTIC",
        "source_sha256": digest(Path(__file__).resolve()),
        "goal02_tail_primary_sha256": digest(TAIL),
        "code": {
            "definition": "C8={(o,U(o),...,U^7(o)):o in F_2^64}",
            "length": LENGTH,
            "dimension": DIMENSION,
            "systematic": True,
        },
        "decision_problem": {
            "information_word_nonzero": True,
            "output_weight_upper_bound": REJECTED_WEIGHT,
            "four_node_block_weight_range": [36, 56],
        },
        "global_implication_if_unsat": {
            "eight_node_distance_lower_bound": 93,
            "eight_node_block_count": 4096,
            "remaining_four_node_distance": 36,
            "full_response_distance_lower_bound": 4096 * 93 + 36,
            "low_spectrum_threshold": 377_530,
        },
        "encoding": {
            "solver": "CryptoMiniSat through python-sat",
            "linear_constraints": "native XOR",
            "cardinality_encoding": "totalizers",
            "information_variables": DIMENSION,
            "output_variables": LENGTH,
            "cnf_variables": top_variable,
            "cnf_clauses_before_xors": len(clauses),
            "native_xor_constraints": LENGTH,
            "xor_width_minimum": min(xor_widths),
            "xor_width_maximum": max(xor_widths),
        },
        "timeout_seconds": args.timeout_seconds,
        "elapsed_seconds": elapsed,
        "solver_statistics": statistics,
        "result": result,
        "scope_limitation": (
            "UNSAT proves the exact predicate inside this solver but has no "
            "independently checkable proof trace. UNKNOWN records only a timeout."
        ),
    }
    if model is not None:
        positive = {literal for literal in model if literal > 0}
        information = sum(1 << bit for bit in range(64) if bit + 1 in positive)
        codeword = 0
        for bit, word in enumerate(generators):
            if (information >> bit) & 1:
                codeword ^= word
        weight = codeword.bit_count()
        if information == 0 or weight > REJECTED_WEIGHT:
            raise RuntimeError("Goal 02 C8: SAT model replay failed")
        payload["counterexample"] = {
            "initial_output_hex": hex(information),
            "node_hex": [
                hex((codeword >> (64 * node)) & MASK64) for node in range(NODES)
            ],
            "node_weights": [
                ((codeword >> (64 * node)) & MASK64).bit_count()
                for node in range(NODES)
            ],
            "four_node_weights": [
                ((codeword >> (256 * block)) & ((1 << 256) - 1)).bit_count()
                for block in range(2)
            ],
            "total_weight": weight,
        }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"result={result}")
    print(f"elapsed_seconds={elapsed:.6f}")
    if "counterexample" in payload:
        print(f"counterexample_weight={payload['counterexample']['total_weight']}")
    print(f"output={OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
