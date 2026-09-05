#!/usr/bin/env python3
"""Solve the Goal 05 distance decision with native XOR constraints."""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

from pysat.card import CardEnc, EncType
from pysat.solvers import Solver

from analyze_riffle_dp_g4_autonomous_map import accumulate, build_apply
from analyze_systematic_group_kernel import systematic_state_columns


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / "constructions" / "riffle_packetmul_2lap_g4"
GOAL04 = CANDIDATE / "receipts" / "goal04_four_node_primary.json"
OUTPUT = CANDIDATE / "receipts" / "goal05_distance_primary.json"
DIMENSION = 64
BLOCKS = 4
LENGTH = DIMENSION * BLOCKS
REJECTED_WEIGHT = 41


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def observation_generators() -> tuple[int, ...]:
    apply_p = build_apply(systematic_state_columns())

    def step(value: int) -> int:
        return accumulate(apply_p(value))

    words = []
    for bit in range(DIMENSION):
        state = 1 << bit
        word = 0
        for block in range(BLOCKS):
            word |= state << (DIMENSION * block)
            state = step(state)
        words.append(word)
    return tuple(words)


def replay_goal04_witness(step) -> dict:
    receipt = json.loads(GOAL04.read_text())
    witness = receipt["minimum_witness"]
    if witness["anchor_alignment"] != 1:
        raise RuntimeError("Goal 05 native XOR: Goal 04 witness alignment changed")
    states = [
        int(witness["left2_hex"], 16),
        int(witness["left1_hex"], 16),
        int(witness["anchor_hex"], 16),
        int(witness["right1_hex"], 16),
    ]
    if any(step(states[index]) != states[index + 1] for index in range(3)):
        raise RuntimeError("Goal 05 native XOR: Goal 04 witness recurrence failed")
    weight = sum(state.bit_count() for state in states)
    if weight != 42:
        raise RuntimeError("Goal 05 native XOR: Goal 04 witness weight changed")
    return {
        "initial_output_hex": hex(states[0]),
        "state_hex": [hex(state) for state in states],
        "block_weights": [state.bit_count() for state in states],
        "total_weight": weight,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--timeout-seconds", type=float, default=900.0)
    args = parser.parse_args()
    if args.timeout_seconds <= 0:
        raise SystemExit("Goal 05 native XOR: timeout must be positive")

    apply_p = build_apply(systematic_state_columns())

    def step(value: int) -> int:
        return accumulate(apply_p(value))

    generators = observation_generators()
    if any((generators[bit] & ((1 << DIMENSION) - 1)) != 1 << bit for bit in range(64)):
        raise RuntimeError("Goal 05 native XOR: observation code is not systematic")
    witness = replay_goal04_witness(step)

    information_variables = list(range(1, DIMENSION + 1))
    output_variables = list(range(DIMENSION + 1, DIMENSION + LENGTH + 1))
    global_cardinality = CardEnc.atmost(
        lits=output_variables,
        bound=REJECTED_WEIGHT,
        top_id=output_variables[-1],
        encoding=EncType.seqcounter,
    )
    clauses = list(global_cardinality.clauses)
    top_variable = global_cardinality.nv
    block_cardinalities = []
    for block in range(BLOCKS):
        literals = output_variables[DIMENSION * block : DIMENSION * (block + 1)]
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
        block_cardinalities.append((lower, upper))
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
        "schema": "riffle-packetmul-2lap-g4-goal05-native-xor-v1",
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
        },
        "encoding": {
            "solver": "CryptoMiniSat through python-sat",
            "linear_constraints": "native XOR",
            "cardinality_encoding": "sequential counter",
            "information_variables": DIMENSION,
            "output_variables": LENGTH,
            "cnf_variables": top_variable,
            "cnf_clauses_before_xors": len(clauses),
            "native_xor_constraints": LENGTH,
            "xor_width_minimum": min(xor_widths),
            "xor_width_maximum": max(xor_widths),
            "proved_goal04_block_weight_range_under_decision_event": [9, 14],
        },
        "timeout_seconds": args.timeout_seconds,
        "elapsed_seconds": elapsed,
        "solver_statistics": statistics,
        "result": result,
        "scope_limitation": (
            "UNSAT is an exact solver conclusion but this receipt does not contain an "
            "independently checkable UNSAT proof. Goal 05 requires a separate CNF "
            "reconstruction. UNKNOWN records only a timeout."
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
            raise RuntimeError("Goal 05 native XOR: SAT model replay failed")
        payload["counterexample"] = {
            "information_hex": hex(information),
            "codeword_weight": weight,
            "block_hex": [
                hex((codeword >> (DIMENSION * block)) & ((1 << DIMENSION) - 1))
                for block in range(BLOCKS)
            ],
        }

    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"authenticated_weight_42_state={witness['initial_output_hex']}", flush=True)
    print(f"result={result}", flush=True)
    print(f"elapsed_seconds={elapsed:.6f}", flush=True)
    print(f"output={OUTPUT.relative_to(ROOT)}", flush=True)


if __name__ == "__main__":
    main()
