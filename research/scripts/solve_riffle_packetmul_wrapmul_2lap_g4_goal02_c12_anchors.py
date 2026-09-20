#!/usr/bin/env python3
"""Split the C12 weight-at-most-138 decision by a mandatory low-node anchor."""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

from pysat.card import CardEnc, EncType
from pysat.solvers import Solver

from solve_riffle_packetmul_wrapmul_2lap_g4_goal02_c12_native_xor import (
    DIMENSION,
    LENGTH,
    MASK64,
    NODES,
    REJECTED_WEIGHT,
    observation_generators,
)


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / "constructions" / "riffle_packetmul_wrapmul_2lap_g4"
BASE_SOURCE = (
    ROOT
    / "scripts"
    / "solve_riffle_packetmul_wrapmul_2lap_g4_goal02_c12_native_xor.py"
)
PART1 = CANDIDATE / "receipts" / "goal02_tail_list_independent.json"
OUTPUT = CANDIDATE / "receipts" / "goal02_c12_anchor_split.json"
ANCHOR_WEIGHT_MAXIMUM = REJECTED_WEIGHT // NODES


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_base_clauses(output_variables: list[int]) -> tuple[list[list[int]], int]:
    top_variable = output_variables[-1]
    cardinality = CardEnc.atmost(
        lits=output_variables,
        bound=REJECTED_WEIGHT,
        top_id=top_variable,
        encoding=EncType.seqcounter,
    )
    clauses = list(cardinality.clauses)
    top_variable = cardinality.nv
    for block in range(3):
        literals = output_variables[256 * block : 256 * (block + 1)]
        lower = CardEnc.atleast(
            lits=literals,
            bound=36,
            top_id=top_variable,
            encoding=EncType.seqcounter,
        )
        clauses.extend(lower.clauses)
        top_variable = lower.nv
        upper = CardEnc.atmost(
            lits=literals,
            bound=66,
            top_id=top_variable,
            encoding=EncType.seqcounter,
        )
        clauses.extend(upper.clauses)
        top_variable = upper.nv
    return clauses, top_variable


def replay_model(
    model: list[int], generators: tuple[int, ...]
) -> tuple[int, int, list[int]]:
    positive = {literal for literal in model if literal > 0}
    information = sum(
        1 << bit for bit in range(DIMENSION) if bit + 1 in positive
    )
    codeword = 0
    for bit, word in enumerate(generators):
        if (information >> bit) & 1:
            codeword ^= word
    return (
        information,
        codeword.bit_count(),
        [
            ((codeword >> (64 * node)) & MASK64).bit_count()
            for node in range(NODES)
        ],
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-timeout-seconds", type=float, default=20.0)
    args = parser.parse_args()
    if args.case_timeout_seconds <= 0:
        raise SystemExit("Goal 02 C12 anchors: timeout must be positive")

    generators = observation_generators()
    information_variables = list(range(1, DIMENSION + 1))
    output_variables = list(range(DIMENSION + 1, DIMENSION + LENGTH + 1))
    base_clauses, base_top = build_base_clauses(output_variables)
    xor_rows = []
    for coordinate, output_variable in enumerate(output_variables):
        inputs = [
            bit + 1
            for bit, word in enumerate(generators)
            if (word >> coordinate) & 1
        ]
        xor_rows.append(inputs + [output_variable])

    rows = []
    counterexample = None
    total_started = time.perf_counter()
    for anchor_node in range(NODES):
        anchor_literals = output_variables[
            64 * anchor_node : 64 * (anchor_node + 1)
        ]
        anchor_upper = CardEnc.atmost(
            lits=anchor_literals,
            bound=ANCHOR_WEIGHT_MAXIMUM,
            top_id=base_top,
            encoding=EncType.seqcounter,
        )
        clauses = list(base_clauses)
        clauses.extend(anchor_upper.clauses)
        clauses.append(anchor_literals)
        clauses.append(information_variables)

        started = time.perf_counter()
        with Solver(name="cms", bootstrap_with=clauses) as solver:
            for literals in xor_rows:
                solver.add_xor_clause(literals, value=False)
            solver.solver.time_budget(args.case_timeout_seconds)
            solver.conf_budget((1 << 31) - 1)
            satisfiable = solver.solve_limited()
            elapsed = time.perf_counter() - started
            try:
                statistics = solver.accum_stats()
            except NotImplementedError:
                statistics = {}
            model = solver.get_model() if satisfiable is True else None
        result = (
            "SAT"
            if satisfiable is True
            else "UNSAT"
            if satisfiable is False
            else "UNKNOWN"
        )
        row = {
            "anchor_node": anchor_node,
            "anchor_weight_upper_bound": ANCHOR_WEIGHT_MAXIMUM,
            "result": result,
            "elapsed_seconds": elapsed,
            "solver_statistics": statistics,
            "cnf_variables": anchor_upper.nv,
            "cnf_clauses": len(clauses),
        }
        rows.append(row)
        print(
            f"anchor_node={anchor_node} result={result} "
            f"elapsed_seconds={elapsed:.6f}",
            flush=True,
        )
        if model is not None:
            information, weight, node_weights = replay_model(model, generators)
            if (
                information == 0
                or weight > REJECTED_WEIGHT
                or node_weights[anchor_node] > ANCHOR_WEIGHT_MAXIMUM
            ):
                raise RuntimeError("Goal 02 C12 anchors: SAT replay failed")
            counterexample = {
                "anchor_node": anchor_node,
                "initial_output_hex": hex(information),
                "node_weights": node_weights,
                "four_node_weights": [
                    sum(node_weights[4 * block : 4 * block + 4])
                    for block in range(3)
                ],
                "total_weight": weight,
            }
            break

    if counterexample is not None:
        result = "SAT"
    elif len(rows) == NODES and all(row["result"] == "UNSAT" for row in rows):
        result = "UNSAT"
    else:
        result = "UNKNOWN"
    payload = {
        "schema": "riffle-packetmul-wrapmul-2lap-g4-goal02-c12-anchor-split-v1",
        "candidate": "Riffle PacketMul-WrapMul-2Lap g=4",
        "candidate_id": "riffle_packetmul_wrapmul_2lap_g4",
        "evidence_label": "EXACT_SPLIT_SAT_RESULT_OR_TIMEOUT_DIAGNOSTIC",
        "source_sha256": digest(Path(__file__).resolve()),
        "base_source_sha256": digest(BASE_SOURCE),
        "goal02_part1_independent_sha256": digest(PART1),
        "decision_problem": {
            "code": "C12={(o,U(o),...,U^11(o)):o in F_2^64}",
            "information_nonzero": True,
            "total_weight_upper_bound": REJECTED_WEIGHT,
            "mandatory_anchor_reason": "one of 12 nodes has weight at most floor(138/12)=11",
            "anchor_weight_upper_bound": ANCHOR_WEIGHT_MAXIMUM,
            "anchor_phases": NODES,
        },
        "solver": "CryptoMiniSat through python-sat with native XOR",
        "case_timeout_seconds": args.case_timeout_seconds,
        "total_elapsed_seconds": time.perf_counter() - total_started,
        "rows": rows,
        "result": result,
        "scope_limitation": (
            "The split is exhaustive only if every one of the twelve rows is decided. "
            "UNKNOWN rows leave the exact C12 predicate open."
        ),
    }
    if counterexample is not None:
        payload["counterexample"] = counterexample
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"result={result}", flush=True)
    print(f"output={OUTPUT.relative_to(ROOT)}", flush=True)


if __name__ == "__main__":
    main()
