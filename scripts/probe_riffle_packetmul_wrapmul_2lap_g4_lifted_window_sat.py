#!/usr/bin/env python3
"""Bounded SAT probe for lifted-window minimum distances.

The small-weight enumerator is the theorem-facing source through weight four.
This script extends discovery beyond that boundary.  UNSAT results remain
solver evidence until an independently checkable proof is recorded.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from pysat.card import CardEnc, EncType
from pysat.formula import CNF, IDPool
from pysat.solvers import Solver

from analyze_riffle_packetmul_wrapmul_2lap_g4_lifted_windows import (
    build_apply,
    coordinate_syndromes,
    lifted_step,
    observation_rows,
    parity_checks,
)
from analyze_systematic_group_kernel import systematic_state_columns


ROOT = Path(__file__).resolve().parents[1]
EXTENDED_NODES = 64


def add_xor(cnf: CNF, left: int, right: int, result: int) -> None:
    cnf.append([left, right, -result])
    cnf.append([-left, -right, -result])
    cnf.append([left, -right, result])
    cnf.append([-left, right, result])


def add_even_parity(cnf: CNF, literals: list[int], pool: IDPool) -> None:
    if not literals:
        return
    if len(literals) == 1:
        cnf.append([-literals[0]])
        return
    accumulator = literals[0]
    for literal in literals[1:]:
        result = pool.id()
        add_xor(cnf, accumulator, literal, result)
        accumulator = result
    cnf.append([-accumulator])


def solve_bound(
    base: CNF,
    coordinate_variables: list[int],
    pool_top: int,
    bound: int,
) -> tuple[bool, list[int]]:
    cnf = CNF(from_clauses=base.clauses)
    pool = IDPool(start_from=pool_top + 1)
    cnf.extend(
        CardEnc.atmost(
            lits=coordinate_variables,
            bound=bound,
            vpool=pool,
            encoding=EncType.seqcounter,
        ).clauses
    )
    cnf.append(coordinate_variables)
    with Solver(name="cadical195", bootstrap_with=cnf.clauses) as solver:
        satisfiable = solver.solve()
        if not satisfiable:
            return False, []
        model = set(literal for literal in solver.get_model() if literal > 0)
    support = [
        coordinate - 1
        for coordinate in coordinate_variables
        if coordinate in model
    ]
    return True, support


def recover_initial_state(generator: list[int], prefix: int) -> int:
    mask = (1 << 128) - 1
    basis: list[tuple[int, int] | None] = [None] * 128
    for source, full_row in enumerate(generator):
        row = full_row & mask
        coefficients = 1 << source
        while row:
            pivot = (row & -row).bit_length() - 1
            entry = basis[pivot]
            if entry is None:
                basis[pivot] = (row, coefficients)
                break
            row ^= entry[0]
            coefficients ^= entry[1]
        else:
            raise RuntimeError("lifted SAT probe: first two outputs lost rank")

    value = prefix & mask
    initial = 0
    while value:
        pivot = (value & -value).bit_length() - 1
        entry = basis[pivot]
        if entry is None:
            raise RuntimeError("lifted SAT probe: prefix is outside the code")
        value ^= entry[0]
        initial ^= entry[1]
    return initial


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--nodes", type=int, required=True)
    parser.add_argument("--minimum-weight", type=int, default=1)
    parser.add_argument("--maximum-weight", type=int, default=32)
    args = parser.parse_args()
    if args.nodes < 2 or args.nodes > 32:
        raise SystemExit("lifted SAT probe: nodes must lie in [2,32]")

    parity = build_apply(systematic_state_columns())
    columns = 64 * args.nodes
    generator = observation_rows(args.nodes, parity)
    checks, pivots = parity_checks(generator, columns)
    syndromes = coordinate_syndromes(checks, columns)

    pool = IDPool(start_from=columns + 1)
    base = CNF()
    for syndrome_bit in range(len(checks)):
        literals = [
            coordinate + 1
            for coordinate, syndrome in enumerate(syndromes)
            if (syndrome >> syndrome_bit) & 1
        ]
        add_even_parity(base, literals, pool)

    trials = []
    minimum = None
    witness: list[int] = []
    for bound in range(args.minimum_weight, args.maximum_weight + 1):
        satisfiable, support = solve_bound(
            base,
            list(range(1, columns + 1)),
            pool.top,
            bound,
        )
        trials.append({"bound": bound, "result": "SAT" if satisfiable else "UNSAT"})
        print(f"nodes={args.nodes} bound={bound} result={trials[-1]['result']}")
        if satisfiable:
            minimum = len(support)
            witness = support
            break

    if minimum is not None:
        word = sum(1 << coordinate for coordinate in witness)
        if any((row & word).bit_count() & 1 for row in checks):
            raise RuntimeError("lifted SAT probe: model replay failed")
        initial = recover_initial_state(generator, word)
        a = initial & ((1 << 64) - 1)
        b = initial >> 64
        extended_weights = []
        for _ in range(EXTENDED_NODES):
            a, b, output_word = lifted_step(a, b, parity)
            extended_weights.append(output_word.bit_count())
        if sum(extended_weights[: args.nodes]) != minimum:
            raise RuntimeError("lifted SAT probe: recovered state does not replay")
    else:
        initial = 0
        extended_weights = []

    payload = {
        "schema": "riffle-packetmul-wrapmul-2lap-g4-lifted-window-sat-v1",
        "candidate": "Riffle PacketMul-WrapMul-2Lap g=4",
        "evidence_label": "BOUNDED_SAT_WITH_EXACT_WITNESS_REPLAY",
        "nodes": args.nodes,
        "length": columns,
        "dimension": len(pivots),
        "solver": "cadical195",
        "maximum_weight": args.maximum_weight,
        "trials": trials,
        "minimum_distance_if_found": minimum,
        "initial_a_hex": f"0x{initial & ((1 << 64) - 1):016x}" if minimum else None,
        "initial_b_hex": f"0x{initial >> 64:016x}" if minimum else None,
        "extended_node_weights": extended_weights,
        "extended_weight": sum(extended_weights),
        "witness": [
            {
                "coordinate": coordinate,
                "node": coordinate // 64,
                "bit": coordinate % 64,
            }
            for coordinate in witness
        ],
        "scope_limitation": (
            "SAT witnesses are replayed exactly. UNSAT results are solver evidence "
            "and are not theorem-facing certificates without a checked proof."
        ),
    }
    output = (
        ROOT
        / "constructions"
        / "riffle_packetmul_wrapmul_2lap_g4"
        / "receipts"
        / f"goal04_lifted_window_{args.nodes:02d}_sat.json"
    )
    output.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"output={output}")
    print("status=BOUNDED_LIFTED_WINDOW_SAT_PROBE")


if __name__ == "__main__":
    main()
