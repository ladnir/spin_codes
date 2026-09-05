#!/usr/bin/env python3
"""SAT probe for two low nine-node windows on one lifted orbit."""

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
    observation_rows,
    parity_checks,
)
from analyze_systematic_group_kernel import systematic_state_columns
from probe_riffle_packetmul_wrapmul_2lap_g4_lifted_window_sat import (
    add_even_parity,
    recover_initial_state,
)


ROOT = Path(__file__).resolve().parents[1]
WINDOW = 9
LOW_WEIGHT_MAXIMUM = 52
EXTENDED_NODES = 96
MASK64 = (1 << 64) - 1


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--offset", type=int, required=True)
    args = parser.parse_args()
    if args.offset < 1:
        raise SystemExit("repeated low windows: offset must be positive")

    nodes = args.offset + WINDOW
    parity = build_apply(systematic_state_columns())
    generator = observation_rows(nodes, parity)
    columns = 64 * nodes
    checks, pivots = parity_checks(generator, columns)
    if len(pivots) != 128:
        raise RuntimeError("repeated low windows: observation rank changed")
    syndromes = coordinate_syndromes(checks, columns)

    coordinate_variables = list(range(1, columns + 1))
    pool = IDPool(start_from=columns + 1)
    cnf = CNF()
    for syndrome_bit in range(len(checks)):
        literals = [
            coordinate + 1
            for coordinate, syndrome in enumerate(syndromes)
            if (syndrome >> syndrome_bit) & 1
        ]
        add_even_parity(cnf, literals, pool)

    first_window = coordinate_variables[: 64 * WINDOW]
    second_window = coordinate_variables[
        64 * args.offset : 64 * (args.offset + WINDOW)
    ]
    cnf.extend(
        CardEnc.atmost(
            lits=first_window,
            bound=LOW_WEIGHT_MAXIMUM,
            vpool=pool,
            encoding=EncType.seqcounter,
        ).clauses
    )
    cnf.extend(
        CardEnc.atmost(
            lits=second_window,
            bound=LOW_WEIGHT_MAXIMUM,
            vpool=pool,
            encoding=EncType.seqcounter,
        ).clauses
    )
    cnf.append(coordinate_variables[:128])

    with Solver(name="cadical195", bootstrap_with=cnf.clauses) as solver:
        satisfiable = solver.solve()
        model = set(literal for literal in solver.get_model() if literal > 0) if satisfiable else set()

    witness = None
    if satisfiable:
        support = [
            coordinate - 1
            for coordinate in coordinate_variables
            if coordinate in model
        ]
        word = sum(1 << coordinate for coordinate in support)
        if any((check & word).bit_count() & 1 for check in checks):
            raise RuntimeError("repeated low windows: syndrome replay failed")
        initial = recover_initial_state(generator, word)
        full_generator = observation_rows(EXTENDED_NODES, parity)
        extended_word = 0
        for bit in range(128):
            if (initial >> bit) & 1:
                extended_word ^= full_generator[bit]
        weights = [
            ((extended_word >> (64 * node)) & MASK64).bit_count()
            for node in range(EXTENDED_NODES)
        ]
        first_weight = sum(weights[:WINDOW])
        second_weight = sum(weights[args.offset : args.offset + WINDOW])
        if first_weight > LOW_WEIGHT_MAXIMUM or second_weight > LOW_WEIGHT_MAXIMUM:
            raise RuntimeError("repeated low windows: exact weight replay failed")
        low_starts = [
            start
            for start in range(EXTENDED_NODES - WINDOW + 1)
            if sum(weights[start : start + WINDOW]) <= LOW_WEIGHT_MAXIMUM
        ]
        witness = {
            "initial_a_hex": f"0x{initial & MASK64:016x}",
            "initial_b_hex": f"0x{initial >> 64:016x}",
            "first_window_weight": first_weight,
            "second_window_weight": second_weight,
            "extended_node_weights": weights,
            "low_window_starts": low_starts,
        }

    payload = {
        "schema": "riffle-packetmul-wrapmul-2lap-g4-repeated-low-windows-v1",
        "candidate": "Riffle PacketMul-WrapMul-2Lap g=4",
        "evidence_label": "BOUNDED_SAT_WITH_EXACT_WITNESS_REPLAY",
        "window_nodes": WINDOW,
        "low_weight_maximum": LOW_WEIGHT_MAXIMUM,
        "second_window_offset": args.offset,
        "solver": "cadical195",
        "result": "SAT" if satisfiable else "UNSAT",
        "exact_witness": witness,
        "scope_limitation": (
            "A SAT witness is replayed through the exact recurrence. An UNSAT "
            "result is not theorem-facing without an independently checked proof."
        ),
    }
    output = (
        ROOT
        / "constructions"
        / "riffle_packetmul_wrapmul_2lap_g4"
        / "receipts"
        / f"goal05_repeated_low_window_offset_{args.offset:03d}.json"
    )
    output.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"offset={args.offset}")
    print(f"result={payload['result']}")
    if witness is not None:
        print(f"first_window_weight={witness['first_window_weight']}")
        print(f"second_window_weight={witness['second_window_weight']}")
        print(f"low_window_starts={witness['low_window_starts']}")
    print(f"output={output}")
    print("status=BOUNDED_REPEATED_LOW_WINDOW_PROBE")


if __name__ == "__main__":
    main()
