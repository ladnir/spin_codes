#!/usr/bin/env python3
"""Native-XOR probe for a low C24 word with all even outputs zero."""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

from pysat.card import CardEnc, EncType
from pysat.solvers import Solver

from analyze_riffle_packetmul_wrapmul_2lap_g4_lifted_windows import (
    build_apply,
    observation_rows,
)
from analyze_systematic_group_kernel import systematic_state_columns
from probe_riffle_packetmul_wrapmul_2lap_g4_goal05_zero_anchor import (
    MASK64,
    kernel_basis,
    xor_selected,
)


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / "constructions" / "riffle_packetmul_wrapmul_2lap_g4"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def select_nodes(word: int, nodes: tuple[int, ...]) -> int:
    result = 0
    for output_node, input_node in enumerate(nodes):
        result |= ((word >> (64 * input_node)) & MASK64) << (64 * output_node)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nodes", type=int, default=24)
    parser.add_argument("--bound", type=int, default=143)
    parser.add_argument("--timeout-seconds", type=float, default=180.0)
    args = parser.parse_args()
    if args.nodes < 2:
        raise SystemExit("even-zero probe: nodes must be at least two")
    if args.bound < 0 or args.timeout_seconds <= 0:
        raise SystemExit("even-zero probe: invalid bound or timeout")

    parity = build_apply(systematic_state_columns())
    full_rows = observation_rows(args.nodes, parity)
    even_nodes = tuple(range(0, args.nodes, 2))
    odd_nodes = tuple(range(1, args.nodes, 2))
    even_columns = tuple(select_nodes(row, even_nodes) for row in full_rows)
    state_basis = kernel_basis(even_columns, 64 * len(even_nodes))
    generators = tuple(
        select_nodes(xor_selected(full_rows, state), odd_nodes)
        for state in state_basis
    )
    dimension = len(state_basis)
    length = 64 * len(odd_nodes)

    information_variables = list(range(1, dimension + 1))
    output_variables = list(range(dimension + 1, dimension + length + 1))
    cardinality = CardEnc.atmost(
        lits=output_variables,
        bound=args.bound,
        top_id=output_variables[-1],
        encoding=EncType.seqcounter,
    )
    clauses = list(cardinality.clauses)
    clauses.append(information_variables)

    started = time.perf_counter()
    with Solver(name="cms", bootstrap_with=clauses) as solver:
        for coordinate, output_variable in enumerate(output_variables):
            inputs = [
                bit + 1
                for bit, generator in enumerate(generators)
                if (generator >> coordinate) & 1
            ]
            solver.add_xor_clause(inputs + [output_variable], value=False)
        solver.solver.time_budget(args.timeout_seconds)
        solver.conf_budget((1 << 31) - 1)
        satisfiable = solver.solve_limited()
        elapsed = time.perf_counter() - started
        model = solver.get_model() if satisfiable is True else None

    result = (
        "SAT" if satisfiable is True else "UNSAT" if satisfiable is False else "UNKNOWN"
    )
    payload = {
        "schema": "riffle-packetmul-wrapmul-2lap-g4-goal06-even-zero-v1",
        "candidate": "Riffle PacketMul-WrapMul-2Lap g=4",
        "evidence_label": "EXACT_SAT_RESULT_OR_TIMEOUT_DIAGNOSTIC",
        "source_sha256": digest(Path(__file__).resolve()),
        "nodes": args.nodes,
        "bound": args.bound,
        "forced_zero_nodes": list(even_nodes),
        "retained_nodes": list(odd_nodes),
        "subcode_dimension": dimension,
        "subcode_length": length,
        "solver": "CryptoMiniSat through python-sat with native XOR",
        "timeout_seconds": args.timeout_seconds,
        "elapsed_seconds": elapsed,
        "result": result,
        "scope_limitation": (
            "The decision covers only states whose even-node outputs are zero. "
            "SAT witnesses are replayed exactly. UNSAT has no independent proof "
            "certificate in this receipt."
        ),
    }
    if model is not None:
        positive = {literal for literal in model if literal > 0}
        coordinates = sum(
            1 << bit
            for bit, variable in enumerate(information_variables)
            if variable in positive
        )
        state = xor_selected(state_basis, coordinates)
        full_word = xor_selected(full_rows, state)
        weights = [
            ((full_word >> (64 * node)) & MASK64).bit_count()
            for node in range(args.nodes)
        ]
        if state == 0 or any(weights[node] for node in even_nodes) or sum(weights) > args.bound:
            raise RuntimeError("even-zero probe: SAT witness replay failed")
        payload["counterexample"] = {
            "initial_a_hex": f"0x{state & MASK64:016x}",
            "initial_b_hex": f"0x{state >> 64:016x}",
            "node_weights": weights,
            "total_weight": sum(weights),
        }

    output = (
        CANDIDATE
        / "receipts"
        / f"goal06_even_zero_nodes_{args.nodes:02d}_bound_{args.bound:03d}.json"
    )
    output.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"nodes={args.nodes}")
    print(f"dimension={dimension}")
    print(f"result={result}")
    print(f"elapsed_seconds={elapsed:.6f}")
    if model is not None:
        print(f"counterexample={payload['counterexample']}")
    print(f"output={output}")
    print("status=EVEN_ZERO_DECISION")


if __name__ == "__main__":
    main()
