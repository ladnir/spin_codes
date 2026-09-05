#!/usr/bin/env python3
"""Native-XOR probe for a zero-node anchored C18 shortening."""

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


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / "constructions" / "riffle_packetmul_wrapmul_2lap_g4"
NODES = 18
STATE_DIMENSION = 128
SHORTENED_DIMENSION = 64
SHORTENED_LENGTH = 17 * 64
REJECTED_WEIGHT = 105
MASK64 = (1 << 64) - 1


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def xor_selected(rows: list[int] | tuple[int, ...], selector: int) -> int:
    result = 0
    while selector:
        bit = selector & -selector
        result ^= rows[bit.bit_length() - 1]
        selector ^= bit
    return result


def kernel_basis(columns: tuple[int, ...], output_width: int) -> tuple[int, ...]:
    pivot_values = [0] * output_width
    pivot_representations = [0] * output_width
    result = []
    for column, original in enumerate(columns):
        value = original
        representation = 1 << column
        while value:
            pivot = value.bit_length() - 1
            if pivot_values[pivot]:
                value ^= pivot_values[pivot]
                representation ^= pivot_representations[pivot]
            else:
                pivot_values[pivot] = value
                pivot_representations[pivot] = representation
                break
        if value == 0:
            result.append(representation)
    return tuple(result)


def delete_node(word: int, node: int) -> int:
    low_width = 64 * node
    low_mask = (1 << low_width) - 1
    return (word & low_mask) | ((word >> (low_width + 64)) << low_width)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--anchor-node", type=int, required=True)
    parser.add_argument("--timeout-seconds", type=float, default=300.0)
    args = parser.parse_args()
    if not 0 <= args.anchor_node < NODES:
        raise SystemExit("zero anchor: anchor node must be in [0,17]")
    if args.timeout_seconds <= 0:
        raise SystemExit("zero anchor: timeout must be positive")

    parity = build_apply(systematic_state_columns())
    full_rows = observation_rows(NODES, parity)
    anchor_columns = tuple(
        (row >> (64 * args.anchor_node)) & MASK64 for row in full_rows
    )
    state_basis = kernel_basis(anchor_columns, 64)
    if len(state_basis) != SHORTENED_DIMENSION:
        raise RuntimeError(
            f"zero anchor: kernel dimension {len(state_basis)} is not 64"
        )
    generators = tuple(
        delete_node(xor_selected(full_rows, state), args.anchor_node)
        for state in state_basis
    )
    if any(generator.bit_length() > SHORTENED_LENGTH for generator in generators):
        raise RuntimeError("zero anchor: shortened generator length changed")

    information_variables = list(range(1, SHORTENED_DIMENSION + 1))
    output_variables = list(
        range(SHORTENED_DIMENSION + 1, SHORTENED_DIMENSION + SHORTENED_LENGTH + 1)
    )
    cardinality = CardEnc.atmost(
        lits=output_variables,
        bound=REJECTED_WEIGHT,
        top_id=output_variables[-1],
        encoding=EncType.seqcounter,
    )
    clauses = list(cardinality.clauses)
    clauses.append(information_variables)

    started = time.perf_counter()
    with Solver(name="cms", bootstrap_with=clauses) as solver:
        xor_widths = []
        for coordinate, output_variable in enumerate(output_variables):
            inputs = [
                bit + 1
                for bit, generator in enumerate(generators)
                if (generator >> coordinate) & 1
            ]
            xor_widths.append(len(inputs))
            solver.add_xor_clause(inputs + [output_variable], value=False)
        solver.solver.time_budget(args.timeout_seconds)
        solver.conf_budget((1 << 31) - 1)
        satisfiable = solver.solve_limited()
        elapsed = time.perf_counter() - started
        model = solver.get_model() if satisfiable is True else None
        try:
            statistics = solver.accum_stats()
        except NotImplementedError:
            statistics = {}

    result = (
        "SAT" if satisfiable is True else "UNSAT" if satisfiable is False else "UNKNOWN"
    )
    payload = {
        "schema": "riffle-packetmul-wrapmul-2lap-g4-goal05-zero-anchor-v1",
        "candidate": "Riffle PacketMul-WrapMul-2Lap g=4",
        "evidence_label": "EXACT_SAT_RESULT_OR_TIMEOUT_DIAGNOSTIC",
        "source_sha256": digest(Path(__file__).resolve()),
        "anchor_node": args.anchor_node,
        "shortened_code": {
            "length": SHORTENED_LENGTH,
            "dimension": SHORTENED_DIMENSION,
            "zero_anchor_bits_deleted": 64,
        },
        "decision_problem": {
            "state_nonzero": True,
            "anchor_output_zero": True,
            "remaining_output_weight_upper_bound": REJECTED_WEIGHT,
        },
        "encoding": {
            "solver": "CryptoMiniSat through python-sat",
            "linear_constraints": "native XOR",
            "cardinality_encoding": "sequential counter",
            "information_variables": SHORTENED_DIMENSION,
            "output_variables": SHORTENED_LENGTH,
            "cnf_variables": cardinality.nv,
            "cnf_clauses_before_xors": len(clauses),
            "native_xor_constraints": SHORTENED_LENGTH,
            "xor_width_minimum": min(xor_widths),
            "xor_width_maximum": max(xor_widths),
        },
        "timeout_seconds": args.timeout_seconds,
        "elapsed_seconds": elapsed,
        "solver_statistics": statistics,
        "result": result,
        "scope_limitation": (
            "SAT is replayed exactly. UNSAT is an exact solver conclusion but this "
            "receipt does not contain an independently checkable proof. UNKNOWN "
            "records only a timeout."
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
        node_weights = [
            ((full_word >> (64 * node)) & MASK64).bit_count()
            for node in range(NODES)
        ]
        if (
            state == 0
            or node_weights[args.anchor_node] != 0
            or sum(node_weights) > REJECTED_WEIGHT
        ):
            raise RuntimeError("zero anchor: SAT witness replay failed")
        payload["counterexample"] = {
            "initial_a_hex": f"0x{state & MASK64:016x}",
            "initial_b_hex": f"0x{state >> 64:016x}",
            "node_weights": node_weights,
            "total_18_weight": sum(node_weights),
        }

    output = (
        CANDIDATE
        / "receipts"
        / f"goal05_zero_anchor_node_{args.anchor_node:02d}.json"
    )
    output.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"anchor_node={args.anchor_node}", flush=True)
    print(f"result={result}", flush=True)
    print(f"elapsed_seconds={elapsed:.6f}", flush=True)
    if model is not None:
        print(f"counterexample={payload['counterexample']}", flush=True)
    print(f"output={output}", flush=True)
    print("status=ZERO_ANCHOR_SHORTENED_DECISION", flush=True)


if __name__ == "__main__":
    main()
