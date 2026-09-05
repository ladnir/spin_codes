#!/usr/bin/env python3
"""Native-XOR decision probe for an 18-node return-cost violation."""

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
OUTPUT_STEM = "goal05_return_native_xor"
NODES = 18
LENGTH = 64 * NODES
DIMENSION = 128
REJECTED_WEIGHT = 105
MASK64 = (1 << 64) - 1


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def xor_rows(rows: list[int], state: int) -> int:
    result = 0
    while state:
        bit = state & -state
        result ^= rows[bit.bit_length() - 1]
        state ^= bit
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--timeout-seconds", type=float, default=300.0)
    parser.add_argument(
        "--low-half", choices=("none", "first", "second"), default="none"
    )
    parser.add_argument("--low-node", type=int)
    parser.add_argument("--low-node-maximum", type=int, default=5)
    args = parser.parse_args()
    if args.timeout_seconds <= 0:
        raise SystemExit("return native XOR: timeout must be positive")
    if args.low_node is not None and not 0 <= args.low_node < NODES:
        raise SystemExit("return native XOR: low node must be in [0,17]")
    if not 0 <= args.low_node_maximum <= 64:
        raise SystemExit("return native XOR: low node maximum must be in [0,64]")

    parity = build_apply(systematic_state_columns())
    generators = observation_rows(NODES, parity)
    known_state = int("ba22ae96294ba476", 16) | (
        int("961e658de73963d2", 16) << 64
    )
    known_word = xor_rows(generators, known_state)
    known_node_weights = [
        ((known_word >> (64 * node)) & MASK64).bit_count()
        for node in range(NODES)
    ]
    if sum(known_node_weights[:9]) != 44 or sum(known_node_weights[9:]) != 275:
        raise RuntimeError("return native XOR: known witness replay changed")

    information_variables = list(range(1, DIMENSION + 1))
    output_variables = list(range(DIMENSION + 1, DIMENSION + LENGTH + 1))
    cardinality = CardEnc.atmost(
        lits=output_variables,
        bound=REJECTED_WEIGHT,
        top_id=output_variables[-1],
        encoding=EncType.seqcounter,
    )
    clauses = list(cardinality.clauses)
    top_variable = cardinality.nv
    if args.low_half != "none":
        split = 64 * 9
        low_variables = (
            output_variables[:split]
            if args.low_half == "first"
            else output_variables[split:]
        )
        low_cardinality = CardEnc.atmost(
            lits=low_variables,
            bound=52,
            top_id=top_variable,
            encoding=EncType.seqcounter,
        )
        clauses.extend(low_cardinality.clauses)
        top_variable = low_cardinality.nv
    if args.low_node is not None:
        node_variables = output_variables[
            64 * args.low_node : 64 * (args.low_node + 1)
        ]
        if args.low_node_maximum == 0:
            clauses.extend([-variable] for variable in node_variables)
        else:
            node_cardinality = CardEnc.atmost(
                lits=node_variables,
                bound=args.low_node_maximum,
                top_id=top_variable,
                encoding=EncType.seqcounter,
            )
            clauses.extend(node_cardinality.clauses)
            top_variable = node_cardinality.nv
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
    payload = {
        "schema": "riffle-packetmul-wrapmul-2lap-g4-goal05-return-native-xor-v1",
        "candidate": "Riffle PacketMul-WrapMul-2Lap g=4",
        "evidence_label": "EXACT_SAT_RESULT_OR_TIMEOUT_DIAGNOSTIC",
        "source_sha256": digest(Path(__file__).resolve()),
        "code": {
            "definition": "C18={q_0(s),...,q_17(s):s in F_2^128}",
            "length": LENGTH,
            "dimension": DIMENSION,
        },
        "decision_problem": {
            "initial_lifted_state_nonzero": True,
            "output_weight_upper_bound": REJECTED_WEIGHT,
            "required_half_of_weight_at_most_52": args.low_half,
            "required_node_of_weight_at_most_5": args.low_node,
            "required_node_weight_maximum": (
                args.low_node_maximum if args.low_node is not None else None
            ),
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
        },
        "known_weight_44_prefix_replay": {
            "initial_a_hex": f"0x{known_state & MASK64:016x}",
            "initial_b_hex": f"0x{known_state >> 64:016x}",
            "node_weights": known_node_weights,
            "total_18_node_weight": sum(known_node_weights),
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
        state = sum(
            1 << bit
            for bit in range(DIMENSION)
            if information_variables[bit] in positive
        )
        word = xor_rows(generators, state)
        weight = word.bit_count()
        if state == 0 or weight > REJECTED_WEIGHT:
            raise RuntimeError("return native XOR: SAT model replay failed")
        node_weights = [
            ((word >> (64 * node)) & MASK64).bit_count()
            for node in range(NODES)
        ]
        payload["counterexample"] = {
            "initial_a_hex": f"0x{state & MASK64:016x}",
            "initial_b_hex": f"0x{state >> 64:016x}",
            "node_weights": node_weights,
            "first_9_weight": sum(node_weights[:9]),
            "second_9_weight": sum(node_weights[9:]),
            "total_18_weight": weight,
        }

    suffix = "" if args.low_half == "none" else f"_{args.low_half}_low"
    if args.low_node is not None:
        suffix += (
            f"_node_{args.low_node:02d}_max_{args.low_node_maximum:02d}"
        )
    output = CANDIDATE / "receipts" / f"{OUTPUT_STEM}{suffix}.json"
    output.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"result={result}", flush=True)
    print(f"elapsed_seconds={elapsed:.6f}", flush=True)
    if model is not None:
        print(f"counterexample={payload['counterexample']}", flush=True)
    print(f"output={output}", flush=True)
    print("status=RETURN_NATIVE_XOR_DECISION", flush=True)


if __name__ == "__main__":
    main()
