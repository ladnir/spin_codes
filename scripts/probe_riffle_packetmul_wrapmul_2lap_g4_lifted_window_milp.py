#!/usr/bin/env python3
"""MILP probe for the lifted observation-code minimum distance."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix

from analyze_riffle_packetmul_wrapmul_2lap_g4_lifted_windows import (
    build_apply,
    observation_rows,
)
from analyze_systematic_group_kernel import systematic_state_columns


ROOT = Path(__file__).resolve().parents[1]
STATE_BITS = 128


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--nodes", type=int, required=True)
    parser.add_argument("--time-limit", type=float, default=120.0)
    args = parser.parse_args()
    if args.nodes < 2:
        raise SystemExit("lifted MILP probe: nodes must be at least two")

    parity = build_apply(systematic_state_columns())
    generator = observation_rows(args.nodes, parity)
    output_bits = 64 * args.nodes
    output_offset = STATE_BITS
    carry_offset = output_offset + output_bits
    variable_count = carry_offset + output_bits

    row_indices: list[int] = []
    column_indices: list[int] = []
    coefficients: list[float] = []
    lower: list[float] = []
    upper: list[float] = []

    for output in range(output_bits):
        inputs = [
            state
            for state, row in enumerate(generator)
            if (row >> output) & 1
        ]
        equation = len(lower)
        for state in inputs:
            row_indices.append(equation)
            column_indices.append(state)
            coefficients.append(1.0)
        row_indices.extend((equation, equation))
        column_indices.extend((output_offset + output, carry_offset + output))
        coefficients.extend((-1.0, -2.0))
        lower.append(0.0)
        upper.append(0.0)

    # Exclude the zero initial state.
    equation = len(lower)
    for state in range(STATE_BITS):
        row_indices.append(equation)
        column_indices.append(state)
        coefficients.append(1.0)
    lower.append(1.0)
    upper.append(np.inf)

    matrix = coo_matrix(
        (coefficients, (row_indices, column_indices)),
        shape=(len(lower), variable_count),
    ).tocsc()
    objective = np.zeros(variable_count)
    objective[output_offset:carry_offset] = 1.0
    integrality = np.ones(variable_count, dtype=np.uint8)
    variable_lower = np.zeros(variable_count)
    variable_upper = np.ones(variable_count)
    for output in range(output_bits):
        input_count = sum((row >> output) & 1 for row in generator)
        variable_upper[carry_offset + output] = input_count // 2

    result = milp(
        c=objective,
        integrality=integrality,
        bounds=Bounds(variable_lower, variable_upper),
        constraints=LinearConstraint(
            matrix,
            np.asarray(lower),
            np.asarray(upper),
        ),
        options={
            "time_limit": args.time_limit,
            "mip_rel_gap": 0.0,
            "presolve": True,
        },
    )

    witness = None
    if result.x is not None:
        state = sum(
            int(round(result.x[bit])) << bit for bit in range(STATE_BITS)
        )
        word = 0
        for bit in range(STATE_BITS):
            if (state >> bit) & 1:
                word ^= generator[bit]
        node_weights = [
            ((word >> (64 * node)) & ((1 << 64) - 1)).bit_count()
            for node in range(args.nodes)
        ]
        exact_weight = sum(node_weights)
        if state == 0 or (
            result.fun is not None and abs(exact_weight - result.fun) > 1e-6
        ):
            raise RuntimeError("lifted MILP probe: exact replay failed")
        witness = {
            "initial_a_hex": f"0x{state & ((1 << 64) - 1):016x}",
            "initial_b_hex": f"0x{state >> 64:016x}",
            "node_weights": node_weights,
            "exact_weight": exact_weight,
        }

    payload = {
        "schema": "riffle-packetmul-wrapmul-2lap-g4-lifted-window-milp-v1",
        "candidate": "Riffle PacketMul-WrapMul-2Lap g=4",
        "evidence_label": "DIAGNOSTIC_MILP_EXACT_WITNESS_ONLY",
        "nodes": args.nodes,
        "length": output_bits,
        "dimension": STATE_BITS,
        "time_limit_seconds": args.time_limit,
        "solver_status": int(result.status),
        "solver_message": result.message,
        "solver_objective": None if result.fun is None else float(result.fun),
        "solver_dual_bound": getattr(result, "mip_dual_bound", None),
        "solver_mip_gap": getattr(result, "mip_gap", None),
        "solver_mip_node_count": getattr(result, "mip_node_count", None),
        "exact_witness": witness,
        "scope_limitation": (
            "The exact recurrence replays any witness. HiGHS optimality and dual "
            "bounds are diagnostic until independently certified."
        ),
    }
    output = (
        ROOT
        / "constructions"
        / "riffle_packetmul_wrapmul_2lap_g4"
        / "receipts"
        / f"goal04_lifted_window_{args.nodes:02d}_milp.json"
    )
    output.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"nodes={args.nodes}")
    print(f"solver_status={result.status}")
    print(f"solver_objective={result.fun}")
    print(f"solver_dual_bound={getattr(result, 'mip_dual_bound', None)}")
    if witness is not None:
        print(f"witness_weight={witness['exact_weight']}")
        print(f"node_weights={witness['node_weights']}")
    print(f"output={output}")
    print("status=DIAGNOSTIC_LIFTED_WINDOW_MILP")


if __name__ == "__main__":
    main()
