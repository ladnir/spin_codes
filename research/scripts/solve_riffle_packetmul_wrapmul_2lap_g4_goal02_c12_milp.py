#!/usr/bin/env python3
"""Minimize the 12-node autonomous Hamming weight with an exact-integer MILP."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix

from analyze_systematic_group_kernel import systematic_state_columns
from analyze_riffle_dp_g4_autonomous_map import accumulate, build_apply


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / "constructions" / "riffle_packetmul_wrapmul_2lap_g4"
PART1 = CANDIDATE / "receipts" / "goal02_tail_list_independent.json"
OUTPUT = CANDIDATE / "receipts" / "goal02_c12_milp.json"
STATE_BITS = 64
NODES = 12
TARGET_LOWER_BOUND = 139


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def transpose_columns(columns: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        sum(
            ((columns[input_bit] >> output_bit) & 1) << input_bit
            for input_bit in range(STATE_BITS)
        )
        for output_bit in range(STATE_BITS)
    )


def apply_columns(columns: tuple[int, ...], value: int) -> int:
    result = 0
    while value:
        low = value & -value
        result ^= columns[low.bit_length() - 1]
        value ^= low
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--time-limit", type=float, default=300.0)
    args = parser.parse_args()
    if args.time_limit <= 0:
        raise SystemExit("Goal 02 C12 MILP: time limit must be positive")

    apply_parity = build_apply(systematic_state_columns())

    def step(output: int) -> int:
        return accumulate(apply_parity(output))

    step_columns = tuple(step(1 << bit) for bit in range(STATE_BITS))
    step_rows = transpose_columns(step_columns)
    row_inputs = tuple(
        tuple(bit for bit in range(STATE_BITS) if (row >> bit) & 1)
        for row in step_rows
    )

    state_variables = NODES * STATE_BITS
    parity_variables = (NODES - 1) * STATE_BITS
    variable_count = state_variables + parity_variables

    def state_index(node: int, bit: int) -> int:
        return node * STATE_BITS + bit

    def parity_index(node: int, bit: int) -> int:
        return state_variables + node * STATE_BITS + bit

    row_indices: list[int] = []
    column_indices: list[int] = []
    coefficients: list[float] = []
    lower_bounds: list[float] = []
    upper_bounds: list[float] = []

    def add_constraint(
        entries: list[tuple[int, float]], lower: float, upper: float
    ) -> None:
        row = len(lower_bounds)
        for column, coefficient in entries:
            row_indices.append(row)
            column_indices.append(column)
            coefficients.append(coefficient)
        lower_bounds.append(lower)
        upper_bounds.append(upper)

    for node in range(NODES - 1):
        for output_bit, inputs in enumerate(row_inputs):
            entries = [(state_index(node, bit), 1.0) for bit in inputs]
            entries.append((state_index(node + 1, output_bit), -1.0))
            entries.append((parity_index(node, output_bit), -2.0))
            add_constraint(entries, 0.0, 0.0)

    add_constraint(
        [(state_index(0, bit), 1.0) for bit in range(STATE_BITS)],
        1.0,
        np.inf,
    )
    for block in range(3):
        add_constraint(
            [
                (state_index(node, bit), 1.0)
                for node in range(4 * block, 4 * block + 4)
                for bit in range(STATE_BITS)
            ],
            36.0,
            np.inf,
        )

    matrix = coo_matrix(
        (coefficients, (row_indices, column_indices)),
        shape=(len(lower_bounds), variable_count),
    ).tocsc()
    objective = np.zeros(variable_count, dtype=np.float64)
    objective[:state_variables] = 1.0
    integrality = np.ones(variable_count, dtype=np.uint8)
    variable_lower = np.zeros(variable_count, dtype=np.float64)
    variable_upper = np.ones(variable_count, dtype=np.float64)
    for node in range(NODES - 1):
        for output_bit, inputs in enumerate(row_inputs):
            variable_upper[parity_index(node, output_bit)] = len(inputs) // 2

    result = milp(
        c=objective,
        integrality=integrality,
        bounds=Bounds(variable_lower, variable_upper),
        constraints=LinearConstraint(
            matrix,
            np.asarray(lower_bounds),
            np.asarray(upper_bounds),
        ),
        options={
            "time_limit": args.time_limit,
            "mip_rel_gap": 0.0,
            "presolve": True,
        },
    )

    witness = None
    if result.x is not None:
        rounded = np.rint(result.x[:state_variables]).astype(np.uint8)
        states = [
            sum(
                int(rounded[state_index(node, bit)]) << bit
                for bit in range(STATE_BITS)
            )
            for node in range(NODES)
        ]
        if states[0] == 0:
            raise RuntimeError("Goal 02 C12 MILP: solver returned the zero state")
        if any(step(states[node]) != states[node + 1] for node in range(NODES - 1)):
            raise RuntimeError("Goal 02 C12 MILP: witness violates the recurrence")
        weights = [state.bit_count() for state in states]
        exact_weight = sum(weights)
        if result.fun is not None and abs(exact_weight - float(result.fun)) > 1e-6:
            raise RuntimeError("Goal 02 C12 MILP: objective replay differs")
        witness = {
            "initial_output_hex": hex(states[0]),
            "node_hex": [hex(state) for state in states],
            "node_weights": weights,
            "four_node_weights": [
                sum(weights[4 * block : 4 * block + 4]) for block in range(3)
            ],
            "total_weight": exact_weight,
        }

    lower_bound = getattr(result, "mip_dual_bound", None)
    payload = {
        "schema": "riffle-packetmul-wrapmul-2lap-g4-goal02-c12-milp-v1",
        "candidate": "Riffle PacketMul-WrapMul-2Lap g=4",
        "candidate_id": "riffle_packetmul_wrapmul_2lap_g4",
        "evidence_label": "DIAGNOSTIC_MILP_EXACT_WITNESS_ONLY",
        "source_sha256": digest(Path(__file__).resolve()),
        "goal02_part1_independent_sha256": digest(PART1),
        "code": {
            "definition": "C12={(o,U(o),...,U^11(o)):o in F_2^64}",
            "length": NODES * STATE_BITS,
            "dimension": STATE_BITS,
        },
        "target_minimum_distance": TARGET_LOWER_BOUND,
        "formulation": {
            "binary_state_variables": state_variables,
            "integer_parity_variables": parity_variables,
            "variables": variable_count,
            "constraints": len(lower_bounds),
            "four_node_lower_bounds": [36, 36, 36],
        },
        "time_limit_seconds": args.time_limit,
        "solver_status": int(result.status),
        "solver_success": bool(result.success),
        "solver_message": result.message,
        "solver_objective": None if result.fun is None else float(result.fun),
        "solver_dual_bound": None if lower_bound is None else float(lower_bound),
        "solver_mip_gap": getattr(result, "mip_gap", None),
        "solver_mip_node_count": getattr(result, "mip_node_count", None),
        "exact_witness": witness,
        "scope_limitation": (
            "An exact witness is independently replayable. HiGHS optimality or dual "
            "bounds are solver evidence and are not independently checkable proofs."
        ),
        "result": (
            "OPTIMAL"
            if result.status == 0
            else "INFEASIBLE"
            if result.status == 2
            else "UNKNOWN"
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"solver_status={result.status}")
    print(f"solver_objective={result.fun}")
    print(f"solver_dual_bound={lower_bound}")
    if witness is not None:
        print(f"witness_weight={witness['total_weight']}")
    print(f"output={OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
