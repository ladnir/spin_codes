#!/usr/bin/env python3
"""MILP probe for short-window nibble distance in the transpose recurrence."""

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
from probe_riffle_dp_g4_component_mixing import transpose_columns


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / "constructions" / "riffle_packetmul_2lap_g4"
MANIFEST = CANDIDATE / "manifest.json"
DEFAULT_OUTPUT = CANDIDATE / "receipts" / "goal02_window_symbol_distance_probe.json"
STATE_BITS = 64
NIBBLES = 16


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def nibble_weight(value: int) -> int:
    return sum(((value >> (4 * nibble)) & 15) != 0 for nibble in range(NIBBLES))


def apply_columns(columns: tuple[int, ...], value: int) -> int:
    result = 0
    while value:
        bit = value & -value
        result ^= columns[bit.bit_length() - 1]
        value ^= bit
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--window", type=int, default=2)
    parser.add_argument("--time-limit", type=float, default=300.0)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.window < 1:
        raise SystemExit("window must be positive")

    apply_p = build_apply(systematic_state_columns())

    def step(state: int) -> int:
        return apply_p(accumulate(state))

    t_columns = tuple(step(1 << bit) for bit in range(STATE_BITS))
    s_columns = transpose_columns(t_columns)
    s_rows = tuple(
        tuple(
            input_bit
            for input_bit in range(STATE_BITS)
            if (s_columns[input_bit] >> output_bit) & 1
        )
        for output_bit in range(STATE_BITS)
    )

    state_variables = args.window * STATE_BITS
    symbol_variables = args.window * NIBBLES
    parity_variables = max(0, args.window - 1) * STATE_BITS
    variable_count = state_variables + symbol_variables + parity_variables

    def state_index(time: int, bit: int) -> int:
        return time * STATE_BITS + bit

    def symbol_index(time: int, nibble: int) -> int:
        return state_variables + time * NIBBLES + nibble

    def parity_index(time: int, output_bit: int) -> int:
        return state_variables + symbol_variables + time * STATE_BITS + output_bit

    row_indices: list[int] = []
    column_indices: list[int] = []
    coefficients: list[float] = []
    lower_bounds: list[float] = []
    upper_bounds: list[float] = []

    def add_constraint(entries: list[tuple[int, float]], lower: float, upper: float) -> None:
        row = len(lower_bounds)
        for column, coefficient in entries:
            row_indices.append(row)
            column_indices.append(column)
            coefficients.append(coefficient)
        lower_bounds.append(lower)
        upper_bounds.append(upper)

    # State transitions: sum_i S[r,i] x_i - y_r - 2k_r = 0.
    for time in range(args.window - 1):
        for output_bit, inputs in enumerate(s_rows):
            entries = [(state_index(time, bit), 1.0) for bit in inputs]
            entries.append((state_index(time + 1, output_bit), -1.0))
            entries.append((parity_index(time, output_bit), -2.0))
            add_constraint(entries, 0.0, 0.0)

    # A symbol indicator dominates every bit in its nibble.
    for time in range(args.window):
        for nibble in range(NIBBLES):
            for bit in range(4):
                add_constraint(
                    [
                        (state_index(time, 4 * nibble + bit), 1.0),
                        (symbol_index(time, nibble), -1.0),
                    ],
                    -np.inf,
                    0.0,
                )

    # The first state is nonzero. Invertibility then makes every state nonzero.
    add_constraint(
        [(state_index(0, bit), 1.0) for bit in range(STATE_BITS)],
        1.0,
        np.inf,
    )

    matrix = coo_matrix(
        (coefficients, (row_indices, column_indices)),
        shape=(len(lower_bounds), variable_count),
    ).tocsc()
    objective = np.zeros(variable_count, dtype=np.float64)
    objective[state_variables : state_variables + symbol_variables] = 1.0
    integrality = np.ones(variable_count, dtype=np.uint8)
    variable_lower = np.zeros(variable_count, dtype=np.float64)
    variable_upper = np.ones(variable_count, dtype=np.float64)
    for time in range(args.window - 1):
        for output_bit, inputs in enumerate(s_rows):
            variable_upper[parity_index(time, output_bit)] = len(inputs) // 2

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
        states = []
        for time in range(args.window):
            state = sum(
                int(rounded[state_index(time, bit)]) << bit
                for bit in range(STATE_BITS)
            )
            states.append(state)
        if states[0] == 0:
            raise RuntimeError("MILP returned the excluded zero state")
        for time in range(args.window - 1):
            if apply_columns(s_columns, states[time]) != states[time + 1]:
                raise RuntimeError("MILP witness violates the exact recurrence")
        exact_weight = sum(nibble_weight(state) for state in states)
        if abs(exact_weight - float(result.fun)) > 1e-6:
            raise RuntimeError("MILP objective differs from exact witness weight")
        witness = {
            "initial_character_hex": hex(states[0]),
            "state_hex": [hex(state) for state in states],
            "state_nibble_weights": [nibble_weight(state) for state in states],
            "exact_window_symbol_weight": exact_weight,
        }

    payload = {
        "schema": "riffle-packetmul-2lap-g4-goal02-window-symbol-distance-probe-v1",
        "candidate": "Riffle PacketMul-2Lap g=4",
        "candidate_id": "riffle_packetmul_2lap_g4",
        "active_manifest_sha256": digest(MANIFEST),
        "evidence_label": "DIAGNOSTIC_MILP_EXACT_WITNESS_ONLY",
        "window_states": args.window,
        "target_weight": 6 * args.window,
        "variable_count": variable_count,
        "constraint_count": len(lower_bounds),
        "solver_status": int(result.status),
        "solver_success": bool(result.success),
        "solver_message": result.message,
        "solver_objective": None if result.fun is None else float(result.fun),
        "solver_mip_gap": getattr(result, "mip_gap", None),
        "solver_mip_node_count": getattr(result, "mip_node_count", None),
        "solver_mip_dual_bound": getattr(result, "mip_dual_bound", None),
        "exact_witness": witness,
        "scope_limitation": (
            "The exact replay authenticates any displayed witness and can refute a "
            "proposed local bound. HiGHS optimality status is not an independently "
            "checkable proof of a lower bound."
        ),
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"window={args.window}")
    print(f"solver_status={result.status}")
    print(f"solver_objective={result.fun}")
    if witness is not None:
        print(f"witness={witness['initial_character_hex']}")
        print(f"weights={witness['state_nibble_weights']}")
    print(f"output={args.output}")
    print("status=DIAGNOSTIC_WINDOW_SYMBOL_DISTANCE")


if __name__ == "__main__":
    main()
