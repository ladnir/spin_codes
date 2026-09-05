#!/usr/bin/env python3
"""Exact component-split SAT prototype for one 24-node local-distance case."""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

from pysat.card import CardEnc, EncType
from pysat.solvers import Solver

from analyze_systematic_group_kernel import systematic_state_columns
from analyze_riffle_dp_g4_autonomous_map import (
    accumulate,
    apply_polynomial,
    binary_rank,
    build_apply,
    kernel_basis,
)
from audit_riffle_dp_2lap_g4_component_endpoint_dimension30 import (
    observations,
    support_annihilator,
    support_degrees,
    support_dimension,
)
from probe_riffle_dp_g4_component_mixing import transpose_columns


ROOT = Path(__file__).resolve().parents[1]
LENGTH = 384
WINDOW = 24
REQUIRED_DISTANCE = 97


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def component_basis(transpose_step, mask: int) -> tuple[int, ...]:
    return kernel_basis(
        tuple(
            apply_polynomial(
                transpose_step,
                support_annihilator(mask),
                1 << bit,
            )
            for bit in range(64)
        )
    )


def balanced_split(mask: int) -> tuple[int, int]:
    dimension = support_dimension(mask)
    candidates = []
    submask = (mask - 1) & mask
    while submask:
        complement = mask ^ submask
        if complement and submask < complement:
            left_dimension = support_dimension(submask)
            right_dimension = dimension - left_dimension
            candidates.append(
                (
                    abs(left_dimension - right_dimension),
                    max(left_dimension, right_dimension),
                    submask,
                    complement,
                )
            )
        submask = (submask - 1) & mask
    if not candidates:
        raise RuntimeError("component-split SAT: support cannot be split")
    _, _, left, right = min(candidates)
    if support_dimension(left) < support_dimension(right):
        left, right = right, left
    return left, right


def generator_words(
    states: tuple[int, ...], basis: tuple[int, ...]
) -> tuple[int, ...]:
    return tuple(
        sum(
            (((character & state).bit_count() & 1) << coordinate)
            for coordinate, state in enumerate(states)
        )
        for character in basis
    )


def solve_side(
    *,
    side: str,
    left_words: tuple[int, ...],
    right_words: tuple[int, ...],
    timeout_seconds: float,
) -> dict:
    left_dimension = len(left_words)
    right_dimension = len(right_words)
    character_variables = list(range(1, left_dimension + right_dimension + 1))
    next_variable = character_variables[-1] + 1
    left_output_variables = list(range(next_variable, next_variable + LENGTH))
    next_variable += LENGTH
    right_output_variables = list(range(next_variable, next_variable + LENGTH))
    next_variable += LENGTH
    mismatch_variables = list(range(next_variable, next_variable + LENGTH))
    next_variable += LENGTH

    cardinality_literals = (
        mismatch_variables if side == "weight" else [-variable for variable in mismatch_variables]
    )
    cardinality = CardEnc.atmost(
        lits=cardinality_literals,
        bound=REQUIRED_DISTANCE - 1,
        top_id=next_variable - 1,
        encoding=EncType.seqcounter,
    )
    clauses = list(cardinality.clauses)
    clauses.append(character_variables)

    started = time.perf_counter()
    with Solver(name="cms", bootstrap_with=clauses) as solver:
        for coordinate in range(LENGTH):
            left_literals = [
                index + 1
                for index, word in enumerate(left_words)
                if (word >> coordinate) & 1
            ]
            left_literals.append(left_output_variables[coordinate])
            solver.add_xor_clause(left_literals, value=False)

            right_literals = [
                left_dimension + index + 1
                for index, word in enumerate(right_words)
                if (word >> coordinate) & 1
            ]
            right_literals.append(right_output_variables[coordinate])
            solver.add_xor_clause(right_literals, value=False)

            solver.add_xor_clause(
                [
                    left_output_variables[coordinate],
                    right_output_variables[coordinate],
                    mismatch_variables[coordinate],
                ],
                value=False,
            )

        solver.solver.time_budget(timeout_seconds)
        solver.conf_budget((1 << 31) - 1)
        satisfiable = solver.solve_limited()
        elapsed = time.perf_counter() - started
        try:
            stats = solver.accum_stats()
        except NotImplementedError:
            stats = {}
        model = solver.get_model() if satisfiable is True else None

    row = {
        "side": side,
        "timeout_seconds": timeout_seconds,
        "elapsed_seconds": elapsed,
        "solver": "CryptoMiniSat through python-sat",
        "cardinality_encoding": "sequential_counter",
        "cnf_variable_count": cardinality.nv,
        "cnf_clause_count_before_xors": len(clauses),
        "native_xor_count": 3 * LENGTH,
        "solver_statistics": stats,
        "result": (
            "COUNTEREXAMPLE" if satisfiable is True else "PASS" if satisfiable is False else "UNKNOWN"
        ),
    }
    if model is not None:
        positive = {literal for literal in model if literal > 0}
        left_character = 0
        for index, basis_character in enumerate(left_words):
            if index + 1 in positive:
                left_character ^= basis_character
        right_character = 0
        for index, basis_character in enumerate(right_words):
            if left_dimension + index + 1 in positive:
                right_character ^= basis_character
        row["left_information_hex"] = hex(
            sum((1 << index) for index in range(left_dimension) if index + 1 in positive)
        )
        row["right_information_hex"] = hex(
            sum(
                (1 << index)
                for index in range(right_dimension)
                if left_dimension + index + 1 in positive
            )
        )
        row["left_character_placeholder"] = hex(left_character)
        row["right_character_placeholder"] = hex(right_character)
    return row


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--support", type=lambda value: int(value, 0), required=True)
    parser.add_argument("--packet-value", type=int, default=1)
    parser.add_argument("--timeout-seconds", type=float, default=300.0)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not 1 <= args.support <= 0x7F:
        raise SystemExit("component-split SAT: support must lie in 1..127")
    if not 1 <= args.packet_value <= 14:
        raise SystemExit("component-split SAT: this prototype covers packet values 1..14")
    if args.timeout_seconds <= 0:
        raise SystemExit("component-split SAT: timeout must be positive")

    apply_state = build_apply(systematic_state_columns())

    def step(state: int) -> int:
        return apply_state(accumulate(state))

    transpose_step = build_apply(
        transpose_columns(tuple(step(1 << bit) for bit in range(64)))
    )
    left_mask, right_mask = balanced_split(args.support)
    left_basis = component_basis(transpose_step, left_mask)
    right_basis = component_basis(transpose_step, right_mask)
    total_basis = left_basis + right_basis
    dimension = support_dimension(args.support)
    if (
        len(total_basis) != dimension
        or binary_rank(total_basis) != dimension
        or len(left_basis) != support_dimension(left_mask)
        or len(right_basis) != support_dimension(right_mask)
    ):
        raise RuntimeError("component-split SAT: component basis mismatch")

    states = observations(step, args.packet_value, WINDOW)
    left_code_words = generator_words(states, left_basis)
    right_code_words = generator_words(states, right_basis)
    rows = [
        solve_side(
            side=side,
            left_words=left_code_words,
            right_words=right_code_words,
            timeout_seconds=args.timeout_seconds,
        )
        for side in ("weight", "complement")
    ]

    # Convert the information-coordinate witnesses into actual characters and
    # replay them independently of the SAT output variables.
    for row in rows:
        if row["result"] != "COUNTEREXAMPLE":
            continue
        left_information = int(row.pop("left_information_hex"), 16)
        right_information = int(row.pop("right_information_hex"), 16)
        row.pop("left_character_placeholder")
        row.pop("right_character_placeholder")
        character = 0
        for index, basis_character in enumerate(left_basis):
            if (left_information >> index) & 1:
                character ^= basis_character
        for index, basis_character in enumerate(right_basis):
            if (right_information >> index) & 1:
                character ^= basis_character
        weight = sum((character & state).bit_count() & 1 for state in states)
        side_weight = weight if row["side"] == "weight" else LENGTH - weight
        if not character or side_weight >= REQUIRED_DISTANCE:
            raise RuntimeError("component-split SAT: model replay failed")
        row["character_hex"] = hex(character)
        row["weight"] = weight
        row["complement_weight"] = LENGTH - weight

    result = (
        "COUNTEREXAMPLE"
        if any(row["result"] == "COUNTEREXAMPLE" for row in rows)
        else "PASS"
        if all(row["result"] == "PASS" for row in rows)
        else "UNKNOWN"
    )
    payload = {
        "schema": "riffle-dp-2lap-g4-component-split-sat-v1",
        "candidate": "Riffle DP-2Lap g=4",
        "evidence_label": "EXACT_SAT_RESULT_OR_TIMEOUT_DIAGNOSTIC",
        "source_sha256": digest(Path(__file__).resolve()),
        "support_mask_hex": hex(args.support),
        "support_degrees": list(support_degrees(args.support)),
        "code_dimension": dimension,
        "packet_value": args.packet_value,
        "window_nodes": WINDOW,
        "code_length": LENGTH,
        "required_two_sided_distance": REQUIRED_DISTANCE,
        "left_support_mask_hex": hex(left_mask),
        "left_support_degrees": list(support_degrees(left_mask)),
        "left_dimension": len(left_basis),
        "right_support_mask_hex": hex(right_mask),
        "right_support_degrees": list(support_degrees(right_mask)),
        "right_dimension": len(right_basis),
        "rows": rows,
        "result": result,
        "scope_limitation": (
            "PASS or COUNTEREXAMPLE is an exact solver conclusion with an exact "
            "witness replay for counterexamples. UNKNOWN reports only a timeout. "
            "This prototype does not emit an independently checkable UNSAT proof."
        ),
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"support={hex(args.support)} split={hex(left_mask)}+{hex(right_mask)}")
    for row in rows:
        print(
            f"side={row['side']} result={row['result']} "
            f"elapsed={row['elapsed_seconds']:.6f}"
        )
    print(f"output={args.output}")
    print(f"status={result}")


if __name__ == "__main__":
    main()
