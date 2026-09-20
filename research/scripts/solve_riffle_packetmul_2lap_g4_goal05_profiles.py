#!/usr/bin/env python3
"""Solve Goal 05 by exact four-block weight profiles with native XORs."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import time
from pathlib import Path

from pysat.card import CardEnc, EncType
from pysat.solvers import Solver

from analyze_riffle_dp_g4_autonomous_map import accumulate, build_apply
from analyze_systematic_group_kernel import systematic_state_columns


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / "constructions" / "riffle_packetmul_2lap_g4"
GOAL04 = CANDIDATE / "receipts" / "goal04_four_node_primary.json"
OUTPUT = CANDIDATE / "receipts" / "goal05_distance_primary.json"
DIMENSION = 64
BLOCKS = 4
LENGTH = DIMENSION * BLOCKS


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def step_function():
    apply_p = build_apply(systematic_state_columns())
    return lambda value: accumulate(apply_p(value))


def generators(step) -> tuple[int, ...]:
    result = []
    for bit in range(DIMENSION):
        state = 1 << bit
        word = 0
        for block in range(BLOCKS):
            word |= state << (DIMENSION * block)
            state = step(state)
        result.append(word)
    return tuple(result)


def profiles() -> list[tuple[int, int, int, int]]:
    result = []
    for excess_total in range(6):
        for profile in itertools.product(range(excess_total + 1), repeat=4):
            if sum(profile) == excess_total:
                result.append(tuple(9 + value for value in profile))
    if len(result) != 126:
        raise RuntimeError("Goal 05 profiles: profile count mismatch")
    return result


def replay_goal04(step) -> dict:
    receipt = json.loads(GOAL04.read_text())
    witness = receipt["minimum_witness"]
    states = [
        int(witness["left2_hex"], 16),
        int(witness["left1_hex"], 16),
        int(witness["anchor_hex"], 16),
        int(witness["right1_hex"], 16),
    ]
    if witness["anchor_alignment"] != 1:
        raise RuntimeError("Goal 05 profiles: witness alignment changed")
    if any(step(states[index]) != states[index + 1] for index in range(3)):
        raise RuntimeError("Goal 05 profiles: witness recurrence failed")
    if sum(state.bit_count() for state in states) != 42:
        raise RuntimeError("Goal 05 profiles: witness weight changed")
    return {
        "initial_output_hex": hex(states[0]),
        "block_hex": [hex(state) for state in states],
        "block_weights": [state.bit_count() for state in states],
        "total_weight": 42,
    }


def solve_profile(
    profile: tuple[int, int, int, int],
    code_generators: tuple[int, ...],
    timeout_seconds: float,
) -> tuple[dict, int | None]:
    information_variables = list(range(1, DIMENSION + 1))
    output_variables = list(range(DIMENSION + 1, DIMENSION + LENGTH + 1))
    top_variable = output_variables[-1]
    clauses = [information_variables]
    for block, weight in enumerate(profile):
        literals = output_variables[DIMENSION * block : DIMENSION * (block + 1)]
        equality = CardEnc.equals(
            lits=literals,
            bound=weight,
            top_id=top_variable,
            encoding=EncType.seqcounter,
        )
        clauses.extend(equality.clauses)
        top_variable = equality.nv

    started = time.perf_counter()
    with Solver(name="cms", bootstrap_with=clauses) as solver:
        for coordinate, output_variable in enumerate(output_variables):
            literals = [
                bit + 1
                for bit, word in enumerate(code_generators)
                if (word >> coordinate) & 1
            ]
            solver.add_xor_clause(literals + [output_variable], value=False)
        solver.solver.time_budget(timeout_seconds)
        solver.conf_budget((1 << 31) - 1)
        satisfiable = solver.solve_limited()
        elapsed = time.perf_counter() - started
        try:
            statistics = solver.accum_stats()
        except NotImplementedError:
            statistics = {}
        model = solver.get_model() if satisfiable is True else None

    result = "SAT" if satisfiable is True else "UNSAT" if satisfiable is False else "UNKNOWN"
    row = {
        "block_weight_profile": list(profile),
        "total_weight": sum(profile),
        "timeout_seconds": timeout_seconds,
        "elapsed_seconds": elapsed,
        "result": result,
        "solver_statistics": statistics,
        "cnf_variables": top_variable,
        "cnf_clauses_before_xors": len(clauses),
        "native_xor_constraints": LENGTH,
    }
    if model is None:
        return row, None
    positive = {literal for literal in model if literal > 0}
    information = sum(
        1 << bit for bit in range(DIMENSION) if bit + 1 in positive
    )
    return row, information


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-timeout-seconds", type=float, default=30.0)
    parser.add_argument("--global-timeout-seconds", type=float, default=900.0)
    args = parser.parse_args()
    if args.case_timeout_seconds <= 0 or args.global_timeout_seconds <= 0:
        raise SystemExit("Goal 05 profiles: timeouts must be positive")

    step = step_function()
    code_generators = generators(step)
    witness = replay_goal04(step)
    rows = []
    counterexample = None
    global_started = time.perf_counter()
    for index, profile in enumerate(profiles()):
        remaining = args.global_timeout_seconds - (time.perf_counter() - global_started)
        if remaining <= 0:
            break
        row, information = solve_profile(
            profile,
            code_generators,
            min(args.case_timeout_seconds, remaining),
        )
        rows.append(row)
        print(
            f"case={index + 1}/126 profile={profile} result={row['result']} "
            f"elapsed={row['elapsed_seconds']:.6f}",
            flush=True,
        )
        if information is not None:
            codeword = 0
            for bit, word in enumerate(code_generators):
                if (information >> bit) & 1:
                    codeword ^= word
            block_words = [
                (codeword >> (DIMENSION * block)) & ((1 << DIMENSION) - 1)
                for block in range(BLOCKS)
            ]
            weights = [word.bit_count() for word in block_words]
            if tuple(weights) != profile or sum(weights) > 41 or information == 0:
                raise RuntimeError("Goal 05 profiles: SAT witness replay failed")
            counterexample = {
                "information_hex": hex(information),
                "block_hex": [hex(word) for word in block_words],
                "block_weights": weights,
                "total_weight": sum(weights),
            }
            break

    if counterexample is not None:
        result = "SAT"
    elif len(rows) == 126 and all(row["result"] == "UNSAT" for row in rows):
        result = "UNSAT"
    else:
        result = "UNKNOWN"
    elapsed = time.perf_counter() - global_started
    payload = {
        "schema": "riffle-packetmul-2lap-g4-goal05-profile-native-xor-v1",
        "candidate": "Riffle PacketMul-2Lap g=4",
        "evidence_label": "EXACT_PROFILE_RESULTS_OR_TIMEOUT_DIAGNOSTIC",
        "source_sha256": digest(Path(__file__).resolve()),
        "goal04_primary_sha256": digest(GOAL04),
        "code": {
            "definition": "C4={(o,U(o),U^2(o),U^3(o)):o in F_2^64}",
            "length": LENGTH,
            "dimension": DIMENSION,
            "systematic": True,
        },
        "authenticated_upper_bound_witness": witness,
        "profile_reduction": {
            "reason": (
                "Goal 04 proves that a word of total weight at most 41 cannot "
                "contain a block of weight at most 8. The other three blocks then "
                "force every block weight to lie in 9..14."
            ),
            "profiles": "All ordered four-tuples with every entry at least 9 and sum at most 41.",
            "profile_count": 126,
        },
        "solver": "CryptoMiniSat through python-sat",
        "linear_constraints": "native XOR",
        "cardinality_encoding": "four exact sequential counters per profile",
        "case_timeout_seconds": args.case_timeout_seconds,
        "global_timeout_seconds": args.global_timeout_seconds,
        "elapsed_seconds": elapsed,
        "completed_profile_count": len(rows),
        "unsat_profile_count": sum(row["result"] == "UNSAT" for row in rows),
        "unknown_profile_count": sum(row["result"] == "UNKNOWN" for row in rows),
        "rows": rows,
        "result": result,
        "scope_limitation": (
            "UNSAT rows are exact solver conclusions but this receipt contains no "
            "independently checkable UNSAT proofs. Goal 05 requires an independent "
            "pure-CNF reconstruction. UNKNOWN rows remain open."
        ),
    }
    if counterexample is not None:
        payload["counterexample"] = counterexample
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"result={result}", flush=True)
    print(f"completed_profiles={len(rows)}", flush=True)
    print(f"output={OUTPUT.relative_to(ROOT)}", flush=True)


if __name__ == "__main__":
    main()
