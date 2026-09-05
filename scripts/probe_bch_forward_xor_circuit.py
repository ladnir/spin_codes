#!/usr/bin/env python3
"""Synthesize low-XOR forward circuits for the committed EBCH [128,64] map.

This probe uses the Paar common-subexpression heuristic.  Every live column
is an already-computed linear form.  When two columns occur together in t
outputs, computing their XOR once and substituting it in those outputs saves
t-1 gates.  The resulting straight-line circuit is verified exactly against
the committed cyclic generator matrix.
"""

from __future__ import annotations

import argparse
import heapq
import json
import random
from dataclasses import dataclass
from pathlib import Path


DIMENSION = 64
LENGTH = 128
GENERATOR_POLYNOMIAL = 0xF4845518B9582A1F


@dataclass
class Circuit:
    gates: list[tuple[int, int]]
    output_signals: list[list[int]]
    xor_count: int
    max_depth: int
    max_live_signals: int


def target_outputs() -> list[int]:
    outputs = [0] * LENGTH
    for message_index in range(DIMENSION):
        product = GENERATOR_POLYNOMIAL << message_index
        product |= 1 << 127
        for output_index in range(LENGTH):
            if product >> output_index & 1:
                outputs[output_index] |= 1 << message_index
    return outputs


def baseline_xors(outputs: list[int]) -> int:
    return sum(max(0, output.bit_count() - 1) for output in outputs)


def transpose_terms(circuit: Circuit) -> list[list[tuple[str, int]]]:
    terms: list[list[tuple[str, int]]] = [
        [] for _ in range(DIMENSION + len(circuit.gates))
    ]
    for gate_output, (left, right) in enumerate(circuit.gates, start=DIMENSION):
        terms[left].append(("adjoint", gate_output))
        terms[right].append(("adjoint", gate_output))
    for output_index, signals in enumerate(circuit.output_signals):
        for signal in signals:
            terms[signal].append(("word", output_index))
    return terms


def transpose_xor_count(circuit: Circuit) -> int:
    return sum(max(0, len(terms) - 1) for terms in transpose_terms(circuit))


def descending_transpose_schedule(circuit: Circuit) -> list[int]:
    """Return the historical reverse-creation order."""
    return list(range(DIMENSION + len(circuit.gates) - 1, -1, -1))


def transpose_live_stats(
    circuit: Circuit, schedule: list[int]
) -> tuple[int, int]:
    """Return peak live adjoints and the sum of live adjoints over the schedule."""
    signal_count = DIMENSION + len(circuit.gates)
    if sorted(schedule) != list(range(signal_count)):
        raise ValueError("transpose schedule is not a permutation of the signals")

    dependencies: list[list[int]] = [[] for _ in range(signal_count)]
    for gate_output, (left, right) in enumerate(circuit.gates, start=DIMENSION):
        dependencies[left].append(gate_output)
        dependencies[right].append(gate_output)

    computed = [False] * signal_count
    remaining_uses = [0] * signal_count
    for gate_output in range(DIMENSION, signal_count):
        remaining_uses[gate_output] = 2

    live = 0
    peak = 0
    area = 0
    for signal in schedule:
        if any(not computed[dependency] for dependency in dependencies[signal]):
            raise ValueError("transpose schedule violates a dependency")
        for dependency in dependencies[signal]:
            remaining_uses[dependency] -= 1
            if remaining_uses[dependency] == 0:
                live -= 1
        computed[signal] = True
        if signal >= DIMENSION:
            live += 1
        peak = max(peak, live)
        area += live
    if live != 0:
        raise AssertionError(f"transpose schedule leaves {live} live adjoints")
    return peak, area


def greedy_transpose_schedule(circuit: Circuit, seed: int = 0) -> list[int]:
    """List-schedule the transpose to release live adjoints early.

    A ready computation is preferred when it consumes the final use of one or
    more live adjoints. Ties prefer computations that make their parents ready,
    then final message outputs, with seeded randomness as the last tie-breaker.
    """
    signal_count = DIMENSION + len(circuit.gates)
    dependencies: list[list[int]] = [[] for _ in range(signal_count)]
    for gate_output, (left, right) in enumerate(circuit.gates, start=DIMENSION):
        dependencies[left].append(gate_output)
        dependencies[right].append(gate_output)

    pending = [len(items) for items in dependencies]
    remaining_uses = [0] * signal_count
    for gate_output in range(DIMENSION, signal_count):
        remaining_uses[gate_output] = 2
    ready = {signal for signal, count in enumerate(pending) if count == 0}
    rng = random.Random(seed)
    schedule: list[int] = []

    while ready:
        scored = []
        for signal in ready:
            frees = sum(
                remaining_uses[dependency] == 1
                for dependency in dependencies[signal]
            )
            creates = int(signal >= DIMENSION)
            unlocks = 0
            if signal >= DIMENSION:
                left, right = circuit.gates[signal - DIMENSION]
                unlocks = int(pending[left] == 1) + int(pending[right] == 1)
            scored.append(
                (
                    frees - creates,
                    frees,
                    unlocks,
                    int(signal < DIMENSION),
                    rng.random(),
                    signal,
                )
            )
        signal = max(scored)[-1]
        ready.remove(signal)
        schedule.append(signal)

        for dependency in dependencies[signal]:
            remaining_uses[dependency] -= 1
        if signal >= DIMENSION:
            left, right = circuit.gates[signal - DIMENSION]
            for parent in (left, right):
                pending[parent] -= 1
                if pending[parent] == 0:
                    ready.add(parent)

    if len(schedule) != signal_count:
        raise AssertionError("transpose list scheduler did not emit every signal")
    transpose_live_stats(circuit, schedule)
    return schedule


def synthesize(
    seed: int, min_overlap: int, target: list[int] | None = None
) -> Circuit:
    target = target_outputs() if target is None else target
    # Column masks say which output equations currently contain each signal.
    columns = [0] * DIMENSION
    for output_index, form in enumerate(target):
        for signal in range(DIMENSION):
            if form >> signal & 1:
                columns[signal] |= 1 << output_index

    rng = random.Random(seed)
    versions = [0] * DIMENSION
    depths = [0] * DIMENSION
    gates: list[tuple[int, int]] = []
    heap: list[tuple[int, float, int, int, int, int]] = []

    def push(left: int, right: int) -> None:
        if left == right:
            return
        if left > right:
            left, right = right, left
        overlap = (columns[left] & columns[right]).bit_count()
        heapq.heappush(
            heap,
            (-overlap, rng.random(), left, right, versions[left], versions[right]),
        )

    for left in range(DIMENSION):
        for right in range(left + 1, DIMENSION):
            push(left, right)

    while heap:
        neg_overlap, _, left, right, left_version, right_version = heapq.heappop(heap)
        if left_version != versions[left] or right_version != versions[right]:
            continue
        common = columns[left] & columns[right]
        overlap = common.bit_count()
        if overlap != -neg_overlap:
            push(left, right)
            continue
        if overlap < min_overlap:
            break

        columns[left] ^= common
        columns[right] ^= common
        versions[left] += 1
        versions[right] += 1

        new_signal = len(columns)
        columns.append(common)
        versions.append(0)
        depths.append(max(depths[left], depths[right]) + 1)
        gates.append((left, right))

        for other in range(new_signal):
            if other != left:
                push(left, other)
            if other != right:
                push(right, other)
            push(new_signal, other)

    output_signals = [[] for _ in range(len(target))]
    for signal, usage in enumerate(columns):
        while usage:
            output_index = (usage & -usage).bit_length() - 1
            output_signals[output_index].append(signal)
            usage &= usage - 1

    xor_count = len(gates) + sum(
        max(0, len(signals) - 1) for signals in output_signals
    )
    max_depth = max(
        (
            max((depths[signal] for signal in signals), default=0)
            + max(0, (len(signals) - 1).bit_length())
            for signals in output_signals
        ),
        default=0,
    )
    last_use = list(range(len(columns)))
    for gate_index, (left, right) in enumerate(gates, start=DIMENSION):
        last_use[left] = max(last_use[left], gate_index)
        last_use[right] = max(last_use[right], gate_index)
    for signals in output_signals:
        ready_at = max(signals, default=0)
        for signal in signals:
            last_use[signal] = max(last_use[signal], ready_at)
    max_live_signals = 0
    for point in range(DIMENSION, len(columns)):
        live = sum(
            (0 if signal < DIMENSION else signal) <= point <= last_use[signal]
            for signal in range(len(columns))
        )
        max_live_signals = max(max_live_signals, live)

    circuit = Circuit(gates, output_signals, xor_count, max_depth, max_live_signals)
    verify(circuit, target)
    verify_transpose(circuit, target)
    return circuit


def verify(circuit: Circuit, target: list[int]) -> None:
    forms = [1 << index for index in range(DIMENSION)]
    for gate_index, (left, right) in enumerate(circuit.gates):
        if left >= len(forms) or right >= len(forms):
            raise AssertionError(f"gate {gate_index} has a forward reference")
        forms.append(forms[left] ^ forms[right])

    actual = []
    for signals in circuit.output_signals:
        form = 0
        for signal in signals:
            form ^= forms[signal]
        actual.append(form)
    if actual != target:
        mismatch = next(i for i, pair in enumerate(zip(actual, target)) if pair[0] != pair[1])
        raise AssertionError(f"circuit verification failed at output {mismatch}")


def verify_transpose(circuit: Circuit, target: list[int]) -> None:
    terms = transpose_terms(circuit)
    adjoints = [0] * len(terms)
    for signal in range(len(terms) - 1, -1, -1):
        form = 0
        for kind, index in terms[signal]:
            form ^= 1 << index if kind == "word" else adjoints[index]
        adjoints[signal] = form

    expected = [0] * DIMENSION
    for output_index, form in enumerate(target):
        for message_index in range(DIMENSION):
            if form >> message_index & 1:
                expected[message_index] |= 1 << output_index
    if adjoints[:DIMENSION] != expected:
        mismatch = next(
            i
            for i, pair in enumerate(zip(adjoints[:DIMENSION], expected))
            if pair[0] != pair[1]
        )
        raise AssertionError(f"transposed circuit verification failed at output {mismatch}")


def write_json(path: Path, circuit: Circuit, seed: int, min_overlap: int) -> None:
    payload = {
        "name": "EBCH [128,64] forward Paar XOR circuit",
        "generator_polynomial": hex(GENERATOR_POLYNOMIAL),
        "seed": seed,
        "min_overlap": min_overlap,
        "xor_count": circuit.xor_count,
        "gate_count": len(circuit.gates),
        "output_xor_count": sum(max(0, len(row) - 1) for row in circuit.output_signals),
        "max_depth_upper": circuit.max_depth,
        "max_live_signals": circuit.max_live_signals,
        "transpose_xor_count": transpose_xor_count(circuit),
        "gates": circuit.gates,
        "outputs": circuit.output_signals,
    }
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def write_header(path: Path, circuit: Circuit) -> None:
    def signal_name(signal: int) -> str:
        if signal < DIMENSION:
            return f"message[{signal}]"
        return f"s{signal}"

    lines = [
        "#pragma once",
        "",
        "#include <cryptoTools/Common/Defines.h>",
        "",
        "namespace osuCrypto",
        "{",
        "\ttemplate<typename T>",
        "\t__declspec(noinline) void ebch128x64DirectForward(",
        "\t\tconst T* __restrict message,",
        "\t\tT* __restrict codeword)",
        "\t{",
    ]
    for output_index, form in enumerate(target_outputs()):
        inputs = [
            f"message[{message_index}]"
            for message_index in range(DIMENSION)
            if form >> message_index & 1
        ]
        lines.append(f"\t\tcodeword[{output_index}] = {' ^ '.join(inputs)};")
    lines.extend([
        "\t}",
        "",
        "\ttemplate<typename T>",
        "\t__declspec(noinline) void ebch128x64PaarForward(",
        "\t\tconst T* __restrict message,",
        "\t\tT* __restrict codeword)",
        "\t{",
    ])
    outputs_at: dict[int, list[int]] = {}
    for output_index, signals in enumerate(circuit.output_signals):
        ready_at = max(signals, default=0)
        outputs_at.setdefault(ready_at, []).append(output_index)

    def emit_outputs(ready_at: int) -> None:
        for output_index in outputs_at.get(ready_at, []):
            signals = circuit.output_signals[output_index]
            expression = " ^ ".join(signal_name(signal) for signal in signals)
            lines.append(f"\t\tcodeword[{output_index}] = {expression};")

    for ready_at in range(DIMENSION):
        emit_outputs(ready_at)
    for gate_index, (left, right) in enumerate(circuit.gates, start=DIMENSION):
        lines.append(
            f"\t\tconst auto s{gate_index} = {signal_name(left)} ^ {signal_name(right)};"
        )
        emit_outputs(gate_index)
    lines.extend(["\t}", "}", ""])
    path.write_text("\n".join(lines), encoding="utf-8")


def write_production_header(
    path: Path, circuit: Circuit, seed: int, min_overlap: int
) -> None:
    def signal_name(signal: int) -> str:
        if signal < DIMENSION:
            return f"message[{signal}]"
        return f"s{signal}"

    lines = [
        "#pragma once",
        "",
        "// Generated by scripts/probe_bch_forward_xor_circuit.py.",
        f"// Exact {circuit.xor_count}-XOR circuit; Paar seed {seed}, minimum overlap {min_overlap}.",
        "",
        "#include <cryptoTools/Common/Defines.h>",
        "",
        "namespace osuCrypto::detail",
        "{",
        "#if defined(_MSC_VER)",
        "#define LIBOTE_RIFFLE_NOINLINE __declspec(noinline)",
        "#elif defined(__GNUC__) || defined(__clang__)",
        "#define LIBOTE_RIFFLE_NOINLINE __attribute__((noinline))",
        "#else",
        "#define LIBOTE_RIFFLE_NOINLINE",
        "#endif",
        "\tLIBOTE_RIFFLE_NOINLINE inline void extendedBch128x64ForwardCircuit(",
        "\t\tconst block* __restrict message,",
        "\t\tblock* __restrict codeword)",
        "\t{",
    ]

    outputs_at: dict[int, list[int]] = {}
    for output_index, signals in enumerate(circuit.output_signals):
        outputs_at.setdefault(max(signals, default=0), []).append(output_index)

    def emit_outputs(ready_at: int) -> None:
        for output_index in outputs_at.get(ready_at, []):
            expression = " ^ ".join(
                signal_name(signal) for signal in circuit.output_signals[output_index]
            )
            lines.append(f"\t\tcodeword[{output_index}] = {expression};")

    for ready_at in range(DIMENSION):
        emit_outputs(ready_at)
    for gate_index, (left, right) in enumerate(circuit.gates, start=DIMENSION):
        lines.append(
            f"\t\tconst block s{gate_index} = {signal_name(left)} ^ {signal_name(right)};"
        )
        emit_outputs(gate_index)

    lines.extend(["\t}", "#undef LIBOTE_RIFFLE_NOINLINE", "}", ""])
    path.write_text("\n".join(lines), encoding="utf-8")


def write_production_transpose_header(
    path: Path,
    circuit: Circuit,
    seed: int,
    min_overlap: int,
    function_name: str = "extendedBch128x64TransposeCircuit",
    split_inputs: bool = False,
) -> None:
    terms = transpose_terms(circuit)

    def term_name(term: tuple[str, int]) -> str:
        kind, index = term
        if kind != "word":
            return f"a{index}"
        if split_inputs:
            return f"word{index // DIMENSION}[{index % DIMENSION}]"
        return f"word[{index}]"

    xor_count = transpose_xor_count(circuit)
    lines = [
        "#pragma once",
        "",
        "// Generated by scripts/probe_bch_forward_xor_circuit.py.",
        f"// Exact {xor_count}-XOR transposition of the seed-{seed} forward circuit.",
        f"// Forward synthesis minimum overlap: {min_overlap}.",
        "",
        "#include <cryptoTools/Common/Defines.h>",
        "",
        "namespace osuCrypto::detail",
        "{",
        "#if defined(_MSC_VER)",
        "#define LIBOTE_RIFFLE_NOINLINE __declspec(noinline)",
        "#elif defined(__GNUC__) || defined(__clang__)",
        "#define LIBOTE_RIFFLE_NOINLINE __attribute__((noinline))",
        "#else",
        "#define LIBOTE_RIFFLE_NOINLINE",
        "#endif",
        f"\tLIBOTE_RIFFLE_NOINLINE inline void {function_name}(",
    ]
    if split_inputs:
        lines.extend([
            "\t\tconst block* __restrict word0,",
            "\t\tconst block* __restrict word1,",
            "\t\tblock* __restrict message)",
            "\t{",
        ])
    else:
        lines.extend([
            "\t\tconst block* __restrict word,",
            "\t\tblock* __restrict message)",
            "\t{",
        ])
    for signal in range(len(terms) - 1, -1, -1):
        expression = " ^ ".join(term_name(term) for term in terms[signal])
        if signal < DIMENSION:
            lines.append(
                f"\t\tmessage[{signal}] = {expression if expression else 'ZeroBlock'};"
            )
        else:
            if not expression:
                raise ValueError("transpose circuit contains an unused internal signal")
            lines.append(f"\t\tconst block a{signal} = {expression};")
    lines.extend(["\t}", "#undef LIBOTE_RIFFLE_NOINLINE", "}", ""])
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trials", type=int, default=32)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--min-overlap", type=int, default=2)
    parser.add_argument("--output-json", type=Path)
    parser.add_argument("--output-header", type=Path)
    parser.add_argument("--output-production-header", type=Path)
    parser.add_argument("--output-production-transpose-header", type=Path)
    parser.add_argument(
        "--transpose-function-name",
        default="extendedBch128x64TransposeCircuit",
    )
    parser.add_argument("--transpose-split-inputs", action="store_true")
    args = parser.parse_args()

    outputs = target_outputs()
    baseline = baseline_xors(outputs)
    best: Circuit | None = None
    best_seed = 0
    counts: list[int] = []
    for trial in range(args.trials):
        seed = args.seed + trial
        circuit = synthesize(seed, args.min_overlap)
        counts.append(circuit.xor_count)
        if best is None or circuit.xor_count < best.xor_count:
            best = circuit
            best_seed = seed

    assert best is not None
    print(f"baseline_xors,{baseline}")
    print(f"trials,{args.trials}")
    print(f"best_seed,{best_seed}")
    print(f"best_xors,{best.xor_count}")
    print(f"reduction,{baseline - best.xor_count}")
    print(f"reduction_percent,{100.0 * (baseline - best.xor_count) / baseline:.3f}")
    print(f"shared_gates,{len(best.gates)}")
    print(
        "output_xors,"
        f"{sum(max(0, len(signals) - 1) for signals in best.output_signals)}"
    )
    print(f"max_depth_upper,{best.max_depth}")
    print(f"max_live_signals,{best.max_live_signals}")
    print(f"transpose_xors,{transpose_xor_count(best)}")
    print(f"min_trial_xors,{min(counts)}")
    print(f"max_trial_xors,{max(counts)}")

    if args.output_json:
        write_json(args.output_json, best, best_seed, args.min_overlap)
        print(f"output_json,{args.output_json}")
    if args.output_header:
        write_header(args.output_header, best)
        print(f"output_header,{args.output_header}")
    if args.output_production_header:
        write_production_header(
            args.output_production_header, best, best_seed, args.min_overlap
        )
        print(f"output_production_header,{args.output_production_header}")
    if args.output_production_transpose_header:
        write_production_transpose_header(
            args.output_production_transpose_header,
            best,
            best_seed,
            args.min_overlap,
            args.transpose_function_name,
            args.transpose_split_inputs,
        )
        print(
            "output_production_transpose_header,"
            f"{args.output_production_transpose_header}"
        )


if __name__ == "__main__":
    main()
