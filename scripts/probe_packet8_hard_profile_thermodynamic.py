#!/usr/bin/env python3
"""Finite-size thermodynamic diagnostic for the hard packet-8 profile.

This samples the exact systematic Riffle inner process conditioned on the
nine packet-weight counts of the first non-closing face.  It deliberately
does not model outer-profile attainability.  Packet order, concrete byte
values within each weight class, and every 64-coordinate state permutation
are part of the sampled ensemble.

Parallel tempering estimates ``<W>_beta`` for emitted weight W and integrates

    log Z(beta) = - integral_0^beta <W>_t dt,

where Z is normalized so Z(0)=1.  This gives the Chernoff diagnostic

    log P[W <= D] <= beta D + log Z(beta).

The result is a mixing-sensitive finite-size diagnostic, not a proof bound.
"""

from __future__ import annotations

import argparse
import collections
import json
import math
import random
from dataclasses import dataclass, field
from pathlib import Path

from analyze_packet8_impulse import accumulate
from analyze_systematic_group_kernel import apply, systematic_state_columns


GROUP_WIDTH = 64
PACKET_WIDTH = 8
PACKETS_PER_GROUP = GROUP_WIDTH // PACKET_WIDTH
FULL_PACKETS = 262144
HARD_PROFILE = (0, 12483, 12483, 12483, 0, 34328, 34328, 34329, 121710)

VALUES_BY_WEIGHT = tuple(
    tuple(value for value in range(256) if value.bit_count() == weight)
    for weight in range(9)
)


def encode_weight_tuple(values: list[int]) -> int:
    code = 0
    scale = 1
    for value in values:
        code += value.bit_count() * scale
        scale *= 9
    return code


@dataclass
class TransitionTrace:
    beta_index: int
    beta: float
    reservoir_limit: int
    rng: random.Random
    observations: int = 0
    state_transition: collections.Counter = field(default_factory=collections.Counter)
    drive_transition: collections.Counter = field(default_factory=collections.Counter)
    packet_transition: collections.Counter = field(default_factory=collections.Counter)
    concrete_transition: collections.Counter = field(default_factory=collections.Counter)
    reservoir: list[dict[str, int]] = field(default_factory=list)

    def observe(self, replica: "Replica", state_columns: tuple[int, ...]) -> None:
        for group in range(len(replica.emitted)):
            base = group * PACKETS_PER_GROUP
            packet_values = replica.packets[base : base + PACKETS_PER_GROUP]
            current_input = 0
            for slot, value in enumerate(packet_values):
                current_input |= value << (slot * PACKET_WIDTH)
            state = replica.states[group]
            permuted_state = permute_bits(state, replica.permutations[group])
            drive = current_input ^ permuted_state
            output = accumulate(drive)
            next_state = apply(state_columns, output)
            if output.bit_count() != replica.emitted[group]:
                raise AssertionError("transition trace emitted-weight mismatch")
            if next_state != replica.states[group + 1]:
                raise AssertionError("transition trace next-state mismatch")

            q = state.bit_count()
            input_weight = current_input.bit_count()
            overlap = (current_input & permuted_state).bit_count()
            drive_weight = drive.bit_count()
            emitted_weight = output.bit_count()
            next_weight = next_state.bit_count()
            weight_code = encode_weight_tuple(packet_values)

            self.state_transition[(q, emitted_weight, next_weight)] += 1
            self.drive_transition[
                (q, input_weight, overlap, drive_weight, emitted_weight, next_weight)
            ] += 1
            self.packet_transition[(q, weight_code, emitted_weight, next_weight)] += 1
            self.concrete_transition[(q, current_input, emitted_weight, next_weight)] += 1

            record = {
                "group": group,
                "state": state,
                "permuted_state": permuted_state,
                "input": current_input,
                "drive": drive,
                "output": output,
                "next_state": next_state,
                "state_weight": q,
                "input_weight": input_weight,
                "overlap": overlap,
                "drive_weight": drive_weight,
                "emitted_weight": emitted_weight,
                "next_state_weight": next_weight,
                "packet_weight_code": weight_code,
            }
            self.observations += 1
            if len(self.reservoir) < self.reservoir_limit:
                self.reservoir.append(record)
            else:
                replacement = self.rng.randrange(self.observations)
                if replacement < self.reservoir_limit:
                    self.reservoir[replacement] = record

    @staticmethod
    def counter_rows(counter: collections.Counter, names: tuple[str, ...]) -> list[dict[str, int]]:
        return [
            {**dict(zip(names, key)), "count": count}
            for key, count in sorted(counter.items())
        ]

    @staticmethod
    def top_counter_rows(
        counter: collections.Counter, names: tuple[str, ...], limit: int
    ) -> list[dict[str, int]]:
        return [
            {**dict(zip(names, key)), "count": count}
            for key, count in counter.most_common(limit)
        ]

    def payload(self, top: int) -> dict[str, object]:
        return {
            "beta_index": self.beta_index,
            "beta": self.beta,
            "observations": self.observations,
            "state_transition": self.counter_rows(
                self.state_transition,
                ("state_weight", "emitted_weight", "next_state_weight"),
            ),
            "drive_transition": self.counter_rows(
                self.drive_transition,
                (
                    "state_weight",
                    "input_weight",
                    "overlap",
                    "drive_weight",
                    "emitted_weight",
                    "next_state_weight",
                ),
            ),
            "packet_transition_top": self.top_counter_rows(
                self.packet_transition,
                ("state_weight", "packet_weight_code", "emitted_weight", "next_state_weight"),
                top,
            ),
            "packet_transition_distinct": len(self.packet_transition),
            "concrete_transition_top": self.top_counter_rows(
                self.concrete_transition,
                ("state_weight", "input", "emitted_weight", "next_state_weight"),
                top,
            ),
            "concrete_transition_distinct": len(self.concrete_transition),
            "reservoir": self.reservoir,
        }


def rounded_profile(groups: int) -> tuple[int, ...]:
    packets = groups * PACKETS_PER_GROUP
    scaled = [count * packets / FULL_PACKETS for count in HARD_PROFILE]
    result = [math.floor(value) for value in scaled]
    remaining = packets - sum(result)
    order = sorted(
        range(9), key=lambda weight: (scaled[weight] - result[weight], weight), reverse=True
    )
    for weight in order[:remaining]:
        result[weight] += 1
    if sum(result) != packets:
        raise AssertionError("rounded profile changed packet count")
    if result[0] or result[4]:
        raise AssertionError("rounding populated a forbidden packet class")
    return tuple(result)


def permute_bits(value: int, permutation: list[int]) -> int:
    result = 0
    for output, source in enumerate(permutation):
        result |= ((value >> source) & 1) << output
    return result


@dataclass
class Replica:
    packets: list[int]
    permutations: list[list[int]]
    states: list[int]
    emitted: list[int]
    energy: int
    rng: random.Random
    label: int
    local_attempts: int = 0
    local_accepts: int = 0

    @classmethod
    def sample(
        cls,
        *,
        profile: tuple[int, ...],
        groups: int,
        rng: random.Random,
        state_columns: tuple[int, ...],
        label: int,
    ) -> "Replica":
        packets = [
            rng.choice(VALUES_BY_WEIGHT[weight])
            for weight, count in enumerate(profile)
            for _ in range(count)
        ]
        rng.shuffle(packets)
        permutations: list[list[int]] = []
        for _ in range(groups):
            permutation = list(range(GROUP_WIDTH))
            rng.shuffle(permutation)
            permutations.append(permutation)
        replica = cls(
            packets=packets,
            permutations=permutations,
            states=[0] * (groups + 1),
            emitted=[0] * groups,
            energy=0,
            rng=rng,
            label=label,
        )
        replica.recompute(0, state_columns)
        return replica

    def recompute(self, start: int, state_columns: tuple[int, ...]) -> None:
        old_suffix = sum(self.emitted[start:])
        state = self.states[start]
        for group in range(start, len(self.emitted)):
            value = 0
            base = group * PACKETS_PER_GROUP
            for slot in range(PACKETS_PER_GROUP):
                value |= self.packets[base + slot] << (slot * PACKET_WIDTH)
            drive = value ^ permute_bits(state, self.permutations[group])
            output = accumulate(drive)
            self.emitted[group] = output.bit_count()
            state = apply(state_columns, output)
            self.states[group + 1] = state
        self.energy += sum(self.emitted[start:]) - old_suffix

    def metropolis_step(self, beta: float, state_columns: tuple[int, ...]) -> None:
        move = self.rng.randrange(3)
        groups = len(self.emitted)
        old_energy = self.energy

        if move == 0:
            first = self.rng.randrange(len(self.packets))
            second = self.rng.randrange(len(self.packets) - 1)
            if second >= first:
                second += 1
            self.packets[first], self.packets[second] = (
                self.packets[second],
                self.packets[first],
            )
            start = min(first, second) // PACKETS_PER_GROUP
            undo = ("swap", first, second)
        elif move == 1:
            index = self.rng.randrange(len(self.packets))
            old_value = self.packets[index]
            choices = VALUES_BY_WEIGHT[old_value.bit_count()]
            if len(choices) == 1:
                return
            old_index = choices.index(old_value)
            candidate_index = self.rng.randrange(len(choices) - 1)
            if candidate_index >= old_index:
                candidate_index += 1
            candidate = choices[candidate_index]
            self.packets[index] = candidate
            start = index // PACKETS_PER_GROUP
            undo = ("value", index, old_value)
        else:
            group = self.rng.randrange(groups)
            first = self.rng.randrange(GROUP_WIDTH)
            second = self.rng.randrange(GROUP_WIDTH - 1)
            if second >= first:
                second += 1
            self.permutations[group][first], self.permutations[group][second] = (
                self.permutations[group][second],
                self.permutations[group][first],
            )
            start = group
            undo = ("perm", group, first, second)

        old_states = self.states[start + 1 :]
        old_emitted = self.emitted[start:]
        self.recompute(start, state_columns)
        delta = self.energy - old_energy
        self.local_attempts += 1
        accept = delta <= 0 or self.rng.random() < math.exp(-beta * delta)
        if accept:
            self.local_accepts += 1
            return

        kind = undo[0]
        if kind == "swap":
            _, first, second = undo
            self.packets[first], self.packets[second] = (
                self.packets[second],
                self.packets[first],
            )
        elif kind == "value":
            _, index, old_value = undo
            self.packets[index] = old_value
        else:
            _, group, first, second = undo
            self.permutations[group][first], self.permutations[group][second] = (
                self.permutations[group][second],
                self.permutations[group][first],
            )
        self.states[start + 1 :] = old_states
        self.emitted[start:] = old_emitted
        self.energy = old_energy


def beta_ladder(replicas: int, beta_max: float) -> list[float]:
    if replicas == 1:
        return [0.0]
    # Quadratic spacing gives the high-temperature region more resolution,
    # where the energy changes fastest and replica overlap is usually worst.
    return [beta_max * (index / (replicas - 1)) ** 2 for index in range(replicas)]


def parse_trace_indices(specification: str, replicas: int) -> list[int]:
    if not specification:
        return []
    result: set[int] = set()
    for token in specification.split(","):
        token = token.strip().lower()
        if token == "last":
            index = replicas - 1
        elif token == "middle":
            index = replicas // 2
        else:
            index = int(token)
            if index < 0:
                index += replicas
        if not 0 <= index < replicas:
            raise ValueError(f"trace beta index {token!r} is out of range")
        result.add(index)
    return sorted(result)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--groups", type=int, default=32)
    parser.add_argument("--replicas", type=int, default=20)
    parser.add_argument("--beta-max", type=float, default=0.6)
    parser.add_argument("--burnin", type=int, default=1000)
    parser.add_argument("--sweeps", type=int, default=4000)
    parser.add_argument("--steps-per-sweep", type=int, default=2)
    parser.add_argument("--seed", type=int, default=20260810)
    parser.add_argument(
        "--trace-indices",
        default="",
        help="comma-separated beta indices, with aliases middle and last",
    )
    parser.add_argument("--trace-stride", type=int, default=1)
    parser.add_argument("--trace-reservoir", type=int, default=1000)
    parser.add_argument("--trace-top", type=int, default=1000)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.groups <= 0:
        raise SystemExit("groups must be positive")
    if args.replicas < 2:
        raise SystemExit("at least two replicas are required")
    if args.beta_max <= 0:
        raise SystemExit("beta max must be positive")
    if args.burnin < 0 or args.sweeps <= 0 or args.steps_per_sweep <= 0:
        raise SystemExit("invalid sweep count")
    if args.trace_stride <= 0 or args.trace_reservoir < 0 or args.trace_top < 0:
        raise SystemExit("invalid trace size")

    profile = rounded_profile(args.groups)
    betas = beta_ladder(args.replicas, args.beta_max)
    try:
        trace_indices = parse_trace_indices(args.trace_indices, args.replicas)
    except ValueError as error:
        raise SystemExit(f"invalid transition trace: {error}") from error
    state_columns = systematic_state_columns()
    master = random.Random(args.seed)
    replicas = [
        Replica.sample(
            profile=profile,
            groups=args.groups,
            rng=random.Random(master.getrandbits(64)),
            state_columns=state_columns,
            label=index,
        )
        for index in range(args.replicas)
    ]
    swap_attempts = [0] * (args.replicas - 1)
    swap_accepts = [0] * (args.replicas - 1)
    energy_sums = [0.0] * args.replicas
    energy_squares = [0.0] * args.replicas
    batches = min(8, args.sweeps)
    batch_sums = [[0.0] * args.replicas for _ in range(batches)]
    batch_counts = [0] * batches
    extreme_state = [0] * args.replicas
    round_trips = [0] * args.replicas
    extreme_state[replicas[0].label] = 1
    traces = [
        TransitionTrace(
            beta_index=index,
            beta=betas[index],
            reservoir_limit=args.trace_reservoir,
            rng=random.Random(args.seed ^ 0x5452414345000000 ^ index),
        )
        for index in trace_indices
    ]
    samples = 0

    total_sweeps = args.burnin + args.sweeps
    for sweep in range(total_sweeps):
        for index, beta in enumerate(betas):
            for _ in range(args.steps_per_sweep):
                replicas[index].metropolis_step(beta, state_columns)
        parity = sweep & 1
        for edge in range(parity, args.replicas - 1, 2):
            swap_attempts[edge] += 1
            left = replicas[edge]
            right = replicas[edge + 1]
            exponent = (betas[edge + 1] - betas[edge]) * (right.energy - left.energy)
            if exponent >= 0 or master.random() < math.exp(exponent):
                replicas[edge], replicas[edge + 1] = right, left
                swap_accepts[edge] += 1
        hot_label = replicas[0].label
        cold_label = replicas[-1].label
        if extreme_state[hot_label] == 0:
            extreme_state[hot_label] = 1
        elif extreme_state[hot_label] == 2:
            round_trips[hot_label] += 1
            extreme_state[hot_label] = 1
        if extreme_state[cold_label] == 1:
            extreme_state[cold_label] = 2
        if sweep >= args.burnin:
            samples += 1
            batch = min(
                batches - 1,
                (sweep - args.burnin) * batches // args.sweeps,
            )
            batch_counts[batch] += 1
            for index, replica in enumerate(replicas):
                energy_sums[index] += replica.energy
                energy_squares[index] += replica.energy * replica.energy
                batch_sums[batch][index] += replica.energy
            if (sweep - args.burnin) % args.trace_stride == 0:
                for trace in traces:
                    trace.observe(replicas[trace.beta_index], state_columns)

    means = [value / samples for value in energy_sums]
    standard_errors = [
        math.sqrt(max(0.0, energy_squares[index] / samples - means[index] ** 2) / samples)
        for index in range(args.replicas)
    ]
    batch_means = [
        [batch_sums[batch][index] / batch_counts[batch] for batch in range(batches)]
        for index in range(args.replicas)
    ]
    batch_standard_errors = [
        math.sqrt(
            sum((value - sum(values) / batches) ** 2 for value in values)
            / (batches * (batches - 1))
        )
        if batches > 1
        else 0.0
        for values in batch_means
    ]
    log_z = [0.0] * args.replicas
    for index in range(1, args.replicas):
        width = betas[index] - betas[index - 1]
        log_z[index] = log_z[index - 1] - 0.5 * width * (means[index - 1] + means[index])

    distance = math.floor(0.09 * GROUP_WIDTH * args.groups)
    bounds = [
        (beta * distance + value) / math.log(2.0)
        for beta, value in zip(betas, log_z)
    ]
    best_index = min(range(args.replicas), key=bounds.__getitem__)
    rows = []
    for index, beta in enumerate(betas):
        swap_rate = None
        if index < args.replicas - 1:
            swap_rate = swap_accepts[index] / swap_attempts[index] if swap_attempts[index] else 0.0
        rows.append(
            {
                "index": index,
                "beta": beta,
                "mean_energy": means[index],
                "naive_standard_error": standard_errors[index],
                "batch_standard_error": batch_standard_errors[index],
                "batch_means": batch_means[index],
                "log_z": log_z[index],
                "chernoff_log2": bounds[index],
                "swap_rate_to_next": swap_rate,
            }
        )

    print("packet-8 hard-profile exact-inner thermodynamic probe")
    print(f"groups={args.groups} packets={args.groups * PACKETS_PER_GROUP}")
    print(f"rounded_profile={','.join(map(str, profile))}")
    print(f"replicas={args.replicas} beta_max={args.beta_max}")
    print(f"burnin={args.burnin} sweeps={args.sweeps} steps_per_sweep={args.steps_per_sweep}")
    print(f"distance={distance}")
    print(
        "swap_rates="
        + ",".join(
            f"{accepted / attempted:.4f}" if attempted else "0.0000"
            for accepted, attempted in zip(swap_accepts, swap_attempts)
        )
    )
    print(f"round_trips_total={sum(round_trips)}")
    print("round_trips=" + ",".join(map(str, round_trips)))
    print(
        f"best_beta={betas[best_index]:.9f} "
        f"best_chernoff_log2={bounds[best_index]:.9f} "
        f"per_group={bounds[best_index] / args.groups:.9f}"
    )
    print(
        "energy_means="
        + ",".join(f"{value:.6f}" for value in means)
    )
    for trace in traces:
        print(
            f"trace_beta_index={trace.beta_index} beta={trace.beta:.9f} "
            f"observations={trace.observations} "
            f"state_cells={len(trace.state_transition)} "
            f"drive_cells={len(trace.drive_transition)} "
            f"packet_cells={len(trace.packet_transition)} "
            f"concrete_cells={len(trace.concrete_transition)}"
        )
    print("status=FINITE_SIZE_MIXING_SENSITIVE_DIAGNOSTIC")

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "groups": args.groups,
            "profile": list(profile),
            "distance": distance,
            "seed": args.seed,
            "burnin": args.burnin,
            "sweeps": args.sweeps,
            "steps_per_sweep": args.steps_per_sweep,
            "best_index": best_index,
            "best_beta": betas[best_index],
            "best_chernoff_log2": bounds[best_index],
            "best_chernoff_per_group": bounds[best_index] / args.groups,
            "rows": rows,
            "round_trips": round_trips,
            "transition_traces": [trace.payload(args.trace_top) for trace in traces],
            "local_acceptance": [
                replica.local_accepts / replica.local_attempts if replica.local_attempts else 0.0
                for replica in replicas
            ],
        }
        args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        print(f"output={args.output}")


if __name__ == "__main__":
    main()
