#!/usr/bin/env python3
"""Exact product-kernel diagnostic for the first non-closing packet-8 face.

This is the first stage of a concrete outer-codeword attack, not a distance
certificate.  It samples the frozen 256-tile, three-band outer layout,
including the 128 distinct-tile graph holes.  It then anchors the target
number of weight-eight packets by forcing a message delta to vanish on a
uniform set of non-hole packets.

The resulting homogeneous coordinate constraints split over the 16384 data
EBCH blocks.  Exact GF(2) elimination reports the dimension of every local
kernel and their product.  A large product dimension means that the hard
packet profile can be searched without relaxing EBCH membership; it does not
by itself show that the complete target profile is attainable.
"""

from __future__ import annotations

import argparse
import array
import collections
import json
import random
import statistics
from pathlib import Path

from certify_three_band_tile_map import LANES, TILES, tile_labels
from probe_bch_forward_xor_circuit import target_outputs


LOCAL_LENGTH = 128
LOCAL_DIMENSION = 64
BAND_SIZES = (42, 43, 43)
BAND_OFFSETS = (0, 42, 85)
GROUPS_PER_TILE = LOCAL_LENGTH
PACKET_WIDTH = 8
DATA_BLOCKS = TILES * LANES
GROUPS = TILES * GROUPS_PER_TILE
PACKETS = GROUPS * (LANES // PACKET_WIDTH)
GRAPH_HOLES = 128

HARD_PROFILE = (0, 12483, 12483, 12483, 0, 34328, 34328, 34329, 121710)
ANCHOR_PACKETS = HARD_PROFILE[8]


def gf2_rank(rows: list[int]) -> int:
    """Return the exact rank of rows represented by 64-bit integers."""

    basis = [0] * LOCAL_DIMENSION
    rank = 0
    for row in rows:
        value = row
        while value:
            pivot = value.bit_length() - 1
            if basis[pivot]:
                value ^= basis[pivot]
            else:
                basis[pivot] = value
                rank += 1
                break
    return rank


def gf2_nullspace(rows: list[int]) -> list[int]:
    """Return a basis for the nullspace of a 64-column binary matrix."""

    echelon = [0] * LOCAL_DIMENSION
    for row in rows:
        value = row
        while value:
            pivot = value.bit_length() - 1
            if echelon[pivot]:
                value ^= echelon[pivot]
            else:
                echelon[pivot] = value
                break

    result: list[int] = []
    for free in range(LOCAL_DIMENSION):
        if echelon[free]:
            continue
        vector = 1 << free
        for pivot in range(free + 1, LOCAL_DIMENSION):
            row = echelon[pivot]
            if row and ((row & vector).bit_count() & 1):
                vector |= 1 << pivot
        if any((row & vector).bit_count() & 1 for row in rows):
            raise AssertionError("nullspace construction failed")
        result.append(vector)
    return result


def solve_affine(rows: list[int], right_sides: list[int]) -> int:
    """Solve a consistent 64-variable affine system, setting free bits to zero."""

    augmented = [row | ((rhs & 1) << LOCAL_DIMENSION) for row, rhs in zip(rows, right_sides)]
    pivots: dict[int, int] = {}
    for original in augmented:
        value = original
        while value & ((1 << LOCAL_DIMENSION) - 1):
            variables = value & ((1 << LOCAL_DIMENSION) - 1)
            pivot = variables.bit_length() - 1
            if pivot in pivots:
                value ^= pivots[pivot]
            else:
                pivots[pivot] = value
                break
        else:
            if value >> LOCAL_DIMENSION:
                raise ValueError("inconsistent affine system")

    solution = 0
    for pivot in sorted(pivots):
        equation = pivots[pivot]
        lower = equation & ((1 << pivot) - 1)
        rhs = (equation >> LOCAL_DIMENSION) & 1
        bit = rhs ^ ((lower & solution).bit_count() & 1)
        solution |= bit << pivot
    return solution


def apply_local(message: int, outputs: list[int]) -> int:
    word = 0
    for coordinate, form in enumerate(outputs):
        word |= ((form & message).bit_count() & 1) << coordinate
    return word


def sample_punctures(rng: random.Random) -> dict[tuple[int, int], int]:
    """Map each punctured (data block, coordinate) to its graph-hole index."""

    punctures: dict[tuple[int, int], int] = {}
    for hole, tile in enumerate(rng.sample(range(TILES), GRAPH_HOLES)):
        lane = rng.randrange(LANES)
        coordinate = rng.randrange(BAND_SIZES[0])
        punctures[(tile * LANES + lane, coordinate)] = hole
    return punctures


def sample_packets(
    rng: random.Random, punctures: dict[tuple[int, int], int]
) -> tuple[list[list[int]], list[int]]:
    """Sample the exact outer coordinate layout and balanced packet grouping.

    Nonnegative entries encode ``block * 128 + coordinate``.  Negative entries
    encode graph holes.  The graph-coordinate-to-hole bijection is immaterial
    for the homogeneous data-kernel constraints and is therefore not sampled.
    """

    groups: list[list[int]] = [[] for _ in range(GROUPS)]
    for base in range(TILES):
        for lane in range(LANES):
            block = base * LANES + lane
            labels = tile_labels(base, lane)
            for band, size in enumerate(BAND_SIZES):
                columns = list(range(BAND_OFFSETS[band], BAND_OFFSETS[band] + size))
                rng.shuffle(columns)
                for local, column in enumerate(columns):
                    coordinate = BAND_OFFSETS[band] + local
                    hole = punctures.get((block, coordinate))
                    reference = -(hole + 1) if hole is not None else block * LOCAL_LENGTH + coordinate
                    groups[labels[band] * GROUPS_PER_TILE + column].append(reference)

    packets: list[list[int]] = []
    hole_packets: list[int] = []
    for group in groups:
        if len(group) != LANES:
            raise AssertionError("outer layout lost balanced group occupancy")
        rng.shuffle(group)
        for begin in range(0, LANES, PACKET_WIDTH):
            packet = group[begin : begin + PACKET_WIDTH]
            packet_index = len(packets)
            packets.append(packet)
            if any(reference < 0 for reference in packet):
                hole_packets.append(packet_index)

    if len(packets) != PACKETS:
        raise AssertionError("packet count changed")
    if len(hole_packets) != GRAPH_HOLES:
        # Distinct band-zero tiles imply distinct physical groups, hence distinct
        # packets, even after independently shuffling the 64 lanes in each group.
        raise AssertionError("graph holes did not occupy distinct packets")
    return packets, hole_packets


def xor_columns(message: int, columns: list[int]) -> int:
    value = 0
    while message:
        low = message & -message
        value ^= columns[low.bit_length() - 1]
        message ^= low
    return value


def profile_objective(profile: list[int]) -> float:
    """Smooth search score; exact equality is always checked separately."""

    return sum(
        (observed - target) ** 2 / (target + 1024.0)
        for observed, target in zip(profile, HARD_PROFILE)
    )


def profile_l1(profile: list[int]) -> int:
    return sum(abs(value - target) for value, target in zip(profile, HARD_PROFILE))


def move_packet_masks(
    block: int,
    output_delta: int,
    graph_delta: int,
    outputs: list[int],
    data_locations: array.array,
    hole_by_graph_coordinate: list[int],
    hole_locations: list[int],
) -> dict[int, int]:
    """Return packet-index -> byte mask toggled by one local kernel direction."""

    masks: dict[int, int] = {}
    value = output_delta
    while value:
        low = value & -value
        coordinate = low.bit_length() - 1
        location = data_locations[block * LOCAL_LENGTH + coordinate]
        if location != 0xFFFFFFFF:
            packet, bit = divmod(location, PACKET_WIDTH)
            masks[packet] = masks.get(packet, 0) ^ (1 << bit)
        value ^= low

    graph_word = apply_local(graph_delta, outputs)
    while graph_word:
        low = graph_word & -graph_word
        coordinate = low.bit_length() - 1
        hole = hole_by_graph_coordinate[coordinate]
        location = hole_locations[hole]
        packet, bit = divmod(location, PACKET_WIDTH)
        masks[packet] = masks.get(packet, 0) ^ (1 << bit)
        graph_word ^= low
    return {packet: mask for packet, mask in masks.items() if mask}


def greedy_search(
    *,
    rng: random.Random,
    packets: list[list[int]],
    coordinate_sets: list[int],
    outputs: list[int],
    all_one_message: int,
    sweeps: int,
    domain_sweeps: int,
    domain_max_bits: int,
    pair_steps: int,
) -> dict[str, object]:
    """Exact graph-coupled search inside the anchored product kernel.

    The first phase toggles individual nullspace basis vectors.  The optional
    second phase enumerates every vector in local kernels no larger than
    ``domain_max_bits`` and replaces the block by its best complete local
    state.  Both phases update the encoded graph word, so every incumbent is a
    legal full outer word.
    """

    data_locations = array.array("I", [0xFFFFFFFF]) * (DATA_BLOCKS * LOCAL_LENGTH)
    hole_locations = [0] * GRAPH_HOLES
    for packet_index, packet in enumerate(packets):
        for bit, reference in enumerate(packet):
            location = packet_index * PACKET_WIDTH + bit
            if reference < 0:
                hole_locations[-reference - 1] = location
            else:
                data_locations[reference] = location

    graph_coordinates = list(range(LOCAL_LENGTH))
    rng.shuffle(graph_coordinates)
    hole_by_graph_coordinate = [0] * LOCAL_LENGTH
    for hole, coordinate in enumerate(graph_coordinates):
        hole_by_graph_coordinate[coordinate] = hole

    # A random 24 x K graph map, generated blockwise to avoid storing one
    # million Python integers.  Every 24-bit column is independent uniform.
    graph_rng = random.Random(rng.getrandbits(64))
    baseline_graph = 0
    moves: list[tuple[int, int, int]] = []
    block_bases: list[list[tuple[int, int]]] = []
    for block, coordinate_set in enumerate(coordinate_sets):
        rows = [
            outputs[coordinate]
            for coordinate in range(LOCAL_LENGTH)
            if coordinate_set >> coordinate & 1
        ]
        local_basis = gf2_nullspace(rows)
        graph_columns = [graph_rng.getrandbits(24) for _ in range(LOCAL_DIMENSION)]
        baseline_graph ^= xor_columns(all_one_message, graph_columns)
        encoded_basis: list[tuple[int, int]] = []
        for message_delta in local_basis:
            direction = (
                apply_local(message_delta, outputs),
                xor_columns(message_delta, graph_columns),
            )
            encoded_basis.append(direction)
            moves.append((block, direction[0], direction[1]))
        block_bases.append(encoded_basis)

    baseline_graph_word = apply_local(baseline_graph, outputs)
    packet_values = [0] * PACKETS
    for packet_index, packet in enumerate(packets):
        value = 0
        for bit, reference in enumerate(packet):
            if reference >= 0:
                value |= 1 << bit
            else:
                hole = -reference - 1
                coordinate = graph_coordinates[hole]
                value |= ((baseline_graph_word >> coordinate) & 1) << bit
        packet_values[packet_index] = value

    profile = [0] * 9
    for value in packet_values:
        profile[value.bit_count()] += 1
    score = profile_objective(profile)
    initial_profile = tuple(profile)
    initial_score = score
    accepted_total = 0
    sweep_rows: list[dict[str, object]] = []
    current_output = [0] * DATA_BLOCKS
    current_graph = [0] * DATA_BLOCKS

    order = list(range(len(moves)))
    for sweep in range(sweeps):
        rng.shuffle(order)
        accepted = 0
        for move_index in order:
            block, output_delta, graph_delta = moves[move_index]
            masks = move_packet_masks(
                block,
                output_delta,
                graph_delta,
                outputs,
                data_locations,
                hole_by_graph_coordinate,
                hole_locations,
            )
            changes = [0] * 9
            for packet, mask in masks.items():
                old_weight = packet_values[packet].bit_count()
                new_weight = (packet_values[packet] ^ mask).bit_count()
                changes[old_weight] -= 1
                changes[new_weight] += 1
            candidate_profile = [value + delta for value, delta in zip(profile, changes)]
            candidate_score = profile_objective(candidate_profile)
            if candidate_score + 1e-12 < score:
                for packet, mask in masks.items():
                    packet_values[packet] ^= mask
                profile = candidate_profile
                score = candidate_score
                current_output[block] ^= output_delta
                current_graph[block] ^= graph_delta
                accepted += 1
        accepted_total += accepted
        row = {
            "sweep": sweep + 1,
            "accepted": accepted,
            "score": score,
            "l1": profile_l1(profile),
            "profile": list(profile),
        }
        sweep_rows.append(row)
        print(
            f"search_sweep={sweep + 1} accepted={accepted} score={score:.9f} "
            f"l1={row['l1']} profile={','.join(map(str, profile))}",
            flush=True,
        )
        if not accepted or tuple(profile) == HARD_PROFILE:
            break

    domain_rows: list[dict[str, object]] = []
    eligible_blocks = [
        block
        for block, basis in enumerate(block_bases)
        if 0 < len(basis) <= domain_max_bits
    ]
    for domain_sweep in range(domain_sweeps):
        rng.shuffle(eligible_blocks)
        accepted_blocks = 0
        enumerated_states = 0
        for block in eligible_blocks:
            basis = block_bases[block]
            basis_masks = [
                move_packet_masks(
                    block,
                    output_delta,
                    graph_delta,
                    outputs,
                    data_locations,
                    hole_by_graph_coordinate,
                    hole_locations,
                )
                for output_delta, graph_delta in basis
            ]
            current_masks = move_packet_masks(
                block,
                current_output[block],
                current_graph[block],
                outputs,
                data_locations,
                hole_by_graph_coordinate,
                hole_locations,
            )
            affected = sorted(
                set(current_masks).union(
                    *(set(masks) for masks in basis_masks)
                )
            )
            local_index = {packet: index for index, packet in enumerate(affected)}
            base_values = [
                packet_values[packet] ^ current_masks.get(packet, 0)
                for packet in affected
            ]
            candidate_values = list(base_values)
            outside_profile = list(profile)
            for packet in affected:
                outside_profile[packet_values[packet].bit_count()] -= 1
            local_profile = [0] * 9
            for value in candidate_values:
                local_profile[value.bit_count()] += 1

            best_gray = 0
            best_profile = [
                outside + local
                for outside, local in zip(outside_profile, local_profile)
            ]
            best_score = profile_objective(best_profile)
            previous_gray = 0
            for state in range(1, 1 << len(basis)):
                gray = state ^ (state >> 1)
                changed = gray ^ previous_gray
                direction = (changed & -changed).bit_length() - 1
                for packet, mask in basis_masks[direction].items():
                    index = local_index[packet]
                    old_weight = candidate_values[index].bit_count()
                    candidate_values[index] ^= mask
                    new_weight = candidate_values[index].bit_count()
                    local_profile[old_weight] -= 1
                    local_profile[new_weight] += 1
                candidate_profile = [
                    outside + local
                    for outside, local in zip(outside_profile, local_profile)
                ]
                candidate_score = profile_objective(candidate_profile)
                if candidate_score + 1e-12 < best_score:
                    best_score = candidate_score
                    best_gray = gray
                    best_profile = candidate_profile
                previous_gray = gray
            enumerated_states += 1 << len(basis)

            if best_score + 1e-12 < score:
                best_output = 0
                best_graph = 0
                best_masks: dict[int, int] = {}
                for direction, (output_delta, graph_delta) in enumerate(basis):
                    if best_gray >> direction & 1:
                        best_output ^= output_delta
                        best_graph ^= graph_delta
                        for packet, mask in basis_masks[direction].items():
                            best_masks[packet] = best_masks.get(packet, 0) ^ mask
                for packet, base_value in zip(affected, base_values):
                    packet_values[packet] = base_value ^ best_masks.get(packet, 0)
                current_output[block] = best_output
                current_graph[block] = best_graph
                profile = best_profile
                score = best_score
                accepted_blocks += 1

        row = {
            "sweep": domain_sweep + 1,
            "max_bits": domain_max_bits,
            "eligible_blocks": len(eligible_blocks),
            "enumerated_states": enumerated_states,
            "accepted_blocks": accepted_blocks,
            "score": score,
            "l1": profile_l1(profile),
            "profile": list(profile),
        }
        domain_rows.append(row)
        print(
            f"domain_sweep={domain_sweep + 1} max_bits={domain_max_bits} "
            f"eligible_blocks={len(eligible_blocks)} states={enumerated_states} "
            f"accepted={accepted_blocks} score={score:.9f} l1={row['l1']} "
            f"profile={','.join(map(str, profile))}",
            flush=True,
        )
        if not accepted_blocks or tuple(profile) == HARD_PROFILE:
            break

    pair_accepted = 0
    pair_start_score = score
    pair_start_l1 = profile_l1(profile)
    for _step in range(pair_steps):
        first_index = rng.randrange(len(moves))
        second_index = rng.randrange(len(moves) - 1)
        if second_index >= first_index:
            second_index += 1
        first_block, first_output, first_graph = moves[first_index]
        second_block, second_output, second_graph = moves[second_index]
        if first_block == second_block:
            continue
        masks = move_packet_masks(
            first_block,
            first_output,
            first_graph,
            outputs,
            data_locations,
            hole_by_graph_coordinate,
            hole_locations,
        )
        second_masks = move_packet_masks(
            second_block,
            second_output,
            second_graph,
            outputs,
            data_locations,
            hole_by_graph_coordinate,
            hole_locations,
        )
        for packet, mask in second_masks.items():
            masks[packet] = masks.get(packet, 0) ^ mask
            if not masks[packet]:
                del masks[packet]
        changes = [0] * 9
        for packet, mask in masks.items():
            old_weight = packet_values[packet].bit_count()
            new_weight = (packet_values[packet] ^ mask).bit_count()
            changes[old_weight] -= 1
            changes[new_weight] += 1
        candidate_profile = [value + delta for value, delta in zip(profile, changes)]
        candidate_score = profile_objective(candidate_profile)
        if candidate_score + 1e-12 < score:
            for packet, mask in masks.items():
                packet_values[packet] ^= mask
            current_output[first_block] ^= first_output
            current_graph[first_block] ^= first_graph
            current_output[second_block] ^= second_output
            current_graph[second_block] ^= second_graph
            profile = candidate_profile
            score = candidate_score
            pair_accepted += 1

    pair_row = {
        "steps": pair_steps,
        "accepted": pair_accepted,
        "start_score": pair_start_score,
        "start_l1": pair_start_l1,
        "score": score,
        "l1": profile_l1(profile),
        "profile": list(profile),
    }
    if pair_steps:
        print(
            f"pair_steps={pair_steps} accepted={pair_accepted} "
            f"score={score:.9f} l1={pair_row['l1']} "
            f"profile={','.join(map(str, profile))}",
            flush=True,
        )

    return {
        "moves": len(moves),
        "initial_profile": list(initial_profile),
        "initial_score": initial_score,
        "final_profile": list(profile),
        "final_score": score,
        "final_l1": sum(abs(value - target) for value, target in zip(profile, HARD_PROFILE)),
        "accepted_moves": accepted_total,
        "sweeps": sweep_rows,
        "domain_max_bits": domain_max_bits,
        "domain_sweeps": domain_rows,
        "pair_search": pair_row,
        "exact_profile": tuple(profile) == HARD_PROFILE,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=20260810)
    parser.add_argument("--anchors", type=int, default=ANCHOR_PACKETS)
    parser.add_argument("--search-sweeps", type=int, default=0)
    parser.add_argument("--domain-sweeps", type=int, default=0)
    parser.add_argument("--domain-max-bits", type=int, default=6)
    parser.add_argument("--pair-steps", type=int, default=0)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if not 0 <= args.anchors <= PACKETS - GRAPH_HOLES:
        raise SystemExit("anchor count exceeds the non-hole packet population")
    if args.search_sweeps < 0:
        raise SystemExit("search sweeps must be nonnegative")
    if args.domain_sweeps < 0:
        raise SystemExit("domain sweeps must be nonnegative")
    if not 1 <= args.domain_max_bits <= 16:
        raise SystemExit("domain max bits must lie in [1,16]")
    if args.pair_steps < 0:
        raise SystemExit("pair steps must be nonnegative")

    outputs = target_outputs()
    all_one_message = solve_affine(outputs, [1] * LOCAL_LENGTH)
    if apply_local(all_one_message, outputs) != (1 << LOCAL_LENGTH) - 1:
        raise AssertionError("failed to recover the all-one EBCH message")

    rng = random.Random(args.seed)
    punctures = sample_punctures(rng)
    packets, hole_packets = sample_packets(rng, punctures)
    forbidden = set(hole_packets)
    candidates = [index for index in range(PACKETS) if index not in forbidden]
    selected = rng.sample(candidates, args.anchors)

    coordinate_sets = [0] * DATA_BLOCKS
    for packet_index in selected:
        packet = packets[packet_index]
        if any(reference < 0 for reference in packet):
            raise AssertionError("selected anchor contains a graph hole")
        for reference in packet:
            block, coordinate = divmod(reference, LOCAL_LENGTH)
            coordinate_sets[block] |= 1 << coordinate

    ranks: list[int] = []
    constraints: list[int] = []
    nullities: list[int] = []
    for coordinate_set in coordinate_sets:
        coordinates: list[int] = []
        value = coordinate_set
        while value:
            low = value & -value
            coordinates.append(low.bit_length() - 1)
            value ^= low
        rank = gf2_rank([outputs[coordinate] for coordinate in coordinates])
        ranks.append(rank)
        constraints.append(len(coordinates))
        nullities.append(LOCAL_DIMENSION - rank)

    constraint_histogram = collections.Counter(constraints)
    nullity_histogram = collections.Counter(nullities)
    total_constraints = sum(constraints)
    total_nullity = sum(nullities)

    print("packet-8 hard-profile product-kernel probe")
    print(f"seed={args.seed}")
    print(f"layout=FROZEN_256_TILE_THREE_BAND_WITH_DISTINCT_GRAPH_HOLES")
    print(f"packets={PACKETS} graph_hole_packets={len(hole_packets)}")
    print(f"target_profile={','.join(map(str, HARD_PROFILE))}")
    print(f"anchored_weight8_packets={args.anchors}")
    print(f"anchored_coordinates={total_constraints}")
    print(f"all_one_local_message=0x{all_one_message:016x}")
    print(
        "constraints_per_block="
        f"min:{min(constraints)},median:{statistics.median(constraints):g},"
        f"mean:{statistics.fmean(constraints):.9f},max:{max(constraints)}"
    )
    print(
        "rank_per_block="
        f"min:{min(ranks)},median:{statistics.median(ranks):g},"
        f"mean:{statistics.fmean(ranks):.9f},max:{max(ranks)}"
    )
    print(
        "nullity_per_block="
        f"min:{min(nullities)},median:{statistics.median(nullities):g},"
        f"mean:{statistics.fmean(nullities):.9f},max:{max(nullities)}"
    )
    print(f"product_kernel_dimension={total_nullity}")
    print(f"rigid_blocks={sum(value == 0 for value in nullities)}")
    print(
        "constraint_histogram="
        + ",".join(f"{key}:{constraint_histogram[key]}" for key in sorted(constraint_histogram))
    )
    print(
        "nullity_histogram="
        + ",".join(f"{key}:{nullity_histogram[key]}" for key in sorted(nullity_histogram))
    )
    search_result: dict[str, object] | None = None
    if args.search_sweeps or args.domain_sweeps or args.pair_steps:
        search_result = greedy_search(
            rng=rng,
            packets=packets,
            coordinate_sets=coordinate_sets,
            outputs=outputs,
            all_one_message=all_one_message,
            sweeps=args.search_sweeps,
            domain_sweeps=args.domain_sweeps,
            domain_max_bits=args.domain_max_bits,
            pair_steps=args.pair_steps,
        )
        print(f"search_moves={search_result['moves']}")
        print(
            "search_initial_profile="
            + ",".join(map(str, search_result["initial_profile"]))
        )
        print(
            "search_final_profile="
            + ",".join(map(str, search_result["final_profile"]))
        )
        print(f"search_final_l1={search_result['final_l1']}")
        print(f"search_exact_profile={search_result['exact_profile']}")
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "seed": args.seed,
            "target_profile": list(HARD_PROFILE),
            "anchors": args.anchors,
            "product_kernel_dimension": total_nullity,
            "rigid_blocks": sum(value == 0 for value in nullities),
            "constraint_histogram": dict(sorted(constraint_histogram.items())),
            "nullity_histogram": dict(sorted(nullity_histogram.items())),
            "search": search_result,
        }
        args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        print(f"output={args.output}")
    print("status=EXACT_INSTANCE_KERNEL_DIMENSION_DIAGNOSTIC")


if __name__ == "__main__":
    main()
