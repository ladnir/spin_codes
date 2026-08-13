#!/usr/bin/env python3
"""Exact first-activation impulse table for one nonzero eight-bit packet.

Both local codes remain EBCH [128,64,22].  Starting from zero recursive state,
place one nonzero packet in a uniformly random one of the eight packet slots,
apply the fixed length-64 accumulator, and then the systematic EBCH map
``v -> (v,Pv)``.  The table is exact over every packet value and slot.

This is a local proof input, not a complete distance certificate: later zero
steps and interactions between multiple occupied packets still require the
packet-profile transfer ledger.
"""

from __future__ import annotations

import hashlib
import math
from collections import Counter

from analyze_systematic_group_kernel import B, MASK, apply, systematic_state_columns


PACKET_BITS = 8
PACKET_SLOTS = B // PACKET_BITS


def accumulate(value: int) -> int:
    parity = 0
    result = 0
    for index in range(B):
        parity ^= (value >> index) & 1
        result |= parity << index
    return result


def main() -> None:
    columns = systematic_state_columns()
    rows: dict[int, Counter[tuple[int, int]]] = {
        weight: Counter() for weight in range(1, PACKET_BITS + 1)
    }
    for slot in range(PACKET_SLOTS):
        for packet in range(1, 1 << PACKET_BITS):
            drive = packet << (PACKET_BITS * slot)
            emitted = accumulate(drive)
            state = apply(columns, emitted)
            emitted_weight = emitted.bit_count()
            state_weight = state.bit_count()
            if emitted_weight + state_weight < 22:
                raise SystemExit("packet8 impulse: EBCH distance check failed")
            rows[packet.bit_count()][(emitted_weight, state_weight)] += 1

    digest = hashlib.sha256()
    print("exact one-packet activation impulse")
    print("inner_local_code=[128,64,22]")
    print(f"packet_bits={PACKET_BITS} packet_slots={PACKET_SLOTS}")
    for weight, histogram in rows.items():
        population = PACKET_SLOTS * math.comb(PACKET_BITS, weight)
        if sum(histogram.values()) != population:
            raise SystemExit(f"packet8 impulse: row {weight} has wrong mass")
        for (emitted, state), count in sorted(histogram.items()):
            digest.update(bytes((weight, emitted, state)))
            digest.update(count.to_bytes(4, "little"))
        mean_emitted = sum(y * count for (y, _q), count in histogram.items()) / population
        mean_state = sum(q * count for (_y, q), count in histogram.items()) / population
        print(
            f"packet_weight={weight} population={population} cells={len(histogram)} "
            f"emitted_min={min(y for y, _q in histogram)} "
            f"emitted_max={max(y for y, _q in histogram)} "
            f"state_min={min(q for _y, q in histogram)} "
            f"state_max={max(q for _y, q in histogram)} "
            f"total_min={min(y+q for y, q in histogram)} "
            f"mean_emitted={mean_emitted:.9f} mean_state={mean_state:.9f}"
        )
    all_ones = accumulate(MASK)
    if all_ones.bit_count() != 32:
        # Prefix accumulation of 64 ones alternates 1,0 and has weight 32.
        raise SystemExit("packet8 impulse: accumulator orientation regression")
    print(f"impulse_table_sha256={digest.hexdigest()}")
    print("all_rows_exact_mass=PASS")
    print("all_nonzero_packets_total_impulse_at_least_22=PASS")
    print("status=EXACT_PACKET8_FIRST_IMPULSE_TABLE")


if __name__ == "__main__":
    main()
