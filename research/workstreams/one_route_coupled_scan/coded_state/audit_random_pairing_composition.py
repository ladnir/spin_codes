"""Bounded exact checks for RANDOM_PAIRING_CONDITIONAL_AUDIT.md.

No search, setup sampling, timing benchmark, or files written. The full-size
calculation counts a specified family; it does not certify the code.
"""

import itertools
import math
import os
import shutil
import threading
import time
from pathlib import Path


def rank_binary(rows):
    basis = {}
    for row in rows:
        while row:
            pivot = row.bit_length() - 1
            if pivot not in basis:
                basis[pivot] = row
                break
            row ^= basis[pivot]
    return len(basis)


def check_identities():
    # All width-four, two-update support histories and matching schedules.
    # Four toy [2,1] outer rows: row i occupies coordinate i at each update.
    matchings = (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2)))
    geometry_checks = exponent_checks = 0
    for hist in itertools.product(range(16), repeat=2):
        states = (0,) + hist
        activity = sum(h.bit_count() for h in hist)
        deaths = sum((a & ~b).bit_count() for a, b in zip(states, hist))
        births = sum((b & ~a).bit_count() for a, b in zip(states, hist))
        overlap = sum((a & b).bit_count() for a, b in zip(states, hist))
        assert births == deaths + hist[-1].bit_count()
        assert overlap == activity - births
        for schedule in itertools.product(matchings, repeat=2):
            live_total = pins = outside = 0
            counts = [0] * 4
            outside_rows = 0
            for before, after, matching in zip(states, hist, schedule):
                live = 0
                for a, b in matching:
                    pair = (1 << a) | (1 << b)
                    if before & pair:
                        live |= pair
                extra = live & ~(before | after)
                cover = live | after
                live_total += live.bit_count()
                pins += (after & ~live).bit_count()
                outside += extra.bit_count()
                outside_rows |= extra
                for i in range(4):
                    counts[i] += (cover >> i) & 1
            assert live_total + pins == activity + deaths + outside
            geometry_checks += 1
            for active in range(1, 16):
                if outside_rows & ~active:
                    continue
                if any(counts[i] < 2 for i in range(4) if active >> i & 1):
                    continue
                inactive_count = 4 - active.bit_count()
                target_cover = sum(counts[i] for i in range(4) if not active >> i & 1)
                dimension = sum(counts[i] - 1 for i in range(4) if active >> i & 1)
                for weight in range(9):
                    zeros = max(0, activity + deaths - weight)
                    if zeros > overlap:
                        assert births + deaths > weight
                        continue
                    lhs = dimension - outside - deaths - zeros
                    assert lhs == activity - 4 + inactive_count - target_cover - zeros
                    if activity + deaths > weight:
                        assert lhs == weight - 4 + inactive_count - target_cover - deaths
                    exponent_checks += 1

    # Zero-output edges and empty-pair pins have disjoint path-component
    # roots. Rescaling along each coordinate path makes edge weights one.
    pin_checks = 0
    for hist in itertools.product(range(4), repeat=3):
        states = (0,) + hist
        vertices = [(t, i) for t, h in enumerate(hist) for i in range(2) if h >> i & 1]
        index = {v: j for j, v in enumerate(vertices)}
        edges, pins = [], []
        for t, (before, after) in enumerate(zip(states, hist)):
            for i in range(2):
                if before & after & (1 << i):
                    edges.append((1 << index[t - 1, i]) | (1 << index[t, i]))
                if before == 0 and after & (1 << i):
                    pins.append(1 << index[t, i])
        for selected in range(1 << len(edges)):
            equations = [e for j, e in enumerate(edges) if selected >> j & 1]
            assert rank_binary(equations) == len(equations)
            assert rank_binary(equations + pins) == len(equations) + len(pins)
            pin_checks += 1
    return geometry_checks, exponent_checks, pin_checks


def count_example():
    n, updates, rows, target_rows = 512, 4096, 8192, 1974
    weight = 46 * (1 << 21) // 100
    activity = updates * 256
    deaths = (updates - 1) * 32
    overlap = (updates - 1) * 256 - deaths
    zeros = activity + deaths - weight
    threshold = 128 * target_rows - ((1 << 20) - weight)
    histories = math.comb(n, 256) * pow(math.comb(256, 32), 2 * (updates - 1))
    choices = math.comb(overlap, zeros)
    row_sets = math.comb(rows, target_rows)
    for name, count in (("histories", histories), ("zero_choices", choices), ("row_sets", row_sets)):
        bits = count.bit_length()
        print(f"log2_{name}: [{bits - 1}, {bits})")
    # 2 is the accepted-map density loss. Target 2^-41 leaves an equal
    # 2^-41 budget for geometry failure when aiming at 40 total bits.
    budgeted_count = histories * choices * row_sets * (1 << 42)
    ceil_log = (budgeted_count - 1).bit_length()
    reserve = (ceil_log + 31) // 32
    assert (1 << (32 * reserve)) >= budgeted_count
    assert (1 << (32 * (reserve - 1))) < budgeted_count
    required_coverage = threshold + reserve
    loss_allowance = activity - 3 * n - 4 * required_coverage
    print(f"W={weight}, A={activity}, D={deaths}, I={overlap}, z={zeros}")
    print(f"threshold={threshold}, sufficient_extra_field_symbols={reserve}")
    print(f"sufficient_coverage={required_coverage}, causal_loss_allowance={loss_allowance}")
    print(f"loss_fraction_diagnostic={loss_allowance / ((updates - 1) * n):.9f}")
    print(f"good_geometry_log2_upper={ceil_log - 41 - 32 * reserve}")


if __name__ == "__main__":
    # Fixed combinatorial sizes keep working storage well below 256 MiB.
    # Independent watchdog limits the entire arithmetic check to 60 seconds.
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free >= 20 * (1 << 30)
    watchdog = threading.Timer(60, lambda: os._exit(124))
    watchdog.daemon = True
    watchdog.start()
    started = time.monotonic()
    try:
        print("geometry, exponent, pin checks:", check_identities())
        count_example()
        print(f"elapsed_seconds={time.monotonic() - started:.3f}")
    finally:
        watchdog.cancel()
