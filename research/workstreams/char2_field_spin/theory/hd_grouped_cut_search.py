"""Targeted grouped-L4 interval-cut diagnostics; never a distance certificate.

Counts and row-set changes are exact. The objective uses the floating
stationary interval pressure from hd_affine_cut. Search is deliberately
bounded and seeded by the previous closest Fourier-threshold row sets.
"""
import argparse
from functools import lru_cache
import json

import numpy as np

from hd_affine_cut import grouped_safe_counts, screen_profile, stationary_pressure
from hd_two_state_regional_strides import ROWS, STRIDES


def local_interior_options(first, last):
    """Best within-group edge count for each number of interior bits."""
    return tuple(max(first*x+x*y+y*last
                     for x in (0, 1) for y in (0, 1) if x+y == count)
                 for count in range(3))


def flip_delta(selected, row, strides=STRIDES):
    """Exact change in every regional safe count when one source row flips."""
    blocks = np.asarray(selected, dtype=bool).reshape((-1, 4))
    groups = len(blocks)
    g, lane = divmod(int(row), 4)
    sign = 1-2*int(blocks[g, lane])
    within = ((int(blocks[g, lane-1]) if lane else 0)
              + (int(blocks[g, lane+1]) if lane < 3 else 0))
    delta = np.full(len(strides), sign*within, dtype=np.int64)
    shift = np.asarray(strides, dtype=np.int64)
    if lane == 0:
        delta += sign*blocks[(g-shift) % groups, 3]
        if g == 0:
            delta[-1] -= sign*int(blocks[-strides[-1] % groups, 3])
    if lane == 3:
        delta += sign*blocks[(g+shift) % groups, 0]
        if g == -strides[-1] % groups:
            delta[-1] += sign*int(not blocks[0, 0])
    return delta


def change_delta(selected, rows, strides=STRIDES):
    """Exact multi-flip delta; restore the caller's row set before returning."""
    delta = np.zeros(len(strides), dtype=np.int64)
    for row in rows:
        delta += flip_delta(selected, row, strides)
        selected[row] = ~selected[row]
    for row in reversed(rows):
        selected[row] = ~selected[row]
    return delta


def weighted_degrees(selected, weights, strides=STRIDES):
    """Weighted marginal edge count, with the final singleton correction."""
    blocks = np.asarray(selected, dtype=bool).reshape((-1, 4))
    result = np.zeros(blocks.shape)
    total = sum(weights)
    result[:, 1:] += total*blocks[:, :-1]
    result[:, :-1] += total*blocks[:, 1:]
    for weight, stride in zip(weights, strides):
        result[:, 0] += weight*np.roll(blocks[:, 3], stride)
        result[:, 3] += weight*np.roll(blocks[:, 0], -stride)
    tail = -strides[-1] % len(blocks)
    result[0, 0] -= weights[-1]*blocks[tail, 3]
    result[tail, 3] += weights[-1]*(1-int(blocks[0, 0]))
    return result.ravel()


@lru_cache(maxsize=None)
def pressure_derivative(safe, numerator, denominator):
    low, high = max(0, safe-1), min(ROWS, safe+1)
    return ((stationary_pressure(high, numerator, denominator)
             - stationary_pressure(low, numerator, denominator))/(high-low))


def pressure(counts, mark):
    return sum(stationary_pressure(int(a), *mark) for a in counts)/len(counts)


def dual_gap(selected, counts, mark):
    eta = mark[0]/mark[1]
    alpha = np.count_nonzero(selected)/len(selected)
    return pressure(counts, mark)+eta*.4-(1+eta)*alpha/2


def ascend(selected, mark=(1, 12), unit=4, max_moves=128, shortlist=8):
    """Exact-count, fixed-cardinality ascent using a pressure tangent.

    Unit4 moves one complete group to an empty group. Unit1 allows arbitrary
    heterogeneous patterns. Only top/bottom `shortlist` tangent candidates
    are tested; every accepted move improves the actual floating objective.
    """
    selected = np.array(selected, dtype=bool, copy=True)
    counts = np.array(grouped_safe_counts(selected, 4), dtype=np.int64)
    score = pressure(counts, mark)
    trace = [score]
    for step in range(max_moves):
        weights = np.array([pressure_derivative(int(a), *mark) for a in counts])
        degree = weighted_degrees(selected, weights)
        if unit == 4:
            # Common within-group terms cancel in a full-to-empty swap.
            eligible_remove = np.all(selected.reshape((-1, 4)), axis=1)
            eligible_add = ~np.any(selected.reshape((-1, 4)), axis=1)
            degree = degree.reshape((-1, 4)).sum(axis=1)
            degree[eligible_remove] -= 3*sum(weights)
            degree[eligible_add] += 3*sum(weights)
        elif unit == 1:
            eligible_remove, eligible_add = selected, ~selected
        else:
            raise ValueError("unit must be 1 or 4")
        removes = np.flatnonzero(eligible_remove)
        adds = np.flatnonzero(eligible_add)
        removes = removes[np.argsort(degree[removes], kind="stable")[:shortlist]]
        adds = adds[np.argsort(-degree[adds], kind="stable")[:shortlist]]
        best = None
        for remove in removes:
            for add in adds:
                rows = (*range(unit*int(remove), unit*(int(remove)+1)),
                        *range(unit*int(add), unit*(int(add)+1)))
                delta = change_delta(selected, rows)
                gain = float(np.dot(weights, delta))
                if best is None or gain > best[0]:
                    best = gain, rows, delta
        if best is None or best[0] <= 1e-15:
            break
        _, rows, delta = best
        next_score = pressure(counts+delta, mark)
        if next_score <= score+1e-14:
            break
        selected[list(rows)] = ~selected[list(rows)]
        counts += delta
        score = next_score
        trace.append(score)
    if tuple(counts) != grouped_safe_counts(selected, 4):
        raise AssertionError("incremental safe-count drift")
    return selected, counts, trace


def seed(mode=1599, groups=819):
    size = ROWS//4
    order = np.argsort(-np.cos(2*np.pi*mode*np.arange(size)/size), kind="stable")
    selected = np.zeros(size, dtype=bool)
    selected[order[:groups]] = True
    return np.repeat(selected, 4)


def lifted_seed(mode, rows=3276, phase_only=True):
    """Threshold the leading four-lane Fourier eigenvector, with fixed size.

    The normalized symbol is one half of the Hermitian path matrix with
    cross-edge entry mean_c exp(2 pi i mode*stride_c/groups). Its phases
    allow first/last lane sets to differ, unlike whole-group seeds.
    """
    groups = ROWS//4
    z = np.exp(2j*np.pi*mode*np.asarray(STRIDES)/groups).mean()
    symbol = np.zeros((4, 4), dtype=complex)
    for lane in range(3):
        symbol[lane, lane+1] = symbol[lane+1, lane] = 1
    symbol[3, 0], symbol[0, 3] = z, z.conjugate()
    _, vectors = np.linalg.eigh(symbol)
    vector = vectors[:, -1]
    vector *= np.exp(-1j*np.angle(vector[0]))
    if phase_only:
        vector /= abs(vector)
    scores = (np.exp(2j*np.pi*mode*np.arange(groups)[:, None]/groups)
              * vector[None, :]).real
    selected = np.zeros(ROWS, dtype=bool)
    selected[np.argsort(-scores.ravel(), kind="stable")[:rows]] = True
    return selected


def periodic_profiles():
    """All 2^16 four-group patterns, collapsed only by exact safe profile."""
    words = np.arange(1 << 16, dtype=np.uint32)
    bits = ((words[:, None] >> np.arange(16)) & 1).reshape((-1, 4, 4))
    repeat = ROWS//16
    cardinality = bits.sum(axis=(1, 2))*repeat
    within = (bits[:, :, :-1] & bits[:, :, 1:]).sum(axis=(1, 2))
    first, last = bits[:, :, 0], bits[:, :, 3]
    a1 = (within+(last & np.roll(first, -1, axis=1)).sum(axis=1))*repeat
    a3 = (within+(last & np.roll(first, -3, axis=1)).sum(axis=1))*repeat
    terminal = last[:, -STRIDES[-1] % 4] & (1-first[:, 0])
    profiles, multiplicity = np.unique(
        np.column_stack((cardinality, a1, a3, terminal)), axis=0, return_counts=True)
    records = []
    for (a, one, three, correction), count in zip(profiles, multiplicity):
        counts = [int(one if stride % 4 == 1 else three) for stride in STRIDES]
        counts[-1] += int(correction)
        record = screen_profile(int(a), counts, extra_marks=((1, 12), (1, 24), (1, 48)))
        record.pop("counts")
        records.append(dict(pattern_count=int(count), **record))
    return records


def parity_endpoint_seed(rows):
    """Opposite endpoint fibers exploit that every frozen stride is odd.

    Between one and three selected rows per group on average, fill inward
    from the selected endpoint. All cross-group safe edges are then present.
    """
    groups = ROWS//4
    if not groups <= rows <= 3*groups:
        raise ValueError("row count must be between R/4 and 3R/4")
    blocks = np.zeros((groups, 4), dtype=bool)
    blocks[::2, 0], blocks[1::2, 3] = True, True
    remaining = rows-groups
    for depth in (1, 2):
        order = np.column_stack((4*np.arange(0, groups, 2)+depth,
                                 4*np.arange(1, groups, 2)+3-depth)).ravel()
        blocks.ravel()[order[:remaining]] = True
        remaining = max(0, remaining-len(order))
    return blocks.ravel()


def finite_score(bits, numerator, denominator):
    """Integer maximum of denominator*G-numerator*C over interval unions."""
    opening = 2*(numerator+denominator)
    outside, difference = 0, -opening
    for bit in bits:
        outside += max(0, difference)
        difference = (denominator if bit else -numerator)+min(0, max(-opening, difference))
    return outside+max(0, difference)


def screen(selected, counts, denominators=range(8, 19)):
    gaps = [(dual_gap(selected, counts, (1, denominator)), denominator)
            for denominator in denominators]
    gap, denominator = min(gaps)
    return dict(rows=int(np.count_nonzero(selected)),
                safe_min=int(min(counts)), safe_max=int(max(counts)),
                safe_total=int(sum(counts)), gap=gap, mark=f"1/{denominator}",
                partial_groups=int(np.count_nonzero(
                    (selected.reshape((-1, 4)).sum(axis=1) > 0)
                    & (selected.reshape((-1, 4)).sum(axis=1) < 4))),
                stationary_only=True, finite_failure_certificate=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--moves", type=int, default=128)
    parser.add_argument("--mode", type=int, default=1599)
    parser.add_argument("--groups", type=int, default=819)
    parser.add_argument("--denominator", type=int, default=12)
    parser.add_argument("--followup", action="store_true")
    args = parser.parse_args()
    if args.followup:
        profiles = periodic_profiles()
        screened = [r for r in profiles if "dual_gap" in r]
        print(json.dumps(dict(stage="all-period-four-groups", patterns=1 << 16,
                              profiles=len(profiles),
                              positive=sum(r["dual_gap"] > 0 for r in screened),
                              largest=sorted(screened, key=lambda r:r["dual_gap"],
                                             reverse=True)[:3])), flush=True)
        for numerator in (3, 4, 5, 6, 7):
            selected = parity_endpoint_seed(numerator*ROWS//10)
            record = screen_profile(int(selected.sum()), grouped_safe_counts(selected, 4),
                                    extra_marks=((1, 12),))
            record.pop("counts")
            print(json.dumps(dict(stage="opposite-parity-endpoints", **record)), flush=True)
        for mode in (87, 619, 1599):
            selected = lifted_seed(mode)
            counts = grouped_safe_counts(selected, 4)
            print(json.dumps(dict(stage="heterogeneous-seed", mode=mode,
                                  **screen(selected, counts))), flush=True)
            for unit in (1, 4, 1):
                selected, counts, trace = ascend(selected, (1, 12), unit, args.moves)
                print(json.dumps(dict(stage=f"heterogeneous-unit-{unit}", mode=mode,
                                      moves=len(trace)-1,
                                      **screen(selected, counts))), flush=True)
        return
    selected = seed(args.mode, args.groups)
    counts = grouped_safe_counts(selected, 4)
    print(json.dumps(dict(stage="seed", **screen(selected, counts))), flush=True)
    for unit in (4, 1):
        selected, counts, trace = ascend(selected, (1, args.denominator), unit, args.moves)
        print(json.dumps(dict(stage=f"unit-{unit}", moves=len(trace)-1,
                              objective_start=trace[0], objective_end=trace[-1],
                              **screen(selected, counts))), flush=True)


if __name__ == "__main__":
    main()
