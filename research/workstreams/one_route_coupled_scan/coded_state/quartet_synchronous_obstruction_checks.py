"""Exact small checks for QUARTET_SYNCHRONOUS_OBSTRUCTION.md.

No production matrix, field arithmetic, optimization, or floating-point
decision is used.  The proof's negative-association theorem and analytic
logarithm bounds remain mathematical inputs; toy laws and all finite
constant comparisons are checked independently here.
"""

from fractions import Fraction
from itertools import combinations, product
from math import comb
import os
import shutil
import threading
import time


SYNC_WIDTH, SYNC_T, SYNC_M, SYNC_W = 512, 4096, 8192, 943718


def _sync_constants():
    a = Fraction(3, 2)
    upper_log_a = sum((Fraction(1 if j % 2 else -1, j*2**j) for j in range(1, 8)), Fraction(0))
    checks = (
        upper_log_a < Fraction(13, 32),
        Fraction(13, 32)**2/(8*Fraction(2, 3)) == Fraction(507, 16384),
        Fraction(507, 16384) < Fraction(1, 32),
        9**8 > 2*8**8,
        Fraction(9, 8)*Fraction(2, 3)**12 == Fraction(512, 59049),
        Fraction(512, 59049) < Fraction(1, 112),
        SYNC_M*Fraction(3, 224) == Fraction(768, 7),
        Fraction(768, 7) < 128,
        a**12 > 2**7,
        256*Fraction(2, 3)**12 < 2,
        3**5 < 2**8,
        57**2 > 2**7,
        75**8 > 2**37,
        Fraction(14)-Fraction(434, 3) < -130,
        Fraction(14)-Fraction(266, 3) < -74,
        Fraction(SYNC_T)+SYNC_M*Fraction(8, 5)-(512*58+1)*Fraction(7, 12) == Fraction(-2401, 20),
        Fraction(-2401, 20) < -120,
        Fraction(SYNC_T)+SYNC_M*Fraction(8, 5)-(416*76+1)*Fraction(7, 12) == Fraction(-24801, 20),
        Fraction(-24801, 20) < -1200,
        Fraction(1, 2**130)+Fraction(1, 2**120) < Fraction(1, 2**119),
        Fraction(1, 2**74)+Fraction(1, 2**1200) < Fraction(1, 2**73),
        SYNC_W < 116*SYNC_M,
        SYNC_W < 128*SYNC_M,
        512-6*16 == 416,
        128+26-Fraction(512*7, 12) == -Fraction(434, 3),
        128+26-Fraction(416*7, 12) == -Fraction(266, 3),
    )
    assert all(checks)
    # This is the concavity consequence used for both small-r ranges.
    for r in range(1, 76):
        assert r**16 >= 2**(r-1)
    return len(checks), 75


def _sync_histories():
    histories = partial_regions = pascal = suffixes = 0
    for updates in (4, 6, 8, 12):
        histogram = {}
        for bits in product((0, 1), repeat=updates):
            r = sum(bits[t] and not bits[t+1] for t in range(updates-1))
            footprint = tuple(int(bits[t] or (t and bits[t-1])) for t in range(updates))
            assert sum(footprint) == sum(bits)+r
            histogram[r] = histogram.get(r, 0)+1
            histories += 1
            for region_length in (2, 3, 4):
                if updates % region_length:
                    continue
                counts = [sum(footprint[start:start+region_length])
                          for start in range(0, updates, region_length)]
                partial = sum(0 < value < region_length for value in counts)
                assert partial <= min(len(counts), 2*r+2)
                partial_regions += 1
                if r == 0:
                    assert footprint == bits
                    assert partial <= 1
                    # With one partial region, the rank is independent of
                    # which rows occupy that region's covered slots.
                    width, rows = 4, 4*region_length
                    cap = len(counts)//2
                    whole = sum(value == region_length for value in counts)
                    fringe = sum(value*width for value in counts if 0 < value < region_length)
                    rank = rows*min(cap, whole)+(fringe if whole < cap else 0)
                    volume = sum(footprint)*width
                    assert rank == min(volume, cap*rows)
                    suffixes += 1
        for r in range(updates//2+1):
            assert histogram.get(r, 0) == comb(updates+1, 2*r+1)
            pascal += 1
        assert sum(histogram.values()) == 2**updates
    return histories, partial_regions, pascal, suffixes


def _sync_stats(support, partitions, row_maps, cap):
    width = len(row_maps[0])
    degrees = [0]*width
    previous = volume = deaths = footprint_volume = 0
    footprints = []
    for state, partition, row_map in zip(support, partitions, row_maps):
        live = 0
        for block in partition:
            if block & previous:
                live |= block
        footprint = live | state
        footprints.append(footprint)
        volume += state.bit_count()
        deaths += (previous & ~state).bit_count()
        footprint_volume += footprint.bit_count()
        remaining = footprint
        while remaining:
            bit = remaining & -remaining
            degrees[row_map[bit.bit_length()-1]] += 1
            remaining ^= bit
        previous = state
    rank = sum(min(cap, degree) for degree in degrees)
    return volume, deaths, footprint_volume, rank, tuple(footprints)


def _sync_lipschitz():
    flips = threshold_comparisons = 0
    fixtures = (
        (4, ((15,),)*3, (0, 1, 2), tuple(range(16)), 2),
        (8, ((15, 240), (23, 232), (113, 142)), (0, 1, 3), (0, 1, 3, 15, 23, 113, 255), 2),
    )
    for width, partitions, shifts, alphabet, cap in fixtures:
        row_maps = tuple(tuple((coordinate+shift) % width for coordinate in range(width)) for shift in shifts)
        for support in product(alphabet, repeat=len(partitions)):
            L, D, F, R, V = _sync_stats(support, partitions, row_maps, cap)
            for t in range(len(partitions)):
                for coordinate in range(width):
                    changed = list(support)
                    changed[t] ^= 1 << coordinate
                    new_L, new_D, new_F, new_R, new_V = _sync_stats(changed, partitions, row_maps, cap)
                    assert abs(new_L-L) <= 1
                    assert abs(new_D-D) <= 1
                    assert abs(new_R-R) <= 5
                    assert abs(new_F-F) <= 5
                    assert sum((left ^ right).bit_count() for left, right in zip(V, new_V)) <= 5
                    flips += 1
                    if new_L > new_R:
                        cutoff = new_R+new_D+1
                        assert R <= L-1+6
                        assert R <= cutoff-D-1+6
                        threshold_comparisons += 1
    return flips, threshold_comparisons


def _sync_subsets(rows, size):
    return tuple(sum(1 << i for i in choice) for choice in combinations(range(rows), size))


def _sync_na_controls():
    """Enumerate complete small dependent row laws, not independent-row surrogates."""
    a = Fraction(3, 2)
    profiles = ((0, 1, 2, 4), (1, 1, 1, 1), (2, 2, 2, 2),
                (3, 3, 3, 3), (1, 2, 3, 4), (0, 0, 4, 4))
    controls = outcomes = centered = 0
    rows, cap = 4, 2
    for volumes in profiles:
        choices = tuple(_sync_subsets(rows, size) for size in volumes)
        joint = Fraction(0)
        total = 0
        for covered in product(*choices):
            degrees = [sum(bool(region & (1 << row)) for region in covered) for row in range(rows)]
            joint += a**sum(max(0, degree-cap) for degree in degrees)
            total += 1
        joint /= total
        distribution = [Fraction(1)]
        for size in volumes:
            p = Fraction(size, rows)
            following = [Fraction(0)]*(len(distribution)+1)
            for degree, probability in enumerate(distribution):
                following[degree] += probability*(1-p)
                following[degree+1] += probability*p
            distribution = following
        row_moment = sum(probability*a**max(0, degree-cap)
                         for degree, probability in enumerate(distribution))
        assert joint <= row_moment**rows
        # Raise to32*rows to keep the centered Hoeffding comparison rational.
        ordinary = sum(probability*a**degree for degree, probability in enumerate(distribution))
        partial = sum(0 < size < rows for size in volumes)
        assert ordinary**(32*rows) <= a**(32*sum(volumes))*2**(partial*rows)
        controls += 1
        outcomes += total
        centered += 1
    return controls, outcomes, centered


def synchronous_obstruction_checks():
    constants, small_ranges = _sync_constants()
    histories, partial, pascal, suffixes = _sync_histories()
    flips, thresholds = _sync_lipschitz()
    na_controls, na_outcomes, centered = _sync_na_controls()
    print(dict(exact_constant_comparisons=constants, exact_small_range_comparisons=small_ranges,
               synchronous_histories=histories, partial_region_comparisons=partial,
               exact_history_count_identities=pascal, deterministic_suffix_rank_checks=suffixes,
               exact_one_bit_lipschitz_checks=flips, obstruction_threshold_comparisons=thresholds,
               negative_association_toy_controls=na_controls, dependent_route_outcomes=na_outcomes,
               centered_moment_toy_controls=centered,
               synchronous_obstruction_bound_bits=-119, radius_16r_obstruction_bound_bits=-73,
               no_perturbation_union=True, full_code_certificate=False), flush=True)


def synchronous_obstruction_main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--checks', action='store_true', help='Run the small exact checks (default action)')
    parser.parse_args()
    import resource
    if shutil.disk_usage('/tmp').free < 20*(1 << 30):
        raise RuntimeError('20GiB /tmp free-space guard')
    resource.setrlimit(resource.RLIMIT_AS, (2*(1 << 30), 2*(1 << 30)))
    timer = threading.Timer(175, lambda: os._exit(124))
    timer.daemon = True
    timer.start()
    started = time.monotonic()
    try:
        synchronous_obstruction_checks()
        print(dict(seconds=round(time.monotonic()-started, 3)), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    synchronous_obstruction_main()
