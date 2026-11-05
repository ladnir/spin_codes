"""Support counts for two alternating pairings of a circular group array.

This is a combinatorial diagnostic, not a distance certificate. All counts
are exact integers; logarithms printed by the CLI are floating diagnostics.
"""
from collections import defaultdict
from math import comb, log2
import argparse


def exact_profiles(pairs):
    """Count indexed bit strings by (weight, even-edge11, odd-edge11).

    A source array has 2*pairs groups. Pairing 0 uses (0,1),(2,3),...;
    pairing 1 uses (1,2),(3,4),...,(2*pairs-1,0). Strings are not identified
    under rotation. Intended only for small exact validation instances.
    """
    if pairs < 1:
        raise ValueError("positive pair count required")
    result = defaultdict(int)
    for first in (0, 1):
        states = {(first, 0, 0, first): 1}
        for edge in range(2 * pairs - 1):
            following = defaultdict(int)
            for (r, a, b, last), multiplicity in states.items():
                for bit in (0, 1):
                    product = last * bit
                    key = (r + bit, a + (product if edge % 2 == 0 else 0),
                           b + (product if edge % 2 else 0), bit)
                    following[key] += multiplicity
            states = following
        for (r, a, b, last), multiplicity in states.items():
            result[r, a, b + last * first] += multiplicity
    return dict(result)


def one_pairing_count(pairs, weight, full_pairs):
    singleton_pairs = weight - 2 * full_pairs
    if (full_pairs < 0 or singleton_pairs < 0
            or full_pairs + singleton_pairs > pairs):
        return 0
    return (comb(pairs, full_pairs)
            * comb(pairs - full_pairs, singleton_pairs)
            * (1 << singleton_pairs))


def circular_run_count(length, weight, runs):
    if weight == 0:
        return int(runs == 0)
    if weight == length:
        return int(runs == 0)
    if not 1 <= runs <= min(weight, length - weight):
        return 0
    # Mark a run start, choose positive run/gap compositions, then remove
    # the runs possible markings. Rotation is not factored out.
    numerator = length * comb(weight - 1, runs - 1) * comb(length - weight - 1, runs - 1)
    assert numerator % runs == 0
    return numerator // runs


def profile_count_upper(pairs, weight, even_full, odd_full):
    length = 2 * pairs
    if not (0 <= weight <= length and 0 <= even_full <= pairs
            and 0 <= odd_full <= pairs):
        return 0
    if weight == 0:
        return int(even_full == odd_full == 0)
    if weight == length:
        return int(even_full == odd_full == pairs)
    # Every nonconstant cyclic bit string has runs = weight - total11.
    runs = weight - even_full - odd_full
    # With no singleton pairs in either matching, contract that matching.
    # This is an exact count of a weight-(weight/2) cyclic string, not just
    # the looser length-2*pairs run bound. It covers the original obstructing
    # class in which all active packets were full.
    if weight == 2 * even_full or weight == 2 * odd_full:
        return circular_run_count(pairs, weight // 2, runs)
    return min(one_pairing_count(pairs, weight, even_full),
               one_pairing_count(pairs, weight, odd_full),
               circular_run_count(length, weight, runs))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pairs', type=int, default=4096)
    parser.add_argument('--weight', type=int, default=1024)
    parser.add_argument('--runs', type=int, nargs='+', default=[1, 8, 32, 64, 128, 256, 448, 512])
    args = parser.parse_args()
    a = args.weight // 2
    print('weight,even_full,odd_full,runs,log2_count_upper,log2_even_count,support_bits_saved')
    old = one_pairing_count(args.pairs, args.weight, a)
    for runs in args.runs:
        b = args.weight - a - runs
        upper = profile_count_upper(args.pairs, args.weight, a, b)
        if upper:
            print(f'{args.weight},{a},{b},{runs},{log2(upper):.9f},'
                  f'{log2(old):.9f},{log2(old)-log2(upper):.9f}')


if __name__ == '__main__':
    main()
