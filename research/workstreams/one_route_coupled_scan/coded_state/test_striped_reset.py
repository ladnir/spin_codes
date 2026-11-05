"""Exact finite checks of reset rank and activation-support geometry, not a cert."""
import itertools
import random
import time
from check_maps import Field
from cheap_reset_maps import BinaryField, novel_values


def main():
    start = time.monotonic()
    field = BinaryField(4, 0x13)
    rank_field = Field(4, 0x13)
    # Two independent d=2 groups, over three batches. Each group has 6 columns.
    columns = []
    labels = []
    for t in range(3):
        for group in range(2):
            for j in range(2):
                column = [0]*4
                column[2*group:2*group+2] = novel_values(field, 1, 2*t+j)
                columns.append(column)
                labels.append(group)
    ranks = 0
    for mask in range(1 << 12):
        support = [j for j in range(12) if mask >> j & 1]
        weights = [sum(labels[j] == g for j in support) for g in range(2)]
        matrix = [[columns[j][i] for j in support] for i in range(4)]
        assert rank_field.rank(matrix) == sum(min(2, w) for w in weights)
        if all(w == 0 or w >= 3 for w in weights):
            assert rank_field.rank(matrix) == 2*sum(w > 0 for w in weights)
        ranks += 1
    # A three-symbol kernel in one group violates the full-state MDS interface.
    support = [0, 4, 8]
    assert rank_field.rank([[columns[j][i] for j in support] for i in range(4)]) == 2
    # Count every group pattern, independently of the binomial upper bound.
    from math import comb
    patterns = 0
    for groups in (2, 3):
        for closed in (1, 2, 3):
            counts = {}
            for choices in itertools.product(range(1, 1 << groups), repeat=closed):
                active = sum(x.bit_count() for x in choices)
                counts[active] = counts.get(active, 0)+1
            for active, count in counts.items():
                assert count <= comb(groups*closed, active)
                patterns += 1
    # Formation of S: unrestricted noninitial batches, selected first groups,
    # free final start. Holes remain inside excursion endpoint regions.
    rng = random.Random(20261022)
    geometry = 0
    b, d, n, per_region = 8, 2, 32, 4
    for _ in range(1000):
        c = rng.randrange(1, 5)
        cuts = sorted(rng.sample(range(n+1), 2*c))
        intervals = list(zip(cuts[::2], cuts[1::2]))
        if any(end-begin < (2 if i < c-1 else 1) for i, (begin, end) in enumerate(intervals)):
            continue
        chosen = [rng.randrange(1, 1 << (b//d)) for _ in range(c-1)]
        active = sum(x.bit_count() for x in chosen)
        v = sum(end-begin for begin, end in intervals)
        S = set()
        allowed_output = 0
        for i, (begin, end) in enumerate(intervals):
            groups = b//d if i == c-1 else chosen[i].bit_count()
            allowed_output += d*groups*(end-begin)
            for t in range(begin, end):
                for j in range(b):
                    if t != begin or i == c-1 or chosen[i] >> (j//d) & 1:
                        S.add(t*b+j)
        assert len(S) == b*(v-c+1)+d*active
        occupancy = [sum(x//(b*per_region) == region for x in S) for region in range(n//per_region)]
        assert sum(0 < x < b*per_region for x in occupancy) <= 2*c
        assert allowed_output <= b*(v-2*(c-1))+2*d*active
        geometry += 1
    assert time.monotonic()-start < 15
    print(dict(kind='exact_local_interfaces_NOT_certificate', ranks=ranks,
               pattern_counts=patterns, activation_geometries=geometry,
               seconds=round(time.monotonic()-start, 3)), flush=True)


if __name__ == '__main__':
    main()
