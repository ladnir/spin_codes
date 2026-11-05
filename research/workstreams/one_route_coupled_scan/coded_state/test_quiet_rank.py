"""Exact small-field rank/moment tests; no output files or K20 certificate."""
import itertools
import random
import time
from fractions import Fraction
from check_maps import Field
from quiet_rank import structural_rank, nullity_moment, binomial_tail, optimistic_mixed_screen


def product(field, left, scales, right):
    return [[field.xor(field.table[field.table[left[i][k]][scales[k]]][right[k][j]]
                       for k in range(len(scales)))
             for j in range(len(right[0]))] for i in range(len(left))]


def main():
    started = time.monotonic()
    field = Field(2,0x7)
    rng = random.Random(94715)
    def cap():
        if time.monotonic()-started > 50:
            raise RuntimeError('50-second test cap')
    comparisons = 0
    # Every 2x2 constant perturbation, plus independently sampled 3x3 ones.
    for size, constants in ((2,itertools.product(range(4),repeat=4)),
                            (3,([rng.randrange(4) for _ in range(9)] for _ in range(64)))):
        for entries in constants:
            cap()
            base = [list(entries[i*size:(i+1)*size]) for i in range(size)]
            deficits = []
            for diagonal in itertools.product(range(1,4),repeat=size):
                matrix = [row[:] for row in base]
                for i, value in enumerate(diagonal):
                    matrix[i][i] ^= value
                deficits.append(size-field.rank(matrix))
            for threshold in range(1,size+1):
                actual = Fraction(sum(d>=threshold for d in deficits),len(deficits))
                assert actual <= binomial_tail(size,threshold,3)
                comparisons += 1
            for tau in (Fraction(3,2),Fraction(4),Fraction(16)):
                actual = sum((tau**d for d in deficits),Fraction(0))/len(deficits)
                assert actual <= nullity_moment(size,size,tau,4)
    print(f'Random diagonal plus constant: {comparisons} exact tail comparisons PASS',flush=True)

    products = 0
    for _ in range(64):
        cap()
        left = [[rng.randrange(4) for _ in range(4)] for _ in range(2)]
        right = [[rng.randrange(4) for _ in range(3)] for _ in range(4)]
        common = 0
        for size in range(1,3):
            for selected in itertools.combinations(range(4),size):
                ls = [[row[k] for k in selected] for row in left]
                rs = [right[k] for k in selected]
                if field.rank(ls)==size and field.rank(rs)==size:
                    common = size
        moments = [Fraction(0)]*3
        markers = (Fraction(3,2),Fraction(4),Fraction(16))
        for scales in itertools.product(range(1,4),repeat=4):
            deficit = 3-field.rank(product(field,left,scales,right))
            for i,tau in enumerate(markers):
                moments[i] += tau**deficit/81
            products += 1
        for actual,tau in zip(moments,markers):
            assert actual <= nullity_moment(3,common,tau,4)
    print(f'B diag(R) C: {products} exact matrix ranks; all moments PASS',flush=True)

    # Two RS[4,2] blocks; four regions with two labels in either order.
    # Each selected quiet batch is one region and has compression rank one.
    # Shortening leaves max(0,t_i-2) dimensions on source block i.
    flow_cases = 0
    for swaps in itertools.product(range(2),repeat=4):
        labels = [i^swaps[j] for j in range(4) for i in range(2)]
        # Coordinate-wise allowed U makes irregular outer occupancies too.
        for mask in range(256):
            allowed = [k for k in range(8) if mask>>k&1]
            counts = [sum(labels[k]==i for k in allowed) for i in range(2)]
            dimensions = [max(0,t-2) for t in counts]
            edges = [[sum(labels[k]==i and k//2==j for k in allowed)
                      for j in range(4)] for i in range(2)]
            capacities = [1,0,1,1]
            flow = structural_rank(edges,dimensions,capacities)
            best = 0
            for submask in range(1<<len(allowed)):
                selected = [allowed[z] for z in range(len(allowed)) if submask>>z&1]
                if all(sum(labels[k]==i for k in selected)<=dimensions[i] for i in range(2)) and all(
                        sum(k//2==j for k in selected)<=capacities[j] for j in range(4)):
                    best = max(best,len(selected))
            assert flow==best
            flow_cases +=1
        cap()
    print(f'MDS-capacity max flow versus exhaustive coordinate sets: {flow_cases} PASS',flush=True)
    optimistic_mixed_screen()
    print(f'Elapsed {time.monotonic()-started:.3f}s; no whole-code certificate claimed.',flush=True)


if __name__ == '__main__':
    main()
