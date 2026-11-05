"""Exact local and small counting checks. Not a whole-code certificate."""
from fractions import Fraction
import itertools
import math
import random
import time
from check_maps import Field
from cheap_reset_maps import BinaryField, novel_values
from group32_joint import weighted_lengths, variance_cap, old
from group_reset_cover import linear_max


def main():
    start = time.monotonic()
    # Exact coefficients of the joint mask/length generating function.
    checks = 0
    for G, d, r, extra in ((2, 1, 2, 3), (3, 2, 2, 2), (2, 2, 3, 2)):
        tau = Fraction(5, 4)
        exact = {}
        for masks in itertools.product(range(1, 1 << G), repeat=r):
            sizes = [x.bit_count() for x in masks]
            active = sum(sizes)
            for lengths in itertools.product(range(extra+1), repeat=r):
                if sum(lengths) > extra:
                    continue  # remainder belongs to the final free interval
                H = sum(d*(G-a)*e for a, e in zip(sizes, lengths))
                exact[active] = exact.get(active, Fraction(0))+tau**(-H)
        for active, value in exact.items():
            bound = weighted_lengths(G, d, r, active, extra, math.log(float(tau)))
            assert math.log(float(value)) <= bound+1e-10
            checks += 1

    # Actual region probabilities, not sampled routes, check the variance cap.
    rng = random.Random(20261022)
    polytope_checks = 0
    for _ in range(200):
        G = rng.randrange(2,6)
        ranges = [(rng.randrange(1,8),rng.randrange(8,16)),
                  (rng.randrange(0,3),rng.randrange(3,7)),
                  (rng.randrange(0,4),rng.randrange(4,12))]
        coeffs = [rng.randrange(-10,11) for _ in range(3)]
        vertices = [(v,r,a) for v in range(ranges[0][0],ranges[0][1]+1)
                    for r in range(ranges[1][0],ranges[1][1]+1)
                    for a in range(ranges[2][0],ranges[2][1]+1)
                    if r<=a<=G*r and v>=2*r+1]
        exact = max((sum(x*y for x,y in zip(coeffs, vertex)) for vertex in vertices),
                    default=-math.inf)
        assert linear_max(coeffs,ranges,G) == exact
        polytope_checks += 1
    geometries = 0
    for B, b in ((256, 512), (512, 256)):
        n, G, d, m = 2**21//b, b//32, 32, 16
        for _ in range(120):
            c = rng.randrange(2, 12)
            cuts = sorted(rng.sample(range(n+1), 2*c))
            intervals = list(zip(cuts[::2], cuts[1::2]))
            if any(end-begin < 2 for begin,end in intervals[:-1]):
                continue
            sizes = [rng.randrange(1,G+1) for _ in range(c-1)]
            H = sum((b-d*a)*(end-begin-2)
                    for (begin,end),a in zip(intervals[:-1],sizes))
            for forward in (False, True):
                masses = [Fraction(0) for _ in range(B)]
                for i,(begin,end) in enumerate(intervals):
                    for t in range(begin,end):
                        weight = Fraction(1)
                        if t == begin:
                            weight = (Fraction(sizes[i],G) if i<c-1 else Fraction(1)) if forward else Fraction(0)
                        masses[t//m] += weight/m
                variance = sum(p*(1-p) for p in masses)
                cap = variance_cap(B,b,d,c,sum(sizes),H,forward)
                assert float(variance) <= cap+1e-10
                geometries += 1

    # Entropy tangents must remain above the full-space branch below W too.
    tangents = 0
    for E, W in ((100, 70), (100, 99), (10000, 9000)):
        for frac in (1, .5, .25, .125):
            anchor = W+(E-W)*frac
            slope = math.log(anchor/(anchor-W))/(128*old.LN2)
            constant = W+old.entropy_logchoose(anchor,W)/(128*old.LN2)-slope*anchor
            for remaining in range(0,E+1,max(1,E//100)):
                reference = float(remaining)
                if remaining > W:
                    reference = min(reference,W+old.entropy_logchoose(remaining,W)/(128*old.LN2))
                assert reference <= constant+slope*remaining+1e-9
                tangents += 1

    # Independently merging two reset groups restores rank on the old
    # weight-(d+1) shortcut; all four-column subsets satisfy the toy MDS law.
    f, poly = Field(4, 0x13), BinaryField(4, 0x13)
    columns = [novel_values(poly, 2, t) for t in range(16)]
    ranks = 0
    for support in itertools.combinations(range(16), 4):
        assert f.rank([[columns[j][i] for j in support] for i in range(4)]) == 4
        ranks += 1
    assert f.rank([[columns[j][i] for j in (0,4,8)] for i in range(4)]) == 3

    # For a fixed nonzero state, all output vectors are equally likely under
    # uniform feedback. This checks both interior holes and the terminal
    # reset charge without double-counting terminal inactive coordinates.
    f = Field(2, 7)
    laws = 0
    matrices = list(itertools.product(range(4), repeat=4))
    for q in itertools.product(range(4), repeat=2):
        if q == (0,0):
            continue
        counts = {}
        for m in matrices:
            y = tuple(f.apply([list(m[:2]),list(m[2:])],q))
            counts[y] = counts.get(y,0)+1
        assert set(counts.values()) == {16}
        assert sum(count for y,count in counts.items() if y[0] == 0) == 64
        assert counts[(0,0)] == 16
        laws += 1
    # Zero coset maximizes the weight enumerator for z in [0,1]. Check
    # every syndrome of a small [4,2] MDS check matrix using exact fractions.
    cosets = {}
    for word in itertools.product(range(4),repeat=4):
        syndrome = tuple(f.apply([[1]*4,list(range(4))],word))
        polynomial = cosets.setdefault(syndrome,[0]*5)
        polynomial[sum(x!=0 for x in word)] += 1
    coset_checks = 0
    for z in (Fraction(1,4),Fraction(1,2),Fraction(3,4),Fraction(1)):
        zero = sum(count*z**i for i,count in enumerate(cosets[(0,0)]))
        for polynomial in cosets.values():
            assert sum(count*z**i for i,count in enumerate(polynomial)) <= zero
            coset_checks += 1
        upper = 1+sum(math.comb(4,w)*4**(w-2)*z**w for w in range(3,5))
        assert zero <= upper
    assert time.monotonic()-start < 30
    print(dict(kind='LOCAL_CHECKS_NOT_CERTIFICATE', weighted_coefficients=checks,
               geometry_checks=geometries, tangent_checks=tangents,
               paired_reset_ranks=ranks, uniform_feedback_laws=laws,
               exact_polytope_checks=polytope_checks,
               exact_coset_checks=coset_checks,
               seconds=round(time.monotonic()-start,3)),flush=True)


if __name__ == '__main__':
    main()
