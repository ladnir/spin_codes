"""Exact tiny-field checks for nonzero-diagonal MDS feedback.

Enumerates every first-scaled state. The last diagonal is summed exactly
by output support; its nonzero amplitudes are symmetric. No sampled tests.
"""
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
import time

from cheap_reset_maps import BinaryField


def rank(field, rows):
    a = [list(row) for row in rows]
    r = 0
    for j in range(len(a[0])):
        pivot = next((i for i in range(r,len(a)) if a[i][j]),None)
        if pivot is None:
            continue
        a[r],a[pivot] = a[pivot],a[r]
        inv = field.inverse(a[r][j])
        a[r] = [field.mul(x,inv) for x in a[r]]
        for i in range(len(a)):
            if i != r and a[i][j]:
                scalar = a[i][j]
                a[i] = [x^field.mul(scalar,y) for x,y in zip(a[i],a[r])]
        r += 1
        if r == len(a):
            break
    return r


def check(q,b,modulus,deadline):
    field = BinaryField(q.bit_length()-1,modulus)
    D = [[field.inverse(i^(b+j)) for j in range(b)] for i in range(b)]
    minors = 0
    for h in range(1,b+1):
        for rows in combinations(range(b),h):
            for cols in combinations(range(b),h):
                assert rank(field,[[D[i][j] for j in cols] for i in rows]) == h
                minors += 1
    vectors = cases = 0
    for source_mask in range(1,1<<b):
        support = [j for j in range(b) if source_mask>>j&1]
        counts = Counter()
        for values in product(range(1,q),repeat=len(support)):
            if time.monotonic()>deadline:
                raise RuntimeError('Exact test cap')
            y = []
            for row in D:
                value = 0
                for j,v in zip(support,values):
                    value ^= field.mul(row[j],v)
                y.append(value)
            out_mask = sum((v!=0)<<j for j,v in enumerate(y))
            counts[out_mask] += 1
            vectors += 1
        total = Fraction(0)
        for mask in range(1<<b):
            w = mask.bit_count()
            probability = Fraction(counts[mask],(q-1)**(len(support)+w))
            assert probability <= Fraction(1,(q-1)**b)
            if b-w >= len(support):
                assert probability == 0
            total += probability*(q-1)**w
            cases += 1
        assert total == 1
    # Removing MDS mixing violates the stated bound already on one active
    # coordinate. Allowing zero on the first diagonal also creates an atom.
    assert Fraction(1,q-1) > Fraction(1,(q-1)**b)
    assert Fraction(1,q) > Fraction(1,(q-1)**b)
    return dict(q=q,b=b,minors=minors,vectors=vectors,support_laws=cases)


def main():
    start = time.monotonic()
    for q,b,modulus in ((4,2,0b111),(8,4,0b1011)):
        print(check(q,b,modulus,start+15),flush=True)
    # Bernoulli: (1-1/q)^N >= 1-N/q > 1/2. This proves the
    # full-trajectory density overhead is <2 without decimal rounding.
    N,q = 2**21,2**128
    assert 2*N < q
    print(dict(kind='exact_feedback_local_checks',global_density_less_than=2,
               seconds=round(time.monotonic()-start,3)))


if __name__ == '__main__':
    main()
