"""Independent small checks of the grouped transfer/conditioning machinery."""

from itertools import combinations
from math import comb,log
import numpy as np
from scipy.special import logsumexp
from flint import arb,ctx

import group_moment as model
import group_rank_one_verify as exact


def conditioning():
    rng = np.random.default_rng(42)
    for n in (2,3,7):
        z = rng.random((n,n))/n
        z[0] = 0
        z[0,0] = 1
        a = rng.random((n,n))/n
        for length in range(1,9):
            direct = model.log_moments(z,a,length)
            conditional = model.first_support_moments(z,a,length)
            for w in range(1,length+1):
                mixed = logsumexp([log(comb(l,w-1))-log(comb(length,w))+conditional[l,w]
                                   for l in range(w-1,length)])
                assert abs(direct[w]-mixed) < 1e-12
            # Independently sum every occupied-region placement.
            expected = np.zeros(length+1)
            for w in range(length+1):
                for places in combinations(range(length),w):
                    current = np.zeros(n)
                    current[0] = 1
                    for j in range(length):
                        current = current @ (a if j in places else z)
                    expected[w] += current.sum()/comb(length,w)
            assert np.allclose(np.exp(direct),expected,rtol=1e-12,atol=1e-14)
    print('Conditional/unconditional identities and exhaustive placement checks passed')


def state_law():
    for s in range(2,6):
        m = (1 << s)-1
        denominator = m*(1 << (s-1))
        for q in range(1,1 << s):
            counts = [0]*(1 << s)
            for u in range(1,1 << s):
                for v in range(1 << s):
                    if (u&v).bit_count()%2 == 0:
                        counts[q ^ (u if (q&v).bit_count()%2 else 0)] += 1
            assert counts[0] == 0
            for r in range(1,1 << s):
                assert 2*m*counts[r] == denominator*(1+m*(q==r))
    print('Exact lazy/uniform transvection law passed exhaustive s=2..5 checks')


def transfers():
    data = model.census()
    ctx.prec = 192
    for weight in (1,2,8,16):
        for tilt in ('0','0.0004','0.0025'):
            approximate = model.transfers(data,float(tilt),weight)
            verified = exact.transfers(data,tilt,weight)
            for array,matrix in zip(approximate,verified):
                expected = np.array([[float(matrix[i,j]) for j in range(matrix.ncols())]
                                     for i in range(matrix.nrows())])
                assert np.allclose(array,expected,rtol=1e-12,atol=1e-15)
            # Empty input never creates or kills a live state. At zero tilt,
            # its transfer conserves mass, unlike the active upper envelope.
            if tilt == '0':
                empty = verified[0]
                for i in range(empty.nrows()):
                    assert exact.up(sum((empty[i,j] for j in range(empty.ncols())),arb(0))) >= 1
                    assert float(sum((empty[i,j] for j in range(empty.ncols())),arb(0))) < 1+1e-14
    print('Binary64/outward transfers agree; zero-input invariants passed')


if __name__ == '__main__':
    state_law()
    conditioning()
    transfers()
