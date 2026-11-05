"""Exact toy route/rank checks, directed-arithmetic checks; no output files."""
from collections import defaultdict
from decimal import localcontext
from fractions import Fraction
import itertools
import math
import random
import time
from feedback_support import one_block_moment, binomial_moment, fresh_dimension
from feedback_support_certificate import Interval, rational, certify
from check_maps import Field


def main():
    started = time.monotonic()
    def cap():
        if time.monotonic()-started>45:
            raise RuntimeError('45-second test cap')
    groups = defaultdict(list)
    for columns in itertools.product(range(8),repeat=4):
        counts = tuple(x.bit_count() for x in columns)
        weights = [sum((x>>i)&1 for x in columns) for i in range(3)]
        d = sum(max(0,w-2) for w in weights)
        complement = sum(max(0,2-w) for w in weights)
        assert d==sum(counts)-6+complement
        groups[counts].append(d)
    comparisons = 0
    for counts,dimensions in groups.items():
        for tau in (2,3,4):
            actual = sum(Fraction(tau)**d for d in dimensions)/len(dimensions)
            direct = one_block_moment(counts,3,2,tau)**3
            complement = Fraction(tau)**(sum(counts)-6)*one_block_moment(tuple(3-a for a in counts),3,2,tau)**3
            homogeneous = binomial_moment(4,2,sum(counts),12,tau)**3
            assert actual<=direct<=homogeneous and actual<=complement
            comparisons +=1
        cap()
    print(f'4096 incidence matrices; {comparisons} exact route-moment comparisons PASS.',flush=True)

    # Fresh block constraints across successive persistent states. Two
    # alternating observations of each state share one root set, not two states.
    f = Field(4,0x13)
    s,b,steps = 2,4,3
    B = f.interpolation_matrix(f.points[:b],f.points[b:b+s])
    A = f.interpolation_matrix(f.points[b:b+s],f.points[:2*b])
    rng = random.Random(221022)
    for case in range(512):
        supports = [rng.sample(range(b),rng.randrange(b+1)) for _ in range(steps)]
        holes = [rng.sample(range(2*b),rng.randrange(s)) for _ in range(steps)]
        total = sum(map(len,supports))
        forms = [[0]*total for _ in range(s)]
        constraints,offset,rank_floor = [],0,0
        for support,H in zip(supports,holes):
            for j,col in enumerate(support):
                scale = rng.randrange(1,16)
                for r in range(s):
                    forms[r][offset+j] ^= f.table[B[r][col]][scale]
            for h in H:
                constraints.append([f.xor(f.table[A[h][r]][forms[r][j]] for r in range(s))
                                    for j in range(total)])
            rank_floor += len(support)-fresh_dimension(s,len(H),len(support))
            offset +=len(support)
        actual = f.rank(constraints) if constraints and total else 0
        assert actual>=rank_floor
        cap()
    for s in (1,2,4,8,16,32,64):
        for h in range(s):
            for w in range(2*s+1):
                dim = fresh_dimension(s,h,w)
                assert (2*s-1)*(2*s-h-dim)>=s*(2*s-w)
    print('512 multitime fresh-column ranks and all scalar tradeoffs through state64 PASS.',flush=True)

    # Check the output volume argument on its full small integer domain.
    for u in range(1,13):
        q = 16
        for w in range(u+1):
            for ceiling in range(u+1):
                v = min(w,ceiling)
                actual = sum(math.comb(u,j)*q**min(j,ceiling) for j in range(w+1))
                bound = (u+1)*math.comb(u,max(v,u//2))*q**v
                assert actual<=bound
    # Directed rational operations contain exact Fraction results.
    with localcontext() as context:
        context.prec = 40
        for _ in range(512):
            a,b = (Fraction(rng.randrange(-100,101),rng.randrange(1,100)) for _ in range(2))
            x,y = rational(a),rational(b)
            for enclosure,value in ((x+y,a+b),(x-y,a-b),(x*y,a*b)):
                assert Fraction(enclosure.lo)<=value<=Fraction(enclosure.hi)
        # Independent rational atanh series checks the logarithm enclosure.
        for value in (Fraction(1,2),Fraction(2),Fraction(3,2),Fraction(7,8),Fraction(17,16)):
            invert = value<1
            target = 1/value if invert else value
            z = (target-1)/(target+1)
            lower = 2*sum(z**(2*j+1)/(2*j+1) for j in range(100))
            upper = lower+2*z**201/(201*(1-z*z))
            if invert:
                lower,upper = -upper,-lower
            enclosure = rational(value).log()
            assert Fraction(enclosure.lo)<=lower<=upper<=Fraction(enclosure.hi)
    try:
        certify(40,0)
    except RuntimeError as error:
        assert 'time cap' in str(error)
    else:
        raise AssertionError('Timeout must fail closed')
    print('Output-volume and directed rational checks PASS; timeout rejects without certificate.',flush=True)
    print(f'Elapsed {time.monotonic()-started:.3f}s. Exact local tests, not a whole-code certificate.')


if __name__=='__main__':
    main()
