"""Exact small-field tests of changing-pair fresh-variable accounting."""
from decimal import Decimal, localcontext
from fractions import Fraction
import itertools
import math
import random
import time
from check_maps import Field
from changing_pair_probe import pair_rank_lower
from changing_pair_certificate import certify, box, mask_cap
from feedback_support_certificate import box as old_box


def main():
    started = time.monotonic()
    def cap():
        if time.monotonic()-started>45:
            raise RuntimeError('45-second exact-test cap')
    rng = random.Random(20261022)
    f = Field(4,0x13)
    improved = 0
    for s in (2,4):
        b = 2*s
        B = f.interpolation_matrix(f.points[:b],f.points[b:b+s])
        A = f.interpolation_matrix(f.points[b:b+s],f.points[:2*b])
        # Every trajectory has three disjoint fresh variable pairs. Earlier
        # states are retained as linear forms, not sampled again at a pair.
        for case in range(1024):
            supports = [rng.sample(range(b),rng.randrange(b+1)) for _ in range(6)]
            scales = [[rng.randrange(1,16) for _ in W] for W in supports]
            holes = [rng.sample(range(b),rng.randrange(s)) for _ in supports]
            total = sum(map(len,supports))
            forms = [[0]*total for _ in range(s)]
            constraints,offset,lower = [],0,0
            for step,(W,R,H) in enumerate(zip(supports,scales,holes)):
                for j,(col,scale) in enumerate(zip(W,R)):
                    for r in range(s):
                        forms[r][offset+j] ^= f.table[B[r][col]][scale]
                for h in H:
                    # First observation may use either alternating phase.
                    row = ((step+case)%2)*b+h
                    constraints.append([f.xor(f.table[A[row][r]][forms[r][j]]
                                               for r in range(s)) for j in range(total)])
                offset += len(W)
                if step%2:
                    h0,h1 = len(holes[step-1]),len(H)
                    u,v = len(supports[step-1]),len(W)
                    rank = pair_rank_lower(s,h0,h1,u,v)
                    separate = max(0,h0+min(s,u)-s)+max(0,h1+min(s,v)-s)
                    improved += rank>separate
                    if case%2:
                        rank = max(rank,min((h0+1)//2,u,s)+max(0,h1+min(s,v)-s))
                    lower += rank
            actual = f.rank(constraints) if constraints and total else 0
            assert actual>=lower, (s,case,actual,lower)
            cap()
    print(f'2048 three-pair exact GF16 trajectories PASS; '
          f'joint floor improves {improved} local pairs over separate floors.',flush=True)

    # The two individual A1 evaluation cosets each give a full-spark
    # restriction of A1 B, even for W crossing the two input halves.
    s,b = 4,8
    B = f.interpolation_matrix(f.points[:b],f.points[b:b+s])
    A = f.interpolation_matrix(f.points[b:b+s],f.points[b:2*b])
    AB = [[f.xor(f.table[row[r]][B[r][c]] for r in range(s))
           for c in range(b)] for row in A]
    minors = 0
    for half in (0,1):
        for size in range(1,s+1):
            for rows in itertools.combinations(range(half*s,(half+1)*s),size):
                for cols in itertools.combinations(range(b),size):
                    assert f.rank([[AB[r][c] for c in cols] for r in rows])==size
                    minors += 1
    assert minors==988
    print('988 exact within-coset A1 B minors PASS (cross-half columns included).',flush=True)

    # Exhaust all local parameters at small state, and boundary/stride cases
    # at state64. Verify the scalar dimension used by the certificate.
    for s in (1,2,4,8,16,64):
        hs = range(s) if s<=16 else (0,1,15,31,32,47,62,63)
        ws = range(2*s+1) if s<=16 else (0,1,2,31,32,63,64,65,96,127,128)
        for h0,h1,u,v in itertools.product(hs,hs,ws,ws):
            r = pair_rank_lower(s,h0,h1,u,v)
            assert 0<=r<=min(h0+h1,u+v)
            assert (2*s-1)*(4*s-h0-h1-(u+v-r))>=s*(4*s-u-v)
        cap()

    # Independent small polynomial enumeration: all local hole-mask choices
    # are counted, including when an allowed changing pair does not change.
    def multiply(a,b):
        c = [0]*(len(a)+len(b)-1)
        for i,x in enumerate(a):
            for j,y in enumerate(b):
                c[i+j] += x*y
        return c
    for s in (1,2,3):
        b = 2*s
        same = [math.comb(2*b,h) for h in range(s)]
        single = [math.comb(b,h) for h in range(s)]
        changed = multiply(single,single)
        for d in range(4):
            coeff = single
            for factor in [same]*(3-d)+[changed]*d:
                coeff = multiply(coeff,factor)
            assert len(coeff)==(s-1)*(4+d)+1
            assert coeff[-1]==math.comb(2*b,s-1)**(3-d)*math.comb(b,s-1)**(2*d+1)
            for x in (Fraction(1,2),Fraction(1),Fraction(2)):
                value = sum(c*x**h for h,c in enumerate(coeff))
                assert all(c<=value/x**h for h,c in enumerate(coeff))
    # Exact extremizers of the coarse two-state envelope. The sparse
    # extremizer cannot repeat across high-hole boundaries with live states.
    scores = [(127*(u+v-pair_rank_lower(64,63,63,u,v))-64*(u+v),u,v)
              for u in range(129) for v in range(129)]
    peak = max(a for a,_,_ in scores)
    assert peak==126
    assert [(u,v) for a,u,v in scores if a==peak]==[(1,128),(128,128)]
    with localcontext() as ctx:
        ctx.prec = 80
        for a,b in ((0,0),(500000,500128),(516096,516096)):
            # Adding the position factor log binom(P,0) introduces one
            # extra outward-rounding operation, so compare to precision.
            assert abs(box(0,a,b,10)-old_box(a,b,10))<Decimal('1e-65')
        assert math.isfinite(float(mask_cap(32,518112,518112)))
    try:
        certify(32,40,0)
    except RuntimeError as error:
        assert 'time cap' in str(error)
    else:
        raise AssertionError('Immediate timeout must fail closed')
    print('Local tradeoffs, mask coefficients, d=0 regression and timeout PASS.',flush=True)
    print(f'Elapsed {time.monotonic()-started:.3f}s; no whole-code claim.')


if __name__=='__main__':
    main()
