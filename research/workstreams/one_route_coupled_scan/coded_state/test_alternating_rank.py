"""Bounded exact two-observation checks; not a whole-code certificate."""
import itertools
import random
import time
from collections import Counter
from check_maps import Field
from alternating_maps import AlternatingMaps


def rank(f, a):
    return f.rank(a) if a and a[0] else 0


def product(f, a, b):
    return [[f.xor(f.table[x][y] for x,y in zip(row,col))
             for col in zip(*b)] for row in a]


def main():
    started = time.monotonic()
    f = Field(4,0x13)
    maps = AlternatingMaps(f,4)
    aa = [maps.a_matrix, f.interpolation_matrix(f.points[8:12],f.points[8:16])]
    ab = [product(f,a,maps.b_matrix) for a in aa]
    hist = [Counter(),Counter()]
    example = [None,None]
    checks = 0
    for phase in (0,1):
        for size in range(1,5):
            for rows in itertools.combinations(range(8),size):
                for cols in itertools.combinations(range(8),size):
                    r = rank(f,[[ab[phase][i][j] for j in cols] for i in rows])
                    hist[phase][size,r] +=1
                    if r<size and example[phase] is None:
                        example[phase] = (rows,cols,r)
                    checks +=1
        if time.monotonic()-started>40:
            raise RuntimeError('40-second cap')
    print(f'{checks} square AB minors: phase0={dict(hist[0])}; phase1={dict(hist[1])}',flush=True)
    print(f'First deficient minors: {example}; all-full-rank is not assumed.',flush=True)

    # Incoming state and fresh supported output remain one linear system.
    rng = random.Random(211984)
    pair_checks = strict = 0
    for _ in range(4096):
        phase = rng.randrange(2)
        h0,h1,w = (rng.randrange(4),rng.randrange(4),rng.randrange(9))
        rows0 = rng.sample(range(8),h0)
        rows1 = rng.sample(range(8),h1)
        cols = rng.sample(range(8),w)
        scales = [rng.randrange(1,16) for _ in cols]
        # A0_H Q=0 and A1_H (Q+B R y)=0.
        matrix = [aa[phase][i]+[0]*w for i in rows0]
        matrix += [aa[1-phase][i]+[f.table[ab[1-phase][i][j]][r]
                   for j,r in zip(cols,scales)] for i in rows1]
        actual = rank(f,matrix)
        h = h0+h1
        lower = min(4,h)+max(0,h-4+min(4,w)-4)
        assert actual>=lower,(phase,rows0,rows1,cols,actual,lower)
        if w==0:
            assert actual==min(4,h)
        if w>=4:
            assert actual==h
        strict +=actual>lower
        pair_checks +=1

    # Eliminate incoming state, rather than charge its rank as fresh output rank.
    # The consistency rank on y is exactly joint rank minus rank on Q.
    eliminated = 0
    for phase in (0,1):
        rows0,rows1,cols = [0,1,2],[0,1,2],[0,1]
        matrix = [aa[phase][i]+[0]*2 for i in rows0]
        matrix += [aa[1-phase][i]+[ab[1-phase][i][j] for j in cols] for i in rows1]
        r = rank(f,matrix)
        by_output = Counter()
        for q in itertools.product(range(16),repeat=4):
            if any(f.apply([aa[phase][i] for i in rows0],q)):
                continue
            for y in itertools.product(range(16),repeat=2):
                if not any(f.apply(matrix,list(q)+list(y))):
                    by_output[y] +=1
        assert len(by_output)==16**(2-(r-4))
        assert set(by_output.values())=={1}
        eliminated +=1
    print(f'{pair_checks} two-observation rank checks PASS ({strict} strict); '
          f'{eliminated} exact eliminated-state fiber counts PASS.',flush=True)

    # Does a one-symbol update allow both expansions to have maximal holes?
    # Five variables (four state, one output), six specified holes.
    one_update = Counter()
    live_update = Counter()
    witness = None
    for phase in (0,1):
        for rows0 in itertools.combinations(range(8),3):
            for rows1 in itertools.combinations(range(8),3):
                for col in range(8):
                    matrix = [aa[phase][i]+[0] for i in rows0]
                    matrix += [aa[1-phase][i]+[ab[1-phase][i][col]] for i in rows1]
                    r = rank(f,matrix)
                    one_update[phase,r] +=1
                    if r<5 and witness is None:
                        witness = phase,rows0,rows1,col
                    if r<5:
                        old_zero = rank(f,matrix+[[int(i==j) for j in range(4)]+[0]
                                                  for i in range(4)])==r
                        new_zero = rank(f,matrix+[[int(i==j) for j in range(4)]+[maps.b_matrix[i][col]]
                                                  for i in range(4)])==r
                        live_update[phase,old_zero,new_zero] +=1
        if time.monotonic()-started>50:
            raise RuntimeError('50-second cap')
    print(f'Maximal holes separated by one output: {dict(one_update)}; '
          f'first deficient fixture={witness}',flush=True)
    print(f'Deficient-fixture state-zero classification: {dict(live_update)}',flush=True)
    assert one_update==Counter({(0,5):24640,(0,4):448,(1,5):24640,(1,4):448})
    assert live_update==Counter({(0,False,True):448,(1,True,False):448})
    # Every nonzero one-coordinate output can be projectively scaled to 1.
    # Exhaust the remaining state, excluding activation/reset endpoints.
    minima = [16,16]
    witnesses = [None,None]
    pairs = 0
    for q in itertools.product(range(16),repeat=4):
        if not any(q):
            continue
        aq = [f.apply(a,q) for a in aa]
        for col in range(8):
            if all(q[i]==maps.b_matrix[i][col] for i in range(4)):
                continue
            for phase in (0,1):
                w = sum(x!=0 for x in aq[phase]) + sum(
                    (x^ab[1-phase][i][col])!=0 for i,x in enumerate(aq[1-phase]))
                if w<minima[phase]:
                    minima[phase] = w
                    witnesses[phase] = q,col
                pairs +=1
        if pairs%16384<16 and time.monotonic()-started>50:
            raise RuntimeError('50-second cap')
    print(f'{pairs} projective one-update live pairs: minima={minima}, '
          f'witnesses={witnesses}. Toy evidence only.',flush=True)
    assert pairs==1048544 and minima==[11,11]
    print(f'Elapsed {time.monotonic()-started:.3f}s. No K20 certificate.',flush=True)


if __name__=='__main__':
    main()
