"""Exact local checks for the alternating-map candidate; no full certificate."""
import itertools
import random
import time
from check_maps import Field,weight
from alternating_maps import AlternatingMaps,run_hole_bound


def main():
    started = time.monotonic()
    f = Field(4,0x13)
    maps = AlternatingMaps(f,4)
    minimum,states = 16,0
    for q in itertools.product(range(16),repeat=4):
        if not any(q):
            continue
        outputs = [maps.a_at(list(q),phase) for phase in (0,1)]
        weights = [weight(y) for y in outputs]
        assert min(weights)>=5
        assert sum(weights)>=13
        for length in range(1,8):
            for phase in (0,1):
                holes = sum(8-weights[(phase+t)%2] for t in range(length))
                assert holes<=run_hole_bound(4,length)
        minimum = min(minimum,sum(weights))
        states +=1
        if states%1024==0 and time.monotonic()-started>50:
            raise RuntimeError('50-second cap')
    assert minimum==13
    print(f'{states} GF16 nonzero states: per-batch>=5, pair minimum13; '
          'both run parities through length7 PASS',flush=True)

    f = Field(2,0x7)
    maps = AlternatingMaps(f,1)
    fixtures = 0
    for phase in (0,1):
        for scales in itertools.product(range(1,4),repeat=2):
            for q in range(1,4):
                aq = maps.a_at([q],phase)
                for x in itertools.product(range(4),repeat=2):
                    counts = [0]*4
                    for d in itertools.product(range(4),repeat=2):
                        y = [a^f.table[c][z] for a,c,z in zip(x,d,aq)]
                        next_q = q^maps.bmap([f.table[r][z] for r,z in zip(scales,y)])[0]
                        counts[next_q] +=1
                    assert counts==[4]*4
                    fixtures +=1
    print(f'{fixtures} exact GF4 uniform-refresh fixtures across both phases PASS',flush=True)

    f = Field(8,0x11b)
    maps = AlternatingMaps(f,64)
    alternate = f.interpolation_matrix(f.points[128:192],f.points[128:256])
    rng = random.Random(210347)
    for _ in range(16):
        q = [rng.randrange(256) for _ in range(64)]
        assert maps.a_at(q,0)==maps.a(q)
        assert maps.a_at(q,1)==f.apply(alternate,q)
        assert maps.a_at(q,1)[:64]==q
    supports = [[],[]]
    for column in zip(*maps.b_matrix):
        for phase in (0,1):
            supports[phase].append(weight(maps.a_at(list(column),phase)))
    assert supports[0]==[65]*128 and supports[1]==[128]*128
    print('Actual GF256: FFT/Lagrange and copied state half PASS; '
          'A0 B column weights all65, A1 B column weights all128.',flush=True)
    state = [0]*64
    records = []
    for t in range(16):
        x = [rng.randrange(256) for _ in range(128)]
        d = [rng.randrange(256) for _ in range(128)]
        r = [rng.randrange(1,256) for _ in range(128)]
        aq = maps.a_at(state,t%2)
        y = [a^f.table[c][z] for a,c,z in zip(x,d,aq)]
        state = [a^b for a,b in zip(state,maps.bmap([f.table[c][z] for c,z in zip(r,y)]))]
        records.append((x,d,r,y))
    state = [0]*64
    for t,(x,d,r,y) in enumerate(records):
        aq = maps.a_at(state,t%2)
        assert x==[a^f.table[c][z] for a,c,z in zip(y,d,aq)]
        state = [a^b for a,b in zip(state,maps.bmap([f.table[c][z] for c,z in zip(r,y)]))]
    print(f'16 actual-size recurrence/inverse steps PASS; elapsed {time.monotonic()-started:.3f}s. '
          'No native timing or full distance certificate.',flush=True)


if __name__=='__main__':
    main()
