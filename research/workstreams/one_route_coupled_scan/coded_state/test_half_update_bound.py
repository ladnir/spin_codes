"""Exact map/root-bound checks. No output files or whole-code claim."""
import itertools
import random
import time
from check_maps import Field,weight
from alternating_maps import AlternatingMaps
from half_update_bound import second_holes_bound,pair_weight_bound,live_nonzero_update_count,dimension_cap


def root_values(f,roots,points):
    result = [1]*len(points)
    for root in roots:
        result = [f.table[a][x^root] for a,x in zip(result,points)]
    return result


def check_pair(f,maps,q,y):
    delta = maps.bmap(y)
    nxt = [a^b for a,b in zip(q,delta)]
    if not any(q) or not any(nxt):
        return None
    support = {j for j,x in enumerate(y) if x}
    s,w = maps.s,len(support)
    half = min(support)//s
    assert all(j//s==half for j in support)
    h0set = {j for j,x in enumerate(maps.a_at(q,0)) if not x}
    h1 = 2*s-weight(maps.a_at(nxt,1))
    common = set(range(half*s,(half+1)*s))
    r = len(h0set & (common-support))
    i = len(h0set & support)
    upper = second_holes_bound(s,len(h0set),r,i,w)
    assert h1<=upper,(s,w,len(h0set),r,i,h1,upper)
    # This root set is the structural reason for the first half of the bound.
    ad = maps.a_at(delta,0)
    assert all(ad[j]==0 for j in common-support)
    return 4*s-len(h0set)-h1


def main():
    started = time.monotonic()
    f = Field(4,0x13)
    maps = AlternatingMaps(f,4)
    count,minimum = 0,16
    # Projective normalization covers every amplitude for a one-symbol update.
    for q in itertools.product(range(16),repeat=4):
        if not any(q):
            continue
        for j in range(8):
            y = [int(k==j) for k in range(8)]
            result = check_pair(f,maps,list(q),y)
            if result is not None:
                minimum = min(minimum,result)
                count +=1
        if count%4096<8 and time.monotonic()-started>45:
            raise RuntimeError('45-second cap')
    assert count==524272 and minimum==11
    print(f'{count} exhaustive GF16 one-update pairs satisfy the root bound; exact minimum11.',flush=True)

    # Independently count both live endpoints and a nonzero update, including
    # high-hole fixtures whose whole kernel consists only of endpoint zeros.
    maps = AlternatingMaps(f,2)
    aa = [maps.a_matrix,f.interpolation_matrix(f.points[4:6],f.points[4:8])]
    checked = 0
    for support,holes0,holes1 in (([0],[1],[0]),([0],[0],[1]),
                                  ([0,1],[0],[0]),([2,3],[0],[1]),
                                  ([0,1],[],[0]),([2,3],[0],[])):
        w = len(support)
        ab1 = [[f.xor(f.table[a][maps.b_matrix[k][j]] for k,a in enumerate(row))
                for j in support] for row in aa[1]]
        matrix = [aa[0][i]+[0]*w for i in holes0]
        matrix += [aa[1][i]+ab1[i] for i in holes1]
        d = 2+w-f.rank(matrix)
        half = support[0]//2
        repeated = len(set(holes0)&(set(range(2*half,2*half+2))-set(support)))
        a,b,c = max(0,w-len(holes1)),max(0,w-len(holes0)+repeated),max(0,2-len(holes0)-len(holes1))
        expected = live_nonzero_update_count(16,d,a,b,c)
        cap = dimension_cap(2,w,len(holes0),len(holes1),repeated)
        assert d<=cap
        assert expected<=live_nonzero_update_count(16,cap,a,b,c)
        actual = 0
        for q in itertools.product(range(16),repeat=2):
            if not any(q):
                continue
            for active in itertools.product(range(16),repeat=w):
                if not any(active):
                    continue
                if any(f.apply(matrix,list(q)+list(active))):
                    continue
                y = [0]*4
                for j,value in zip(support,active):
                    y[j] = value
                nxt = [u^v for u,v in zip(q,maps.bmap(y))]
                actual +=any(nxt)
        assert actual==expected,(support,holes0,holes1,actual,expected)
        checked +=1
    print(f'{checked} exact live-endpoint/nonzero-update fiber counts PASS.',flush=True)

    f = Field(8,0x11b)
    rng = random.Random(211083)
    fixtures = 0
    for s in (4,8,16,32,64):
        maps = AlternatingMaps(f,s)
        for half in (0,1):
            for w in sorted({1,2,s//2,s}):
                for _ in range(8):
                    support = rng.sample(range(half*s,(half+1)*s),w)
                    y = [rng.randrange(1,256) if j in support else 0 for j in range(2*s)]
                    # Construct many holes explicitly rather than expect random states to have them.
                    h = rng.randrange(s//2,s)
                    roots = rng.sample(f.points[:2*s],h)
                    q = root_values(f,roots,f.points[2*s:3*s])
                    result = check_pair(f,maps,q,y)
                    assert result is not None
                    fixtures +=1
        if time.monotonic()-started>55:
            raise RuntimeError('55-second cap')
    print(f'{fixtures} GF256 structured-hole fixtures at state4,8,16,32,64 PASS.',flush=True)
    values = [pair_weight_bound(64,w)[0] for w in (1,2,4,8,12,16,32,64)]
    assert values[0]==150 and values[-1]==130
    assert all(a>=b for a,b in zip(values,values[1:]))
    print(f'K20 map-size local pair floors for w=1,2,4,8,12,16,32,64: {values}',flush=True)
    print(f'Elapsed {time.monotonic()-started:.3f}s. No full-code certificate.',flush=True)


if __name__=='__main__':
    main()
