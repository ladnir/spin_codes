"""Two-group diagnostic retaining activation-state distributions within a region.

The extra coordinates are discarded conservatively at region boundaries.
Same-epoch collision births still enter the coarse nonzero coordinate.
No production changes; binary64 exploration, not an outward certificate.
"""
import numpy as np

from random_group_verify import epoch_transfers
from two_group_moment import collision_transfers
from memory_moment import memory_data


def epoch_matrices(single_data,pair_data,tilt):
    spectrum,allowed,moments,atoms,cancel = single_data
    shapes,choices,nonzero,minimum,pair_probability = memory_data(single_data,windows=16)
    exact_shapes,old = epoch_transfers(single_data,str(tilt),windows=16)
    assert shapes == exact_shapes
    n = 7+len(shapes)
    empty = np.zeros((n,n))
    empty[:7,:7] = [[float(old[0][i,j]) for j in range(7)] for i in range(7)]
    levels = sorted(spectrum)
    m = (1 << 19)-1
    powers = np.exp(-tilt*np.arange(145))
    refresh = np.array([0.,0.]+[spectrum[v]/(2*m) for v in levels])
    for a,s in enumerate(shapes):
        f = sum(sum(row.values())*powers[v] for v,row in cancel[s].items())/nonzero[a]
        empty[7+a,:7] = f*refresh
        empty[7+a,7+a] = powers[minimum[a]]/2
    single = np.zeros((len(shapes),n,n))
    for a,s in enumerate(shapes):
        single[a,:7,:7] = [[float(old[a+1][i,j]) for j in range(7)] for i in range(7)]
        single[a,0,7+a] = single[a,0,1]
        single[a,0,1] = 0
        weight = sum(s)
        least = min(w for row in cancel[s].values() for w in row)
        for b in range(len(shapes)):
            f = powers[minimum[b]-weight]
            c = powers[max(minimum[b]-weight,least)]*pair_probability[a,b]
            single[a,7+b,:7] = f*refresh
            single[a,7+b,0] = c/2+f/(2*m)
            single[a,7+b,1] = f/2
    pair = {}
    for (s,t),old_pair in collision_transfers(pair_data,tilt).items():
        matrix = np.zeros((n,n))
        matrix[:7,:7] = old_pair
        choices2,zero,maximum,counts = pair_data[1][s,t]
        weight = sum(s)+sum(t)
        least = min(w for row in counts.values() for w in row)
        for b in range(len(shapes)):
            f = powers[minimum[b]-weight]
            c = powers[max(minimum[b]-weight,least)]*maximum/choices2
            matrix[7+b,:7] = f*refresh
            matrix[7+b,0] = c/2+f/(2*m)
            matrix[7+b,1] = f/2
        pair[s,t] = matrix
    return shapes,empty,single,pair


def collapse(matrix):
    """Discard memory labels, sending their total mass to arbitrary nonzero."""
    result = matrix[...,:7,:7].copy()
    result[...,1] += matrix[...,:7,7:].sum(axis=-1)
    return result


def raw_regions(single_data,pair_data,tilt,epochs=128,check=False):
    shapes,empty,active,pair = epoch_matrices(single_data,pair_data,tilt)
    ia = np.repeat(np.arange(len(shapes)),len(shapes))
    ib = np.tile(np.arange(len(shapes)),len(shapes))
    collision = np.array([pair[shapes[a],shapes[b]] for a,b in zip(ia,ib)])
    one = np.zeros_like(active)
    two = np.zeros_like(collision)
    same = np.zeros_like(collision)
    power = np.eye(len(empty))
    for _ in range(epochs):
        two = two@empty + one[ia]@active[ib] + one[ib]@active[ia]
        same = same@empty + power@collision
        one = one@empty + power@active
        power = power@empty
    pairs = (16*two+15*same)/(epochs*(epochs*16-1))
    if check:
        # Enumerate ordered distinct physical slots, not just epoch pairs.
        for a,b in ((0,0),(0,13),(13,13)):
            direct = np.zeros_like(empty)
            for x in range(epochs*16):
                for y in range(epochs*16):
                    if x == y:
                        continue
                    product = np.eye(len(empty))
                    for e in range(epochs):
                        step = (pair[shapes[a],shapes[b]] if x//16 == e == y//16
                                else active[a] if x//16 == e
                                else active[b] if y//16 == e else empty)
                        product = product@step
                    direct += product
            direct /= epochs*16*(epochs*16-1)
            assert np.allclose(direct,pairs[a*len(shapes)+b],rtol=2e-12,atol=1e-16)
    return collapse(power),collapse(one)/epochs,collapse(pairs),shapes,ia,ib


def self_test(single_data,pair_data):
    for tilt in (0.,.0004,.0025):
        shapes,empty,single,pair = epoch_matrices(single_data,pair_data,tilt)
        if tilt == 0:
            assert np.all(empty.sum(axis=1)>=1-1e-14)
            assert np.all(single.sum(axis=2)>=1-1e-14)
            assert all(np.all(matrix.sum(axis=1)>=1-1e-14) for matrix in pair.values())
        for epochs in (1,2):
            raw_regions(single_data,pair_data,tilt,epochs,check=True)
    # A collapsed state assigns every nonzero memory mass to D, so any
    # pointwise D bound must also dominate the original normalized law.
    test = np.arange(21*21).reshape(21,21)
    assert np.array_equal(collapse(test)[:,1],test[:7,1]+test[:7,7:].sum(axis=1))
    print('Memory transfers: zero-tilt, ordered-slot enumeration, and collapse tests passed',flush=True)


if __name__=='__main__':
    from two_column_moment import census
    from two_group_moment import collision_census
    single,pair = census(),collision_census()
    self_test(single,pair)
