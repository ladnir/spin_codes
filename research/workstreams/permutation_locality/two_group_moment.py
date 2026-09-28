"""Exact same-epoch census and conservative transfers for two active groups.

The c=1 or c=2 distribution places groups in distinct four- or eight-bit
windows. This module counts two windows sharing an IMT epoch. Floating-point
transfer matrices are diagnostic, not outward certificates.
"""
from collections import Counter, defaultdict
from itertools import combinations, combinations_with_replacement
from fractions import Fraction

import numpy as np

from group_moment import maps
from joint_support import gf2_rank
from rank_two_moment import region_transfers


def collision_census(columns_per_bundle=2):
    assert columns_per_bundle in (1,2)
    width=4*columns_per_bundle
    window_count=128//width
    mask_count=1<<width
    images, columns, spectrum = maps()
    allowed = defaultdict(list)
    for mask in range(1, mask_count):
        allowed[tuple(sorted(((mask>>(4*j))&15).bit_count() for j in range(columns_per_bundle)))].append(mask)
    syndromes = np.zeros((window_count, mask_count), dtype=np.int32)
    for j in range(window_count):
        for mask in range(1, mask_count):
            bit = mask & -mask
            syndromes[j, mask] = syndromes[j, mask ^ bit] ^ columns[width*j + bit.bit_length()-1]
        assert len(set(map(int, syndromes[j]))) == mask_count
    ranks = Counter(gf2_rank(columns[width*j:width*(j+1)] + columns[width*k:width*(k+1)])
                    for j, k in combinations(range(window_count), 2))
    weights = np.array([v.bit_count() for v in images], dtype=np.int16)
    windows = np.array([[(v >> (width*j)) & (mask_count-1) for v in images] for j in range(window_count)], dtype=np.uint8)
    result = {}
    for a, b in combinations_with_replacement(sorted(allowed), 2):
        ma = np.array(allowed[a], dtype=np.int32)[:, None]
        mb = np.array(allowed[b], dtype=np.int32)[None, :]
        atoms = np.zeros(1 << 19, dtype=np.int64)
        cancel = np.zeros((129, 145), dtype=np.int64)
        for j in range(window_count):
            for k in range(window_count):
                if j == k:
                    continue
                q = syndromes[j, ma] ^ syndromes[k, mb]
                v = weights[q]
                overlap = np.bitwise_count(windows[j, q] & ma) + np.bitwise_count(windows[k, q] & mb)
                out = v + sum(a) + sum(b) - 2*overlap.astype(np.int16)
                assert np.all(out >= 0)
                np.add.at(atoms, q.ravel(), 1)
                nonzero = q != 0
                np.add.at(cancel, (v[nonzero], out[nonzero]), 1)
        choices = window_count*(window_count-1)*len(allowed[a])*len(allowed[b])
        assert int(atoms.sum()) == choices
        assert int(cancel.sum()) == choices-int(atoms[0])
        row = (choices, int(atoms[0]), int(atoms[1:].max()),
               {v: Counter({w: int(cancel[v,w]) for w in np.flatnonzero(cancel[v])})
                for v in spectrum})
        result[a,b] = result[b,a] = row
    # A rank-r map on 2*width coordinates has 2^(2*width-r)-1 nonzero kernel
    # vectors. Individual windows are injective, so both halves are nonzero.
    expected = 2*sum(count*((1 << (2*width-rank))-1) for rank,count in ranks.items())
    assert sum(row[1] for row in result.values()) == expected
    worst = max((Fraction(row[1],row[0]), key) for key,row in result.items())
    print('Exact two-window rank histogram:', dict(sorted(ranks.items())), flush=True)
    print('Same-epoch zero feedback: total', expected, 'worst shape probability', worst, flush=True)
    return spectrum, result


def collision_transfers(data, tilt):
    """Positive seven-coordinate envelope, with triangle bounds for moments."""
    spectrum, counts = data
    levels = sorted(spectrum)
    m = (1 << 19)-1
    mass = np.array([spectrum[v]/m for v in levels])
    result = {}
    for (a,b), (choices,zero,maxatom,cancel) in counts.items():
        weight = sum(a)+sum(b)
        matrix = np.zeros((7,7))
        matrix[0,0] = np.exp(-tilt*weight)*zero/choices
        matrix[0,1] = np.exp(-tilt*weight)*(choices-zero)/choices
        least = min(w for row in cancel.values() for w in row)
        for i in range(1,7):
            v = 48 if i == 1 else levels[i-2]
            f = np.exp(-tilt*(v-weight))
            if i == 1:
                c = min(f, np.exp(-tilt*least)*maxatom/choices)
            else:
                c = sum(count*np.exp(-tilt*w) for w,count in cancel[v].items())/(choices*spectrum[v])
            matrix[i,0] = c/2+f/(2*m)
            matrix[i,1] = f/2
            matrix[i,2:] = f*mass/2
        result[a,b] = matrix
    return result


def placement_regions(empty, single, collision, epochs=128, windows=16):
    """Return operators indexed by two occupied-column counts (0,1,2).

    Different-epoch terms count ordered group placements. Same-epoch terms
    average distinct window pairs. No independence replacement is made.
    """
    result = {(0,0): np.linalg.matrix_power(empty, epochs)}
    one = [np.zeros_like(empty) for _ in single]
    occupancies=range(1,len(single)+1)
    two = {(a,b): np.zeros_like(empty) for a in occupancies for b in occupancies}
    same = {key: np.zeros_like(empty) for key in two}
    power = np.eye(len(empty))
    for _ in range(epochs):
        for a,b in two:
            two[a,b] = two[a,b]@empty + one[a-1]@single[b-1] + one[b-1]@single[a-1]
            same[a,b] = same[a,b]@empty + power@collision[a,b]
        one = [value@empty + power@active for value,active in zip(one,single)]
        power = power@empty
    for a in occupancies:
        result[a,0] = result[0,a] = one[a-1]/epochs
    # Total ordered slots = E*W*(E*W-1); different epochs have W^2
    # choices per ordered epoch pair, same epochs W*(W-1) per epoch.
    denominator = epochs*(epochs*windows-1)
    for key in two:
        result[key] = (windows*two[key] + (windows-1)*same[key])/denominator
    return result


def worst_regions(single_data, pair_data, tilt, early_max=False, return_shapes=False):
    from random_group_verify import epoch_transfers
    bundle_width=len(next(iter(single_data[1])))
    window_count=32//bundle_width
    epochs=64*bundle_width
    shapes, exact = epoch_transfers(single_data, str(tilt), windows=window_count)
    arrays = [np.array([[float(row[i,j]) for j in range(7)] for i in range(7)]) for row in exact]
    pair = collision_transfers(pair_data,tilt)
    if not early_max:
        # Keep each group's shape fixed throughout the region averaging.
        # Taking an entrywise maximum at the epoch level allows incompatible
        # shape choices at different possible placements and is much looser.
        active = np.array(arrays[1:])
        ia = np.repeat(np.arange(len(shapes)),len(shapes))
        ib = np.tile(np.arange(len(shapes)),len(shapes))
        collision = np.array([pair[shapes[a],shapes[b]] for a,b in zip(ia,ib)])
        one = np.zeros_like(active)
        two = np.zeros_like(collision)
        same = np.zeros_like(collision)
        power = np.eye(7)
        for _ in range(epochs):
            two = two@arrays[0] + one[ia]@active[ib] + one[ib]@active[ia]
            same = same@arrays[0] + power@collision
            one = one@arrays[0] + power@active
            power = power@arrays[0]
        result = {(0,0): power}
        occupancy = np.array([sum(w!=0 for w in s) for s in shapes])
        for b in range(1,bundle_width+1):
            result[b,0] = result[0,b] = np.max(one[occupancy==b],axis=0)/epochs
        pairs = (window_count*two+(window_count-1)*same)/(epochs*2047)
        if return_shapes:
            return power,one/epochs,pairs,shapes,ia,ib
        for a in range(1,bundle_width+1):
            for b in range(1,bundle_width+1):
                result[a,b] = np.max(pairs[(occupancy[ia]==a)&(occupancy[ib]==b)],axis=0)
        return result
    assert not return_shapes and bundle_width==2
    singles = [np.maximum.reduce([arrays[k+1] for k,s in enumerate(shapes)
                                 if sum(w!=0 for w in s)==b]) for b in (1,2)]
    collisions = {(a,b): np.maximum.reduce([matrix for (s,t),matrix in pair.items()
                                           if sum(w!=0 for w in s)==a and sum(w!=0 for w in t)==b])
                  for a in (1,2) for b in (1,2)}
    return placement_regions(arrays[0],singles,collisions)


def self_test():
    # Direct enumeration checks ordered placement, collisions, and ordering
    # of noncommuting matrices. Integer-valued matrices make equality exact.
    empty = np.array([[1.,1.],[0.,1.]])
    single = [np.array([[1.,0.],[1.,1.]]), np.array([[2.,1.],[0.,1.]])]
    pair = {(a,b): np.array([[a+b,1.],[1.,a*b]]) for a in (1,2) for b in (1,2)}
    for epochs in (1,2,3):
        for windows in (2,3):
            actual = placement_regions(empty,single,pair,epochs,windows)
            for a,b in pair:
                total = np.zeros((2,2))
                for x in range(epochs*windows):
                    for y in range(epochs*windows):
                        if x == y:
                            continue
                        product = np.eye(2)
                        for e in range(epochs):
                            if x//windows == e == y//windows:
                                step = pair[a,b]
                            elif x//windows == e:
                                step = single[a-1]
                            elif y//windows == e:
                                step = single[b-1]
                            else:
                                step = empty
                            product = product@step
                        total += product
                assert np.array_equal(actual[a,b],total/(epochs*windows*(epochs*windows-1)))
    print('Two-group placement recurrence passed exhaustive small noncommuting tests',flush=True)


if __name__ == '__main__':
    self_test()
    data = collision_census()
    zero = collision_transfers(data,0)
    assert all(np.all(matrix.sum(axis=1)>=1-1e-14) for matrix in zero.values())
    print('All 196 ordered shape pairs checked; zero-tilt mass checks passed',flush=True)
    print('Local component only: no two-group output-weight certificate yet.')
