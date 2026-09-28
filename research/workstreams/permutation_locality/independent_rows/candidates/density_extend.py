"""Extend the exact-feedback convolution bound with upward integer counts.

For a census n of total D, let n'[s]=ceil(n[s]/2^r). Then n/D is
pointwise dominated by 2^r*n'/D. The smaller integer census permits the
existing guarded exact XOR convolutions. Its normalized bounds are scaled
by sum(n')*2^r/D; rounding is never interpreted as a probability law for
the actual encoder. Only positive density contributions are changed.
"""
from fractions import Fraction as Q
from itertools import combinations_with_replacement

import numpy as np
from flint import arb,ctx

from mass_density_screen import baseline
from exact_feedback import inverse_counts
from feedback_character_census import character_polynomials
from feedback_density import class_transforms,shape_bounds
from group_moment import maps
from group_rank_one_verify import up
from occupancy_memory import C,U
from pair_tail import hadamard


def rounded_counts(counts,denominator,*,census_bits=40):
    counts=np.asarray(counts)
    if (counts.ndim!=1 or counts.dtype!=np.dtype(np.int64) or not len(counts)
            or len(counts)&(len(counts)-1) or int(counts.min())<0
            or not isinstance(denominator,int) or isinstance(denominator,bool)
            or not 0<denominator<1<<63 or sum(map(int,counts))!=denominator
            or not isinstance(census_bits,int) or isinstance(census_bits,bool)
            or not 1<=census_bits<=44):
        raise ValueError('nonnegative int64 counting measure and exact total required')
    shift=max(0,denominator.bit_length()-census_bits)
    unit=1<<shift
    if int(counts.max())+unit-1>=1<<63:
        raise ValueError('upward count rounding exceeds its integer guard')
    rounded=(counts+(unit-1))>>shift
    total=sum(map(int,rounded))
    if len(counts)*total>=1<<63:
        raise ValueError('rounded census exceeds the exact convolution guard')
    assert denominator<=unit*total<denominator+unit*len(counts)
    return rounded,total,unit


def bounds(counts,denominator,weights,expansion,prepared,tilts,*,rounds=2,census_bits=40):
    small,total,unit=rounded_counts(counts,denominator,census_bits=census_bits)
    values=shape_bounds(hadamard(small),total,weights,expansion,prepared,tilts,rounds=rounds)
    factor=arb(total*unit)/denominator
    for record in values.values():
        record['density']=up(record['density']*factor)
        record['uniform']={v:up(x*factor) for v,x in record['uniform'].items()}
        record['rounding_factor']=Q(total*unit,denominator)
    return values


def build(maximum,tilts,*,minimum=7,precision=192,rounds=2):
    if (not isinstance(minimum,int) or not isinstance(maximum,int)
            or not 1<=minimum<=maximum<=10 or precision<128):
        raise ValueError('bounded complete shape range and at least 128-bit arithmetic required')
    ctx.prec=precision
    shapes,coefficients,denominators,_,inverse=character_polynomials(maximum)
    images,_,spectrum=maps()
    expansion=np.array([image.bit_count() for image in images],dtype=np.int64)
    prepared=class_transforms(expansion)
    assert prepared[1]==dict(spectrum)
    tilts=tuple(dict.fromkeys(map(str,tilts)))
    result={tilt:{} for tilt in tilts}; checked=0; inflation=Q(1)
    for shape,row in zip(shapes,coefficients):
        if sum(shape)<minimum: continue
        weights=tuple(b for b in range(1,5) for _ in range(shape[b-1]))
        denominator=int(denominators[shape])
        counts=inverse_counts(row[inverse],denominator)
        local=bounds(counts,denominator,weights,expansion,prepared,tilts,rounds=rounds)
        for tilt,record in local.items():
            result[tilt][weights]=record
            inflation=max(inflation,record['rounding_factor'])
        checked+=1
        if checked%50==0:
            print('EXTENDED DENSITY exact rounded-census convolutions',checked,'shapes',flush=True)
    required={shape for j in range(minimum,maximum+1) for shape in combinations_with_replacement(range(1,5),j)}
    assert checked==len(required) and all(set(rows)==required for rows in result.values())
    print('EXTENDED DENSITY complete shapes',checked,'maximum mass inflation',inflation,flush=True)
    return result


def refine(base,records,spectrum,penalty,*,minimum=7,maximum=10):
    required={shape for j in range(minimum,maximum+1) for shape in combinations_with_replacement(range(1,5),j)}
    if set(records)!=required or not 0<Q(penalty)<=1 or maximum>=len(base):
        raise ValueError('complete shape range and matching operators required')
    values={j:[arb(0) for _ in range(6)] for j in range(minimum,maximum+1)}
    rho=Q(penalty)
    for shape,record in records.items():
        if set(record['uniform'])!=set(spectrum):
            raise ValueError('every nonzero expansion class is required')
        scale=rho**shape.count(4)
        factor=arb(scale.numerator)/scale.denominator
        row=[record['density'],*[record['uniform'][v] for v in sorted(spectrum)]]
        if any(not x.is_finite() or not x>=0 for x in row):
            raise ValueError('finite nonnegative outward coefficients required')
        for i,value in enumerate(row):
            values[len(shape)][i]=max(values[len(shape)][i],up(value*factor))
    result=[t*1 for t in base]
    changes=0
    for j,row in values.items():
        for source,value in zip((C,*range(U,U+5)),row):
            changes+=int(value<result[j][source,C])
            result[j][source,C]=min(result[j][source,C],value)
    print('EXTENDED DENSITY tightened scalar entries',changes,flush=True)
    return result
