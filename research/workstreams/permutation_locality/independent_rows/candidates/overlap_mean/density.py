"""Guarded integer evaluation of the translated overlap-chord bound.

The old normalized feedback law is not substituted for the true law.
Rounded counts give pointwise class-mass uppers; the overlap budget comes
from the original input distribution. See ../../OVERLAP_MEAN.md.
"""
from fractions import Fraction as Q
from itertools import combinations_with_replacement

import numpy as np
from flint import arb,ctx

from mean import build as overlap_build,baseline
from density_extend import rounded_counts,refine
from exact_feedback import inverse_counts
from feedback_character_census import character_polynomials
from feedback_density import class_transforms,class_convolution,dyadic_precision
from group_moment import maps
from pair_tail import hadamard


LIMIT=(1<<63)-1


def powers(levels,weight,tilt,bits):
    tilt=Q(tilt)
    if tilt<=0 or min(levels)<weight or not 1<=bits<=62:
        raise ValueError('positive tilt and nonnegative chord output exponents required')
    lam=arb(tilt.numerator)/tilt.denominator
    scale=1<<bits
    result={}
    for v in levels:
        pair=[]
        for w in (v+weight,v-weight):
            upper=baseline.up((-lam*w).exp()*scale).fmpq()
            n,d=int(upper.numerator),int(upper.denominator)
            pair.append(min(scale,(n+d-1)//d))
        result[v]=tuple(pair)
    return result


def shape_bounds(counts,denominator,weights,expansion,prepared,mean,tilts,
                 *,rounds=2,census_bits=40,maximum_bits=44):
    """Upper C->C and class U_v->C bounds, including one lazy factor."""
    if (type(rounds) is not int or rounds<1 or not weights
            or any(type(w) is not int or not 1<=w<=4 for w in weights)
            or not isinstance(mean,Q) or not 0<=mean<=sum(weights)):
        raise ValueError('positive packet weights, updates and an exact overlap cap required')
    levels,sizes,transforms=prepared
    W=sum(weights)
    if min(levels)<W:
        raise ValueError('this bounded implementation requires W <= the minimum expansion weight')
    small,total,unit=rounded_counts(counts,denominator,census_bits=census_bits)
    if W*total>LIMIT:
        raise ValueError('weighted count capacity exceeds the signed-int64 guard')
    bits=dyadic_precision(W*total,maximum_bits)
    tilts=tuple(dict.fromkeys(map(str,tilts)))
    if not tilts:
        raise ValueError('at least one output tilt required')
    coefficients={tilt:powers(levels,W,tilt,bits) for tilt in tilts}
    budget=mean*denominator/unit
    budget=min(W*total,(budget.numerator+budget.denominator-1)//budget.denominator)
    assert 0<=budget<=W*total<=LIMIT
    remaining=np.full(len(small),budget,dtype=np.int64)
    weighted={tilt:np.zeros(len(small),dtype=np.int64) for tilt in tilts}
    uniforms={tilt:{} for tilt in tilts}
    transformed=hadamard(small)
    combined=np.zeros(len(small),dtype=np.int64)
    scalar_denominator=denominator*W*(1<<bits)*(1<<rounds)
    for v in levels:
        hits=class_convolution(transformed,small,expansion,transforms[v],total,v)
        combined+=hits
        capacity=hits*W
        allocated=np.minimum(remaining,capacity)
        remaining-=allocated
        # A class-specific bound may spend the entire overlap budget in
        # that class; it must not inherit the joint allocation above.
        individual=np.minimum(capacity,budget)
        for tilt in tilts:
            low,high=coefficients[tilt][v]
            weighted[tilt]+=(capacity-allocated)*low+allocated*high
            local=(capacity-individual)*low+individual*high
            assert 0<=int(local.min())<=int(local.max())<=W*total*(1<<bits)<=LIMIT
            peak=int(local[1:].max())
            uniforms[tilt][v]=baseline.up(arb(peak*unit)/(scalar_denominator*sizes[v]))
    assert np.array_equal(combined,total-small)
    result={}
    for tilt,values in weighted.items():
        assert 0<=int(values.min())<=int(values.max())<=W*total*(1<<bits)<=LIMIT
        result[tilt]=dict(density=baseline.up(arb(int(values[1:].max())*unit)/scalar_denominator),
                          uniform=uniforms[tilt],rounding_factor=Q(total*unit,denominator),
                          overlap_budget=budget,dyadic_bits=bits)
    return result


def build(maximum,tilts,*,minimum=5,precision=192,overlap=None):
    if (type(minimum) is not int or type(maximum) is not int
            or not 1<=minimum<=maximum<=10 or precision<128):
        raise ValueError('complete shape range within one through ten required')
    ctx.prec=precision
    overlap=overlap_build(maximum) if overlap is None else overlap
    required={shape for j in range(minimum,maximum+1)
              for shape in combinations_with_replacement(range(1,5),j)}
    if not required<=set(overlap):
        raise ValueError('complete original-law overlap bounds required')
    shapes,coefficients,denominators,_,inverse=character_polynomials(maximum)
    images,_,spectrum=maps()
    expansion=np.array([x.bit_count() for x in images],dtype=np.int64)
    prepared=class_transforms(expansion)
    assert prepared[1]==dict(spectrum)
    tilts=tuple(dict.fromkeys(map(str,tilts)))
    result={tilt:{} for tilt in tilts}
    feedback={};checked=0;inflation=Q(1)
    for shape,row in zip(shapes,coefficients):
        if sum(shape)<minimum: continue
        weights=tuple(b for b in range(1,5) for _ in range(shape[b-1]))
        denominator=int(denominators[shape])
        if denominator!=overlap[weights][0]:
            raise ValueError('feedback and overlap counting measures disagree')
        counts=inverse_counts(row[inverse],denominator)
        local=shape_bounds(counts,denominator,weights,expansion,prepared,overlap[weights][2],tilts)
        for tilt,record in local.items():
            result[tilt][weights]=record
            inflation=max(inflation,record['rounding_factor'])
        # Exact original counts, not rounded counts, for the zero-return
        # refinement. This avoids repeating a feedback census afterward.
        classes={v:int(counts[expansion==v].sum()) for v in prepared[0]}
        feedback[weights]=int(counts[0]),int(counts[1:].max()),denominator,classes
        assert sum(classes.values())+int(counts[0])==denominator
        checked+=1
        if checked%50==0:
            print('OVERLAP DENSITY complete rounded-count convolutions',checked,'shapes',flush=True)
    assert set(feedback)==required and all(set(rows)==required for rows in result.values())
    print('OVERLAP DENSITY complete shapes',checked,'largest count inflation',inflation,flush=True)
    return result,feedback
