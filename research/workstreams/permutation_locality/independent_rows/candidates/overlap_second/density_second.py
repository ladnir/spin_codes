"""Translated quadratic overlap bounds with positive integer convolutions.

The second-moment budget is for the original input law. Rounded class
counts dominate only positive intercept terms; they never redefine that
law. Each scalar bound is valid for every nonzero target separately.
"""
from fractions import Fraction as Q
from itertools import combinations_with_replacement

import numpy as np
from flint import arb,ctx

from quadratic import witness,upper_intercepts,baseline
from density_extend import rounded_counts,refine
from exact_feedback import inverse_counts
from feedback_character_census import character_polynomials
from feedback_density import class_transforms,class_convolution,dyadic_precision
from group_moment import maps
from pair_tail import hadamard


LIMIT=(1<<63)-1


def dyadic_intercepts(levels,weight,tilt,c,bits):
    if type(bits) is not int or not 1<=bits<=62:
        raise ValueError('positive bounded dyadic precision required')
    scale=1<<bits;result={}
    for v,a in upper_intercepts(levels,weight,tilt,c).items():
        exact=baseline.up(a*scale).fmpq()
        n,d=int(exact.numerator),int(exact.denominator)
        # v >= W implies exp(-tilt*(v-d)) <= 1. Thus a_v <= 1
        # is independently valid even if endpoint rounding exceeds one.
        result[v]=min(scale,(n+d-1)//d)
        assert 0<=result[v]<=scale
    return result


def slopes(classes,weight,second,tilt):
    """Only proposes witnesses; all their inequalities are checked later."""
    total=sum(classes.values(),Q(0))
    if not total:return (Q(0),)
    # Normalization is only for this numerical proposal, never for a bound.
    _,c,_=witness({v:p/total for v,p in classes.items()},weight,second,tilt)
    return tuple(sorted({Q(0),c/2,c,c*2}))


def shape_bounds(counts,denominator,weights,expansion,prepared,second,tilts,
                 *,rounds=2,census_bits=40,maximum_bits=44):
    if (type(rounds) is not int or rounds<1 or not weights
            or any(type(w) is not int or not 1<=w<=4 for w in weights)
            or not isinstance(second,Q) or not 0<=second<=sum(weights)**2):
        raise ValueError('exact nonnegative second cap and valid packet geometry required')
    levels,sizes,transforms=prepared;W=sum(weights)
    if min(levels)<W:
        raise ValueError('requires W no larger than the minimum expansion weight')
    small,total,unit=rounded_counts(counts,denominator,census_bits=census_bits)
    bits=dyadic_precision(total,maximum_bits)
    tilts=tuple(dict.fromkeys(map(str,tilts)))
    if not tilts:raise ValueError('at least one positive output tilt required')
    transformed=hadamard(small);hits={};combined=np.zeros(len(small),dtype=np.int64)
    for v in levels:
        hits[v]=class_convolution(transformed,small,expansion,transforms[v],total,v)
        combined+=hits[v]
    assert np.array_equal(combined,total-small)
    # The max-class envelope is used to propose slopes, not as a proof law.
    proposal={v:Q(int(hits[v][1:].max())*unit,denominator) for v in levels}
    scalar_denominator=denominator*(1<<bits)

    def value(tilt,c,selected_levels):
        powers=dyadic_intercepts(selected_levels,W,tilt,c,bits)
        weighted=np.zeros(len(small),dtype=np.int64)
        for v in selected_levels:weighted+=hits[v]*powers[v]
        assert 0<=int(weighted.min())<=int(weighted.max())<=total*(1<<bits)<=LIMIT
        intercept=arb(int(weighted[1:].max())*unit)/scalar_denominator
        product=c*second
        return baseline.up((intercept+arb(product.numerator)/product.denominator)/(1<<rounds))

    result={}
    for tilt in tilts:
        overall=min(value(tilt,c,levels) for c in slopes(proposal,W,second,tilt))
        uniform={}
        for v in levels:
            # Every class may use the whole second-moment budget. These
            # separate upper bounds must not be combined as allocations.
            candidates=slopes({v:proposal[v]},W,second,tilt)
            uniform[v]=baseline.up(min(value(tilt,c,(v,)) for c in candidates)/sizes[v])
        result[tilt]=dict(density=overall,uniform=uniform,
                          rounding_factor=Q(total*unit,denominator),dyadic_bits=bits)
    return result


def build(maximum,tilts,second,*,minimum=5,precision=192):
    if (type(minimum) is not int or type(maximum) is not int
            or not 1<=minimum<=maximum<=10 or precision<128):
        raise ValueError('complete bounded shape range and outward precision required')
    ctx.prec=precision
    required={shape for j in range(minimum,maximum+1)
              for shape in combinations_with_replacement(range(1,5),j)}
    if not required<=set(second):raise ValueError('complete original-law second census required')
    shapes,coefficients,denominators,_,inverse=character_polynomials(maximum)
    images,_,spectrum=maps()
    expansion=np.array([x.bit_count() for x in images],dtype=np.int64)
    prepared=class_transforms(expansion);assert prepared[1]==dict(spectrum)
    tilts=tuple(dict.fromkeys(map(str,tilts)))
    result={tilt:{} for tilt in tilts};feedback={};inflation=Q(1)
    for shape,row in zip(shapes,coefficients):
        if sum(shape)<minimum:continue
        weights=tuple(b for b in range(1,5) for _ in range(shape[b-1]))
        D=int(denominators[shape])
        if D!=second[weights][0]:raise ValueError('inconsistent counting denominators')
        counts=inverse_counts(row[inverse],D)
        records=shape_bounds(counts,D,weights,expansion,prepared,second[weights][2],tilts)
        for tilt,record in records.items():
            result[tilt][weights]=record;inflation=max(inflation,record['rounding_factor'])
        classes={v:int(counts[expansion==v].sum()) for v in prepared[0]}
        feedback[weights]=int(counts[0]),int(counts[1:].max()),D,classes
        assert sum(classes.values())+int(counts[0])==D
        if len(feedback)%50==0:
            print('QUADRATIC DENSITY complete shapes',len(feedback),flush=True)
    assert set(feedback)==required and all(set(records)==required for records in result.values())
    print('QUADRATIC DENSITY census complete',len(feedback),'count inflation',inflation,flush=True)
    return result,feedback
