"""Exact feedback convolutions for lazy-state density bounds.

Fix a packet-weight shape, with total input weight W. Let n[r] count its
inputs with feedback r, and let D=sum(n). For a nonzero source state s,
the output E(s)+x has weight at least |wt(E(s))-W|. Therefore the lazy
contribution to the maximum density at a nonzero target t is at most

    alpha * sum_{s != 0} n[t+s] exp(-lambda*|wt(E(s))-W|) / D,

times the incoming maximum density. Here alpha=2**(-rounds); input and
source are independent before the output tilt. No independence is assumed
after that tilt. Uniform-class mass U_v assigns density U_v/A_v to each
of the A_v states with expansion weight v. Its analogous bound uses only
that class and divides by A_v.

For each nonzero expansion class, compute exact XOR convolutions

    G_v(t) = sum_s n[t+s] 1_{wt(E(s))=v}
           = FWHT(FWHT(n)*FWHT(1_v))[t] / N,  N=2**19.

The five class transforms are shared across shapes and tilts. Full
feedback character polynomials already provide FWHT(n). The convolution
uses signed int64 arithmetic only after proving N*D < 2**63.

Integer-width argument: a partial inverse Walsh transform sums over a
character subspace V. Character orthogonality restricts each pair of
primal indices to one coset and supplies a factor |V|. For each first
index, at most N/|V| second indices satisfy that restriction. Because
n is nonnegative with total D and the class indicator is at most one,
every partial butterfly has absolute value at most N*D. This includes
the initial pointwise products. A naive N**2*D bound is unnecessary.

For the final maximum, round each exponential upward to Q_v/2**b. Since
sum_v G_v(t)=D-n[t] <= D and 0 <= Q_v <= 2**b, every partial weighted sum
is at most D*2**b. Choosing that product below 2**63 makes the maximum
exact too. Only the final rational coefficient is rounded outward in Arb.

This is an optional local refinement, not a whole-code certificate.
No existing operator is changed by this module.
"""
import argparse
from fractions import Fraction as Q
from math import comb, factorial, prod
from pathlib import Path
import sys

import numpy as np
from flint import arb, ctx

sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from pair_tail import hadamard


INT64_MAX = (1 << 63)-1


def _integer_array(values):
    """Reject truncating casts before entering exact integer arithmetic."""
    values = np.asarray(values)
    if values.dtype.kind not in 'iu':
        raise ValueError('exact integer array required; floating casts are not supported')
    if values.dtype.kind == 'u' and values.size and int(values.max()) > INT64_MAX:
        raise ValueError('unsigned integer input exceeds the signed int64 range')
    return values.astype(np.int64,copy=False)


def _rational(value):
    value = Q(value)
    return arb(value.numerator)/value.denominator


def _up(value):
    if not value.is_finite():
        raise ArithmeticError('finite outward density coefficient required')
    return arb(value.upper())


def _dimension(size):
    if not isinstance(size,int) or size < 2 or size & (size-1):
        raise ValueError('a power-of-two state count of at least two is required')


def recover_counts(transformed, denominator):
    """Recover and authenticate a nonnegative integer feedback histogram."""
    transformed = _integer_array(transformed)
    if transformed.ndim != 1:
        raise ValueError('one-dimensional Walsh coefficients required')
    size = len(transformed)
    _dimension(size)
    if not isinstance(denominator,int) or denominator < 1 or size*denominator > INT64_MAX:
        raise ValueError('the exact N*D int64 guard failed')
    if int(transformed.min()) < -denominator or int(transformed.max()) > denominator:
        raise ValueError('feedback character exceeds its exact counting denominator')
    # The elementary l1 bound N*D suffices for this first inverse transform.
    counts = hadamard(transformed)
    if np.any(counts % size):
        raise ValueError('feedback inverse transform is not integral')
    counts //= size
    if int(counts.min()) < 0 or int(counts.sum()) != denominator:
        raise ValueError('feedback histogram is not a nonnegative counting measure')
    return counts


def class_transforms(expansion):
    """Exact transforms of every nonzero expansion class."""
    expansion = _integer_array(expansion)
    if expansion.ndim != 1:
        raise ValueError('one-dimensional expansion weights required')
    _dimension(len(expansion))
    if expansion[0] != 0 or np.any(expansion[1:] <= 0):
        raise ValueError('the expansion map must have only one zero word')
    levels = tuple(int(v) for v in np.unique(expansion[1:]))
    sizes = {v:int(np.count_nonzero(expansion == v)) for v in levels}
    transforms = {v:hadamard((expansion == v).astype(np.int64)) for v in levels}
    assert sum(sizes.values()) == len(expansion)-1
    return levels,sizes,transforms


def class_convolution(transformed, counts, expansion, indicator, denominator, level):
    """One exact translated-class count, including independent spot checks.

    transformed/counts must come from recover_counts; indicator must come
    from class_transforms for this expansion map and level.
    """
    size = len(counts)
    if size*denominator > INT64_MAX:
        raise ValueError('the exact N*D convolution guard failed')
    hits = hadamard(transformed*indicator)
    assert not np.any(hits % size)
    hits //= size
    assert 0 <= int(hits.min()) <= int(hits.max()) <= denominator
    cardinality = int(np.count_nonzero(expansion == level))
    assert int(hits.sum()) == denominator*cardinality
    indices = np.arange(size,dtype=np.int64)
    for target in sorted({0,1,min(17,size-1),size-1,int(hits.argmax())}):
        direct = int(counts[expansion[indices ^ target] == level].sum())
        assert int(hits[target]) == direct
    return hits


def dyadic_precision(denominator, maximum=44):
    if (not isinstance(denominator,int) or denominator < 1 or denominator > INT64_MAX
            or not isinstance(maximum,int) or maximum < 1):
        raise ValueError('positive denominator and maximum precision required')
    bits = min(maximum,(INT64_MAX//denominator).bit_length()-1)
    if bits < 1:
        raise ValueError('insufficient int64 headroom for a dyadic coefficient')
    assert denominator*(1 << bits) <= INT64_MAX
    return bits


def dyadic_powers(levels, input_weight, tilt, bits):
    """Certified ceil(2**bits exp(-lambda |v-W|)), capped by 2**bits."""
    tilt = Q(tilt)
    if tilt <= 0 or not isinstance(input_weight,int) or input_weight < 0:
        raise ValueError('positive tilt and nonnegative integer input weight required')
    if not isinstance(bits,int) or not 1 <= bits <= 62:
        raise ValueError('dyadic precision must lie between 1 and 62')
    scale, result = 1 << bits,{}
    for v in levels:
        upper = _up((-_rational(tilt)*abs(v-input_weight)).exp()*scale).fmpq()
        numerator,denominator = int(upper.numerator),int(upper.denominator)
        result[v] = min(scale,(numerator+denominator-1)//denominator)
        assert 0 < result[v] <= scale
    return result


def shape_bounds(transformed, denominator, weights, expansion, prepared,
                 tilts, *, rounds=2, maximum_bits=44):
    """Local C->C and U_v->C coefficients for one exact feedback law."""
    if not isinstance(rounds,int) or rounds < 1:
        raise ValueError('a positive integer transvection count is required')
    weights = tuple(weights)
    if not weights or any(not isinstance(w,int) or not 1 <= w <= 4 for w in weights):
        raise ValueError('nonempty four-bit packet weights required')
    expansion = _integer_array(expansion)
    transformed = _integer_array(transformed)
    if transformed.shape != expansion.shape:
        raise ValueError('feedback and expansion state dimensions differ')
    levels,sizes,indicators = prepared
    counts = recover_counts(transformed,denominator)
    bits = dyadic_precision(denominator,maximum_bits)
    tilts = tuple(dict.fromkeys(str(tilt) for tilt in tilts))
    if not tilts:
        raise ValueError('at least one positive tilt is required')
    powers = {tilt:dyadic_powers(levels,sum(weights),tilt,bits) for tilt in tilts}
    weighted = {tilt:np.zeros(len(counts),dtype=np.int64) for tilt in tilts}
    uniforms = {tilt:{} for tilt in tilts}
    total = np.zeros(len(counts),dtype=np.int64)
    scalar_denominator = denominator*(1 << bits)*(1 << rounds)
    for v in levels:
        hits = class_convolution(transformed,counts,expansion,indicators[v],denominator,v)
        total += hits
        peak = int(hits[1:].max())
        for tilt in tilts:
            coefficient = powers[tilt][v]
            weighted[tilt] += hits*coefficient
            uniforms[tilt][v] = _up(arb(peak*coefficient)/(scalar_denominator*sizes[v]))
    assert np.array_equal(total,denominator-counts)
    result = {}
    for tilt in tilts:
        assert int(weighted[tilt].min()) >= 0
        assert int(weighted[tilt].max()) <= denominator*(1 << bits)
        peak = int(weighted[tilt][1:].max())
        result[tilt] = {
            'density':_up(arb(peak)/scalar_denominator),
            'uniform':uniforms[tilt],
            'denominator':denominator,
            'dyadic_bits':bits,
        }
    return result


def build(maximum=4, tilts=('.052',), precision=192, rounds=2):
    """Build exact-census refinements; callers integrate them explicitly."""
    if not isinstance(maximum,int) or not 1 <= maximum <= 6:
        raise ValueError('this bounded int64 implementation supports one through six windows')
    if not isinstance(precision,int) or precision < 128:
        raise ValueError('at least 128 bits of outward precision required')
    from feedback_character_census import character_polynomials
    from group_moment import maps
    shapes,coefficients,denominators,_,inverse = character_polynomials(maximum)
    images,_,spectrum = maps()
    ctx.prec = precision
    expansion = np.array([image.bit_count() for image in images],dtype=np.int64)
    prepared = class_transforms(expansion)
    assert prepared[1] == dict(spectrum)
    tilts = tuple(dict.fromkeys(str(tilt) for tilt in tilts))
    result = {tilt:{} for tilt in tilts}
    maximum_product,checked = 0,0
    for shape,row in zip(shapes,coefficients):
        if not sum(shape):
            continue
        weights = tuple(b for b in range(1,5) for _ in range(shape[b-1]))
        denominator = denominators[shape]
        expected = comb(32,len(weights))*factorial(len(weights))//prod(factorial(n) for n in shape)
        expected *= prod(comb(4,b)**shape[b-1] for b in range(1,5))
        assert denominator == expected
        maximum_product = max(maximum_product,len(expansion)*denominator)
        local = shape_bounds(row[inverse],denominator,weights,expansion,prepared,tilts,rounds=rounds)
        for tilt,record in local.items():
            result[tilt][weights] = record
        checked += 1
        print('Exact feedback density shape/denominator/bits:',weights,denominator,
              next(iter(local.values()))['dyadic_bits'],flush=True)
    return {'bounds':result,'provenance':{
        'map':'Map128S19; maps() authenticates expansion and feedback columns',
        'feedback':'exact integer character_polynomials and nonnegative inverse Walsh counts',
        'method':'exact nonzero expansion-class XOR convolutions; dyadic upward output factors',
        'maximum':maximum,'tilts':tilts,'precision':precision,'rounds':rounds,
        'state_count':len(expansion),'spectrum':dict(spectrum),
        'checked_shapes':checked,'maximum_N_times_D':maximum_product,
        'integer_limit':INT64_MAX,
    }}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--maximum',type=int,choices=range(1,7),default=4)
    parser.add_argument('--tilts',nargs='+',default=['.052'])
    parser.add_argument('--precision',type=int,default=192)
    parser.add_argument('--rounds',type=int,default=2)
    args = parser.parse_args()
    result = build(args.maximum,args.tilts,args.precision,args.rounds)
    print('PROVENANCE',result['provenance'],flush=True)
    for tilt,records in result['bounds'].items():
        for weights,record in records.items():
            print('DENSITY',tilt,weights,record,flush=True)
