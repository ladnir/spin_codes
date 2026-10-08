"""Separate directed K16/K18 replays for the unchanged GRS128/87 ensemble.

The K20 files are imported without modifying their constants. Geometry is
explicit. Middle conditioning uses the exact central binomial mode bound;
dense membership is rechecked at exponent1.30 with the frozen rational tilts.
No optimizer, native encoder, benchmark, or saved numerical receipt is used.
"""
from __future__ import annotations

import argparse
from decimal import Context, Decimal, ROUND_CEILING
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb

from flint import arb, ctx

from field128_dense_certificate import TILTS, interval_base
from field128_density_certificate import (WITNESSES, cell_upper, integer_span,
    retarget_witness, transfer, weight_geometry)
from field128_sparse_certificate import certify_zone, complete_cover
from field256_whole_certificate import margin_bits
from joint_outer_screen import WITNESSES_128_87
from return_cover_certificate import interval, terminal_power
from sparse_active_interval_screen import anchor_zones, interpolated_anchors
from state8_density_certificate import refined_cells
from state8_sparse_certificate import integer_bits
from subfield_certificate import ZERO, integer_margin


MIDDLE_FIRST = {16:F(21,400),18:F(3,100)}
MIDDLE_LAST = F(6,25)


def geometry(k_log):
    if k_log not in MIDDLE_FIRST:
        raise ValueError('this checkpoint covers only K16 and K18')
    g = weight_geometry(k_log)
    g['middle_first'],g['middle_last'] = integer_span(g['n'],MIDDLE_FIRST[k_log],MIDDLE_LAST)
    g['sparse_last'] = (g['middle_first']-1)//42
    g['dense_first'] = g['middle_last']+1
    return g


def component_metadata(k_log):
    """Identify the common ensemble and all finite geometry used by a receipt."""
    g = geometry(k_log)
    return dict(k_log=k_log,outer=[128,87],q=256,width=16,state=8,
                construction='symmetric-cauchy-w16-s8-eta2-regional-nonzero-diagonals',
                requested_k=2**k_log,padded_k=g['k'],rows=g['rows'],n=g['n'],
                output_threshold=g['threshold'])


def modal_reciprocal(rows):
    """Exact uniform conditioning penalty for all Bin(L,h/L) modal counts.

    For even L, min_h C(L,h)(h/L)^h(1-h/L)^(L-h)=C(L,L/2)/2^L.
    This is NOT a bound for arbitrary nonmodal h under a fixed parameter.
    """
    if rows<=0 or rows%2:
        raise ValueError('positive even regional length required')
    return F(2**rows,comb(rows,rows//2))


def sparse_certificate(k_log, precision=4096):
    g = geometry(k_log)
    anchors = interpolated_anchors(interpolated_anchors(WITNESSES_128_87))
    zones = anchor_zones(1,g['sparse_last'],anchors)
    previous = ctx.prec
    try:
        ctx.prec = 256
        entries,total = [],arb(0)
        for lo,hi,anchor in zones:
            upper,receipt = certify_zone(lo,hi,anchor,anchors[anchor],g,
                                         cell_bits=64,precision=precision)
            total += upper
            entries.extend(receipt['entries'])
        complete_cover(entries,1,g['sparse_last'])
        upper = total.upper()
        rational = upper.fmpq()
        bits = integer_bits(upper)
        if bits<40:
            raise ArithmeticError('complete sparse class is below forty bits')
        digest = sha256(''.join(f"{e['lower']}:{e['upper']}:{e['anchor']}\n"
                               for e in entries).encode('ascii')).hexdigest()
        return dict(**component_metadata(k_log),active_range=[1,g['sparse_last']],
                    interval_count=len(entries),partition_sha256=digest,
                    failure_upper_dyadic=[str(rational.p),str(rational.q)],
                    margin_bits_lower=bits,complete=True,rigorous=True)
    finally:
        ctx.prec = previous


def middle_prefactor(witness,g):
    return (interval(modal_reciprocal(g['rows']))**128
            *terminal_power(transfer(witness),g['n']//16)
            /(interval(witness.tau)*interval(witness.y)**g['threshold']))


def middle_certificate(k_log, initial_grid=1000, cell_bits=64):
    if initial_grid<1 or cell_bits<40:
        raise ValueError('invalid directed middle replay parameters')
    g = geometry(k_log)
    prefactors = tuple(middle_prefactor(w,g) for _,w in WITNESSES)

    def evaluate(lo,hi):
        midpoint = F(lo+hi,2*g['n'])
        options = [(i,cell_upper(lo,hi,g['n'],g['rows'],
                    retarget_witness(anchor,w,midpoint),prefactors[i]))
                   for i,(anchor,w) in enumerate(WITNESSES)]
        return min(options,key=lambda item:item[1].hi)

    span = (g['n']+initial_grid-1)//initial_grid
    threshold = interval(F(1,2**cell_bits)).lo
    next_weight,total,cells = g['middle_first'],ZERO,0
    digest = sha256()
    for lo,hi,index,upper in refined_cells(g['middle_first'],g['middle_last'],span,evaluate,threshold):
        if lo!=next_weight or hi<lo:
            raise ArithmeticError('middle cover has a gap or overlap')
        next_weight=hi+1
        total+=upper
        cells+=1
        digest.update(f'{lo}:{hi}:{index}\n'.encode('ascii'))
    if next_weight!=g['middle_last']+1:
        raise ArithmeticError('middle cover has the wrong endpoint')
    bits=integer_margin(total.hi)
    if bits<40:
        raise ArithmeticError('complete middle class is below forty bits')
    return dict(**component_metadata(k_log),weights=[g['middle_first'],g['middle_last']],
                covered_weights=g['middle_last']-g['middle_first']+1,
                checked_cells=cells,partition_sha256=digest.hexdigest(),
                failure_upper=str(total.hi),margin_bits_lower=bits,
                conditioning='exact reciprocal of central binomial mode',
                complete=True,rigorous=True)


def dense_interval_passes(index,tilt):
    # exponent13/10: raising by GRID200 requires2^(260*128).
    if not 48<=index<200 or tilt<1:
        raise ValueError('invalid dense rational witness')
    return interval_base(128,87,index,tilt)**200*2**(260*128)<=tilt**(128*index)


def dense_certificate(k_log):
    g=geometry(k_log)
    if len(TILTS)!=152 or not all(dense_interval_passes(i,t) for i,t in enumerate(TILTS,48)):
        raise ArithmeticError('incomplete or failed dense membership witnesses')
    if not F(10,9)**100*2295**10<2**127:
        raise ArithmeticError('failed radius-.1 output ball bound')
    # Retain the old (L+1)^128 conditioning loss; the stronger membership
    # exponent alone gives .03N before that loss. No central lemma used here.
    bits=(3*g['n'])//100-128*g['rows'].bit_length()
    if not g['rows']+1<=2**g['rows'].bit_length() or bits<40:
        raise ArithmeticError('dense finite bound is below forty bits')
    return dict(**component_metadata(k_log),weights=[g['dense_first'],g['n']],
                membership_exponent='13/10',ball_exponent='127/100',
                rational_cells=152,margin_bits_lower=bits,
                failure_upper_dyadic_exponent=bits,complete=True,rigorous=True)


def combine(k_log,sparse,middle,dense):
    g=geometry(k_log)
    metadata=component_metadata(k_log)
    if any(any(r.get(key)!=value for key,value in metadata.items())
           or not r.get('complete') or not r.get('rigorous')
           for r in (sparse,middle,dense)):
        raise ValueError('component size or rigorous coverage mismatch')
    if (sparse['active_range']!=[1,g['sparse_last']]
        or middle['weights']!=[g['middle_first'],g['middle_last']]
        or dense['weights']!=[g['dense_first'],g['n']]
        or 42*(g['sparse_last']+1)<g['middle_first']
        or middle['weights'][1]+1!=dense['weights'][0]):
        raise ValueError('component coverage has a gap')
    p,q=map(int,sparse['failure_upper_dyadic'])
    bounds=(F(p,q),F(middle['failure_upper']),F(1,2**dense['failure_upper_dyadic_exponent']))
    if any(not 0<b<1 for b in bounds):
        raise ArithmeticError('invalid component failure bound')
    total=sum(bounds,F(0))
    bits=margin_bits(total)
    if bits<40:
        raise ArithmeticError('whole-code sum is below forty bits')
    display=Context(prec=80,rounding=ROUND_CEILING).divide(Decimal(total.numerator),Decimal(total.denominator))
    return dict(**metadata,minimum_distance=g['threshold']+1,
                design_rate='87/128',effective_requested_rate=str(F(2**k_log,g['n'])),
                relative_distance='1/10',whole_code_margin_bits_lower=bits,
                failure_upper=str(display),component_margin_bits_lower=[margin_bits(b) for b in bounds],
                sparse=sparse,middle=middle,dense=dense,rigorous=True,whole_code=True,
                scope='ideal independent regional permutations and GF256 nonzero diagonals')


def certify(k_log,progress=None):
    sparse=sparse_certificate(k_log)
    if progress:progress('sparse',sparse)
    middle=middle_certificate(k_log)
    if progress:progress('middle',middle)
    dense=dense_certificate(k_log)
    if progress:progress('dense',dense)
    return combine(k_log,sparse,middle,dense)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--k-log',type=int,nargs='+',default=(16,18))
    parser.add_argument('--component',choices=('all','sparse','middle','dense'),default='all')
    args=parser.parse_args()
    for size in args.k_log:
        if args.component=='all':
            result=certify(size,lambda c,r:print(json.dumps(dict(component=c,**r)),flush=True))
        else:
            result={'sparse':sparse_certificate,'middle':middle_certificate,'dense':dense_certificate}[args.component](size)
        print(json.dumps(result),flush=True)
