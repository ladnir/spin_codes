"""Exhaustive random-window averages for arbitrary states and density caps.

All 2^19 states and all 480 nonzero four-bit-window inputs are enumerated.
Exponential powers are rounded upward to dyadic integers; sums/maxima are
exact uint64 operations. No random sampling or floating maximum is used.
"""
from itertools import product
from math import comb
import numpy as np
from flint import arb,arb_mat

from group_moment import maps
from group_rank_one_verify import up
from occupancy_memory import epoch_operators,rounded,Z,F,M,C,U


def prepare_inputs():
    images,columns,spectrum=maps()
    mask=(1<<64)-1
    low=np.array([x&mask for x in images],dtype=np.uint64)
    high=np.array([x>>64 for x in images],dtype=np.uint64)
    inputs={b:[] for b in range(1,5)}
    for start in range(0,128,4):
        for bits in range(1,16):
            syndrome=0
            for bit in range(4):
                if bits>>bit&1:
                    syndrome^=columns[start+bit]
            word=bits<<start
            inputs[bits.bit_count()].append((syndrome,word&mask,word>>64))
    assert all(len(words)==32*comb(4,b) for b,words in inputs.items())
    return low,high,inputs,spectrum


def dyadic_powers(tilt,bits=44):
    scale=1<<bits
    values=[]
    for w in range(129):
        power=(-arb(tilt)*w).exp()
        integer=scale if w==0 else int((power*scale).upper().ceil().unique_fmpz())
        assert power<=arb(integer)/scale and integer<=scale
        values.append(integer)
    # At most 192 terms enter any row. The exact sum fits uint64.
    assert 192*scale<1<<64
    return np.array(values,dtype=np.uint64),scale


def averages(data,tilt,bits=44):
    low,high,inputs,spectrum=data
    powers,scale=dyadic_powers(tilt,bits)
    states=np.arange(len(low),dtype=np.uint64)
    result={}
    for b,words in inputs.items():
        mass=np.zeros(len(low),dtype=np.uint64)
        density=np.zeros(len(low),dtype=np.uint64)
        for syndrome,lo,hi in words:
            weights=np.bitwise_count(low^np.uint64(lo))+np.bitwise_count(high^np.uint64(hi))
            values=powers[weights]
            mass+=values
            density+=values[states^np.uint64(syndrome)]
        # Incoming and outgoing mature components exclude zero. The density
        # bound may conservatively include the source-zero terms; only the
        # outgoing zero coordinate uses its separate exact average.
        denominator=len(words)*scale
        result[b]=(int(mass[1:].max()),int(density[1:].max()),int(density[0]),denominator)
        assert all(0<v<=denominator for v in result[b][:3])
        # Independently replay both maximizing states and the zero row with
        # Python integers. This also checks the uint64 gather convention.
        winners=(int(mass[1:].argmax())+1,int(density[1:].argmax())+1,0)
        for kind,state in enumerate(winners):
            reference=0
            for syndrome,lo,hi in words:
                index=state if kind==0 else state^syndrome
                weight=(int(low[index])^lo).bit_count()+(int(high[index])^hi).bit_count()
                reference+=int(powers[weight])
            assert reference==result[b][kind]
        print('Exact dyadic window average',tilt,'weight',b,
              'mass/density/zero',[round(v/denominator,10) for v in result[b][:3]],flush=True)
    return result


def refine(prepared,values,tilt,groups,full_penalty=1,base=None,input_penalty=1,odd_penalty=1):
    """Refine only the one-input operator, before taking shape maxima."""
    if base is None:
        base=epoch_operators(prepared,tilt,maximum_groups=groups,pair_conditioned=True,full_penalty=full_penalty,input_penalty=input_penalty,odd_penalty=odd_penalty)
    exact=epoch_operators(prepared,tilt,detailed=True,full_penalty=full_penalty,input_penalty=input_penalty,odd_penalty=odd_penalty)
    spectrum=prepared[0][0][0]
    states=(1<<19)-1
    candidates=[]
    for shape,t in exact[1].items():
        b=shape[0]
        mass,density,zero,denominator=values[b]
        rho=arb(full_penalty)**int(b==4)*arb(input_penalty)**b*arb(odd_penalty)**(b%2)
        mass=up(arb(mass)/denominator*rho)
        density=up(arb(density)/denominator*rho)
        zero=up(arb(zero)/denominator*rho)
        t[M,M]=min(t[M,M],up(mass/2))
        t[M,Z]=min(t[M,Z],up(mass/(2*states)))
        t[C,Z]=min(t[C,Z],up(zero/2))
        t[C,C]=min(t[C,C],up(density/2))
        for k,v in enumerate(sorted(spectrum)):
            t[M,U+k]=min(t[M,U+k],up(mass*spectrum[v]/(2*states)))
            t[U+k,C]=min(t[U+k,C],up(density/(2*spectrum[v])))
        candidates.append(t)
    # Both base and the newly aggregated matrix are coordinatewise valid
    # for the same invariant. Their entrywise minimum remains valid here:
    # each affected entry bounds the same component contribution.
    for i,j in product(range(9),repeat=2):
        base[1][i,j]=min(base[1][i,j],max(t[i,j] for t in candidates))
    return base


def self_test(data,values,tilt):
    low,high,inputs,spectrum=data
    # Independent Python integers replay selected states. These include zero
    # for the separate cancellation average and endpoints for the two maxima.
    for b,words in inputs.items():
        mass,density,zero,denominator=values[b]
        for state in (0,1,17,131071,262144,524287):
            sums=[arb(0),arb(0)]
            for syndrome,lo,hi in words:
                for k,index in enumerate((state,state^syndrome)):
                    weight=(int(low[index])^lo).bit_count()+(int(high[index])^hi).bit_count()
                    sums[k]+=(-arb(tilt)*weight).exp()/len(words)
            if state:
                assert sums[0]<=arb(mass)/denominator
                assert sums[1]<=arb(density)/denominator
            else:
                assert sums[1]<=arb(zero)/denominator
    print('Independent arbitrary-state and cancellation-average spot checks passed',flush=True)
