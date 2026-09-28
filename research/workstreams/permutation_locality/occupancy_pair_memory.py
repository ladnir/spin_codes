"""Retain the nonzero part of two independent feedback inputs as a cone.

Adds one coefficient H to the nine-coordinate envelope. H is a mixture
coefficient for pair-convolution laws restricted to nonzero states, not a
claim that the conditioned encoder state is distributed uniformly.
"""
from collections import Counter
from itertools import combinations_with_replacement
from fractions import Fraction
from math import comb
import numpy as np
from flint import arb,arb_mat

from group_moment import maps
from group_rank_one_verify import up
from occupancy_memory import epoch_operators,rounded,shape_classes,Z,F,M,C,U,TERMINAL

H=9
PAIR_TERMINAL=np.append(TERMINAL,1.)


def census(prepared):
    images,columns,spectrum=maps()
    atoms=prepared[0][0][3]
    inputs={b:[] for b in range(1,5)}
    for start in range(0,128,4):
        for mask in range(1,16):
            syndrome=0
            for bit in range(4):
                if mask>>bit&1:
                    syndrome^=columns[start+bit]
            inputs[mask.bit_count()].append((syndrome,mask<<start))
    records=[]
    for a,b in combinations_with_replacement(range(1,5),2):
        convolution=Counter()
        for x,cx in atoms[a,].items():
            for y,cy in atoms[b,].items():
                convolution[x^y]+=cx*cy
        denominator=sum(atoms[a,].values())*sum(atoms[b,].values())
        assert sum(convolution.values())==denominator
        windows=Counter()
        levels=Counter()
        for state,count in convolution.items():
            if not state or not count:
                continue
            image=images[state]
            weight=image.bit_count()
            levels[weight]+=count
            for start in range(0,128,4):
                windows[weight,((image>>start)&15).bit_count()]+=count
        moments={}
        cancellations={}
        for w,words in inputs.items():
            moment=Counter()
            cancel=Counter()
            for (v,local),count in windows.items():
                for overlap in range(max(0,w+local-4),min(w,local)+1):
                    moment[v+w-2*overlap]+=count*comb(local,overlap)*comb(4-local,w-overlap)
            for state,word in words:
                cancel[(images[state]^word).bit_count()]+=convolution[state]
            assert sum(moment.values())==(denominator-convolution[0])*len(words)
            assert sum(cancel.values())<=sum(moment.values())
            moments[w]=moment
            cancellations[w]=cancel
        density=Fraction(max(count for state,count in convolution.items() if state),denominator)
        records.append(((a,b),denominator,levels,moments,cancellations,density))
    print('Pair-memory census: ten exact convolution laws, moments and cancellation histograms; max nonzero atom',
          max(r[-1] for r in records),flush=True)
    return spectrum,records


def operators(prepared,data,tilt,groups,full_penalty=1):
    spectrum,records=data
    base=epoch_operators(prepared,tilt,maximum_groups=groups,pair_conditioned=True,full_penalty=full_penalty)
    powers=[up((-arb(tilt)*w).exp()) for w in range(129)]
    rho=arb(full_penalty)
    states=(1<<19)-1
    density=max(arb(r[-1].numerator)/r[-1].denominator for r in records)
    def average(hist,denominator,weight=0):
        return up(sum((count*powers[max(0,v-weight)] for v,count in hist.items()),arb(0))/denominator)
    triangle={w:max(average(r[2],r[1],w) for r in records) for w in range(4*groups+1)}
    result=[]
    for j,old in enumerate(base):
        t=arb_mat([[old[a,b] if a<9 and b<9 else arb(0) for b in range(10)] for a in range(10)])
        if j==1:
            # The pointwise lazy multiplier, not its average, preserves
            # domination by a mixture of pair-convolution distributions.
            t[F,H]=t[F,M]
            t[F,M]=arb(0)
        candidates=[]
        shapes=[(0,0,0)] if j==0 else sorted(shape_classes(j,True))
        for weight,choices,full in shapes:
            f=powers[max(0,48-weight)]
            row=[arb(0)]*10
            if j==0:
                row[H]=powers[48]/2
                moment=triangle[0]
            elif j==1:
                denominator=32*comb(4,weight)
                moment=max(average(r[3][weight],r[1]*denominator) for r in records)
                cancellation=max(average(r[4][weight],r[1]*denominator) for r in records)
                row[Z]=cancellation/2+moment/(2*states)
                row[M]=moment/2
                row[C]=f*density/2
            else:
                moment=triangle[weight]
                atom=min(density,arb(1)/((33-j)*choices))
                row[Z]=f*atom/2+moment/(2*states)
                row[M]=moment/2
                row[C]=f*atom/2
            for k,v in enumerate(sorted(spectrum)):
                row[U+k]=moment*spectrum[v]/(2*states)
            candidates.append([up(value*rho**full) for value in row])
        for k in range(10):
            t[H,k]=max(row[k] for row in candidates)
        result.append(rounded(t))
    return result


def self_test(prepared,data):
    spectrum,records=data
    for tilt in ('0','.0064'):
        for penalty in ('1','.5'):
            ops=operators(prepared,data,tilt,4,penalty)
            assert ops[1][F,M]==0 and ops[1][F,H]>0
            for _,denominator,levels,moments,cancellations,_ in records:
                def moment(hist,den):
                    return sum((count*(-arb(tilt)*w).exp() for w,count in hist.items()),arb(0))/den
                empty=moment(levels,denominator)
                for k,v in enumerate(sorted(spectrum)):
                    assert empty*spectrum[v]/(2*((1<<19)-1))<=ops[0][H,U+k]
                for b in range(1,5):
                    scale=arb(penalty)**int(b==4)
                    den=denominator*32*comb(4,b)
                    mass=moment(moments[b],den)*scale
                    zero=moment(cancellations[b],den)*scale
                    assert mass/2<=ops[1][H,M]
                    assert zero/2+mass/(2*((1<<19)-1))<=ops[1][H,Z]
                    for k,v in enumerate(sorted(spectrum)):
                        assert mass*spectrum[v]/(2*((1<<19)-1))<=ops[1][H,U+k]
    print('Pair-memory actual-map average/zero/refresh inequalities passed for both penalties',flush=True)
