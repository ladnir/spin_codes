"""Exact distinct-window output moments using the expansion's 20 histograms.

Integer coefficient extraction and upward dyadic factors provide valid
moment uppers. This does not replace the feedback-cancellation bounds.
"""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations,permutations,product
from math import comb,factorial,prod
import argparse
import numpy as np
from flint import arb,ctx
from group_moment import maps
from group_rank_one_verify import up
from occupancy_memory import Z,F,M,U


def coefficients(hist,factors,maximum):
    """Coefficients of product_i(1+sum_b factors[r_i,b]*x_b)."""
    values={(0,0,0,0):1}
    for r,number in enumerate(hist):
        for _ in range(number):
            previous=list(values.items())
            for counts,value in previous:
                if sum(counts)==maximum:continue
                for b,factor in enumerate(factors[r]):
                    new=list(counts);new[b]+=1;key=tuple(new)
                    values[key]=values.get(key,0)+value*factor
    return values


def denominator(windows,counts):
    j=sum(counts)
    return comb(windows,j)*factorial(j)//prod(factorial(n) for n in counts)


def lane_factor(r,b,z):
    return sum((comb(r,h)*comb(4-r,b-h)*z**(b-2*h)
                for h in range(max(0,r+b-4),min(r,b)+1)),0)/comb(4,b)


def census(prepared=None):
    images,columns,spectrum=maps()
    mask=(1<<64)-1
    low=np.array([x&mask for x in images],dtype=np.uint64)
    high=np.array([x>>64 for x in images],dtype=np.uint64)
    hist=np.zeros((len(images),5),dtype=np.uint8);indices=np.arange(len(images))
    for part in (low,high):
        for shift in range(0,64,4):
            weights=np.bitwise_count((part>>np.uint64(shift))&np.uint64(15))
            np.add.at(hist,(indices,weights),1)
    shapes,inverse,numbers=np.unique(hist,axis=0,return_inverse=True,return_counts=True)
    shapes=[tuple(map(int,row)) for row in shapes]
    assert len(shapes)==20 and sum(map(int,numbers))==1<<19
    sizes=Counter()
    for row,count in zip(shapes,numbers):sizes[sum(r*n for r,n in enumerate(row))]+=int(count)
    assert sizes==Counter({0:1,**spectrum})
    fresh=[]
    # Each fresh distribution is the feedback of one nonzero window.
    for b in range(1,5):
        atoms=Counter()
        for window in range(32):
            for bits in range(1,16):
                if bits.bit_count()!=b:continue
                syndrome=0
                for bit in range(4):
                    if bits>>bit&1:syndrome^=columns[4*window+bit]
                atoms[syndrome]+=1
        assert not atoms[0] and sum(atoms.values())==32*comb(4,b)
        if prepared is not None:assert atoms==prepared[0][0][3][(b,)]
        fresh.append(Counter({h:sum(n for s,n in atoms.items() if int(inverse[s])==h)
                              for h in range(len(shapes))}))
    for state in (0,1,17,(1<<19)-1):
        exact=Counter(((images[state]>>(4*w))&15).bit_count() for w in range(32))
        assert shapes[int(inverse[state])]==tuple(exact[r] for r in range(5))
    print('Expansion window census: 20 exact histograms, all states and fresh distributions checked',flush=True)
    return shapes,list(map(int,numbers)),fresh,spectrum


def moments(data,tilt,maximum=8,bits=48):
    shapes,numbers,fresh,spectrum=data
    z=(-arb(tilt)).exp();scale=1<<bits
    factors=[[int((lane_factor(r,b,z)*scale).upper().ceil().unique_fmpz())
              for b in range(1,5)] for r in range(5)]
    assert all(lane_factor(r,b,z)<=arb(factors[r][b-1])/scale for r in range(5) for b in range(1,5))
    result={}
    for h,hist in enumerate(shapes):
        weight=sum(r*n for r,n in enumerate(hist))
        if not weight:continue
        for counts,value in coefficients(hist,factors,maximum).items():
            j=sum(counts)
            bound=up(z**weight*value/(denominator(32,counts)*scale**j))
            result.setdefault(counts,{})[h]=min(arb(1),bound)
    print('Outward distinct-window moments:',len(result),'shapes, degree',maximum,'dyadic bits',bits,flush=True)
    return result


def refine(base,data,moments,penalty,rounds=2,input_penalty=1,odd_penalty=1):
    histograms,numbers,fresh,spectrum=data
    levels=sorted(spectrum);m=(1<<19)-1
    alpha=arb(2)**(-rounds);beta=1-alpha
    classes={v:[h for h,row in enumerate(histograms) if sum(r*n for r,n in enumerate(row))==v] for v in levels}
    maxima={}
    for counts,by_hist in moments.items():
        j=sum(counts)
        if j==0 or j>=len(base):continue
        scale=arb(penalty)**counts[3]*arb(input_penalty)**sum((b+1)*n for b,n in enumerate(counts))*arb(odd_penalty)**(counts[0]+counts[2])
        arbitrary=max(by_hist.values())
        fresh_bound=max(up(sum((n*by_hist[h] for h,n in dist.items() if n),arb(0))/sum(dist.values())) for dist in fresh)
        uniform={v:up(sum((numbers[h]*by_hist[h] for h in classes[v]),arb(0))/spectrum[v]) for v in levels}
        entries={(M,M):alpha*arbitrary,(M,Z):beta*arbitrary/m,(F,M):alpha*fresh_bound}
        for k,v in enumerate(levels):
            entries[M,U+k]=beta*arbitrary*spectrum[v]/m
            entries[F,U+k]=beta*fresh_bound*spectrum[v]/m
            entries[U+k,M]=alpha*uniform[v]
            for ell,w in enumerate(levels):entries[U+k,U+ell]=beta*uniform[v]*spectrum[w]/m
        for (source,target),value in entries.items():
            key=j,source,target
            maxima[key]=max(maxima.get(key,arb(0)),up(value*scale))
    changed=0
    for (j,source,target),value in maxima.items():
        if value<base[j][source,target]:changed+=1
        base[j][source,target]=min(base[j][source,target],value)
    print('Distinct-window moment refinement:',changed,'local coefficients tightened',flush=True)
    return base


def self_test():
    checks=0
    for pop in ((0,4),(1,3),(0,2,4),(1,2,3),(2,2,2)):
        hist=[pop.count(r) for r in range(5)]
        z=Q(7,8);factors=[[lane_factor(r,b,z) for b in range(1,5)] for r in range(5)]
        poly=coefficients(hist,factors,len(pop))
        for counts,coefficient in poly.items():
            weights=[b for b,n in enumerate(counts,1) for _ in range(n)]
            j=len(weights)
            actual=sum((prod(factors[pop[i]][b-1] for i,b in zip(choice,weights))
                        for choice in permutations(range(len(pop)),j)),Q(0))
            actual/=factorial(len(pop))//factorial(len(pop)-j)
            assert coefficient/denominator(len(pop),counts)==actual
            checks+=1
        for r in range(5):
            base=(1<<r)-1
            for b in range(1,5):
                masks=[x for x in range(16) if x.bit_count()==b]
                direct=sum((z**((base^x).bit_count()) for x in masks),Q(0))/len(masks)
                assert z**r*factors[r][b-1]==direct
                checks+=1
    print('Distinct-window coefficient and lane-mask tests:',checks,'exact checks passed',flush=True)


def production_test(data,bounds,tilt):
    """Direct mask enumeration, independent of histogram coefficients."""
    images,_,_=maps();checks=0
    for state in (1,17,1<<18):
        image=images[state]
        histogram=tuple(sum(((image>>(4*i))&15).bit_count()==r for i in range(32)) for r in range(5))
        h=data[0].index(histogram)
        for weights in ((1,),(4,),(1,2),(2,2),(4,4)):
            counts=tuple(weights.count(b) for b in range(1,5))
            if counts not in bounds:continue
            actual=Counter()
            for windows in combinations(range(32),len(weights)):
                for assignment in set(permutations(weights)):
                    choices=[[x for x in range(1,16) if x.bit_count()==b] for b in assignment]
                    for masks in product(*choices):
                        word=sum(mask<<(4*i) for i,mask in zip(windows,masks))
                        actual[(image^word).bit_count()]+=1
            expected=denominator(32,counts)*prod(comb(4,b) for b in weights)
            assert sum(actual.values())==expected
            saved=ctx.prec
            try:
                ctx.prec=saved+128
                reference=sum((n*(-arb(tilt)*w).exp() for w,n in actual.items()),arb(0))/expected
            finally:ctx.prec=saved
            assert reference<=bounds[counts][h]
            checks+=1
    print('Production direct-mask moments:',checks,'independent outward inequalities passed',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--maximum',type=int,choices=range(1,13),default=8)
    parser.add_argument('--tilt',default='.032')
    parser.add_argument('--precision',type=int,default=192)
    args=parser.parse_args();self_test();ctx.prec=args.precision
    data=census();result=moments(data,args.tilt,args.maximum)
    production_test(data,result,args.tilt)
    for counts in ((2,0,0,0),(0,2,0,0),(0,0,0,2),(1,1,1,0)):
        if counts in result:print(counts,'max moment',max(result[counts].values()),flush=True)
