"""Exact output moments for mature lazy cancellation with 2/3 windows."""
from collections import Counter
from itertools import combinations,combinations_with_replacement,product
from math import comb,factorial,prod
import numpy as np
from flint import arb
from group_moment import maps
from group_rank_one_verify import up
from occupancy_memory import C,Z


def census(maximum=3):
    assert maximum in (2,3)
    images,columns,_=maps()
    mask=(1<<64)-1
    low=np.zeros((32,15),dtype=np.uint64);high=low.copy()
    for window in range(32):
        for bits in range(1,16):
            feedback=0
            for bit in range(4):
                if bits>>bit&1:feedback^=columns[4*window+bit]
            word=images[feedback]^(bits<<(4*window))
            low[window,bits-1]=word&mask;high[window,bits-1]=word>>64
    result={}
    for j in range(2,maximum+1):
        shapes=list(combinations_with_replacement(range(1,5),j))
        indices=np.indices((15,)*j).reshape(j,-1)
        weights=np.array([bits.bit_count() for bits in range(1,16)])
        shape_ids=np.array([shapes.index(tuple(sorted(weights[indices[:,i]]))) for i in range(indices.shape[1])])
        hist=np.zeros(len(shapes)*129,dtype=np.int64)
        for windows in combinations(range(32),j):
            lo=np.zeros(indices.shape[1],dtype=np.uint64);hi=lo.copy()
            for slot,window in enumerate(windows):
                lo^=low[window,indices[slot]];hi^=high[window,indices[slot]]
            output_weight=np.bitwise_count(lo)+np.bitwise_count(hi)
            local=np.bincount(shape_ids*129+output_weight,minlength=len(hist))
            if windows in (tuple(range(j)),tuple(range(32-j,32))):
                reference=Counter()
                for bits in product(range(1,16),repeat=j):
                    word=0;feedback=0
                    for window,mask_bits in zip(windows,bits):
                        word^=mask_bits<<(4*window)
                        for bit in range(4):
                            if mask_bits>>bit&1:feedback^=columns[4*window+bit]
                    weight=(images[feedback]^word).bit_count()
                    shape=tuple(sorted(b.bit_count() for b in bits))
                    reference[shapes.index(shape)*129+weight]+=1
                assert all(int(n)==reference[i] for i,n in enumerate(local))
            hist+=local
        rows={}
        for i,shape in enumerate(shapes):
            multiplicity=prod(factorial(n) for n in Counter(shape).values())
            counts=[int(v)*multiplicity for v in hist[129*i:129*(i+1)]]
            denominator=prod(range(32-j+1,33))*prod(comb(4,a) for a in shape)
            assert sum(counts)==denominator
            rows[shape]=counts,denominator
        result[j]=rows
        print('Exact cancellation-output census:',j,'windows;',len(shapes),'weight shapes; totals checked',flush=True)
    return result


def refine(base,counts,tilt,penalty,input_penalty=1,odd_penalty=1):
    powers=[(-arb(tilt)*w).exp() for w in range(129)]
    for j,rows in counts.items():
        values=[]
        for shape,(hist,denominator) in rows.items():
            value=sum((n*p for n,p in zip(hist,powers)),arb(0))/denominator
            values.append(up(value*arb(penalty)**shape.count(4)*arb(input_penalty)**sum(shape)*arb(odd_penalty)**sum(w%2 for w in shape)/2))
        base[j][C,Z]=min(base[j][C,Z],max(values))
    return base


if __name__=='__main__':census()
