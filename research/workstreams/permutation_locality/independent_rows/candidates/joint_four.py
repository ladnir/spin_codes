"""Exact joint feedback/output census through four distinct windows.

The bounded candidate keeps the existing three-window routine unchanged.
Pairs of windows are precomputed, then whole mask products are combined
with NumPy integer kernels. No sampling or floating histogram enters the
counts. All fresh-state distributions still use all 32 packet windows.
"""
import argparse
from collections import Counter
from itertools import combinations,combinations_with_replacement,product
from math import comb,factorial,prod

import numpy as np

from mass_density_screen import baseline
from group_moment import maps


def atoms(images,columns):
    feedback=np.zeros((32,15),dtype=np.uint32)
    low=np.zeros((32,15),dtype=np.uint64); high=low.copy()
    fresh=np.zeros((4,len(images)),dtype=np.uint8)
    for window in range(32):
        for mask in range(1,16):
            syndrome=0
            for bit in range(4):
                if mask>>bit&1: syndrome^=columns[4*window+bit]
            word=images[syndrome]^(mask<<(4*window))
            feedback[window,mask-1]=syndrome
            low[window,mask-1]=word&((1<<64)-1)
            high[window,mask-1]=word>>64
            fresh[mask.bit_count()-1,syndrome]+=1
    if int(fresh.max())!=1 or fresh[:,0].any():
        raise ValueError('candidate requires the checked fresh-state indicator property')
    assert list(map(int,fresh.sum(axis=1)))==[32*comb(4,b) for b in range(1,5)]
    membership=sum(fresh[a]*(1<<a) for a in range(4))
    return (feedback,low,high),membership


def pair_cache(data,positions):
    return {windows:tuple((row[windows[0],:,None]^row[windows[1],None,:]).ravel()
                          for row in data)
            for windows in combinations(positions,2)}


def batch(windows,data,pairs):
    if len(windows)==1:
        return tuple(row[windows[0]] for row in data)
    if len(windows)==2:
        return pairs[windows]
    left=pairs[windows[:2]]
    right=tuple(row[windows[2]] for row in data) if len(windows)==3 else pairs[windows[2:]]
    return tuple((a[:,None]^b[None,:]).ravel() for a,b in zip(left,right))


def shape_ids(j):
    shapes=list(combinations_with_replacement(range(1,5),j))
    lookup=np.zeros(1+4*5**3,dtype=np.int64)
    for index,shape in enumerate(shapes):
        lookup[sum(5**(b-1) for b in shape)]=index
    mask_codes=np.array([5**(mask.bit_count()-1) for mask in range(1,16)],dtype=np.int64)
    indices=np.indices((15,)*j).reshape(j,-1)
    codes=mask_codes[indices].sum(axis=0)
    return shapes,lookup[codes]


def direct(windows,images,columns,levels,shapes,membership):
    joint=np.zeros((len(shapes),len(levels),129),dtype=np.int64)
    fresh=np.zeros((4,len(shapes),129),dtype=np.int64)
    shape_index={shape:i for i,shape in enumerate(shapes)}
    level_index={v:i for i,v in enumerate(levels)}
    for masks in product(range(1,16),repeat=len(windows)):
        word=0; syndrome=0
        for window,mask in zip(windows,masks):
            word^=mask<<(4*window)
            for bit in range(4):
                if mask>>bit&1: syndrome^=columns[4*window+bit]
        shape=shape_index[tuple(sorted(mask.bit_count() for mask in masks))]
        weight=(images[syndrome]^word).bit_count()
        joint[shape,level_index[images[syndrome].bit_count()],weight]+=1
        for a in range(4):
            if int(membership[syndrome])>>a&1: fresh[a,shape,weight]+=1
    return joint,fresh


def census(j=4,positions=tuple(range(32)),inputs=None,check_samples=True):
    positions=tuple(positions)
    if (not isinstance(j,int) or isinstance(j,bool) or not 1<=j<=4
            or len(positions)<j or tuple(sorted(set(positions)))!=positions
            or any(not isinstance(w,int) or isinstance(w,bool) or not 0<=w<32 for w in positions)):
        raise ValueError('one through four windows in an ordered distinct subset of 0..31 required')
    images,columns,spectrum=maps() if inputs is None else inputs
    levels=[0]+sorted(spectrum)
    classes=np.searchsorted(levels,np.array([image.bit_count() for image in images]))
    data,membership=atoms(images,columns)
    pairs=pair_cache(data,positions)
    shapes,ids=shape_ids(j)
    joint=np.zeros((len(shapes),len(levels),129),dtype=np.int64)
    fresh=np.zeros((4,len(shapes),129),dtype=np.int64)
    total_subsets=comb(len(positions),j)
    assert total_subsets*15**j<1<<63
    selected={positions[:j],positions[-j:]}
    for number,windows in enumerate(combinations(positions,j),1):
        syndrome,low,high=batch(windows,data,pairs)
        weights=np.bitwise_count(low)+np.bitwise_count(high)
        indices=(ids*len(levels)+classes[syndrome])*129+weights
        local=np.bincount(indices,minlength=joint.size).reshape(joint.shape)
        joint+=local
        fresh_indices=ids*129+weights
        members=membership[syndrome]
        local_fresh=np.array([np.bincount(fresh_indices[(members&(1<<a))!=0],
                                         minlength=fresh.shape[1]*129).reshape(fresh.shape[1:])
                              for a in range(4)])
        fresh+=local_fresh
        if check_samples and windows in selected:
            direct_joint,direct_fresh=direct(windows,images,columns,levels,shapes,membership)
            assert np.array_equal(local,direct_joint) and np.array_equal(local_fresh,direct_fresh)
        if number%4000==0:
            print('JOINT FOUR exact window subsets',number,'of',total_subsets,flush=True)
    result={}
    for i,shape in enumerate(shapes):
        denominator=total_subsets*factorial(j)//prod(factorial(n) for n in Counter(shape).values())
        denominator*=prod(comb(4,b) for b in shape)
        rows={v:list(map(int,joint[i,k])) for k,v in enumerate(levels)}
        ff=[list(map(int,fresh[a,i])) for a in range(4)]
        assert sum(map(sum,rows.values()))==denominator
        assert all(sum(row)<=denominator for row in ff)
        assert all(not n or abs(v-sum(shape))<=w<=v+sum(shape)
                   for v,row in rows.items() for w,n in enumerate(row))
        result[shape]=rows,ff,denominator
    print('JOINT FOUR census complete:',j,'windows;',len(shapes),'shapes;',
          len(positions),'positions; exact totals and selected scalar comparisons pass',flush=True)
    return result,spectrum


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--windows',type=int,choices=(1,2,3,4),default=4)
    args=parser.parse_args()
    result=census(args.windows)
    from full_feedback_census import census as feedback
    from cancellation_joint import check_feedback
    check_feedback(result,feedback(max(2,args.windows)))
