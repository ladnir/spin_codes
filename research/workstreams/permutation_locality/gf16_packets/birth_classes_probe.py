"""Diagnostic only: retain the expansion weight of newly activated states.

The character census is exact integer arithmetic. Kernel evaluations and
frontier scores in this module are floating proposals, not certificates.
"""
import argparse
from fractions import Fraction as Q
import json
from math import comb
from pathlib import Path
import sys
import numpy as np
import birth_refresh as base
import scalar_cover as sc


def walsh(values):
    """Exact unnormalized Walsh transform, with bounded int64 intermediates."""
    result=np.asarray(values,dtype=np.int64).copy()
    n=len(result)
    if n==0 or n&(n-1):raise ValueError('power-of-two vector required')
    if int(np.abs(result).max())*n>(1<<62):raise ValueError('unsafe integer transform')
    half=1
    while half<n:
        blocks=result.reshape(-1,2*half)
        left=blocks[:,:half].copy();right=blocks[:,half:].copy()
        blocks[:,:half]=left+right;blocks[:,half:]=left-right
        half*=2
    return result


def attach(data,images):
    S=1<<data['bits'];W=data['windows'];radix=W+1
    if len(images)!=S or data['bits']>24:raise ValueError('bounded full state census required')
    powers=np.array([radix**j for j in range(5)],dtype=np.uint64)
    chars=np.arange(S,dtype=np.uint32);codes=np.zeros(S,dtype=np.uint64)
    for offset in range(0,4*W,4):
        weights=sum(np.bitwise_count(chars&c)&1 for c in data['columns'][offset:offset+4])
        codes+=powers[weights]
    order=np.argsort(codes);unique,starts=np.unique(codes[order],return_index=True)
    record_codes=np.array([sum(int(n)*radix**j for j,n in enumerate(row)) for row in data['records']],dtype=np.uint64)
    remap=np.searchsorted(unique,record_codes)
    if not np.array_equal(unique[remap],record_codes):raise ArithmeticError('character profiles disagree')
    image_weights=np.array([int(x).bit_count() for x in images],dtype=np.int64)
    levels=np.unique(image_weights[1:]);counts=[];census=[]
    for level in levels:
        indicator=(image_weights==level).astype(np.int64);indicator[0]=0
        counts.append(int(indicator.sum()))
        spectrum=walsh(indicator)
        row=np.add.reduceat(spectrum[order],starts)[remap]
        if sum(map(int,row))!=0:raise ArithmeticError('zero state entered a nonzero class')
        census.append(row)
    census=np.array(census,dtype=np.int64)
    expected=-np.array(data['multiplicities'],dtype=np.int64)
    expected[int(np.flatnonzero(record_codes==codes[0])[0])]+=S
    if not np.array_equal(census.sum(axis=0),expected):raise ArithmeticError('nonzero class partition failed')
    return dict(data,birth_class_levels=levels,birth_class_counts=np.array(counts,dtype=np.int64),
                birth_class_census=census,birth_class_float=census.astype(float),
                image_histogram_weights=data['histograms']@np.arange(5))


def prepare(images,columns,bits,updates=2):
    return attach(base.prepare(images,columns,bits,updates),images)


def actual(updates=2):
    images,_,_=sc.kernel.maps()
    data=attach(base.actual(updates),images)
    print('BIRTH CLASS CENSUS',list(map(int,data['birth_class_levels'])),flush=True)
    return data


def class_masses(data,p,z):
    values=np.array([((1+z)**(4-r)*(1-z)**r-1)/15 for r in range(5)])
    fourier=np.prod((1-p+p*values)[None,:]**data['records'],axis=1)
    return data['birth_class_float']@fourier/(1<<data['bits'])


def candidates(data,probabilities,tilt,*,class_returns=False):
    original=base.floating(data,probabilities,tilt)
    p=float(1-probabilities[0]);z=np.exp(-tilt);W=data['windows']
    levels=data['birth_class_levels'];alpha=2.**-data['updates'];beta=1-alpha
    matrix=np.zeros((3+len(levels),3+len(levels)));matrix[:3,:3]=original
    matrix[0,1:3]=0
    matrix[0,3:]=np.maximum(0,class_masses(data,p,z))
    factors=(1-16*p/15)*z**np.arange(5)+(p/15)*(1+z)**4
    moments=np.prod(factors[None,:]**data['histograms'],axis=1)
    if class_returns:
        powers=z**np.arange(4*W+1)
        fibers=np.minimum((data['rank_float']@powers)*(1+z)**(4*np.arange(W+1)-data['bits'])[None,:],
                          data['trimmed_float']@powers)
        occupancies=np.array([comb(W,j)*p**j*(1-p)**(W-j) for j in range(W+1)])
    empty=(1-p)**W
    for i,level in enumerate(levels,3):
        selected=moments[data['image_histogram_weights']==level]
        total=float(selected.max());quiet=empty*z**int(level)
        matrix[i,0]=min(original[1,0],total)
        if class_returns:
            selected=fibers[data['image_histogram_weights']==level]
            matrix[i,0]=min(matrix[i,0],alpha*float(occupancies@selected.max(axis=0))+beta*total/((1<<data['bits'])-1))
        matrix[i,1]=alpha*max(0,total-quiet)
        matrix[i,2]=beta*total
        matrix[i,i]=alpha*quiet
    result={'all_classes':matrix,'original':original}
    values=np.array([((1+z)**(4-r)*(1-z)**r-1)/15 for r in range(5)])
    coefficients=base.proposal.elementary(data['birth_records_float'],values,data['birth_degree'])
    by_occupancy=(data['birth_class_float']@coefficients)/(1<<data['bits'])
    low=np.zeros(len(levels))
    tails=base.proposal.candidates(data,p,z,original)
    for cutoff in range(1,data['birth_degree']+2):
        j=cutoff-1
        low+=p**j*(1-p)**(W-j)*by_occupancy[:,j]
        candidate=matrix.copy()
        candidate[0,2]=tails[f'birth_from_j{cutoff}'][0,2]
        candidate[0,3:]=np.maximum(0,low)
        result[f'classes_below_{cutoff}']=candidate
    return result


def floating(data,probabilities,tilt):
    options=candidates(data,probabilities,tilt,class_returns=data.get('class_returns',False))
    return min(options.values(),key=sc.log_power)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--points',nargs='+',default=['.03','.035','.04','.045','.05'])
    parser.add_argument('--distance',default='9/100')
    parser.add_argument('--minimum-groups',type=int,default=129)
    parser.add_argument('--variance-shuffle',action='store_true')
    parser.add_argument('--class-returns',action='store_true')
    parser.add_argument('--output',type=Path)
    args=parser.parse_args();data=actual();data['class_returns']=args.class_returns
    model=sc.Model(sc.actual_components(129,Q(21,50),Q(1001,1000),True),data,
                   int(Q(args.distance)*sc.N),args.minimum_groups,tilt=Q(1,8),inner=sys.modules[__name__],
                   variance_shuffle=args.variance_shuffle)
    rows=[]
    for x in map(Q,args.points):
        score,witness=model.proposal((x,x));row=dict(x=str(x),proposal=score,witness=witness)
        rows.append(row);print(json.dumps(row),flush=True)
    if args.output:args.output.write_text(json.dumps(dict(diagnostic_only=True,row_parity=True,distance=args.distance,
        minimum_groups=args.minimum_groups,variance_shuffle=args.variance_shuffle,class_returns=args.class_returns,rows=rows),indent=2)+'\n')


if __name__=='__main__':main()
