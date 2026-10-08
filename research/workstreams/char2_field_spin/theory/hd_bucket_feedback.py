"""Exact local moments for MDS emission and block-diagonal RS feedback.

No whole-code claim is made by this module. The feedback rank depends on
the counts in each block, unlike the global-MDS reference construction.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import comb,prod,factorial
from itertools import combinations_with_replacement
from hd_systematic_mds import maps as mds_maps,novel,moments as mds_moments


@lru_cache(None)
def maps(state=16,width=64,local_rank=8,global_nodes=False):
    if local_rank<1 or state%local_rank or width%(state//local_rank):
        raise ValueError('Equal feedback block dimensions required')
    blocks=state//local_rank;size=width//blocks
    a,_=mds_maps(state,width)
    c=tuple(tuple(novel(i if global_nodes else i%size,j%local_rank) if i//size==j//local_rank else 0
                  for i in range(width)) for j in range(state))
    return a,c


def moments(occupancies,z,state=16,width=64,local_rank=8,q=65536):
    z=F(z);size=width//(state//local_rank)
    if len(occupancies)!=state//local_rank or any(not 0<=j<=size for j in occupancies):
        raise ValueError('Invalid feedback-block occupancies')
    j=sum(occupancies)
    k=prod(sum(F(comb(l,c))*z**(l-c)*(1-z)**c/q**(c+min(local_rank,l-c))
               for c in range(l+1)) for l in occupancies)
    _,t,v=mds_moments(j,z,state,width,q)
    return k,t,v


def matrix(occupancies,z,state=16,width=64,local_rank=8,q=65536):
    k,t,v=moments(occupancies,z,state,width,local_rank,q)
    j=sum(occupancies);rank=sum(min(local_rank,l) for l in occupancies)
    gamma=min(k,F(z)**max(0,width+1-j-state)/q**rank,(q**state-1)*v)
    return k,t-k,v/(q**state-1) if j else F(0),(v if j else F(q-1,q)*v)+gamma/q


def digest(local_rank,global_nodes=False):
    a,c=maps(local_rank=local_rank,global_nodes=global_nodes)
    f=tuple(tuple(int(i==j) for j in range(64)) for i in range(64));value=0
    for m in (f,a,c):
        for row in m:
            for x in row:value=((value^x)*0x100000001b3)&((1<<64)-1)
    return value


@lru_cache(None)
def profiles(state=16,width=64,local_rank=8):
    if local_rank<1 or state%local_rank or width%(state//local_rank):
        raise ValueError('Equal feedback block dimensions required')
    blocks=state//local_rank;size=width//blocks
    result=[]
    for row in combinations_with_replacement(range(size+1),blocks):
        count=factorial(blocks)//prod(factorial(row.count(j)) for j in set(row))
        count*=prod(comb(size,j) for j in row)
        result.append((row,count))
    for j in range(width+1):assert sum(n for row,n in result if sum(row)==j)==comb(width,j)
    return tuple(result)


def grouped_packets(z,state=16,width=64,local_rank=8,bits=12,maximum=None,cutoff=None):
    from hd_sumrank_paired_certificate import round_up
    z=F(z);q=65536;dim=q**state-1;size=width//(state//local_rank)
    if maximum is None:maximum=width
    if not 0<=maximum<=width:raise ValueError('Invalid maximum degree')
    local=[sum(F(comb(l,c))*z**(l-c)*(1-z)**c/q**(c+min(local_rank,l-c))
               for c in range(l+1)) for l in range(size+1)]
    tv=[mds_moments(j,z,state,width)[1:] for j in range(maximum+1)]
    groups={};tails={}
    for row,count in profiles(state,width,local_rank):
        j=sum(row)
        if j>maximum:continue
        k=prod(local[l] for l in row);t,v=tv[j];rank=sum(min(local_rank,l) for l in row)
        gamma=min(k,z**max(0,width+1-j-state)/q**rank,dim*v)
        entries=(k,t-k,v/dim if j else F(0),(v if j else F(q-1,q)*v)+gamma/q)
        if cutoff is not None and j>cutoff:
            key=(j,rank)
            if key in tails:
                old_count,old_entries=tails[key]
                tails[key]=(old_count+count,tuple(max(x,y) for x,y in zip(entries,old_entries)))
            else:tails[key]=(count,entries)
            continue
        key=(j,tuple(round_up(x,bits) for x in entries))
        groups[key]=groups.get(key,0)+count
    for (j,_),(count,entries) in tails.items():
        key=(j,tuple(round_up(x,bits) for x in entries))
        groups[key]=groups.get(key,0)+count
    return tuple((j,count,entries) for (j,entries),count in groups.items())
