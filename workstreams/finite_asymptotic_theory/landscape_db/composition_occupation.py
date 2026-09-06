"""Sparse bounds retaining each outer-weight-band composition separately.

Auxiliary Bernoulli laws dominate deterministic counting measures. The
composition is held fixed across every region and every state transition.
All numerical output is nearest-binary64 diagnostic evidence.
"""
import itertools
import math

import numpy as np

import activation_occupation as general
from activation_q1 import matrix_product


def sparse_bands(counts, block, tail_bands=8, singleton_prefix=4):
    if tail_bands < 1 or singleton_prefix < 0:
        raise ValueError('invalid sparse partition')
    weights=sorted(w for w in counts if w!=block)
    prefix=weights[:singleton_prefix]
    remaining=weights[singleton_prefix:]
    groups=[[w] for w in prefix]
    if remaining:
        groups.extend(list(map(int,a)) for a in np.array_split(remaining,min(tail_bands,len(remaining))))
    if block in counts: groups.append([block])
    return groups


def compositions(groups, occupation):
    if groups < 1 or occupation < 1:
        raise ValueError('positive group count and occupation required')
    return np.array(list(itertools.combinations_with_replacement(range(groups),occupation)),dtype=np.int32)


def log_multiplicities(indices):
    q=indices.shape[1]
    result=np.full(len(indices),math.lgamma(q+1))
    # q is small. Runs count repeated bands without per-composition dictionaries.
    run=np.ones(len(indices),dtype=np.int32)
    for j in range(1,q):
        same=indices[:,j]==indices[:,j-1]
        result[~same]-=np.array([math.lgamma(int(v)+1) for v in run[~same]])
        run=np.where(same,run+1,1)
    result-=np.array([math.lgamma(int(v)+1) for v in run])
    return result


def probability_logs(indices, active, inactive):
    q=indices.shape[1]
    probabilities=np.full((len(indices),q+1),-np.inf)
    probabilities[:,0]=0.
    for j in range(q):
        p=active[indices[:,j],None]
        not_p=inactive[indices[:,j],None]
        old=probabilities[:,:j+1].copy()
        probabilities[:,:j+1]=old+not_p
        probabilities[:,1:j+2]=np.logaddexp(probabilities[:,1:j+2],old+p)
    return probabilities


def mixture_logs(regions, probabilities):
    if len(regions)<probabilities.shape[1]:
        raise ValueError('insufficient region coefficients')
    matrices=np.full((len(probabilities),3,3),-np.inf)
    for j in range(probabilities.shape[1]):
        np.logaddexp(matrices,probabilities[:,j,None,None]+regions[j],out=matrices)
    return matrices


def terminal_logs(matrices, block):
    if block < 1:
        raise ValueError('positive block length required')
    current=np.full_like(matrices,-np.inf)
    current[:,0,0]=current[:,1,1]=current[:,2,2]=0.
    while block:
        if block&1: current=matrix_product(current,matrices)
        block>>=1
        if block: matrices=matrix_product(matrices,matrices)
    return np.logaddexp.reduce(current[:,0],axis=1)


class SparseComposition:
    def __init__(self, counts, block, occupation, tail_bands=8, singleton_prefix=4):
        if any(not isinstance(n,int) or n<=0 for n in counts.values()):
            raise ValueError('deterministic integer shell counts or simultaneous caps required')
        self.counts=counts
        self.block=block
        self.occupation=occupation
        self.bands=sparse_bands(counts,block,tail_bands,singleton_prefix)
        self.indices=compositions(len(self.bands),occupation)
        self.log_mult=log_multiplicities(self.indices)
        band_mass=np.array([math.log(sum(counts[w] for w in band)) for band in self.bands])
        self.trivial=np.sum(band_mass[self.indices],axis=1)
        self.auxiliary_cache={}

    def components(self, regions, cutoff, lam, shift):
        if shift not in self.auxiliary_cache:
            roots,active,inactive=general.density_roots(self.counts,self.block,self.bands,shift)
            self.auxiliary_cache[shift]=(probability_logs(self.indices,active,inactive),
                                         self.block*np.sum(roots[self.indices],axis=1))
        probabilities,cost=self.auxiliary_cache[shift]
        matrices=mixture_logs(regions,probabilities)
        return np.minimum(self.trivial,terminal_logs(matrices,self.block)+cost+cutoff*lam)

    def aggregate(self, components, length):
        if len(components)!=len(self.indices) or length<self.occupation:
            raise ValueError('incorrect composition coverage or row count')
        return math.log(math.comb(length,self.occupation))+float(np.logaddexp.reduce(components+self.log_mult))
