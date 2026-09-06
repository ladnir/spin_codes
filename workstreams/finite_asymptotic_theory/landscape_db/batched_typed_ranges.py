"""Batched finite covers of all dense row-type counts.

The cover is deterministic and includes every occupation from minimum to L.
Every box uses the typed positive-coefficient inequality, with its own fixed
tilt, shell probabilities, and coefficient proposal. No performance ranking
or outward-rounding claim follows from these binary64 evaluations.
"""
import math

import numpy as np
from scipy.special import gammaln

import dense_occupation_ranges as dense
import typed_dense_boxes as typed
from composition_occupation import terminal_logs


def three_bands(counts, block):
    bands = [[],[],[]]
    for w in sorted(counts):
        bands[0 if 16*w<3*block else 1 if 16*w<=13*block else 2].append(w)
    return [b for b in bands if b]


def partition(length, groups, minimum, depth):
    if not 1<=minimum<=length or groups<2 or not 0<=depth<=16:
        raise ValueError('invalid type partition')
    lo=np.zeros(groups,dtype=np.int64);hi=np.full(groups,length,dtype=np.int64)
    hi[0]=length-minimum
    current=[(lo,hi)]
    for _ in range(depth):
        following=[]
        for lo,hi in current:
            corners=typed.vertices(lo,hi,length)
            children=typed.split_box(lo,hi,corners)
            if not children:following.append((lo,hi))
            else:
                for a,b in children:
                    vertices=typed.vertices(a,b,length)
                    if len(vertices):following.append((vertices.min(axis=0),vertices.max(axis=0)))
        current=following
    return current


class TypedPartition:
    def __init__(self, counts, block, length, minimum=65, depth=7,
                 scales=(.5,.75,1.), eligibility_shifts=(-.5,0.,.5,1.)):
        self.block=block;self.length=length;self.cutoff=block*length//10
        self.bands,ps,costs=typed.category_parameters(counts,block,three_bands(counts,block))
        self.minimum=minimum;self.boxes=partition(length,len(ps),minimum,depth)
        corners=[typed.vertices(lo,hi,length) for lo,hi in self.boxes]
        self.starts=np.concatenate(([0],np.cumsum([len(c) for c in corners])))[:-1]
        self.vertices=np.concatenate(corners)
        self.owners=np.repeat(np.arange(len(corners)),[len(c) for c in corners])
        self.log_mult=gammaln(length+1)-gammaln(self.vertices+1).sum(axis=1)
        self.log_counts=np.array([typed.lattice_log_count(lo,hi) for lo,hi in self.boxes])
        proposals=np.stack([typed.proposal_for(c,length) for c in corners])
        self.witnesses=[];self.proposals=[];self.probabilities=[];self.costs=[]
        for scale in scales:
            if not math.isfinite(scale) or scale<0:raise ValueError('invalid shell scale')
            p=ps.copy();gamma=costs.copy()
            for g,band in enumerate(self.bands,start=1):
                if band==[block]:continue
                eta=scale*(math.log(ps[g])-math.log1p(-ps[g]))
                lp=-float(np.logaddexp(0.,-eta));ln=-float(np.logaddexp(0.,eta))
                p[g]=math.exp(lp)
                gamma[g]=max(math.log(counts[w])-math.log(math.comb(block,w))-w*lp-(block-w)*ln for w in band)
            for shift in eligibility_shifts:
                if not math.isfinite(shift) or abs(shift)>20:raise ValueError('invalid proposal shift')
                proposal=proposals.copy();proposal[:,1:]*=math.exp(shift)
                proposal/=proposal.sum(axis=1,keepdims=True)
                self.witnesses.append(dict(probability_scale=float(scale),eligibility_shift=float(shift)))
                self.proposals.append(proposal);self.probabilities.append(p);self.costs.append(gamma)
        if not self.witnesses:raise ValueError('empty witness grid')
        self.proposals=np.array(self.proposals);self.probabilities=np.array(self.probabilities);self.costs=np.array(self.costs)
        self.thetas=np.einsum('wbg,wg->wb',self.proposals,self.probabilities)
        self.best=np.full(len(self.boxes),np.inf)
        self.selected=[None]*len(self.boxes)

    def observe(self, epoch, t, lam, label):
        matrices=dense.epoch_mixture_logs(epoch,t,self.thetas.reshape(-1))
        moments=terminal_logs(matrices,self.block*self.length//t).reshape(self.thetas.shape)
        for index,(proposal,cost,moment) in enumerate(zip(self.proposals,self.costs,moments)):
            terms=(self.vertices*(cost-self.block*np.log(proposal[self.owners]))).sum(axis=1)
            terms+=moment[self.owners]-(self.block-1)*self.log_mult+self.cutoff*lam
            bounds=np.maximum.reduceat(terms,self.starts)+self.log_counts
            improve=np.flatnonzero(bounds<self.best)
            self.best[improve]=bounds[improve]
            for j in improve:self.selected[int(j)]=(label,index)
        if not np.isfinite(self.best).all():raise ArithmeticError('nonfinite typed range bound')

    def result(self):
        if any(w is None for w in self.selected):raise ValueError('no evaluated witness for every box')
        rows=[]
        for index,((lo,hi),(label,witness)) in enumerate(zip(self.boxes,self.selected)):
            rows.append(dict(lower=lo.tolist(),upper=hi.tolist(),log_union_upper=float(self.best[index]),
                             log_surprisal=label,**self.witnesses[witness],
                             proposal=self.proposals[witness,index].tolist(),
                             probabilities=self.probabilities[witness].tolist(),
                             log_density_costs=self.costs[witness].tolist()))
        return dict(occupation_min=self.minimum,occupation_max=self.length,
                    log_union_upper=float(np.logaddexp.reduce(self.best)),boxes=rows,bands=self.bands,
                    arithmetic='nearest binary64 diagnostic',all_integers_covered=True)
