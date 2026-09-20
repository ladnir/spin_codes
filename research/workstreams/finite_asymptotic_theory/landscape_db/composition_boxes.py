"""Fixed-occupation type boxes using exact region support coefficients.

Every box fixes lower counts of each outer band and relaxes only its
remaining row types. Density factors stay attached to the same counting
measures across every region. These are binary64 diagnostic bounds.
"""
import heapq
import math

import numpy as np
from scipy.special import gammaln

import balanced_occupation as balanced
import composition_occupation as sparse
import typed_dense_boxes as typed


def balanced_counts(lower, upper, total):
    """Integer vector maximizing the multinomial coefficient in a box."""
    lower = np.asarray(lower,dtype=np.int64)
    upper = np.asarray(upper,dtype=np.int64)
    if (lower.shape != upper.shape or lower.ndim != 1 or np.any(lower<0)
            or np.any(lower>upper) or lower.sum()>total or upper.sum()<total):
        raise ValueError('infeasible count box')
    lo, hi = int(lower.min()), int(upper.max())
    while lo < hi:
        mid = (lo+hi+1)//2
        if np.clip(mid,lower,upper).sum() <= total: lo = mid
        else: hi = mid-1
    counts = np.clip(lo,lower,upper)
    remaining = total-int(counts.sum())
    for j in np.argsort(counts,kind='stable'):
        if remaining and counts[j]<upper[j]: counts[j]+=1;remaining-=1
    if remaining: raise ArithmeticError('integer water filling failed')
    return counts


def log_assignment_count_upper(lower, upper, total):
    mode = balanced_counts(lower,upper,total)
    return (math.lgamma(total+1)-float(gammaln(mode+1).sum())
            +typed.lattice_log_count(lower,upper))


def binomial_logs(count, probability):
    if count < 0 or not 0 <= probability <= 1: raise ValueError('invalid binomial law')
    values = np.full(count+1,-np.inf)
    if probability == 0: values[0]=0.
    elif probability == 1: values[-1]=0.
    else:
        values = np.array([math.log(math.comb(count,j))+j*math.log(probability)
                           +(count-j)*math.log1p(-probability) for j in range(count+1)])
    return values


def positive_log_convolve(left, right):
    """Scale each positive polynomial; retain a log fallback for rare tails."""
    a=float(max(left)); b=float(max(right))
    raw=np.convolve(np.exp(left-a),np.exp(right-b))
    with np.errstate(divide='ignore'): result=np.log(raw)+a+b
    # A zero coefficient may come from a deterministic endpoint or underflow.
    for degree in np.flatnonzero(~np.isfinite(result)):
        lo=max(0,int(degree)-len(right)+1);hi=min(int(degree),len(left)-1)
        result[degree]=np.logaddexp.reduce(left[lo:hi+1]+right[degree-hi:degree-lo+1][::-1])
    return result


def shifted_mixture_logs(regions, distribution, remaining):
    """Sum pi_i R_(i+j) for j=0..remaining with columnwise scaling."""
    if len(regions)<len(distribution)+remaining: raise ValueError('missing region coefficients')
    windows=np.lib.stride_tricks.sliding_window_view(regions,len(distribution),axis=0)[:remaining+1]
    terms=windows+distribution
    # The small array avoids losing rare binomial tails when the region
    # coefficient attaining the unweighted maximum has zero probability.
    return np.logaddexp.reduce(terms,axis=-1)


class CompositionBoxes:
    def __init__(self, counts, block, bands, probabilities, maximum):
        if (maximum<1 or sorted(w for band in bands for w in band)!=sorted(counts)
                or len(probabilities)!=len(bands) or len(bands)<2):
            raise ValueError('complete band partition and positive occupation ceiling required')
        self.block=block;self.bands=bands;self.maximum=maximum
        self.probabilities=np.asarray(probabilities,dtype=float)
        self.log_gamma=[]
        for band,p in zip(bands,probabilities):
            if p==1 and band==[block]:
                self.log_gamma.append(math.log(counts[block]));continue
            if not 0<p<1: raise ValueError('interior fixed band probabilities required')
            self.log_gamma.append(max(math.log(counts[w])-math.log(math.comb(block,w))
                                      -w*math.log(p)-(block-w)*math.log1p(-p) for w in band))
        self.log_gamma=np.array(self.log_gamma)
        with np.errstate(divide='ignore'):
            self.envelope=balanced.Envelope(self.log_gamma/block,np.log(self.probabilities),np.log1p(-self.probabilities))
        self.binomials=[[binomial_logs(q,float(p)) for q in range(maximum+1)] for p in probabilities]
        self.distribution_cache={}

    def matrix(self, regions, lower, total):
        lower=np.asarray(lower,dtype=np.int64)
        if len(lower)!=len(self.bands) or np.any(lower<0) or lower.sum()>total or total>self.maximum:
            raise ValueError('invalid lower type counts')
        remaining=total-int(lower.sum());key=tuple(map(int,lower))
        distribution=self.distribution_cache.get(key)
        if distribution is None:
            distribution=np.array([0.])
            for g,count in enumerate(lower):
                distribution=positive_log_convolve(distribution,self.binomials[g][int(count)])
            if len(self.distribution_cache)>=2048:self.distribution_cache.clear()
            self.distribution_cache[key]=distribution
        current=shifted_mixture_logs(regions,distribution,remaining)
        current+=float(lower@self.log_gamma)/self.block
        for _ in range(remaining):current=self.envelope.apply(current[:-1],current[1:])
        return current[0]

    def bound(self, regions, lower, upper, total, length, cutoff, lam):
        matrix=self.matrix(regions,lower,total)
        return (math.log(math.comb(length,total))+log_assignment_count_upper(lower,upper,total)
                +cutoff*lam+float(sparse.terminal_logs(matrix[None,:,:],self.block)[0]))


def search(witnesses, occupation, length, cutoff, maximum_nodes=255, target_bits=80):
    """Partition all band counts summing to one fixed occupation.

    Witness tuples are (prepared model, region coefficients, lambda, label).
    A box may choose its witness; every witness is fixed inside that box.
    """
    if not witnesses or maximum_nodes<1 or not 1<=occupation<=length:
        raise ValueError('invalid composition search')
    groups=len(witnesses[0][0].bands)
    if any(w[0].bands!=witnesses[0][0].bands for w in witnesses):raise ValueError('mixed band partitions')
    def evaluate(lower,upper):
        corners=typed.vertices(lower,upper,occupation)
        if not len(corners):return None
        lower=corners.min(axis=0);upper=corners.max(axis=0)
        values=[model.bound(regions,lower,upper,occupation,length,cutoff,lam) for model,regions,lam,label in witnesses]
        index=int(np.argmin(values))
        return dict(lower=lower,upper=upper,corners=corners,own=values[index],best=values[index],
                    witness=index,parent=None,children=[])
    root=evaluate(np.zeros(groups,dtype=np.int64),np.full(groups,occupation,dtype=np.int64))
    nodes=[root];heap=[(-root['own'],0)]
    while heap and len(nodes)+2<=maximum_nodes and nodes[0]['best']>-target_bits*math.log(2):
        _,index=heapq.heappop(heap);node=nodes[index]
        for lo,hi in typed.split_box(node['lower'],node['upper'],node['corners']):
            child=evaluate(lo,hi)
            if child is None:continue
            child['parent']=index;node['children'].append(len(nodes))
            heapq.heappush(heap,(-child['own'],len(nodes)));nodes.append(child)
        cursor=index
        while cursor is not None:
            ancestor=nodes[cursor]
            if ancestor['children']:
                ancestor['best']=min(ancestor['own'],float(np.logaddexp.reduce([nodes[j]['best'] for j in ancestor['children']])))
            cursor=ancestor['parent']
    selected=[]
    def collect(index):
        node=nodes[index]
        if not node['children'] or node['own']<=node['best']:
            selected.append(dict(lower=node['lower'].tolist(),upper=node['upper'].tolist(),
                                 log_bound=node['own'],witness=witnesses[node['witness']][3]))
        else:
            for child in node['children']:collect(child)
    collect(0)
    return dict(occupation=occupation,log_union_upper=root['best'],nodes_evaluated=len(nodes),boxes=selected)
