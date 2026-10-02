"""Analytic derivatives for sparse witness proposals, never certificates.

The caller evaluates quantized proposals with its original objective and
replays retained witnesses with outward arithmetic.
"""
from collections import Counter
from math import log
import numpy as np
from scipy.special import expit,gammaln,logsumexp


def binomial_polynomial(count,p):
    result=np.array([1.]);factor=np.array([1-p,p])
    while count:
        if count&1:result=np.convolve(result,factor)
        count>>=1
        if count:factor=np.convolve(factor,factor)
    return result


def count_distribution(probabilities,multiplicities):
    """Count probabilities and their derivatives with respect to logits."""
    factors=[binomial_polynomial(n,p) for n,p in zip(multiplicities,probabilities)]
    prefix=[np.array([1.])]
    for factor in factors:prefix.append(np.convolve(prefix[-1],factor))
    suffix=[None]*len(factors)+[np.array([1.])]
    for i in range(len(factors)-1,-1,-1):suffix[i]=np.convolve(factors[i],suffix[i+1])
    derivatives=[]
    for i,(n,p) in enumerate(zip(multiplicities,probabilities)):
        lower=binomial_polynomial(n-1,p)
        slope=n*p*(1-p)*(np.pad(lower,(1,0))-np.pad(lower,(0,1)))
        derivatives.append(np.convolve(np.convolve(prefix[i],slope),suffix[i+1]))
    return prefix[-1],np.asarray(derivatives)


def log_power_gradient(matrix,derivatives,length,terminal):
    """Scaled matrix powering and its directional derivatives."""
    if type(length) is not int or length<1:raise ValueError('positive integer length required')
    peak=float(matrix.max())
    if not np.isfinite(peak) or peak<=0:raise ArithmeticError('positive finite moment required')
    value=matrix/peak;dvalue=derivatives/peak;scale=log(peak)
    row=np.zeros(len(matrix));row[0]=1.;drow=np.zeros((len(derivatives),len(matrix)))
    total=0.;remaining=length
    # Scaling factors are held constant while differentiating at this
    # evaluation point. Dividing values and derivatives by the same
    # factors therefore preserves the derivative of the original power.
    while remaining:
        if remaining&1:
            drow=drow@value+np.einsum('i,kij->kj',row,dvalue)
            row=row@value;peak=float(row.max())
            if not np.isfinite(peak) or peak<=0:raise ArithmeticError('positive finite row moment required')
            row/=peak;drow/=peak;total+=scale+log(peak)
        remaining>>=1
        if remaining:
            dvalue=dvalue@value+value@dvalue
            value=value@value;peak=float(value.max())
            if not np.isfinite(peak) or peak<=0:raise ArithmeticError('positive finite powered moment required')
            value/=peak;dvalue/=peak;scale=2*scale+log(peak)
    moment=float(row@terminal)
    if not np.isfinite(moment) or moment<=0:raise ArithmeticError('positive finite terminal moment required')
    return total+log(moment),(drow@terminal)/moment


class JointObjective:
    """Grouped Bernoulli objective using the same CDF suffix maximum."""
    def __init__(self,region,box,counts,terminal,offset,length=256):
        groups=Counter(box);self.intervals=sorted(groups)
        self.multiplicities=[groups[interval] for interval in self.intervals]
        self.size=len(terminal);self.terminal=np.asarray(terminal);self.offset=offset;self.length=length
        self.region=np.asarray(region[:len(box)+1]).reshape(len(box)+1,-1)
        self.folds=[]
        for lo,hi in self.intervals:
            supports=np.arange(lo,hi+1)
            choose=gammaln(257)-gammaln(supports+1)-gammaln(257-supports)
            increments=[counts[lo]]+[counts[u]-counts[u-1] for u in range(lo+1,hi+1)]
            logs=np.array([log(c) if c else -np.inf for c in increments])
            self.folds.append((supports,choose,logs))

    def __call__(self,logits):
        probabilities=expit(logits)
        masses,derivatives=count_distribution(probabilities,self.multiplicities)
        matrix=(masses@self.region).reshape(self.size,self.size)
        slopes=(derivatives@self.region).reshape(len(probabilities),self.size,self.size)
        value,gradient=log_power_gradient(matrix,slopes,self.length,self.terminal)
        for i,(p,n,(supports,choose,logs)) in enumerate(zip(probabilities,self.multiplicities,self.folds)):
            weights=-choose-supports*np.log(p)-(256-supports)*np.log1p(-p)
            # Track which support attains each suffix maximum. At a tie
            # this chooses a valid one-sided derivative of the envelope.
            selected=np.empty(len(weights),dtype=int);winner=len(weights)-1
            for j in range(len(weights)-1,-1,-1):
                if weights[j]>=weights[winner]:winner=j
                selected[j]=winner
            terms=logs+weights[selected];fold=float(logsumexp(terms))
            value+=n*fold
            gradient[i]+=n*float(np.exp(terms-fold)@(256*p-supports[selected]))
        return value+self.offset,gradient
