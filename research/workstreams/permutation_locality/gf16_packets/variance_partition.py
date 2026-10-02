"""Outer counting and shuffle bounds on a partition of average variance.

Proposals use floating optimization. Outward evaluation reconstructs the
entire variance partition and checks its rational duals independently.
"""
from fractions import Fraction as Q
import numpy as np
from scipy.optimize import brentq
from scipy.special import logsumexp
from flint import arb,ctx
import scalar_cover as sc
import shuffle_variance as sv


def outer_witness(logs,features,active,mean,occupancy,interval):
    """Propose mean, occupancy, and variance duals for two closed intervals."""
    mean_cell=(float(mean),float(mean)) if np.isscalar(mean) else tuple(map(float,mean))
    occupancy=float(occupancy)
    features=np.asarray(features,dtype=float);active=np.asarray(active,dtype=float)
    variances=features*(1-features);logs=np.asarray(logs,dtype=float)
    lo,hi=map(float,interval);best=None
    for positive in (False,True):
        target=lo if positive else hi
        def evaluate(gamma):
            candidates=[]
            for sign in (False,True):
                mean_target=mean_cell[0 if sign else 1]
                eta,mu=sc.outer_dual(logs+gamma*variances,features,active,mean_target,occupancy,sign)
                values=logs+eta*features+mu*active+gamma*variances;norm=logsumexp(values)
                candidates.append((norm-eta*mean_target-mu*occupancy-gamma*target,
                    float(np.exp(values-norm)@variances-target),eta,mu))
            return min(candidates,key=lambda row:row[0])
        left,right=(0.,20000.) if positive else (-20000.,0.)
        if evaluate(left)[1]>=0:gamma=left
        elif evaluate(right)[1]<=0:gamma=right
        else:gamma=brentq(lambda g:evaluate(g)[1],left,right,xtol=1e-9)
        _,_,eta,mu=evaluate(gamma)
        dual=[Q(round(float(v)*10**8),10**8) for v in (eta,mu,gamma)]
        eta,mu,gamma=map(float,dual)
        value=float(logsumexp(logs+eta*features+mu*active+gamma*variances)
            -min(eta*x for x in mean_cell)-mu*occupancy-min(gamma*lo,gamma*hi))
        if best is None or value<best[0]:best=value,dual
    return best


def maximum_variance(cell):
    lo,hi=map(Q,cell)
    if not 0<lo<=hi<1:raise ValueError('interior mean cell required')
    return min(hi,1-lo,Q(1,4))


def factor(model,cell,lower,variance_dual):
    """Rational upper on the density factor throughout a variance part."""
    lower=Q(lower)
    if not 0<=lower<=maximum_variance(cell):raise ValueError('valid variance lower endpoint required')
    key='variance-partition',cell,lower,tuple(map(Q,variance_dual)),ctx.prec
    if key not in model.loss_cache:
        original=model.comparison_loss(cell,model.tilt,variance_dual)
        value=sv.density_interval_upper(sc.G,[sc.G*x for x in cell],sc.G*lower)
        m,e=value.upper().man_exp()
        model.loss_cache[key]=min(original,Q(int(m))*Q(2)**int(e))
    return model.loss_cache[key]


def propose(model,cell,base):
    if not model.variance_shuffle or model.variance_bins<1:raise ValueError('enabled variance partition required')
    score,witness=base
    if Q(witness['tilt'])!=model.tilt:raise ValueError('base-tilt witness required')
    _,_,logs=model.family(model.tilt);features=np.array(list(map(float,model.features)));active=np.array(model.active)
    eta,mu=map(lambda v:float(Q(v)),witness['parameters'][1:])
    outside=(sc.G*logsumexp(logs+eta*features+mu*active)
        -sc.G*min(eta*float(x) for x in cell)-mu*model.q_min)/np.log(2)
    old_loss=model.comparison_loss(cell,model.tilt,witness['variance_dual'])
    offset=score-outside-sc.REGIONS*sc.logq(old_loss)/np.log(2)
    upper=maximum_variance(cell);parts=[];scores=[]
    for j in range(model.variance_bins):
        interval=(upper*j/model.variance_bins,upper*(j+1)/model.variance_bins)
        count,dual=outer_witness(logs,features,active,cell,model.q_min/sc.G,interval)
        loss=factor(model,cell,interval[0],witness['variance_dual'])
        scores.append(sc.G*count+sc.REGIONS*sc.logq(loss))
        parts.append(dict(interval=list(map(str,interval)),dual=list(map(str,dual))))
    value=float(offset+logsumexp(scores)/np.log(2))
    return value,dict(witness,variance_partition=parts)


def validate(cell,parts):
    upper=maximum_variance(cell)
    if not isinstance(parts,list) or not 1<=len(parts)<=1024:raise ValueError('bounded nonempty variance partition required')
    previous=Q(0);result=[]
    for part in parts:
        if len(part['interval'])!=2 or len(part['dual'])!=3:raise ValueError('variance interval and three duals required')
        lo,hi=map(Q,part['interval']);eta,mu,gamma=map(Q,part['dual'])
        if lo!=previous or not lo<hi<=upper or mu<0:raise ValueError('complete ordered partition and nonnegative occupancy dual required')
        result.append(((lo,hi),(eta,mu,gamma)));previous=hi
    if previous!=upper:raise ValueError('variance partition does not cover its full range')
    return result


def outward_outer(model,cell,witness):
    """Outer comparison mass times regional density factors, with rounding."""
    if not model.variance_shuffle or not model.variance_bins or Q(witness['tilt'])!=model.tilt:
        raise ValueError('variance partition only applies to enabled base tilt')
    parts=validate(cell,witness['variance_partition']);cs,_,_=model.family(model.tilt)
    result=arb(0)
    for (lo,hi),(eta,mu,gamma) in parts:
        total=sum((sc.kernel.aq(c)*sc.kernel.aq(eta*f+mu*a+gamma*f*(1-f)).exp()
            for c,f,a in zip(cs,model.features,model.active)),arb(0))
        exponent=sc.G*(min(eta*x for x in cell)+min(gamma*lo,gamma*hi))+mu*model.q_min
        loss=factor(model,cell,lo,witness['variance_dual'])
        value=(sc.G*total.log()-sc.kernel.aq(exponent)+sc.REGIONS*sc.kernel.aq(loss).log()).exp()
        result=sc.kernel.up(result+value)
    return result
