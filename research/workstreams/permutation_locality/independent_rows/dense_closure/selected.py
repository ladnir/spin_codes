"""Selected-composition diagnostics and outward witnesses, not a cover."""
import argparse
from fractions import Fraction as Q
from math import lgamma,factorial
import numpy as np
from flint import arb
from scipy.optimize import minimize
import kernel
from mixture import actual_components,G,REGIONS,PACKETS,EPOCHS,N,capped_density_loss,logq,log_power
import density_tangent
import shuffle_comparison


def selected(data,components,counts,threshold,*,density_anchors=None,exact_shuffle=False,posterior_shuffle=False):
    if sum((density_anchors is not None,exact_shuffle,posterior_shuffle))>1:raise ValueError('choose one shuffle comparison')
    geometry=shuffle_comparison.structure(components,counts) if exact_shuffle else None
    posterior=shuffle_comparison.posterior_structure(components,counts) if posterior_shuffle else None
    ids=np.flatnonzero(counts);counts=np.asarray(counts)[ids]
    laws=np.array([list(map(float,components[i][2])) for i in ids])
    logcoeff=np.array([logq(components[i][1]) for i in ids])
    outer=float(counts@logcoeff)+lgamma(G+1)-sum(lgamma(int(c)+1) for c in counts)
    caps=tuple(int(sum(c for c,p in zip(counts,laws[:,j]) if p>0)) for j in range(5))
    if density_anchors is None:
        loss=REGIONS*logq(capped_density_loss(G,caps));tau=np.ones(5)
    else:
        if len(density_anchors)!=5:raise ValueError('five density anchors required')
        logloss,logtau=density_tangent.floating(G,density_anchors)
        loss=REGIONS*logloss;tau=np.exp(logtau)
    def objective(w):
        lam=w[0];tilt=np.exp(np.r_[0.,w[1:]])
        normalizers=laws@tilt
        weights=((counts/G/normalizers)@laws)*tau;scale=weights.sum()
        matrix=kernel.floating(data,weights/scale,lam)
        current_loss=REGIONS*shuffle_comparison.floating(geometry,tilt) if geometry is not None else loss
        if posterior is not None:current_loss=min(loss,REGIONS*shuffle_comparison.posterior_floating(posterior,tilt))
        return (outer+current_loss+REGIONS*(counts@np.log(normalizers))+PACKETS*np.log(scale)
                +log_power(matrix)+lam*threshold)/np.log(2)
    best=None
    for start in ((.12,-1,-1,-1,-1),(.5,-1,-2,-3,-4)):
        fit=minimize(objective,start,method='L-BFGS-B',bounds=[(.000001,4.)]+[(-6.,4.)]*4,
                     options={'maxiter':120,'ftol':1e-11})
        if best is None or fit.fun<best.fun:best=fit
    return dict(log2_upper=float(best.fun),tilt=list(map(float,best.x)))


def outward(data,components,counts,threshold,parameters,*,density_anchors=None,exact_shuffle=False,posterior_shuffle=False):
    if (len(counts)!=len(components) or any(type(c) is not int or c<0 for c in counts)
            or sum(counts)!=G or len(parameters)!=5):raise ValueError('complete composition and five tilts required')
    if type(threshold) is not int or not 0<=threshold<N:raise ValueError('integer output cutoff in [0,N) required')
    lam,*tilt=map(Q,parameters);tilt=[Q(1),*tilt]
    if lam<=0 or min(tilt)<=0:raise ValueError('positive tilts required')
    ids=[i for i,c in enumerate(counts) if c]
    z={i:sum(p*t for p,t in zip(components[i][2],tilt)) for i in ids}
    weights=[sum(Q(counts[i],G)*components[i][2][j]/z[i] for i in ids) for j in range(5)]
    caps=tuple(sum(counts[i] for i in ids if components[i][2][j]>0) for j in range(5))
    if sum((density_anchors is not None,exact_shuffle,posterior_shuffle))>1:raise ValueError('choose one shuffle comparison')
    if exact_shuffle:
        geometry=shuffle_comparison.structure(components,counts)
        logloss=kernel.aq(shuffle_comparison.exact(geometry,tilt)).log()
    elif posterior_shuffle:
        geometry=shuffle_comparison.posterior_structure(components,counts)
        loss=min(capped_density_loss(G,caps),shuffle_comparison.posterior_exact(geometry,tilt))
        logloss=kernel.aq(loss).log()
    elif density_anchors is None:logloss=kernel.aq(capped_density_loss(G,caps)).log()
    else:
        if len(density_anchors)!=5:raise ValueError('five density anchors required')
        logloss,logtau=density_tangent.outward(G,density_anchors)
        uppers=[kernel.up(kernel.aq(w)*t.exp()).man_exp() for w,t in zip(weights,logtau)]
        weights=[Q(int(m))*Q(2)**int(e) for m,e in uppers]
    scale=sum(weights);probabilities=[w/scale for w in weights]
    matrix=kernel.outward(data,probabilities,lam)**EPOCHS
    moment=matrix[0,0]+matrix[0,1]
    multiplicity=factorial(G)
    for c in counts:multiplicity//=factorial(c)
    # Dividing successively is integral: the remaining factors describe
    # disjoint blocks of a multinomial coefficient.
    value=(moment.log()+PACKETS*kernel.aq(scale).log()+kernel.aq(lam)*threshold
           +REGIONS*logloss+arb(multiplicity).log())
    for i in ids:value+=counts[i]*(kernel.aq(components[i][1]).log()+REGIONS*kernel.aq(z[i]).log())
    return kernel.up(value.exp())


def reference_weights(components,counts,parameters):
    """Floating reference weights for choosing witnesses only."""
    laws=np.array([list(map(float,row[2])) for row in components])
    tilt=np.array([1.,*map(float,parameters[1:])])
    return (np.array(counts)/G/(laws@tilt))@laws


def posterior_anchors(data,components,counts,parameters,anchors):
    """Integer tangent anchors proposed by numerical log derivatives.

    These derivatives are never used in verification. Every bounded
    integer anchor defines a valid density majorant independently.
    """
    _,logtau=density_tangent.floating(G,anchors)
    weights=reference_weights(components,counts,parameters)*np.exp(logtau)
    lam=float(parameters[0])
    def moment(ws):
        scale=ws.sum()
        return PACKETS*np.log(scale)+log_power(kernel.floating(data,ws/scale,lam))
    result=[];step=1e-5
    for j in range(5):
        if weights[j]==0:result.append(0);continue
        plus=weights.copy();minus=weights.copy()
        plus[j]*=np.exp(step);minus[j]*=np.exp(-step)
        estimate=(moment(plus)-moment(minus))/(2*step*REGIONS)
        if not np.isfinite(estimate):raise ArithmeticError('nonfinite density anchor proposal')
        result.append(max(0,min(G,int(round(estimate)))))
    return tuple(result)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--distance',type=float,default=.05)
    parser.add_argument('--groups',type=int,nargs='+',default=[59,128,256,512,1024,2048])
    parser.add_argument('--types',nargs='+',default=['0001','0002','0003','0022','0222','2222','3333','1111'])
    args=parser.parse_args()
    data=kernel.actual();components=actual_components();names=[c[0] for c in components]
    for q in args.groups:
        for name in args.types:
            counts=np.zeros(len(components),dtype=int);counts[0]=G-q;counts[names.index(name)]=q
            result=selected(data,components,counts,int(args.distance*N))
            print(q,name,result,flush=True)


if __name__=='__main__':main()
