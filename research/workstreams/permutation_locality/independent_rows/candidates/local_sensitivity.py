"""Adjoint sensitivity of the unrestricted bound, not a certificate.

Differentiate the actual finite region product and its normalized
without-replacement placement recurrence. The elasticities describe the
positive envelope calculation, not frequencies in the encoder itself.
No operator coefficient is changed by this diagnostic.
"""
import argparse
from math import comb, log

import numpy as np
from scipy.special import gammaln

from mass_density_screen import epochs, as_array, baseline, score
from mass_verify import alternative


def placement_tape(operators, degree, epochs=64, windows=32):
    operators=np.asarray(operators,dtype=float)
    if (operators.ndim!=3 or operators.shape[1]!=operators.shape[2]
            or not np.isfinite(operators).all() or np.any(operators<0)
            or not 0<=degree<=epochs*windows or len(operators)<=min(degree,windows)):
        raise ValueError('complete nonnegative square local operators required')
    tape=[np.eye(operators.shape[1])[None,:,:]]
    for epoch in range(1,epochs+1):
        previous=tape[-1]
        limit=min(degree,epoch*windows)
        current=np.zeros((limit+1,*operators.shape[1:]))
        for k,rs,old,weights in placement_weights(epoch,len(previous),limit,len(operators),windows):
            current[rs]+=weights[:,None,None]*(previous[old]@operators[k])
        tape.append(current)
    return tape


def placement_weights(epoch, previous_length, limit, local_length, windows):
    total=epoch*windows
    for k in range(min(windows,local_length-1,limit)+1):
        rs=np.arange(k,min(limit,k+previous_length-1)+1)
        old=rs-k
        exponent=(log(comb(windows,k))+gammaln(total-windows+1)-gammaln(old+1)
                  -gammaln(total-windows-old+1)-gammaln(total+1)
                  +gammaln(rs+1)+gammaln(total-rs+1))
        yield k,rs,old,np.exp(exponent)


def moment_gradient(matrix, steps, terminal):
    """Natural log(e_0 M**steps terminal) and its gradient in M."""
    matrix=np.asarray(matrix,dtype=float)
    terminal=np.asarray(terminal,dtype=float)
    if (steps<1 or matrix.shape!=(len(terminal),len(terminal))
            or not np.isfinite(matrix).all() or np.any(matrix<0)
            or not np.isfinite(terminal).all() or np.any(terminal<0)):
        raise ValueError('positive length and finite nonnegative matrix/terminal required')
    forward=[np.eye(len(terminal))[0]]
    value=0.
    for _ in range(steps):
        row=forward[-1]@matrix
        scale=row.sum()
        if not scale>0:
            raise ValueError('positive finite moment required')
        value+=log(scale)
        forward.append(row/scale)
    final=forward[-1]@terminal
    if not final>0:
        raise ValueError('positive terminal moment required')
    value+=log(final)
    gradient=np.zeros_like(matrix)
    backward=terminal.copy()
    for index in range(steps-1,-1,-1):
        denominator=forward[index]@matrix@backward
        if not denominator>0:
            raise ArithmeticError('positive adjoint normalization required')
        gradient+=np.outer(forward[index],backward)/denominator
        backward=matrix@backward
        backward/=backward.max()
    assert np.isfinite(gradient).all() and np.all(gradient>=0)
    return value,gradient


def sensitivity(operators, degree, p, terminal, *, epochs=64, windows=32, regions=256,
                return_gradient=False):
    """Return log moment and d(log moment)/d(log T[j,s,t])."""
    if not 0<p<1:
        raise ValueError('probability in (0,1) required')
    tape=placement_tape(operators,degree,epochs,windows)
    probabilities=np.array([comb(degree,j)*p**j*(1-p)**(degree-j) for j in range(degree+1)])
    matrix=np.einsum('j,jab->ab',probabilities,tape[-1])
    value,gradient=moment_gradient(matrix,regions,terminal)
    assert np.isclose(np.sum(matrix*gradient),regions,rtol=1e-9)
    adjoint=probabilities[:,None,None]*gradient
    local=np.zeros_like(operators)
    for epoch in range(epochs,0,-1):
        previous=tape[epoch-1]
        older=np.zeros_like(previous)
        for k,rs,old,weights in placement_weights(epoch,len(previous),len(adjoint)-1,len(operators),windows):
            weighted=weights[:,None,None]*adjoint[rs]
            local[k]+=np.sum(previous[old].transpose(0,2,1)@weighted,axis=0)
            older[old]+=weighted@operators[k].T
        adjoint=older
    elasticities=local*operators
    assert np.isfinite(elasticities).all() and np.all(elasticities>=0)
    assert np.isclose(elasticities.sum(),regions*epochs,rtol=1e-8)
    return (value,elasticities,local) if return_gradient else (value,elasticities)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--groups',type=int,default=80)
    parser.add_argument('--support',type=int,default=200)
    parser.add_argument('--tilt',default='.064')
    parser.add_argument('--penalty',default='.9')
    args=parser.parse_args()
    if (not 1<=args.groups<=2048 or not 38<=args.support<256
            or baseline.Q(args.tilt)<=0 or not 0<baseline.Q(args.penalty)<=1):
        parser.error('invalid support, occupancy, or tilt')
    print('FINITE-PRODUCT ADJOINT DIAGNOSTIC; no probability certificate',flush=True)
    base,census=epochs(args.tilt,args.penalty)
    caps=baseline.authenticated_caps()
    cdf=baseline.integer_cdf(baseline.weighted_cdf_upper(caps,1<<128,full_weight=1/baseline.Q(args.penalty)))
    shells=baseline.weighted_union_shells(caps,full_weight=1/baseline.Q(args.penalty))
    count=min(baseline.Q(cdf[args.support]),shells[args.support])
    for name,ops in (('baseline',base),('mass-ZC-6-12',alternative(base,census,args.tilt,args.penalty))):
        array=as_array(ops)
        region=placement_tape(array,args.groups)[-1]
        point,p=score(region,args.groups,args.support,count,args.tilt)
        moment,elasticities=sensitivity(array,args.groups,p,baseline.TAIL_TERMINAL)
        reference=baseline.log_power_moment(baseline.matrix_for_probabilities(region,[p]*args.groups),256,baseline.TAIL_TERMINAL)
        assert abs(moment-reference)<1e-7
        print('DIAGNOSTIC',name,'q/u',args.groups,args.support,'log2 score',point,'p',p,flush=True)
        print('EULER CHECK',elasticities.sum(),'expected',64*256,flush=True)
        by_transition=elasticities.sum(axis=0)
        labels=('Z','F','M','C','U48','U56','U64','U72','U80','L48','L56')
        indices=np.dstack(np.unravel_index(np.argsort(by_transition.ravel())[::-1][:15],by_transition.shape))[0]
        print('TOP TRANSITIONS natural-log elasticities',
              [(labels[s],labels[t],float(by_transition[s,t])) for s,t in indices],flush=True)
        selected=((3,0),(3,3),(1,0),(1,3))
        print('CANCELLATION/DENSITY BY LOCAL OCCUPANCY',
              [(j,*[float(elasticities[j,s,t]) for s,t in selected]) for j in range(1,33)],flush=True)
        print('UNIFORM RETURNS BY LOCAL OCCUPANCY',
              [(j,float(elasticities[j,4:9,0].sum())) for j in range(1,33)],flush=True)


if __name__=='__main__':
    main()
