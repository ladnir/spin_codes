"""Larger-field diagonal screens for the full-Cauchy packet construction.

The fixed GF256 matrices are embedded in GF(2**bits). Each input and
state diagonal entry is independently uniform over the entire nonzero
larger field. This is not a claim for GF256-only or restricted-polynomial
randomizers. Floating screens are not complete distance certificates.
"""
from fractions import Fraction as F
from functools import lru_cache
import argparse
import json
import math

import numpy as np
from scipy.optimize import minimize

from cauchy_firstmoment_screen import log_choose,log_terminal_power
from diagonal_refresh_rank import rank_potential_bound
from gf8_punctured_screen import profile
from gf8_rank_history_screen import modal_log_penalty,source_row_log
from higher_distance_screen import output_threshold
from reduced_distance_certificate import log_polynomial_product
from sparse_active_interval_screen import regional_from_coefficients,log_terminal


def alphabet(bits):
    if bits not in (8,16,32,64,128):
        raise ValueError('supported extension degrees are8,16,32,64,128')
    return 1<<bits


def rectangular_cauchy_maps(width=4,state=8):
    """Explicit GF256 maps; embedding them preserves their nonzero minors.

    Normalize all bottom state rows of a full Cauchy matrix by E^-1.
    Selected output rows, with either no state rows or all state rows,
    retain the original Cauchy rank. Arbitrary state-row subsets are not
    claimed to retain that rank merely by this normalization.
    """
    from gf8_packet_maps import MUL,INV,matmul
    if width<1 or state<1 or width+state>128:
        raise ValueError('disjoint GF256 Cauchy node sets required')
    n=width+state
    raw=tuple(tuple(INV[i^(128+j)] for j in range(n)) for i in range(n))
    e=tuple(row[width:] for row in raw[width:])
    work=[list(row)+[int(i==j) for j in range(state)] for i,row in enumerate(e)]
    for j in range(state):
        pivot=next(i for i in range(j,state) if work[i][j])
        work[j],work[pivot]=work[pivot],work[j]
        scale=INV[work[j][j]]
        work[j]=[MUL[scale][v] for v in work[j]]
        for i in range(state):
            if i!=j:
                factor=work[i][j]
                work[i]=[a^MUL[factor][b] for a,b in zip(work[i],work[j])]
    ei=tuple(tuple(row[state:]) for row in work)
    identity=tuple(tuple(int(i==j) for j in range(state)) for i in range(state))
    if matmul(e,ei)!=identity or matmul(ei,e)!=identity:
        raise ArithmeticError('state block inversion failed')
    f=tuple(row[:width] for row in raw[:width])
    a=tuple(row[width:] for row in raw[:width])
    c=matmul(ei,tuple(row[:width] for row in raw[width:]))
    return f,a,c,raw,e,ei


@lru_cache(maxsize=256,typed=True)
def roots(tau,width,state,q):
    # The internal dimension moments may exceed double range even though
    # their final normalized rank potential is at most one.
    exact=F(tau)
    values=tuple(max(rank_potential_bound(r,z,exact,q,state) for r in range(1,state+1))
                 for z in range(min(width,state-1)+1))
    return tuple(map(float,values)) if isinstance(tau,float) else values


def local_matrices(y,tau,width=4,outputs=2,state=4,bits=8):
    q=alphabet(bits)
    if not (width in (4,8,16,32) and 1<=outputs<=width and 1<=state<=32
            and 0<y<=1 and 1<=tau<=q-1):
        raise ValueError('invalid packet geometry or marks')
    zero=y*0+tau*0
    result=[[[zero,zero],[zero,zero]] for _ in range(width+1)]
    result[0][0][0]=y**0
    gamma=roots(tau,width,state,q)
    inverse=1/tau
    if isinstance(tau,(int,F)):
        inverse=F(1,tau)
    powers=[inverse**r for r in range(state+outputs+1)]
    flats=[]
    value=zero
    for k in range(state):
        value+=math.comb(state,k)*powers[k]
        flats.append(value)
    terms=[math.comb(outputs,z)*y**(outputs-z)*(1-y)**z for z in range(outputs+1)]
    for w in range(width+1):
        for z in range(outputs+1):
            pivots,remaining=min(w,z),max(0,z-w)
            if remaining<state:
                result[w][1][1]+=terms[z]*gamma[remaining]*powers[pivots]
            if w and z<w:
                result[w][0][1]+=terms[z]*powers[z]
            if w>state+z:
                result[w][0][0]+=terms[z]*powers[state+z]
        if w:
            minimum=max(0,outputs+1-w)
            result[w][0][1]=min(result[w][0][1],y**minimum)
            pure=y**minimum*powers[state]*flats[min(w-1,state-1)]
            joint=sum(terms[z]*powers[state+z]*flats[min(w-z-1,state-1)]
                      for z in range(min(w,outputs+1)))
            result[w][1][0]=min(pure,joint)
    return result


def bound(point,value,sparse=True,bits=8,distance=F(9,25),width=4,outputs=2,
          state=4,length=64,dimension=16,k_log=20):
    q=alphabet(bits)
    g=profile(width,outputs,length,dimension,k_log)
    threshold=output_threshold(g['output_n'],distance)
    if sparse:
        active=int(value)
        logtilt,logtau=map(float,point[:2])
        loglambda=float(point[2]) if len(point)>2 else 0.
        if loglambda<0:
            raise ValueError('minimum-weight mark must be at least one')
        logy=-math.exp(logtilt)*active/g['output_n']
        tau=min(q-1,math.exp(logtau))
        local=np.asarray(local_matrices(math.exp(logy),tau,width,outputs,state,bits)[:min(width,active)+1],dtype=float)
        with np.errstate(divide='ignore'):
            factor=np.log(local)
        for w in range(len(factor)):
            factor[w]+=math.log(math.comb(width,w))
        result=np.array([[[0.,-math.inf],[-math.inf,0.]]])
        steps=g['rows']//width
        while steps:
            if steps&1:result=log_polynomial_product(result,factor,active)
            steps>>=1
            if steps:factor=log_polynomial_product(factor,factor,active)
        regional=regional_from_coefficients(result,active,g['rows'],math.log(tau)+loglambda)
        return (log_choose(g['rows'],active)-((length-dimension)*active+1)*math.log(tau)
                -(length-dimension+1)*active*loglambda
                -threshold*logy+log_terminal(regional,length))/math.log(2)
    logx,logy,logtau,logtilt=map(float,point)
    if not 0<value<g['input_n']:
        raise ValueError('middle weight must be strictly between zero and raw length')
    tau=min(q-1,math.exp(logtau))
    local=np.asarray(local_matrices(math.exp(logy),tau,width,outputs,state,bits),dtype=float)
    packet=np.einsum('w,wij->ij',[math.comb(width,w)*math.exp(w*logx) for w in range(width+1)],local)
    return (length*modal_log_penalty(g['rows'])
            +g['rows']*source_row_log(value/g['input_n'],math.log(tau),logtilt,length,dimension,True)
            +log_terminal_power(packet,g['input_n']//width)
            -value*(logx+logtilt)-threshold*logy-math.log(tau))/math.log(2)


def search(value,sparse=True,bits=8,distance=F(9,25),width=4,outputs=2,state=4,
           length=64,dimension=16,k_log=20,iterations=30,start=None,minimum_weight=False):
    q=alphabet(bits)
    g=profile(width,outputs,length,dimension,k_log)
    target=int(value) if sparse else max(1,round(float(value)*g['input_n']))
    if not 1<=target<=(g['rows'] if sparse else g['input_n']-1):
        raise ValueError('class outside construction geometry')
    loglimit=math.log(q-1)
    starts=[(5.,math.log(16)),(5.,min(loglimit,math.log(255)))] if sparse else [(-2.,-.2,math.log(16),0.),(0.,-1.,min(loglimit,math.log(255)),-1.)]
    if not sparse and bits>32:
        starts.append((loglimit/4,-.9*loglimit,.9*loglimit,-loglimit/4))
        delta=float(distance)
        alpha=target/g['input_n']
        starts.append((loglimit/4,math.log(delta/(1-delta))-loglimit,loglimit,
                       -dimension/length*loglimit+math.log((1-alpha)/alpha)))
    if sparse and minimum_weight:
        starts=[point+(mark,) for point,mark in zip(starts,(0.,2.))]
    if start is not None:starts.insert(0,start)
    extent=max(30,loglimit+8)
    limits=[(-2,10),(0,loglimit)] if sparse else [(-extent,min(extent,640/width)),
        (-extent,0),(0,loglimit),(-extent,extent)]
    if sparse and minimum_weight:limits.append((0,15))
    objective=lambda p:bound(p,target,sparse,bits,distance,width,outputs,state,length,dimension,k_log)
    choices=[minimize(objective,p,method='L-BFGS-B',bounds=limits,
        options={'maxiter':iterations,'ftol':1e-10,'maxls':20}) for p in starts]
    best=min(choices,key=lambda r:r.fun if math.isfinite(r.fun) else math.inf)
    return dict(outer=[length,dimension],width=width,outputs=outputs,state=state,bits=bits,q=q,
        k_log=k_log,input_n=g['input_n'],output_n=g['output_n'],
        rate=str(F(dimension*width,length*outputs)),distance=str(F(distance)),
        output_threshold=output_threshold(g['output_n'],distance),mode='active' if sparse else 'weight',
        target=target,margin_bits=-float(best.fun),log_parameters=best.x.tolist(),
        optimizer_success=bool(best.success),iterations=int(best.nit),rigorous=False,whole_code=False,
        minimum_source_weight_mark=bool(sparse and minimum_weight),
        diagonal_distribution='independent uniform nonzero elements of the full coefficient field')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bits',type=int,default=8)
    parser.add_argument('--outer',nargs=2,type=int,default=[64,16])
    parser.add_argument('--outputs',type=int,default=2)
    parser.add_argument('--width',type=int,default=4)
    parser.add_argument('--state',type=int,default=4)
    parser.add_argument('--distance',type=F,default=F(9,25))
    parser.add_argument('--active',nargs='*',type=int,default=[1,16,64,96,128])
    parser.add_argument('--density',nargs='*',type=float,default=[])
    parser.add_argument('--iterations',type=int,default=30)
    parser.add_argument('--minimum-weight',action='store_true')
    args=parser.parse_args()
    for sparse,values in ((True,args.active),(False,args.density)):
        previous=None
        for value in values:
            result=search(value,sparse,args.bits,args.distance,args.width,args.outputs,args.state,*args.outer,
                          iterations=args.iterations,start=previous,minimum_weight=args.minimum_weight)
            previous=result['log_parameters']
            print(json.dumps(result),flush=True)
