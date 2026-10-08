"""Bounded floating full-Cauchy screen permitting outputs > input width.

Fixed b=256,k=136,width16,state8,outputs17 is the default next geometry.
The local rank formulas are the existing Cauchy formulas with their output
range extended. Their rectangular extension still requires independent audit.
No inverse-membership/dense certificate is borrowed from a square inner.
"""
import argparse
from fractions import Fraction as F
import json
import math

import numpy as np
from scipy.optimize import minimize

from cauchy_firstmoment_screen import geometry, log_choose, log_terminal_power
from gf8_rank_history_screen import modal_log_penalty, source_row_log
from high_distance_fields import roots
from higher_distance_screen import output_threshold
from reduced_distance_certificate import log_polynomial_product
from sparse_active_interval_screen import regional_from_coefficients, log_terminal


def local_matrices(y,tau,width=16,outputs=17,state=8,q=65536):
    if not (1 <= width+state <= 128 and 1 <= outputs+state <= 128
            and width >= 1 and outputs >= 1 and state >= 1
            and 0 < y <= 1 and 1 <= tau <= q-1):
        raise ValueError('invalid rectangular Cauchy map or marks')
    zero=y*0+tau*0
    result=[[[zero,zero],[zero,zero]] for _ in range(width+1)]
    result[0][0][0]=y**0
    gamma=roots(tau,max(width,outputs),state,q)
    inverse=F(1,tau) if isinstance(tau,(int,F)) else 1/tau
    powers=[inverse**r for r in range(state+outputs+1)]
    flats=[sum(math.comb(state,j)*powers[j] for j in range(k+1)) for k in range(state)]
    terms=[math.comb(outputs,z)*y**(outputs-z)*(1-y)**z for z in range(outputs+1)]
    for w in range(width+1):
        for z in range(outputs+1):
            pivots,remaining=min(w,z),max(0,z-w)
            if remaining < state:
                result[w][1][1] += terms[z]*gamma[remaining]*powers[pivots]
            if w and z < w:
                result[w][0][1] += terms[z]*powers[z]
            if w > state+z:
                result[w][0][0] += terms[z]*powers[state+z]
        if w:
            minimum=max(0,outputs+1-w)
            result[w][0][1]=min(result[w][0][1],y**minimum)
            pure=y**minimum*powers[state]*flats[min(w-1,state-1)]
            joint=sum(terms[z]*powers[state+z]*flats[min(w-z-1,state-1)]
                      for z in range(min(w,outputs+1)))
            result[w][1][0]=min(pure,joint)
    return result


def profile(length=256,dimension=136,width=16,outputs=17,state=8):
    g=geometry(20,length,dimension,math.lcm(64,width))
    g.update(width=width,outputs=outputs,state=state,input_n=g['n'],
             output_n=g['n']//width*outputs)
    g['threshold']=output_threshold(g['output_n'],F(2,5))
    return g


def bound(point,value,sparse=True,length=256,dimension=136,width=16,outputs=17,state=8):
    g=profile(length,dimension,width,outputs,state)
    if sparse:
        active=int(value)
        logc,logtau=map(float,point)
        logy=-math.exp(logc)*active/g['output_n']
    else:
        logx,logy,logtau,logtilt=map(float,point)
    tau=min(65535.,math.exp(logtau))
    local=np.asarray(local_matrices(math.exp(logy),tau,width,outputs,state),dtype=float)
    if sparse:
        with np.errstate(divide='ignore'):
            factor=np.log(local[:min(width,active)+1])
        for w in range(len(factor)):
            factor[w]+=math.log(math.comb(width,w))
        result=np.array([[[0.,-math.inf],[-math.inf,0.]]])
        steps=g['rows']//width
        while steps:
            if steps&1:result=log_polynomial_product(result,factor,active)
            steps>>=1
            if steps:factor=log_polynomial_product(factor,factor,active)
        region=regional_from_coefficients(result,active,g['rows'],logtau)
        return (log_choose(g['rows'],active)-((length-dimension)*active+1)*logtau
                -g['threshold']*logy+log_terminal(region,length))/math.log(2)
    if not 0 < value < g['input_n']:
        raise ValueError('middle bound excludes endpoints')
    packet=np.einsum('w,wij->ij',[math.comb(width,w)*math.exp(w*logx)
                                for w in range(width+1)],local)
    return (length*modal_log_penalty(g['rows'])
            +g['rows']*source_row_log(value/g['input_n'],logtau,logtilt,length,dimension,True)
            +log_terminal_power(packet,g['input_n']//width)
            -value*(logx+logtilt)-g['threshold']*logy-logtau)/math.log(2)


def search(value,sparse=True,iterations=16):
    g=profile()
    target=int(value) if sparse else round(value*g['input_n'])
    starts=[(6.3,2.5),(6.3,4.)] if sparse else [(-2.6,-.15,2.5,2.5),(0.,-1.,math.log(255),-1.)]
    limits=[(-2,10),(0,math.log(65535))] if sparse else [(-30,30),(-30,0),(0,math.log(65535)),(-30,30)]
    results=[minimize(lambda p:bound(p,target,sparse),p,method='L-BFGS-B',bounds=limits,
                      options=dict(maxiter=iterations,ftol=1e-10,maxls=20)) for p in starts]
    best=min(results,key=lambda r:r.fun if math.isfinite(r.fun) else math.inf)
    return dict(outer=[256,136],width=16,state=8,outputs=17,
                raw_expansion='32/17',rate='1/2',distance='2/5',bits=16,
                rows=g['rows'],input_n=g['input_n'],output_n=g['output_n'],
                mode='active' if sparse else 'weight',target=target,
                margin_bits=-float(best.fun),log_parameters=best.x.tolist(),
                optimizer_success=bool(best.success),iterations=int(best.nit),
                rigorous=False,whole_code=False,rectangular_local_audit_complete=False)


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--active',nargs='*',type=int,default=[16])
    p.add_argument('--density',nargs='*',type=float,default=[.03])
    p.add_argument('--iterations',type=int,default=16)
    a=p.parse_args()
    for sparse,values in ((True,a.active),(False,a.density)):
        for value in values:
            print(json.dumps(search(value,sparse,a.iterations)),flush=True)
