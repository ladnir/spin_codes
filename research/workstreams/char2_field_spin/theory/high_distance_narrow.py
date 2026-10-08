"""Floating low-expansion screens with width-one/two full-Cauchy inners.

Uses the same rank/flat cover as high_distance_fields, with explicit geometry
and a width-one/two extension of the local formula. All fixed maps are over
embedded GF256; both diagonals are independently uniform nonzero full-field
scalars. A positive sampled class is neither a certificate nor interval coverage.
The width-four path is a regression control, not a change to frozen helpers.
"""
import argparse
from fractions import Fraction as F
import json
import math

import numpy as np
from scipy.optimize import minimize

from cauchy_firstmoment_screen import geometry, log_choose, log_terminal_power
from gf8_rank_history_screen import modal_log_penalty, source_row_log
from high_distance_fields import alphabet, roots, rectangular_cauchy_maps
from higher_distance_screen import output_threshold
from reduced_distance_certificate import log_polynomial_product
from sparse_active_interval_screen import regional_from_coefficients, log_terminal


def profile(width, outputs, length, dimension, k_log=20):
    if width not in (1,2,4) or not 1 <= outputs <= width or not 1 <= dimension < length <= 256:
        raise ValueError('invalid narrow geometry')
    g=geometry(k_log,length,dimension)
    if g['rows'] % width:
        raise ValueError('width must divide a region')
    g['input_n']=g['n']
    g['output_n']=g['n']//width*outputs
    return g


def local_matrices(y,tau,width=1,outputs=1,state=1,bits=32):
    q=alphabet(bits)
    if width not in (1,2,4) or not 1<=outputs<=width or not 1<=state<=32 or not 0<y<=1 or not 1<=tau<=q-1:
        raise ValueError('invalid local geometry or marks')
    zero=y*0+tau*0
    result=[[[zero,zero],[zero,zero]] for _ in range(width+1)]
    result[0][0][0]=y**0
    gamma=roots(tau,width,state,q)
    inverse=F(1,tau) if isinstance(tau,(int,F)) else 1/tau
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


def bound(point,value,sparse=True,bits=32,distance=F(2,5),width=1,outputs=1,
          state=1,length=64,dimension=32,k_log=20):
    q=alphabet(bits)
    g=profile(width,outputs,length,dimension,k_log)
    threshold=output_threshold(g['output_n'],distance)
    if sparse:
        active=int(value)
        if not 1<=active<=g['rows']:
            raise ValueError('active class outside geometry')
        logtilt,logtau=map(float,point)
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
        regional=regional_from_coefficients(result,active,g['rows'],math.log(tau))
        return (log_choose(g['rows'],active)-((length-dimension)*active+1)*math.log(tau)
                -threshold*logy+log_terminal(regional,length))/math.log(2)
    if not 0<value<g['input_n']:
        raise ValueError('middle density must be strictly between zero and one')
    logx,logy,logtau,logtilt=map(float,point)
    tau=min(q-1,math.exp(logtau))
    local=np.asarray(local_matrices(math.exp(logy),tau,width,outputs,state,bits),dtype=float)
    packet=np.einsum('w,wij->ij',[math.comb(width,w)*math.exp(w*logx) for w in range(width+1)],local)
    return (length*modal_log_penalty(g['rows'])
            +g['rows']*source_row_log(value/g['input_n'],math.log(tau),logtilt,length,dimension,True)
            +log_terminal_power(packet,g['input_n']//width)
            -value*(logx+logtilt)-threshold*logy-math.log(tau))/math.log(2)


def search(value,sparse=True,bits=32,distance=F(2,5),width=1,outputs=1,
           state=1,length=64,dimension=32,k_log=20,iterations=40,start=None):
    q=alphabet(bits)
    g=profile(width,outputs,length,dimension,k_log)
    target=int(value) if sparse else round(float(value)*g['input_n'])
    if not 1<=target<=(g['rows'] if sparse else g['input_n']-1):
        raise ValueError('class outside geometry')
    loglimit=math.log(q-1)
    starts=[(5.,math.log(16)),(5.,min(loglimit,math.log(255)))] if sparse else [(-2.,-.2,math.log(16),0.),(0.,-1.,min(loglimit,math.log(255)),-1.)]
    if start is not None:starts.insert(0,start)
    extent=max(30,loglimit+8)
    limits=[(-2,10),(0,loglimit)] if sparse else [(-extent,min(extent,640/width)),(-extent,0),(0,loglimit),(-extent,extent)]
    objective=lambda p:bound(p,target,sparse,bits,distance,width,outputs,state,length,dimension,k_log)
    choices=[minimize(objective,p,method='L-BFGS-B',bounds=limits,
        options={'maxiter':iterations,'ftol':1e-10,'maxls':20}) for p in starts]
    best=min(choices,key=lambda r:r.fun if math.isfinite(r.fun) else math.inf)
    return dict(outer=[length,dimension],width=width,outputs=outputs,state=state,bits=bits,
        k_log=k_log,input_n=g['input_n'],output_n=g['output_n'],rate=str(F(dimension*width,length*outputs)),
        distance=str(F(distance)),mode='active' if sparse else 'weight',target=target,
        margin_bits=-float(best.fun),log_parameters=best.x.tolist(),optimizer_success=bool(best.success),
        rigorous=False,whole_code=False,diagonal_distribution='independent uniform nonzero full-field scalars')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bits',type=int,default=32)
    p.add_argument('--outer',type=int,nargs=2,default=[64,32])
    p.add_argument('--width',type=int,default=1)
    p.add_argument('--state',type=int,default=1)
    p.add_argument('--outputs',type=int,default=1)
    p.add_argument('--distance',type=F,default=F(2,5))
    p.add_argument('--active',type=int,nargs='*',default=[1,16,64,128])
    p.add_argument('--density',type=float,nargs='*',default=[])
    p.add_argument('--iterations',type=int,default=40)
    a=p.parse_args()
    for sparse,values in ((True,a.active),(False,a.density)):
        previous=None
        for value in values:
            r=search(value,sparse,a.bits,a.distance,a.width,a.outputs,a.state,*a.outer,iterations=a.iterations,start=previous)
            previous=r['log_parameters']
            print(json.dumps(r),flush=True)
