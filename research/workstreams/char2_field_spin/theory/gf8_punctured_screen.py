"""Fixed local puncturing of the full-Cauchy packet, with a stronger outer.

One regional route is retained. The packet consumes `width` symbols and
emits a fixed `outputs`-row restriction of its original output. No second
permutation or random puncturing is assumed. Screens are floating point,
not complete distance certificates.
"""
from functools import lru_cache
import argparse
import json
import math

import numpy as np
from scipy.optimize import minimize

from cauchy_firstmoment_screen import geometry, log_choose, log_terminal_power
from gf8_mds_geometry_screen import roots
from gf8_rank_history_screen import modal_log_penalty, source_row_log
from reduced_distance_certificate import log_polynomial_product
from sparse_active_interval_screen import regional_from_coefficients, log_terminal


def profile(width=8, outputs=7, length=128, dimension=None, k_log=20):
    if dimension is None:
        if 3*length*outputs % (4*width):
            raise ValueError('outer dimension for rate3/4 must be integral')
        dimension = 3*length*outputs//(4*width)
    if not (width in (4,8,16,32) and 1 <= outputs <= width
            and 1 <= dimension < length <= 256):
        raise ValueError('invalid punctured geometry')
    g = geometry(k_log,length,dimension)
    if g['rows'] % width:
        raise ValueError('packets must divide a region')
    g['input_n'] = g['n']
    g['output_n'] = g['n']//width*outputs
    g['threshold'] = (g['output_n']+9)//10-1
    return g


def local_matrices(y,tau,width=8,outputs=7,state=8):
    if not (width in (4,8,16,32) and 1 <= outputs <= width and 1 <= state <= 32
            and 0 < y <= 1 and 1 <= tau <= 255):
        raise ValueError('invalid packet geometry or marks')
    zero=y*0+tau*0
    result=[[[zero,zero],[zero,zero]] for _ in range(width+1)]
    result[0][0][0]=y**0
    gamma=roots(tau,width,state)
    flats=[]
    value=zero
    for k in range(state):
        value+=math.comb(state,k)/tau**k
        flats.append(value)
    terms=[math.comb(outputs,z)*y**(outputs-z)*(1-y)**z for z in range(outputs+1)]
    for w in range(width+1):
        for z in range(outputs+1):
            pivots,remaining=min(w,z),max(0,z-w)
            if remaining<state:
                result[w][1][1]+=terms[z]*gamma[remaining]/tau**pivots
            if w and z<w:
                result[w][0][1]+=terms[z]/tau**z
            if w>state+z:
                result[w][0][0]+=terms[z]/tau**(state+z)
        if w:
            # Rectangular MDS birth/return maps may have nonzero kernels.
            minimum=max(0,outputs+1-w)
            result[w][0][1]=min(result[w][0][1],y**minimum)
            pure=y**minimum/tau**state*flats[min(w-1,state-1)]
            joint=sum(terms[z]/tau**(state+z)*flats[min(w-z-1,state-1)]
                      for z in range(min(w,outputs+1)))
            result[w][1][0]=min(pure,joint)
    return result


def packet(point,width,outputs,state):
    logx,logy,logtau=point[:3]
    local=np.asarray(local_matrices(math.exp(logy),min(255.,math.exp(logtau)),width,outputs,state),dtype=float)
    return np.einsum('w,wij->ij',
        [math.comb(width,w)*math.exp(w*logx) for w in range(width+1)],local)


def middle_bound(point,weight,width=8,outputs=7,state=8,length=128,dimension=None,k_log=20):
    g=profile(width,outputs,length,dimension,k_log)
    logx,logy,logtau,logtilt=map(float,point)
    return (length*modal_log_penalty(g['rows'])
        +g['rows']*source_row_log(weight/g['input_n'],logtau,logtilt,length,g['dimension'],True)
        +log_terminal_power(packet(point,width,outputs,state),g['input_n']//width)
        -weight*(logx+logtilt)-g['threshold']*logy-logtau)/math.log(2)


def sparse_bound(active,point,width=8,outputs=7,state=8,length=128,dimension=None,k_log=20):
    g=profile(width,outputs,length,dimension,k_log)
    logtilt,logtau=map(float,point)
    logy=-math.exp(logtilt)*active/g['output_n']
    tau=min(255.,math.exp(logtau))
    local=np.asarray(local_matrices(math.exp(logy),tau,width,outputs,state)[:min(width,active)+1],dtype=float)
    with np.errstate(divide='ignore'):
        factor=np.log(local)
    for w in range(len(factor)):
        factor[w]+=math.log(math.comb(width,w))
    result=np.array([[[0.,-math.inf],[-math.inf,0.]]])
    steps=g['rows']//width
    while steps:
        if steps&1:
            result=log_polynomial_product(result,factor,active)
        steps>>=1
        if steps:
            factor=log_polynomial_product(factor,factor,active)
    regional=regional_from_coefficients(result,active,g['rows'],logtau)
    return (log_choose(g['rows'],active)-((length-g['dimension'])*active+1)*logtau
        -g['threshold']*logy+log_terminal(regional,length))/math.log(2)


def search(value,sparse=False,width=8,outputs=7,state=8,length=128,dimension=None,
           iterations=60,start=None):
    g=profile(width,outputs,length,dimension)
    if sparse:
        target=int(value)
        objective=lambda p:sparse_bound(target,p,width,outputs,state,length,dimension)
        starts=[(6.5,2.),(7.,3.)]
        bounds=[(-2,10),(0,math.log(255))]
    else:
        target=round(value*g['input_n'])
        objective=lambda p:middle_bound(p,target,width,outputs,state,length,dimension)
        starts=[(-2.,-.2,math.log(16),0.),(0.,-1.,math.log(255),-1.)]
        bounds=[(-30,8),(-20,0),(0,math.log(255)),(-20,20)]
    if start is not None:
        starts.insert(0,start)
    choices=[minimize(objective,p,method='L-BFGS-B',bounds=bounds,
                      options={'maxiter':iterations,'ftol':1e-11,'maxls':25}) for p in starts]
    best=min(choices,key=lambda r:r.fun if math.isfinite(r.fun) else math.inf)
    return dict(outer=[length,g['dimension']],width=width,outputs=outputs,state=state,
        k_log=20,requested_k=2**20,padded_k=g['k'],rows=g['rows'],
        input_n=g['input_n'],output_n=g['output_n'],nominal_rate=g['k']/g['output_n'],
        mode='sparse' if sparse else 'middle',target=target,margin_bits=-float(best.fun),
        log_parameters=best.x.tolist(),optimizer_success=bool(best.success),iterations=int(best.nit),
        rigorous=False,whole_code=False,concrete_matrix_validated=False)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--width',type=int,default=8)
    parser.add_argument('--outputs',type=int,default=7)
    parser.add_argument('--state',type=int,default=8)
    parser.add_argument('--outer',type=int,default=128)
    parser.add_argument('--dimension',type=int)
    parser.add_argument('--sparse',type=int,nargs='*',default=[64,128,256])
    parser.add_argument('--density',type=float,nargs='*',default=[.05,.10,.14])
    parser.add_argument('--iterations',type=int,default=60)
    args=parser.parse_args()
    for sparse,values in ((True,args.sparse),(False,args.density)):
        previous=None
        for value in values:
            result=search(value,sparse,args.width,args.outputs,args.state,args.outer,
                          args.dimension,args.iterations,previous)
            previous=result['log_parameters']
            print(json.dumps(result),flush=True)
