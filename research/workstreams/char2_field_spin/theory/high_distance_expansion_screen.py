"""Lower-expansion floating controls, with chronological puncture periods.

This reuses the full-field local rank cover unchanged. Periodic retention
changes the construction but adds neither a routing stage nor randomness.
Outer lengths above256 require a larger fixed coefficient field; profile
metadata distinguishes it from the embedded GF256 inner coefficients.
All outputs are numerical screening results, never whole-code certificates.
"""
from fractions import Fraction as F
import math

import numpy as np
from scipy.optimize import minimize

from cauchy_firstmoment_screen import geometry,log_choose,log_terminal_power
from gf8_rank_history_screen import modal_log_penalty,source_row_log
from high_distance_fields import alphabet,local_matrices
from higher_distance_screen import output_threshold
from reduced_distance_certificate import log_polynomial_product
from sparse_active_interval_screen import regional_from_coefficients,log_terminal


def profile(length=96,dimension=32,pattern=(3,3,2),width=4,state=4,
            distance=F(3,10),k_log=20,outer_bits=8):
    pattern=tuple(pattern)
    outer_nodes=alphabet(outer_bits)
    if (not pattern or width not in (4,8,16,32) or not 1<=state<=32
            or not 1<=dimension<length<=min(65536,outer_nodes)
            or any(not isinstance(v,int) or not 1<=v<=width for v in pattern)):
        raise ValueError('invalid outer, packet, or retention period')
    if not 0<F(distance)<1:
        raise ValueError('distance must be between zero and one')
    g=geometry(k_log,length,dimension,math.lcm(64,width*len(pattern)))
    output_n=g['n']*sum(pattern)//(width*len(pattern))
    return dict(outer=[length,dimension],width=width,state=state,pattern=list(pattern),
        outer_field_bits=outer_bits,outer_max_nodes=outer_nodes,fixed_inner_field_bits=8,
        k_log=k_log,requested_k=2**k_log,padded_k=g['k'],rows=g['rows'],
        input_n=g['n'],output_n=output_n,rate=str(F(g['k'],output_n)),
        internal_expansion=str(F(length,dimension)),distance=str(F(distance)),
        output_threshold=output_threshold(output_n,distance))


def period_matrix(point,pattern,width,state,bits):
    logx,logy,logtau=map(float,point[:3])
    tau=min(alphabet(bits)-1,math.exp(logtau))
    weights=[math.comb(width,w)*math.exp(w*logx) for w in range(width+1)]
    product,scale=np.eye(2),0.
    for outputs in pattern:
        local=np.asarray(local_matrices(math.exp(logy),tau,width,outputs,state,bits),dtype=float)
        product=product@np.einsum('w,wij->ij',weights,local)
        normalizer=float(product.max())
        product/=normalizer
        scale+=math.log(normalizer)
    return product,scale


def bound(point,value,sparse=True,bits=16,length=96,dimension=32,
          pattern=(3,3,2),width=4,state=4,distance=F(3,10),k_log=20,outer_bits=8):
    alphabet(outer_bits)
    if bits%outer_bits:
        raise ValueError('the outer coefficient field must embed in the diagonal field')
    g=profile(length,dimension,pattern,width,state,distance,k_log,outer_bits)
    period=len(pattern)
    if sparse:
        active=int(value)
        if not 1<=active<=g['rows']:
            raise ValueError('active count outside geometry')
        logc,logtau=map(float,point)
        logy=-math.exp(logc)*active/g['output_n']
        tau=min(alphabet(bits)-1,math.exp(logtau))
        logtau=math.log(tau)
        identity=np.array([[[0.,-math.inf],[-math.inf,0.]]])
        factor=identity
        for outputs in pattern:
            local=np.asarray(local_matrices(math.exp(logy),tau,width,outputs,state,bits)[:min(width,active)+1],dtype=float)
            with np.errstate(divide='ignore'):
                local=np.log(local)
            for w in range(len(local)):
                local[w]+=math.log(math.comb(width,w))
            factor=log_polynomial_product(factor,local,active)
        result=identity
        steps=g['rows']//(width*period)
        while steps:
            if steps&1:result=log_polynomial_product(result,factor,active)
            steps>>=1
            if steps:factor=log_polynomial_product(factor,factor,active)
        regional=regional_from_coefficients(result,active,g['rows'],logtau)
        return (log_choose(g['rows'],active)-((length-dimension)*active+1)*logtau
            -g['output_threshold']*logy+log_terminal(regional,length))/math.log(2)
    if not 0<value<g['input_n']:
        raise ValueError('middle weight must exclude both endpoints')
    logx,logy,logtau,logtilt=map(float,point)
    packet,scale=period_matrix(point,pattern,width,state,bits)
    count=g['input_n']//(width*period)
    return (length*modal_log_penalty(g['rows'])
        +g['rows']*source_row_log(value/g['input_n'],logtau,logtilt,length,dimension,True)
        +log_terminal_power(packet,count)+count*scale
        -value*(logx+logtilt)-g['output_threshold']*logy-logtau)/math.log(2)


def search(value,sparse=True,bits=16,length=96,dimension=32,
           pattern=(3,3,2),width=4,state=4,distance=F(3,10),k_log=20,
           iterations=25,start=None,outer_bits=8):
    alphabet(outer_bits)
    if bits%outer_bits:
        raise ValueError('the outer coefficient field must embed in the diagonal field')
    g=profile(length,dimension,pattern,width,state,distance,k_log,outer_bits)
    target=int(value) if sparse else round(float(value)*g['input_n'])
    loglimit=math.log(alphabet(bits)-1)
    extent=max(30,loglimit+8)
    if sparse:
        starts=[(5.,math.log(16)),(5.,min(loglimit,math.log(255)))]
        bounds=[(-2,10),(0,loglimit)]
    else:
        starts=[(-2.,-.2,math.log(16),0.),(0.,-1.,min(loglimit,math.log(255)),-1.)]
        bounds=[(-extent,min(extent,640/width)),(-extent,0),(0,loglimit),(-extent,extent)]
    if start is not None:starts.insert(0,start)
    objective=lambda p:bound(p,target,sparse,bits,length,dimension,pattern,width,state,distance,k_log,outer_bits)
    candidates=[minimize(objective,p,method='L-BFGS-B',bounds=bounds,
        options={'maxiter':iterations,'ftol':1e-10,'maxls':20}) for p in starts]
    best=min(candidates,key=lambda r:r.fun if math.isfinite(r.fun) else math.inf)
    return dict(**g,bits=bits,mode='active' if sparse else 'weight',target=target,
        margin_bits=-float(best.fun),log_parameters=best.x.tolist(),
        optimizer_success=bool(best.success),iterations=int(best.nit),
        rigorous=False,whole_code=False,
        diagonal_distribution='independent uniform nonzero elements of the full coefficient field')
