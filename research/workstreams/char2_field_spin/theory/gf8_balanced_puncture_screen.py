"""Floating gates for a fixed periodic local output puncture.

Same full-Cauchy8/8 packet and one regional route. Only the number of
retained prefix rows varies with public packet index. These scores are
not whole-code certificates or claims about the optimal rate.
"""
import argparse,json,math
import numpy as np
from scipy.optimize import minimize

from cauchy_firstmoment_screen import geometry,log_choose,log_terminal_power
from gf8_punctured_screen import local_matrices,packet
from gf8_rank_history_screen import modal_log_penalty,source_row_log
from reduced_distance_certificate import log_polynomial_product
from sparse_active_interval_screen import regional_from_coefficients,log_terminal


def profile(dimension=74,pattern=(5,6),k_log=20):
    if not pattern or any(not isinstance(v,int) or not 1<=v<=8 for v in pattern):
        raise ValueError('invalid output pattern')
    g=geometry(k_log,128,dimension)
    if g['rows']%(8*len(pattern)):
        raise ValueError('period must divide every regional packet count')
    g['input_n']=g['n']
    g['output_n']=g['n']//(8*len(pattern))*sum(pattern)
    g['threshold']=(g['output_n']+9)//10-1
    return g


def period_matrix(point,pattern):
    result=np.eye(2)
    for outputs in pattern:result=result@packet(point,8,outputs,8)
    return result


def middle_bound(point,weight,dimension=74,pattern=(5,6)):
    g=profile(dimension,pattern)
    logx,logy,logtau,logtilt=map(float,point)
    return (128*modal_log_penalty(g['rows'])
        +g['rows']*source_row_log(weight/g['input_n'],logtau,logtilt,128,dimension,True)
        +log_terminal_power(period_matrix(point,pattern),g['input_n']//(8*len(pattern)))
        -weight*(logx+logtilt)-g['threshold']*logy-logtau)/math.log(2)


def sparse_bound(active,point,dimension=74,pattern=(5,6)):
    g=profile(dimension,pattern)
    logtilt,logtau=map(float,point)
    logy=-math.exp(logtilt)*active/g['output_n']
    tau=min(255.,math.exp(logtau))
    period=np.array([[[0.,-math.inf],[-math.inf,0.]]])
    for outputs in pattern:
        local=np.asarray(local_matrices(math.exp(logy),tau,8,outputs,8)[:min(8,active)+1],dtype=float)
        with np.errstate(divide='ignore'):factor=np.log(local)
        for w in range(len(factor)):factor[w]+=math.log(math.comb(8,w))
        period=log_polynomial_product(period,factor,active)
    result=np.array([[[0.,-math.inf],[-math.inf,0.]]])
    steps=g['rows']//(8*len(pattern))
    while steps:
        if steps&1:result=log_polynomial_product(result,period,active)
        steps>>=1
        if steps:period=log_polynomial_product(period,period,active)
    regional=regional_from_coefficients(result,active,g['rows'],logtau)
    return (log_choose(g['rows'],active)-((128-dimension)*active+1)*logtau
        -g['threshold']*logy+log_terminal(regional,128))/math.log(2)


def search(value,sparse=False,dimension=74,pattern=(5,6),iterations=85,start=None):
    g=profile(dimension,pattern)
    if sparse:
        target=int(value);objective=lambda p:sparse_bound(target,p,dimension,pattern)
        starts=[(6.6,2.5),(6.5,math.log(255))]
        bounds=[(-2,10),(0,math.log(255))]
    else:
        target=round(value*g['input_n']);objective=lambda p:middle_bound(p,target,dimension,pattern)
        starts=[(-1.,-.5,2.,2.),(0.,-2.,math.log(255),-1.)]
        bounds=[(-30,8),(-20,0),(0,math.log(255)),(-20,20)]
    if start is not None:starts.insert(0,start)
    choices=[minimize(objective,p,method='L-BFGS-B',bounds=bounds,
                      options={'maxiter':iterations,'ftol':1e-11,'maxls':25}) for p in starts]
    best=min(choices,key=lambda r:r.fun if math.isfinite(r.fun) else math.inf)
    return dict(outer=[128,dimension],pattern=list(pattern),state=8,width=8,
        requested_k=2**20,rows=g['rows'],padded_k=g['k'],input_n=g['input_n'],output_n=g['output_n'],
        nominal_rate=g['k']/g['output_n'],mode='sparse' if sparse else 'middle',target=target,
        margin_bits=-float(best.fun),log_parameters=best.x.tolist(),
        optimizer_success=bool(best.success),iterations=int(best.nit),rigorous=False,whole_code=False)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dimension',type=int,default=74)
    parser.add_argument('--pattern',type=int,nargs='+',default=[5,6])
    parser.add_argument('--sparse',type=int,nargs='*',default=[1,64,256])
    parser.add_argument('--density',type=float,nargs='*',default=[.03,.1,.5,.9,.99])
    args=parser.parse_args()
    for sparse,values in ((True,args.sparse),(False,args.density)):
        previous=None
        for value in values:
            result=search(value,sparse,args.dimension,tuple(args.pattern),start=previous)
            previous=result['log_parameters'];print(json.dumps(result),flush=True)
