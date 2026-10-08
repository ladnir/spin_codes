"""Variable-distance screens of the new full-Cauchy GF256 packet maps.

Floating optimization proposes class bounds, never a complete certificate.
The frozen10% helpers and maps are unchanged. Only the forbidden integer
output threshold changes relative to gf8_punctured_screen.
"""
from fractions import Fraction as F
import argparse
import json
import math

from scipy.optimize import minimize
from gf8_punctured_screen import profile,sparse_bound,middle_bound
from higher_distance_screen import output_threshold


def adjusted_bound(point,value,sparse=True,distance=F(3,10),width=4,outputs=4,
                   state=4,length=64,dimension=32,k_log=20):
    g=profile(width,outputs,length,dimension,k_log)
    logy=-math.exp(point[0])*value/g['output_n'] if sparse else point[1]
    reference=sparse_bound if sparse else middle_bound
    args=(value,point) if sparse else (point,value)
    return (reference(*args,width,outputs,state,length,dimension,k_log)
            -(output_threshold(g['output_n'],distance)-g['threshold'])*logy/math.log(2))


def search(value,sparse=True,distance=F(3,10),width=4,outputs=4,state=4,
           length=64,dimension=32,k_log=20,iterations=35,start=None):
    g=profile(width,outputs,length,dimension,k_log)
    target=int(value) if sparse else max(1,round(float(value)*g['input_n']))
    if not 1<=target<=(g['rows'] if sparse else g['input_n']):
        raise ValueError('class outside the construction geometry')
    starts=[(5.,math.log(32)),(5.,math.log(255))] if sparse else [(-2.,-.2,math.log(16),0.),(0.,-1.,math.log(255),-1.)]
    if start is not None:
        starts.insert(0,start)
    bounds=[(-2,10),(0,math.log(255))] if sparse else [(-30,8),(-20,0),(0,math.log(255)),(-20,20)]
    objective=lambda point:adjusted_bound(point,target,sparse,distance,width,outputs,state,length,dimension,k_log)
    candidates=[minimize(objective,point,method='L-BFGS-B',bounds=bounds,
        options={'maxiter':iterations,'ftol':1e-10,'maxls':20}) for point in starts]
    best=min(candidates,key=lambda result:result.fun if math.isfinite(result.fun) else math.inf)
    return dict(outer=[length,dimension],width=width,outputs=outputs,state=state,
        requested_k=2**k_log,padded_k=g['k'],input_n=g['input_n'],output_n=g['output_n'],
        rate=str(F(dimension*width,length*outputs)),distance=str(F(distance)),
        output_threshold=output_threshold(g['output_n'],distance),
        mode='active' if sparse else 'weight',target=target,margin_bits=-float(best.fun),
        log_parameters=best.x.tolist(),optimizer_success=bool(best.success),iterations=int(best.nit),
        rigorous=False,whole_code=False)


def gv_distance(bits,rate=F(1,2)):
    """Random-code asymptotic benchmark, not a general distance ceiling."""
    if not isinstance(bits,int) or bits<1 or not 0<float(rate)<1:
        raise ValueError('invalid alphabet bits or rate')
    logq=bits*math.log(2)
    logqminus1=logq+math.log1p(-math.exp(-logq))
    lo,hi=0.,1-math.exp(-logq)
    for _ in range(70):
        d=(lo+hi)/2
        entropy=(d*logqminus1-d*math.log(d)-(1-d)*math.log1p(-d))/logq
        if entropy<1-float(rate):lo=d
        else:hi=d
    return (lo+hi)/2


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outer',nargs=2,type=int,default=[64,32])
    parser.add_argument('--width',type=int,default=4)
    parser.add_argument('--state',type=int,default=4)
    parser.add_argument('--outputs',type=int,default=4)
    parser.add_argument('--distance',type=F,default=F(3,10))
    parser.add_argument('--active',nargs='*',type=int,default=[1,64])
    parser.add_argument('--density',nargs='*',type=float,default=[.03,.1,.3])
    parser.add_argument('--iterations',type=int,default=35)
    args=parser.parse_args()
    for sparse,values in ((True,args.active),(False,args.density)):
        previous=None
        for value in values:
            result=search(value,sparse,args.distance,args.width,args.outputs,args.state,*args.outer,
                          iterations=args.iterations,start=previous)
            previous=result['log_parameters']
            print(json.dumps(result),flush=True)
