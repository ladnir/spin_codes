"""FLOAT exact-regional coefficient diagnostic for sum-rank maps."""

from fractions import Fraction as F
from functools import lru_cache
from math import comb,log2

import numpy as np

from hd_sumrank_transfer import operator_logs
from hd_randomized_outer_sparse_ladder import float_polynomial_product,float_endpoint


@lru_cache(maxsize=12)
def paths(z,maximum,steps,state=4,input_mark=F(1)):
    if comb(maximum+steps,steps)>150000:
        raise ValueError('Too many short-block paths')
    local=np.zeros((maximum+1,3,3))
    for j,probability,packet in operator_logs(F(z),'rank_one',F(input_mark),state,maximum):
        local[j]+=np.exp2(probability+np.array(packet))
    q,d=65536,65536**state-1
    min_distance=64-8*((state-1)//4)-(state-1)%4
    for j in range(1,min(maximum,3)+1):
        bound=float(input_mark)**j*float(z)**max(0,min_distance-8*j)*q**(-j)*((j/64)/q+(1-1/q)*(q**j-1)/d)
        local[j,2,0]=min(local[j,2,0],bound)
    degrees,weights=np.array([0]),np.array([1.0])
    products=np.eye(3).reshape((1,3,3))
    for _ in range(steps):
        next_degrees,next_weights,next_products=[],[],[]
        for j in range(maximum+1):
            keep=degrees+j<=maximum
            if np.any(keep):
                next_degrees.append(degrees[keep]+j)
                next_weights.append(weights[keep]*comb(64,j))
                next_products.append(products[keep]@local[j])
        degrees=np.concatenate(next_degrees)
        weights=np.concatenate(next_weights)
        products=np.concatenate(next_products)
    return degrees,weights,products


def regions(z,alpha,maximum,steps,state=4,input_mark=F(1)):
    degrees,weights,products=paths(F(z),maximum,steps,state,F(input_mark))
    coefficients=np.zeros((maximum+1,3,3))
    weighted=weights[:,None,None]*products**float(alpha)
    for j in range(maximum+1):
        coefficients[j]=weighted[degrees==j].sum(axis=0)
    result=np.eye(3).reshape((1,3,3))
    exponent=128//steps
    while exponent:
        if exponent&1:
            result=float_polynomial_product(result,coefficients,maximum)
        exponent>>=1
        if exponent:
            coefficients=float_polynomial_product(coefficients,coefficients,maximum)
    return tuple(result[r]/comb(8192,r) for r in range(maximum+1))


if __name__=='__main__':
    import argparse,json
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--state',type=int,default=8)
    parser.add_argument('--maximum',type=int,default=8)
    parser.add_argument('--steps',type=int,default=4)
    parser.add_argument('--input-marks',type=F,nargs='+',default=(F(1),F(2),F(16)))
    args=parser.parse_args()
    best={r:None for r in range(1,args.maximum+1)}
    for log_z in (F(-1,2048),F(-1,1024),F(-1,512),F(-1,256),F(-1,128),F(-1,64),F(-1,32),F(-1,16)):
        z=F(2**float(log_z))
        for v in args.input_marks:
            for alpha in (F(1,32),F(1,16),F(1,8),F(1,4),F(1,2),F(1)):
                all_regions=regions(z,alpha,args.maximum,args.steps,args.state,v)
                for r in best:
                    bound=float_endpoint(r,all_regions[r],log2(z),alpha)-float(alpha)*129*r*log2(v)
                    item=dict(active_rows=r,mark=str(z),alpha=str(alpha),input_mark=str(v),state=args.state,
                              steps=args.steps,log2_upper_screen=bound,directed=False,whole_code_certificate=False)
                    if best[r] is None or bound<best[r]['log2_upper_screen']:
                        best[r]=item
        print(json.dumps(dict(completed_log_z=str(log_z),best=list(best.values()))),flush=True)
