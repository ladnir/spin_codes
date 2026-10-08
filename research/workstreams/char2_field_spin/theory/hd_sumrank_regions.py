"""FLOAT regional matrix-coefficient screens, including short-block grouping.

No omitted occupancy tails and no benchmark. Local profile conditioning is
retained for steps=1. Larger blocks average profiles at each declared packet
count, then multiply matrices before applying the fractional power.
"""

from fractions import Fraction as F
from functools import lru_cache
from math import comb,log2

import numpy as np

from hd_sumrank_transfer import fractional_packets,operator_logs,screen
from hd_randomized_outer_state_profile import power,total


def product(left,right,degree):
    a,sa=left
    b,sb=right
    size=a.shape[1]
    length=min(degree+1,len(a)+len(b)-1)
    result=np.zeros((length,size,size))
    for i in range(size):
        for k in range(size):
            for j in range(size):
                result[:,i,j]+=np.convolve(a[:,i,k],b[:,k,j])[:length]
    peak=np.max(result)
    if not peak>0:
        raise ArithmeticError('FLOAT regional coefficient underflow')
    return result/peak,sa+sb+log2(peak)


@lru_cache(maxsize=12)
def short_paths(z,maximum,steps,state,mixer,input_mark):
    if steps not in (1,2,4) or (min(64,maximum)+1)**steps>20000000:
        raise ValueError('Bounded exact short-path screen')
    size=3 if mixer=='rank_one' else 2
    local=np.full((min(64,maximum)+1,size,size),-np.inf)
    for j,probability,packet in operator_logs(F(z),mixer,F(input_mark),state,min(64,maximum)):
        local[j]=np.logaddexp2(local[j],probability+np.array(packet))
    degrees,weights=np.array([0]),np.array([0.0])
    products=np.full((1,size,size),-np.inf)
    for i in range(size): products[0,i,i]=0
    for _ in range(steps):
        ds,ws,ps=[],[],[]
        for j in range(min(64,maximum)+1):
            keep=degrees+j<=maximum
            if np.any(keep):
                ds.append(degrees[keep]+j)
                ws.append(weights[keep]+log2(comb(64,j)))
                ps.append(np.logaddexp2.reduce(products[keep,:,:,None]+local[j][None,None,:,:],axis=2))
        degrees=np.concatenate(ds)
        weights=np.concatenate(ws)
        products=np.concatenate(ps)
    return degrees,weights,products


def region(rows,z,alpha,marker,*,steps=1,state=4,mixer='feedback_refresh',input_mark=F(1)):
    size=3 if mixer=='rank_one' else 2
    width=64*steps
    degree=min(rows,width)
    if steps==1:
        packets=fractional_packets(F(z),F(alpha),mixer=mixer,state=state,input_mark=F(input_mark))
        coefficients=np.array([np.exp2(log2(comb(64,j))+j*log2(marker)+(64-j)*log2(1-marker)
                                      +np.array(packets[j])) for j in range(degree+1)])
    else:
        degrees,weights,products=short_paths(F(z),degree,steps,state,mixer,F(input_mark))
        weighted=np.exp2(weights[:,None,None]+float(alpha)*products)
        coefficients=np.zeros((degree+1,size,size))
        for j in range(degree+1):
            coefficients[j]=weighted[degrees==j].sum(axis=0)*marker**j*(1-marker)**(width-j)
    result=(np.eye(size).reshape((1,size,size)),0.0)
    base=(coefficients,0.0)
    exponent=128//steps
    while exponent:
        if exponent&1: result=product(result,base,rows)
        exponent>>=1
        if exponent: base=product(base,base,rows)
    coefficient,scale=result
    if not np.all(coefficient[rows]>0):
        return dict(active_rows=rows,log2_upper_screen=float('inf'),
                    float_error='Regional coefficient underflow; no bound accepted',
                    directed=False,whole_code_certificate=False)
    with np.errstate(divide='ignore'):
        logarithms=np.log2(coefficient[rows])+scale
    choose=log2(comb(8192,rows))
    mass=choose+rows*log2(marker)+(8192-rows)*log2(1-marker)
    run=power(tuple(tuple(row) for row in logarithms-mass),256)
    bound=choose+float(alpha)*(rows*(4096-128*log2(65535))-log2(65535)
                              -838860*log2(z)-129*rows*log2(input_mark))+total(run[0])
    return dict(active_rows=rows,mark=str(z),alpha=str(alpha),marker=marker,steps=steps,state=state,
                mixer=mixer,input_mark=str(input_mark),log2_upper_screen=bound,directed=False,
                whole_code_certificate=False)


if __name__=='__main__':
    import argparse,json
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rows',type=int,nargs='+',default=[16,128,512,1024,2048])
    parser.add_argument('--state',type=int,default=8)
    parser.add_argument('--steps',type=int,default=2)
    parser.add_argument('--mixer',choices=('rank_one','uniform_pre','feedback_refresh'),default='feedback_refresh')
    parser.add_argument('--input-mark',type=F,default=F(1))
    args=parser.parse_args()
    for r in args.rows:
        z0=2**(-r/512)
        best=None
        for z in (F(z0),F(z0**0.75),F(z0**1.25)):
            for a in (F(1,32),F(1,16),F(1,12),F(1,8),F(1,4),F(1,2)):
                proposal=screen(r,z,a,state=args.state,mixer=args.mixer,input_mark=args.input_mark)
                marker=1/(1+2**(-proposal['marker_log_odds']))
                item=region(r,z,a,marker,steps=args.steps,state=args.state,mixer=args.mixer,input_mark=args.input_mark)
                if best is None or item['log2_upper_screen']<best['log2_upper_screen']:
                    best=item
        print(json.dumps(best),flush=True)
