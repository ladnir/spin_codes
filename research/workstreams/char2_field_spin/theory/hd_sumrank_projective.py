"""FLOAT projective-state diagnostic for the uniform-pre-refresh control.

Retains sums of state histories before fractional powers. Outward cone
rounding supplies a finite positive operator; acceptance still requires
directed reconstruction of every transition and the Collatz witness.
"""

from fractions import Fraction as F
from functools import lru_cache
from math import comb, log2

import numpy as np
from scipy.sparse import coo_matrix
from scipy.special import expit

from hd_sumrank_transfer import operator_logs


@lru_cache(maxsize=4)
def transitions(mark, alpha, input_mark=F(1), spacing=0.5,mixer='uniform_pre',state=4,adaptive=False):
    """Upper normalized two-vectors, scale^alpha, and profile probabilities."""
    # Store each coordinate independently: 1-expit(x) loses tiny positives.
    if adaptive:
        positive=np.unique(np.concatenate((np.arange(0,4,.05),np.arange(4,16,.25),np.arange(16,64,1),np.arange(64,513,4))))
        odds=np.concatenate((-positive[:0:-1],positive))
    else:
        odds = np.arange(-1024.0,1024.0+spacing/2,spacing)
    zero = np.concatenate(([1.0],expit(-odds*np.log(2)),[0.0]))
    live = np.concatenate(([0.0],expit(odds*np.log(2)),[1.0]))
    finite = np.arange(1,len(odds)+1)
    result=[]
    for j, log_probability, packet in operator_logs(F(mark),mixer,F(input_mark),state):
        p=np.exp2(np.array(packet))
        y0=zero*p[0,0]+live*p[1,0]
        y1=zero*p[0,1]+live*p[1,1]
        with np.errstate(divide='ignore',invalid='ignore'):
            ratio=np.log2(y1)-np.log2(y0)
        # Round toward the balanced ray: this increases the small component.
        # Its outward scaling cost tends to one near the two boundary rays.
        rounded=np.where(ratio<0,np.searchsorted(odds,ratio,side='left'),np.searchsorted(odds,ratio,side='right')-1)
        target=np.clip(rounded,0,len(odds)-1)+1
        target[y1==0]=0
        target[y0==0]=len(zero)-1
        with np.errstate(divide='ignore',invalid='ignore'):
            factor=np.maximum(np.divide(y0,zero[target],out=np.zeros_like(y0),where=zero[target]!=0),
                              np.divide(y1,live[target],out=np.zeros_like(y1),where=live[target]!=0))
        weight=np.power(factor,float(alpha))
        result.append((j,log_probability,target,weight))
    return len(zero),tuple(result)


def screen(rows,mark,alpha,marker,input_mark=F(1),spacing=0.5,iterations=200,mixer='uniform_pre',state=4,adaptive=False):
    size,edges=transitions(F(mark),F(alpha),F(input_mark),spacing,mixer,state,adaptive)
    sources=np.arange(size)
    src,dst,data=[],[],[]
    for j,probability,target,weight in edges:
        probability+=log2(comb(64,j))+j*log2(marker)+(64-j)*log2(1-marker)
        scaled=np.exp2(probability)*weight
        keep=scaled>1e-300
        src.append(sources[keep])
        dst.append(target[keep])
        data.append(scaled[keep])
    matrix=coo_matrix((np.concatenate(data),(np.concatenate(src),np.concatenate(dst))),shape=(size,size)).tocsr()
    vector=np.ones(size)
    best=None
    choose=log2(comb(8192,rows))
    base=(choose+float(alpha)*(rows*(4096-128*log2(65535))-log2(65535)
                             -838860*log2(mark)-129*rows*log2(input_mark))
          -256*(choose+rows*log2(marker)+(8192-rows)*log2(1-marker)))
    for iteration in range(iterations):
        following=matrix@vector
        upper=np.max(following/vector)
        bound=base+32768*log2(upper)+log2(vector[0])-log2(np.min(vector))
        if best is None or bound<best['log2_upper_screen']:
            best=dict(active_rows=rows,mark=str(mark),alpha=str(alpha),marker=marker,
                      input_mark=str(input_mark),spacing=spacing,iteration=iteration,
                      log2_upper_screen=bound,log2_eigenvalue_upper=log2(upper),
                      witness_range_log2=log2(np.max(vector))-log2(np.min(vector)),
                      directed=False,whole_code_certificate=False,mixer=mixer,state=state,adaptive=adaptive)
        vector=following/np.max(following)
        vector=np.maximum(vector,1e-250)
    return best


if __name__=='__main__':
    import argparse,json
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rows',type=int,nargs='+',default=[16,128,512,1024,2048,4096])
    parser.add_argument('--spacing',type=float,default=1.0)
    parser.add_argument('--adaptive',action='store_true')
    parser.add_argument('--state',type=int,default=4)
    parser.add_argument('--mixer',choices=('uniform_pre','feedback_refresh'),default='uniform_pre')
    parser.add_argument('--mark',type=F)
    parser.add_argument('--alphas',type=F,nargs='+',default=(F(1,32),F(1,16),F(1,8),F(1,4),F(1,2)))
    parser.add_argument('--factors',type=float,nargs='+',default=(1,1.25,1.5))
    args=parser.parse_args()
    for rows in args.rows:
        best=None
        z=F(63,64) if rows<=16 else F(7,8) if rows<=128 else F(1,2) if rows<=512 else F(1,8) if rows<=1024 else F(1,64) if rows<=2048 else F(1,256)
        z=args.mark or z
        for alpha in args.alphas:
            for factor in args.factors:
                p=rows/8192*factor
                if p>=1: continue
                candidate=screen(rows,z,alpha,p,spacing=args.spacing,mixer=args.mixer,state=args.state,adaptive=args.adaptive)
                if best is None or candidate['log2_upper_screen']<best['log2_upper_screen']:
                    best=candidate
            print(json.dumps(dict(completed_alpha=str(alpha),best=best)),flush=True)
