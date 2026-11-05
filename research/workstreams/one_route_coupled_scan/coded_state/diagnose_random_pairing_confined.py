"""FLOAT screen of constant-support histories with actual RS row counts.

No certificate. Serial <=30s selected screen, no generated data. The state
is zero before a region-aligned suffix and has one fixed nonempty support
through its free terminal state. Pairings remain independent uniform.
"""
import math
import time
from functools import lru_cache
import numpy as np
from scipy.optimize import minimize
from diagnose_random_pairing import crossing_law, log2_fraction


def add2(a,b):
    top=max(a,b)
    return top+math.log2(2**(a-top)+2**(b-top))


def choose2(n,k):
    return (math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1))/math.log(2)


@lru_cache(maxsize=513)
def crossing_logs(n,w):
    law=crossing_law(n,w)
    return np.array(list(law)),np.array([log2_fraction(x) for x in law.values()])


def score(R,w,s,logit,log_y,n=512,M=8192,field_bits=32,threshold=None):
    """Positive Chernoff upper bound in bits, including row-set union.

    This includes all matchings, not only noncrossing matchings. It includes
    no union over R,w. p is a proof tilt, not a construction parameter.
    """
    assert 128<R<=256 and 1<=w<=n and w<=s<=M and log_y<=0
    if threshold is None:threshold=460*(1<<21)//1000
    steps=R*M//n
    if s==M:lp,lone=0.0,-math.inf
    else:lp=-add2(0,-logit);lone=-add2(0,logit)
    la=lp-field_bits
    lb=add2(lone,la)
    log_pi=choose2(M,s)+s*lp+(M-s)*lone if s<M else 0.0
    ks,probs=crossing_logs(n,w)
    terms=probs-field_bits*ks+(n-w-ks)*lb
    top=float(np.max(terms))
    moment=top+math.log2(float(np.exp2(terms-top).sum()))
    first=w*(field_bits+log_y+la)+(n-w)*lb
    repeat=w*(-field_bits+add2(0,field_bits+log_y))+moment
    result=(1+choose2(n,w)+choose2(M,s)+field_bits*(R-128)*s
            -R*log_pi-threshold*log_y+first+(steps-1)*repeat)
    return result


def selected(R,w,s):
    if s<w:return dict(R=R,w=w,s=s,excluded='initial pins exceed active rows')
    starts=[(0,-4),(24,-4),(-8,-16),(40,-16)]
    answers=[]
    for x in starts:
        result=minimize(lambda a:score(R,w,s,*a),x,method='L-BFGS-B',
                        bounds=((-60,60),(-48,0)),
                        options=dict(maxiter=100,ftol=1e-12))
        answers.append(result)
    best=min(answers,key=lambda a:a.fun)
    return dict(R=R,w=w,s=s,upper_bits=round(float(best.fun),3),
                logit_p=round(float(best.x[0]),5),log2_y=round(float(best.x[1]),5),
                float_only=True)


def main():
    began=time.monotonic()
    for R in (192,256):
        for w,s in ((2,2),(32,32),(128,128),(256,256),(256,4096),
                    (256,8192),(384,4096),(512,8192)):
            print(selected(R,w,s),flush=True)
            if time.monotonic()-began>25:raise TimeoutError('selected screen budget')
    print(dict(seconds=time.monotonic()-began,whole_code_certificate=False))


if __name__=='__main__':main()
