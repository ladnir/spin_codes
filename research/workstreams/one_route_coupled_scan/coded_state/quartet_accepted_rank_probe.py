"""Global rank with accepted-map filters; floating diagnostic, NOT certificate.

Classify rows by covered degree>=129. Keep both superregular branch tests
inside the quartet count. Root-owned file-free180s/2GiB remote run only.
"""
from collections import Counter
from math import comb,lgamma,log,log2,exp,expm1,log1p
import os,resource,shutil,sys,threading,time
import numpy as np

N=1<<21;K=1<<20;M=8192;W=943718;T=4096
LN2=log(2)


def descriptors():
    terms=Counter()
    for u in range(5):
        incoming=(1<<u)-1
        for nxt in range(16):
            v=nxt.bit_count();births=(nxt&~incoming).bit_count();deaths=(incoming&~nxt).bit_count()
            survived=incoming&nxt;mandatory=incoming^nxt
            for slots in range(16):
                k=slots.bit_count()
                if not incoming:
                    if nxt&~slots:continue
                    optional=(1,0,0,0,0)
                else:
                    if u+(nxt|slots).bit_count()<5:continue
                    optional=[0]*5
                    for subset in range(16):
                        if subset&~survived:continue
                        if u+(mandatory|subset|slots).bit_count()>=5:optional[subset.bit_count()]+=1
                    optional=tuple(optional)
                    if not any(optional):continue
                terms[u,v,k,births,deaths,optional]+=comb(4,u)
    rows=list(terms)
    return (np.array([r[:5] for r in rows],dtype=int),
            np.array([r[5] for r in rows],dtype=float),
            np.array([terms[r] for r in rows],dtype=float))


DESC,OPTIONAL,MULT=descriptors()
BINOM512=np.array([(lgamma(513)-lgamma(u+1)-lgamma(513-u))/LN2 for u in range(513)])


def local(theta,logz,p,logcovered=0.,logoutside=0.):
    u,v,k,b,d=DESC.T
    a=2.0**(32*theta);z=2.0**logz;m=a*z;w=2.0**(32+logz/theta)
    optional=(OPTIONAL@np.array([w**i for i in range(5)]))**theta
    probability=(p**k)*((1-p)**(4-k))
    rank=a**np.where(u>0,-(4-k),0)
    covered=np.where(u>0,k,v);outside=np.where(u>0,4-k,0)
    values=MULT*optional*probability*rank*m**b*z**d*np.exp2(logcovered*covered+logoutside*outside)
    out=np.bincount(5*u+v,weights=values,minlength=25).reshape(5,5)
    # This row has a closed form: inactive input permits only H-slot births.
    # In particular the unique zero path has weight exactly one.
    out[0]=[comb(4,v)*(p*m*2**logcovered)**v for v in range(5)]
    return out


def log_expm1(x):
    if x>40:return x+log1p(-exp(-x))
    return log(expm1(x))


def envelope(theta,logz,p,logh,logcovered=0.,logoutside=0.):
    h=2**logh
    coeff=local(theta,logz,p,logcovered,logoutside)@np.array([h**v for v in range(5)])
    if min(coeff)<=0:return None
    terms=np.log(coeff);poly=np.array([0.])
    for _ in range(128):
        new=np.full(len(poly)+4,-np.inf)
        for u in range(5):new[u:u+len(poly)]=np.logaddexp(new[u:u+len(poly)],poly+terms[u])
        poly=new
    ratios=poly/LN2-BINOM512-logh*np.arange(513)
    maximum=float(ratios.max());worst=int(np.argmax(ratios))
    # Pin the first departure from zero, then use the all-source potential.
    logrow0=512*log1p(p*2**(32*theta+logz+logh+logcovered))
    if logrow0<=0 or maximum<=0:return None
    logsum=(log_expm1(logrow0)+log_expm1(T*maximum*LN2)-log_expm1(maximum*LN2))/LN2
    return logsum+max(0.,-512*logh),worst,maximum


def logbinom(n,k):return (lgamma(n+1)-lgamma(k+1)-lgamma(n-k+1))/LN2


def objective(hrows,theta,point):
    logz,logh,logit=point[:3]
    logcovered,logoutside=point[3:] if len(point)>3 else (0.,0.)
    if not (-64<=logz<=0 and -32<=logh<=32 and -48<=logit<=48):return float('inf'),None
    if not (0<=logcovered<=32 and -32<=logoutside<=0):return float('inf'),None
    p=1. if hrows==M else 1/(1+2**(-logit))
    value=envelope(theta,logz,p,logh,logcovered,logoutside)
    if value is None:return float('inf'),None
    logchoose=logbinom(M,hrows)
    logpi=0. if hrows==M else 256*(logchoose+hrows*log2(p)+(M-hrows)*log2(1-p))
    bound=theta-32*theta-W*logz+logchoose-32*theta*128*hrows-logpi+value[0]
    bound-=129*hrows*logcovered+128*(M-hrows)*logoutside
    return bound,(value[1],p,value[2])


def controls():
    # p=1 removes all support restrictions and all local live rank cost.
    theta=.25;zlog=-8;m=1;c=2**theta;z=2**zlog
    p1=local(theta,zlog,1.)
    for u in range(5):
        expected=np.zeros(5)
        for death in range(u+1):
            for birth in range(5-u):
                v=u-death+birth
                expected[v]+=comb(4,u)*comb(u,death)*comb(4-u,birth)*z**death*m**birth*c**(u-death)
        assert np.allclose(p1[u],expected,rtol=1e-13,atol=1e-15)
    p0=local(theta,zlog,0.)
    assert p0[0,0]==1 and np.all(p0[0,1:]==0)
    assert np.all(p0[1,:4]==0)  # Singleton input must spread to all4 state entries.
    assert abs(p0[1,4]-4*z**4)<1e-20  # All4 outputs forced nonzero as well.
    assert p0[4,0]==0  # Invertible A cannot kill a nonzero state without input.
    print(dict(accepted_kernel_controls=True,descriptor_count=len(MULT)),flush=True)


def scan(degrees=False):
    for hrows in (1,8,64,512,2048,4096,6144,8192):
        start=time.monotonic();best=(float('inf'),None)
        true_logit=0. if hrows==M else log2(hrows/(M-hrows))
        for theta in (1/32,1/16,1/8,1/4,1/2,1):
            for dense in (False,True):
                point=[-32*theta*(1 if dense else hrows/M),-2.,true_logit-16*theta]
                if degrees:point.extend([0.,0.])
                value,detail=objective(hrows,theta,point)
                for step in (4.,2.,1.,.5,.25):
                    for _ in range(2):
                        changed=False
                        for axis in range(len(point)):
                            if axis==2 and hrows==M:continue
                            for sign in (-1,1):
                                candidate=point.copy();candidate[axis]+=sign*step
                                trial,meta=objective(hrows,theta,candidate)
                                if trial<value:point=candidate;value=trial;detail=meta;changed=True
                        if not changed:break
                if value<best[0]:best=(value,(theta,point,detail))
        print(dict(eligible_rows=hrows,log2_diagnostic=round(best[0],3),
                   theta=best[1][0],parameters=best[1][1],details=best[1][2],
                   seconds=round(time.monotonic()-start,3),certificate=False),flush=True)


def accepted_update(theta,logz,p,logh,width=512,logcovered=0.,logoutside=0.):
    """Exact block-convolution formula, evaluated without directed rounding."""
    ld=np.longdouble
    matrix=np.asarray(local(theta,logz,p,logcovered,logoutside),dtype=ld)
    for u in range(5):
        for v in range(5):
            matrix[u,v]*=np.exp2(ld(logh)*(v-u))/comb(4,u)
    old=np.ones((1,1),dtype=ld)
    for slots in range(4,width+1,4):
        size=len(old);new=np.zeros((slots+1,slots+1),dtype=ld)
        denominator=ld(comb(slots,4))
        for u in range(5):
            weights=np.array([ld(comb(k,u)*comb(slots-k,4-u))/denominator
                              for k in range(u,u+size)],dtype=ld)
            weighted=old*weights[:,None]
            for v in range(5):
                new[u:u+size,v:v+size]+=weighted*matrix[u,v]
        old=new
    assert old[0,0]==1
    return old


def regional_scan():
    # The loader supplies the separately tested positive/renewal helpers.
    # Exact quota conditioning charges only positive-footprint regions.
    cases=((512,.125,-.25),
           (2048,.125,-1.),(4096,.125,-2.),(6144,.125,-3.))
    for hrows,theta,logz in cases:
        start=time.monotonic();pactual=hrows/M;a=2**(32*theta)
        p=pactual/(a-(a-1)*pactual);logh=-2.
        kernel=accepted_update(theta,logz,p,logh)
        positive,scale=positive_power16(kernel)
        zero_count=int(np.count_nonzero(positive==0))
        # A floating diagnostic only; upward floor makes lost tiny entries
        # pessimistic but is not a substitute for directed rounding.
        positive=np.maximum(positive,np.exp2(np.longdouble(-15000)))
        logchoose=logbinom(M,hrows)
        logpi=logchoose+hrows*log2(p)+(M-hrows)*log2(1-p)
        tail=occupied_tail(positive,scale-logpi)
        bound=theta-32*theta-W*logz+logchoose-32*theta*128*hrows+tail+max(0.,-512*logh)
        print(dict(regional_eligible_rows=hrows,theta=theta,logz=logz,p=p,
                   logpi=logpi,log2_diagnostic=float(bound),floored_entries=zero_count,
                   seconds=round(time.monotonic()-start,3),certificate=False),flush=True)


def tilted_regional_scan():
    for hrows,theta,logz in ((1,.5,-1/512),(512,.125,-.25),(2048,.125,-1.),(4096,.125,-2.)):
        start=time.monotonic();pactual=hrows/M;a=2**(32*theta);logh=-2.
        low=pactual/(a-(a-1)*pactual)
        odds0=low/(1-low);odds1=pactual/(1-pactual)
        logchoose=logbinom(M,hrows);upper_logs=None
        for fraction in (0.,1/3,2/3,1.):
            odds=odds0**(1-fraction)*odds1**fraction;p=odds/(1+odds)
            kernel=accepted_update(theta,logz,p,logh)
            positive,scale=positive_power16(kernel)
            positive=np.maximum(positive,np.exp2(np.longdouble(-15000)))
            logpi=logchoose+hrows*log2(p)+(M-hrows)*log2(1-p)
            logs=np.log2(positive)+scale-logpi
            upper_logs=logs if upper_logs is None else np.minimum(upper_logs,logs)
        scale=upper_logs.max();positive=np.exp2(upper_logs-scale)
        norm=positive.sum(axis=1).max();positive/=norm;scale+=np.log2(norm)
        tail=occupied_tail(positive,scale)
        bound=theta-32*theta-W*logz+logchoose-32*theta*128*hrows+tail+max(0.,-512*logh)
        print(dict(tilted_regional_rows=hrows,log2_diagnostic=float(bound),
                   theta=theta,logz=logz,tilts=4,seconds=round(time.monotonic()-start,3),
                   certificate=False),flush=True)


def main():
    if shutil.disk_usage('/tmp').free<20*(1<<30):raise RuntimeError('20GiB guard')
    resource.setrlimit(resource.RLIMIT_AS,(2*(1<<30),2*(1<<30)))
    timer=threading.Timer(175,lambda:os._exit(124));timer.daemon=True;timer.start()
    controls()
    if '--tilted-regional' in sys.argv:tilted_regional_scan()
    elif '--regional' in sys.argv:regional_scan()
    else:scan('--degrees' in sys.argv)
    timer.cancel()


if __name__=='__main__':main()
