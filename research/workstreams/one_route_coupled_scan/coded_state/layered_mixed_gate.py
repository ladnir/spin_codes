"""Bounded integer diagnostics for mixed layered-inner histories.

Selected screens are not whole-family certificates. No generated files.
"""
import time
from math import comb,isqrt
import numpy as np
from layered_route_gate import graph_envelope
from layered_one_sided_gate import head_envelope


M,b,N=8192,512,2**21
LOG64=[0]+[(a**64).bit_length()-1 for a in range(1,1025)]


def choose_log2(n,k):
    """Integer upper bound via one term of (p+(1-p))^n, dyadic p."""
    if k==0 or k==n:return 0
    a0=1024*k//n
    return min(n,min((640*n-k*LOG64[a]-(n-k)*LOG64[1024-a]+63)//64
                     for a in {max(1,min(1023,a0)),max(1,min(1023,a0+1))}))


def ordered_head_envelope(R):
    # For every cutoff k, sum_v (head_degree(v)-k)_+ concentrates.
    # R*M/2 without-replacement reveals, conditional range <=1 per reveal.
    n=R*M//2;num=7*n*80
    slack=isqrt(num//20);slack+=20*slack*slack<num
    s=np.arange(M+1,dtype=np.int64)
    best=np.minimum(R*s,R*M//2)
    numerator=tail=0;den=1<<R
    for k in range(R-1,-1,-1):
        tail+=comb(R,k+1);numerator+=tail
        mean=(M*numerator+den-1)//den
        best=np.minimum(best,k*s+mean+slack)
    return best


def orientation_floor(R,f,births):
    """Necessary r and l floor from the full nine-bit wiring cycle."""
    T=R*M//b-1;P=T*b//2;remaining=P-f;weighted=0
    for weight in range(10):
        take=min(remaining,T*comb(9,weight))
        weighted+=weight*take;remaining-=take
    assert remaining==0
    return max(0,(weighted-20*(births+b//2)+8)//9)


def occupation_dual(T,a,h,m,bits=9):
    """Exact min of a*r-h*H+m*births, relaxing spatial support constraints.

    One two-state temporal path per coordinate; min-plus period powering.
    Start and terminal hole indicators are free. Integer arithmetic only.
    """
    labels=np.arange(1<<bits,dtype=np.int64)
    def multiply(A,B):
        return np.minimum(A[:,:,:,None]+B[:,None,:,:],
                          np.iinfo(np.int64).max//4).min(axis=2)
    identity=np.full((1<<bits,2,2),10**12,dtype=np.int64)
    identity[:,0,0]=identity[:,1,1]=0
    period=identity.copy();prefix=identity.copy()
    for bit in range(bits):
        cost=a*((labels&(1<<bit))==0)-h
        step=np.zeros_like(identity)
        step[:,0,1]=step[:,1,1]=cost;step[:,1,0]=m
        period=multiply(period,step)
        if bit<T%bits:prefix=multiply(prefix,step)
    result=identity.copy();power=T//bits
    while power:
        if power&1:result=multiply(result,period)
        period=multiply(period,period);power//=2
    result=multiply(result,prefix)
    return int(result.min(axis=(1,2)).sum())


def occupation_markers(R):
    T=R*M//b-1
    return [(16,h,m,occupation_dual(T,16,h,m))
            for h in (4,8,12,16) for m in (4,8,16,32,64)]


def occupation_floor(R,f,births,markers):
    P=(R*M-b)//2;H=P-f
    return max(orientation_floor(R,f,births),
               max((h*H+value-m*births+a-1)//a for a,h,m,value in markers))


def envelopes(R):
    return (np.minimum(np.asarray(head_envelope(M,R),dtype=np.int64),ordered_head_envelope(R)),
            np.asarray(graph_envelope(M,R),dtype=np.int64),
            np.arange(M+1,dtype=np.int64))


def mixed_dimension(R,P,z,r,f,e,env):
    H,E,s=env
    l=P-z-r-f
    D=l+2*r+2*f-2*e
    degree=np.minimum.reduce((np.full(M+1,D,dtype=np.int64),
                              P-z-e+E,H+r+f-e,R*s))
    return min(D//2,max(0,R*M-N//2),max(0,int((degree-128*s).max())))


def no_empty_point(R,f,r,e,env,permille=400,births=None,dimension=None):
    u=R*M;P=(u-b)//2;W=permille*N//1000
    assert 0<=f<=P and 0<=r<=P-f and 0<=e<=f
    l=P-f-r;V=P+f
    d=mixed_dimension(R,P,0,r,f,e,env) if dimension is None else dimension
    pattern=choose_log2(P,f)+choose_log2(P-f,r)+choose_log2(f,e)
    if births is not None:
        positions=u-2*b;changes=min(positions,2*births+b//2)
        change_bits=positions if changes*2>=positions else (changes+1).bit_length()+choose_log2(positions,changes)
        pattern=min(pattern,b+change_bits+choose_log2(f,e))
    C=max(0,V-W,3*u//4-b-min(r,l)-W)
    output=u+P+32*(W-P-r-f+e+d)+1
    state=pattern+V+32*(b+e-r-C+d)+1
    roots=max(e,u//4-b//2-f)
    forward=pattern+32*(max(0,u-N//2)-roots)+1
    if births is not None:
        C=max(C,V+births-b//2-W)
        edge_bits=V-births-b//2
        state=pattern+edge_bits+32*(b+e-r-C+d)+1
        roots_state=max(e,l+r+e-births-b//2)
        forward=min(forward,pattern+32*(max(0,u-N//2)-roots_state)+1)
        D=P+r+f-2*e
        roots_joint=max(roots_state,D+e+births-W-b)
        joint=pattern+edge_bits+32*(max(0,u-N//2)-roots_joint)+1
        forward=min(forward,joint)
    return dict(R=R,f=f,r=r,e=e,dimension=d,rank=C,
                output=output,state=state,forward=forward,
                births=births,best=min(output,state,forward))


def screen():
    started=time.monotonic();worst=None;count=0;failed=0;by_birth={};impossible=0
    for R in (129,160,168,176,192,208,224,240,256):
        P=(R*M-b)//2;env=envelopes(R);markers=occupation_markers(R)
        for f in sorted({P*j//16 for j in range(17)}):
            for r in sorted({(P-f)*j//8 for j in range(9)}):
                for e in sorted({f*j//8 for j in range(9)}):
                    W=400*N//1000
                    dim=mixed_dimension(R,P,0,r,f,e,env)
                    for j,births in enumerate(sorted({min(P-f,W//2)*j//4 for j in range(5)})):
                        floor=occupation_floor(R,f,births,markers)
                        if min(r,P-f-r)<floor:
                            impossible+=1
                            continue
                        value=no_empty_point(R,f,r,e,env,births=births,dimension=dim)
                        count+=1;failed+=value['best']>-100
                        if value['best']>-100:by_birth[j]=by_birth.get(j,0)+1
                        if worst is None or value['best']>worst['best']:worst=value
                    if time.monotonic()-started>30:raise RuntimeError('Screen cap')
    return dict(selected_classes=count,impossible_classes=impossible,failed_classes=failed,failed_birth_grid=by_birth,worst=worst,
                whole_code_certificate=False)


if __name__=='__main__':print(screen())
