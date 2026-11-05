"""Cap conditional bad-history probability before averaging pair geometry.

FLOAT proposals, unchanged construction. Serial <=30s, no data files.
"""
import time
import math
import shutil
import tempfile
from functools import lru_cache


def choose2(n,k):
    return (math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1))/math.log(2)


def sum2(values):
    top=max(values)
    return top+math.log2(sum(2**(x-top) for x in values))


@lru_cache(maxsize=513)
def crossing_logs(n,w):
    def matchings(a):return math.prod(range(1,a,2))
    total=matchings(n)
    ks=list(range(w%2,min(w,n-w)+1,2))
    numerators=[math.comb(w,k)*math.comb(n-w,k)*math.factorial(k)
                *matchings(w-k)*matchings(n-w-k) for k in ks]
    assert sum(numerators)==total
    return ks,[math.log2(x)-math.log2(total) for x in numerators]


def history_count(T,w,B,n=512):
    if not B:return choose2(n,w)
    degrees=range(min(w,n-w)+1)
    coefficients=[choose2(w,j)+choose2(n-w,j) for j in degrees]
    def at(x):
        terms=[a+j*x for j,a in enumerate(coefficients)];normal=sum2(terms)
        mean=sum(j*2**(a-normal) for j,a in enumerate(terms))
        return mean,normal
    low,high=-128.,128.
    for _ in range(64):
        middle=(low+high)/2
        if at(middle)[0]<B/T:low=middle
        else:high=middle
    x=(low+high)/2
    return choose2(n,w)+T*at(x)[1]-B*x


def profile(R,w,B,delta=460):
    assert w==256 and 129<=R<=256 and 0<=B<=(16*R-1)*256
    M=8192;n=512;T=16*R-1;W=delta*(1<<21)//1000
    area=(T+1)*w+B;intersection=T*w-B;zeros=max(0,area-W)
    if zeros>intersection:return dict(R=R,w=w,B=B,excluded=True)
    history=history_count(T,w,B)
    ks,law=crossing_logs(n,w)
    logm=max(float(p)-choose2(n-w,int(k)) for k,p in zip(ks,law))
    density=T*((n-w)+logm)
    # Enlarge the initial union to the whole batch, since the initial
    # zero state supplies no random crossing endpoints.
    padded_area=area+n-w
    cover=(1+padded_area/(R*M))/2
    law=[choose2(R,j)+j*math.log2(cover)+(R-j)*math.log2(1-cover) for j in range(R+1)]
    deficiency=[max(128-j,0) for j in range(R+1)]
    base=padded_area-128*M-B-zeros
    zero_choices=choose2(intersection,zeros)
    overhead=history+density+R*math.log2(M+1)
    def objective(r):
        return overhead+(r/32)*(1+zero_choices+32*base)+M*sum2([a+r*d for a,d in zip(law,deficiency)])
    def derivative(r):
        terms=[a+r*d for a,d in zip(law,deficiency)];normal=sum2(terms)
        return base+(1+zero_choices)/32+M*sum(d*2**(a-normal) for a,d in zip(terms,deficiency))
    low,high=0.,32.
    for _ in range(64):
        middle=(low+high)/2
        if derivative(middle)<0:low=middle
        else:high=middle
    optimum=min((0.,(low+high)/2,32.),key=objective)
    return dict(R=R,w=w,B=B,history_bits=round(history,3),
                pair_comparison_bits=round(density,3),coverage=round(cover,7),
                r=round(optimum,8),upper_bits=round(objective(optimum),3),
                float_only=True)


def minimize_convex(function,low,high,iterations=44):
    """Bounded FLOAT proposal only; no certification of an optimum."""
    left,right=low,high;ratio=(math.sqrt(5)-1)/2
    a=right-ratio*(right-left);b=left+ratio*(right-left)
    fa,fb=function(a),function(b)
    for _ in range(iterations):
        if fa<fb:
            right,b,fb=b,a,fa;a=right-ratio*(right-left);fa=function(a)
        else:
            left,a,fa=a,b,fb;b=left+ratio*(right-left);fb=function(b)
    x=min((low,high,a,b),key=function)
    return x,function(x)


def active_profile(R,w,B,s=None,delta=460):
    """Exact active-row set inside the capped comparison, FLOAT weights."""
    assert w==256 and (s is None or w<=s<=8192)
    M=8192;T=16*R-1;W=delta*(1<<21)//1000
    I=T*w-B;z=max(0,(T+1)*w+B-W)
    if z>I:return dict(R=R,w=w,B=B,s=s,excluded=True)
    A=(T+1)*w+B+512-w;p=A/(R*M)
    h=history_count(T,w,B)
    ks,law=crossing_logs(512,w)
    density=T*(256+max(value-choose2(256,k) for k,value in zip(ks,law)))
    overhead=h+density+R*math.log2(M+1)+(choose2(M,s) if s is not None else 0)
    output=1+choose2(I,z)
    def evaluate(logr,logt,detail=False):
        a=p*2**(logr-logt)+(1-p)/2;b=(1-p)/2
        if R==129:
            active=R*math.log2(a);logmark=math.inf
        elif R*a/(a+b)>=129:
            active=R*math.log2(a+b);logmark=0.
        else:
            mark=129*b/((R-129)*a);logmark=math.log2(mark)
            active=R*math.log2(a*mark+b)-129*logmark
        inactive=R*math.log2(p*2**(-logt)+(1-p)/2)
        if s is None:
            # Sum all row sets, including the harmless overcount S=empty.
            result=(overhead+logr*(output/32-B-z)+A*logt
                    +M*sum2([-128*logr+active,inactive]))
        else:
            result=(overhead+logr*(output/32-128*s-B-z)+A*logt
                    +s*active+(M-s)*inactive)
        return (result,logmark) if detail else result
    def inner(logr):return minimize_convex(lambda logt:evaluate(logr,logt),-64,64)
    logr,value=minimize_convex(lambda lr:inner(lr)[1],0,32)
    logt,value=inner(logr);_,logmark=evaluate(logr,logt,True)
    return dict(R=R,w=w,B=B,s=s,upper_bits=round(value,3),
                log2_r=round(logr,6),log2_t=round(logt,6),log2_marker=round(logmark,6),
                float_only=True)


def main():
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--active',action='store_true')
    parser.add_argument('--active-sum',action='store_true')
    args=parser.parse_args()
    assert shutil.disk_usage(tempfile.gettempdir()).free>=20*(1<<30),'disk guard'
    began=time.monotonic()
    if args.active_sum:
        for R in (192,256):
            for B in (0,8192,16384,32768,65536,131072,262144,450000):
                print(active_profile(R,256,B),flush=True)
                if time.monotonic()-began>25:raise TimeoutError('active sum budget')
        print(dict(seconds=time.monotonic()-began,whole_code_certificate=False))
        return
    if args.active:
        for R in (192,256):
            for B in (32768,131072):
                values=[]
                for s in (256,512,1024,2048,4096,6144,8192):
                    values.append(active_profile(R,256,B,s))
                    if time.monotonic()-began>25:raise TimeoutError('active-row screen budget')
                print(dict(R=R,B=B,selected_rows=values),flush=True)
        print(dict(seconds=time.monotonic()-began,whole_code_certificate=False))
        return
    for R in (192,256):
        for B in (0,256,4096,16384,32768,65536,131072,262144,450000):
            print(profile(R,256,B),flush=True)
            if time.monotonic()-began>25:raise TimeoutError('capped screen budget')
    print(dict(seconds=time.monotonic()-began,whole_code_certificate=False))


if __name__=='__main__':main()
