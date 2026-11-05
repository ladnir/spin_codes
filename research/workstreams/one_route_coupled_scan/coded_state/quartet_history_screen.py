"""Full-width FLOAT lower screens for one history-conditioned potential.

Not a certificate. Tests whether a chosen scalar-envelope expression is
already too large, and separates Maclaurin loss from exact kernel ratios.
Serial,90s,<64MiB,no files. The construction is unchanged.
"""
import math
import itertools
import os
from pathlib import Path
import shutil
import threading
import time
from functools import lru_cache


def logchoose(n,k):
    return (math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1))/math.log(2)


def logsum(values):
    peak=max(values)
    return peak+math.log2(sum(2**(value-peak) for value in values))


@lru_cache(maxsize=64)
def occupation(l):
    width=512;groups=width//4;powers=[[1]]
    for g in range(1,groups+1):
        previous=powers[-1];after=[0]*(len(previous)+3)
        for i,a in enumerate(previous):
            for j,b in enumerate((1,4,6,4)):after[i+j]+=a*b
        powers.append(after)
    result=[]
    for f in range(groups+1):
        degree=l-4*f;poly=powers[groups-f]
        if 0<=degree<len(poly):
            count=math.comb(groups,f)*poly[degree]
            if count:result.append((4*f,count,math.log2(count)))
    assert sum(count for k,count,log in result)==math.comb(width,l)
    return result


def envelope_sample(a,l,logb,slope,tau,h):
    d0=h[1];d1=h[2];d4=h[4]*2**(logb*slope)
    entries=occupation(l);best=(-math.inf,0)
    for j in range(a+1):
        mean=((512-a)*d0+(a-j)*d1+j*d4)/512
        value=logsum([log+k*math.log2(mean/h[1]) for k,count,log in entries])
        value+=(l-a+j)*math.log2(h[1])-j*math.log2(h[4])
        if value>best[0]:best=(value,j)
    return best


def exact_kernel_ratio_float(counts,l,logb,slope,tau,h):
    # Compute the elementary symmetric polynomial as a Poisson-binomial
    # distribution times product(1+d). This avoids intermediate overflow.
    coeff=[1.0]+[0.0]*512;used=0;scale=0.0
    ds=(h[1],h[2],h[3],h[4]*2**(logb*(3*slope-tau)),h[4]*2**(logb*slope))
    for d,count in zip(ds,counts):
        p=d/(1+d);scale+=count*math.log2(1+d)
        for _ in range(count):
            used+=1
            for k in range(used,0,-1):coeff[k]=(1-p)*coeff[k]+p*coeff[k-1]
            coeff[0]*=1-p
    assert abs(sum(coeff)-1)<1e-11
    values=[log+math.log2(coeff[k])+scale-logchoose(512,k)+(l-k)*math.log2(h[1])
            for k,count,log in occupation(l) if coeff[k]>0]
    assert values
    return logsum(values)-sum(n*math.log2(hi) for n,hi in zip(counts,h))


def trace_expression(gamma,logb,slope,tau,h):
    # 44 regions have 15 counts408 followed by407;212 have16 counts408.
    # These are fixed disjoint regional payment counts, totaling111572.
    transitions={(408,408):4007,(408,407):44,(407,408):44}
    value=logchoose(512,408)+408*math.log2(h[1])
    value+=sum(count*gamma[pair] for pair,count in transitions.items())
    value+=44*(logchoose(6527,435)-logchoose(8192,2100))
    value+=212*(logchoose(6528,436)-logchoose(8192,2100))
    rewards=(slope,2*slope,3*slope,slope+tau)
    value+=408*max(logb*reward-math.log2(hi) for reward,hi in zip(rewards,h[1:]))
    return value-logb*165556+logchoose(8192,2100)


def point_mass_upper(variance,refined=False):
    upper=min(1.0,math.sqrt(math.pi/(8*variance)))
    if refined:
        for theta in (.5,.75,1.,1.25):
            c=1-theta*theta/12
            candidate=1/math.sqrt(2*math.pi*variance*c)+math.exp(-variance*theta*theta*c/2)
            upper=min(upper,candidate)
    return upper


def coefficient_cauchy_upper(a,l,logb,slope,tau,h,return_pair=False,incoming=None,point_mass='coarse'):
    """FLOAT evaluation of a valid type-separated one-step upper formula.

    Each chosen positive z bounds the elementary symmetric coefficient.
    Maximize the product over the permitted active types before summing k.
    """
    g=h if incoming is None else incoming
    assert h[0]==g[0]==1
    ds=(h[1],h[2],h[3],h[4]*2**(logb*(3*slope-tau)),h[4]*2**(logb*slope))
    terms=[];concentrated=[]
    for k,count,log in occupation(l):
        lo=-80.;hi=80.
        for _ in range(60):
            mid=(lo+hi)/2;z=2**mid
            i=max(range(1,5),key=lambda i:(1+ds[i]*z)/g[i])
            derivative=(512-a)*ds[0]*z/(1+ds[0]*z)+a*ds[i]*z/(1+ds[i]*z)-k
            if derivative<0:lo=mid
            else:hi=mid
        x=(lo+hi)/2;z=2**x
        value=(512-a)*math.log2(1+ds[0]*z)
        value+=a*max(math.log2((1+ds[i]*z)/g[i]) for i in range(1,5))-k*x
        term=log-logchoose(512,k)+(l-k)*math.log2(h[1])+value
        terms.append(term)
        variances=[d*z/(1+d*z)**2 for d in ds]
        minimum_variance=(512-a)*variances[0]+a*min(variances[1:])
        correction=math.log2(point_mass_upper(minimum_variance,point_mass=='refined'))
        concentrated.append(term+correction)
    result=logsum(terms)
    return (result,logsum(concentrated)) if return_pair else result


def column_perron(logb,slope,tau,active=408):
    p=active/512;k=math.prod((active-i)/(512-i) for i in range(4))
    wa=2**(logb*(3*slope-tau));wl=2**(logb*slope)
    h=[1.]*5
    for _ in range(512):
        common=(1-p)*h[0]+(p-k)*h[1]
        after=[(1-p)*h[0]+p*h[1],common+k*h[2],common+k*h[3],
               common+k*wa*h[4],common+k*wl*h[4]]
        after=[value/after[0] for value in after]
        if max(abs(math.log(x/y)) for x,y in zip(after,h))<1e-13:return after
        h=after
    raise AssertionError('Perron iteration did not converge')


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(90,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:
        slope=707/2560;tau=381/512
        for label,logb in itertools.product(('full_run','column_perron'),(1,2,3,4,6,8)):
            h=[1.,1.,2**(logb*slope),2**(2*logb*slope),2**(logb*tau)] if label=='full_run' else column_perron(logb,slope,tau)
            envelope={};exact={};cauchy={};selected=[]
            for a,l in ((408,408),(408,407),(407,408)):
                value,j=envelope_sample(a,l,logb,slope,tau,h)
                envelope[a,l]=value
                candidates=sorted({0,a//2,a,j})
                histograms=[(512-a,a-k,0,0,k) for k in candidates]
                sampled=max(exact_kernel_ratio_float(n,l,logb,slope,tau,h) for n in histograms)
                assert sampled<=value+1e-6
                histograms.extend(((512-a,0,a,0,0),(512-a,0,0,a,0)))
                exact[a,l]=max(exact_kernel_ratio_float(n,l,logb,slope,tau,h) for n in histograms)
                cauchy[a,l]=coefficient_cauchy_upper(a,l,logb,slope,tau,h)
                assert exact[a,l]<=cauchy[a,l]+1e-6
                selected.append(j)
            print(dict(log2_base=logb,potential=label,
                       float_envelope_expression_lower=trace_expression(envelope,logb,slope,tau,h),
                       float_exact_kernel_sample_expression_lower=trace_expression(exact,logb,slope,tau,h),
                       float_cauchy_envelope_fixed_path_term=trace_expression(cauchy,logb,slope,tau,h),
                       envelope_selected_n4=selected,regional_quotas_exact=True,
                       paid_count_exact=True,input_count_path_fixed=True,
                       partial_updates_only_after_start=True,
                       directed_certificate=False,whole_code_distance=False),flush=True)
        print(dict(seconds=time.monotonic()-start,files_written=0),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
