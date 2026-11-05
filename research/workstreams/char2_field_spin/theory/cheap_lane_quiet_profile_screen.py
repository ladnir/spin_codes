"""Bounded FLOAT screen for all paired profiles at actual support size256.

Each profile is (f0,f1), where fi counts full pairs in matching i and
si=256-2fi. The cyclic run count is256-f0-f1. Source multiplicities use
the pinned joint-profile upper bound, never an artificial repeated-full
profile. A44-packet window is allowed in each region;83 quiet packets
outside it are selected, excluding the last packet of every region.

At y=3/4, one quiet half has the positive polynomial
 P(u,v)=sum_(f+s<=16) C(16,f)C(16-f,s)u**f v**s /2**s
           *sum_j C(s,j)y**(min(4,f+j)+min(4,f+s-j)).
For a matching histogram (f,s), the regional Cauchy bound is
 P(u,v)**166*(1+u+v)**1440 / (u**f*v**s*M(4096;f,s)).
Positive markers are optimized by a bounded convex Newton method.
A zero marker is used exactly only when its required coefficient is zero.
Any finite positive marker gives an upper bound, regardless of convergence.

The smoothed variant multiplies this bound by the largest point mass of
Multinomial(1440;u,v,1)/(1+u+v). After normalization, the quiet polynomial
and this multinomial describe independent integer-valued random vectors.
The point mass of their sum is at most the multinomial's maximum point
mass. The mode is computed over every possible full count and the exact
conditional-binomial mode candidates. No heavy-occupancy tail is omitted.

The source relaxation is2*(4/3)**32768*(256/255)**512. The final sum includes
85**256 window choices and all ordered compatible profiles. Values are
FLOAT diagnostics, not certificates. No encoder or frozen helper changes.
One48-second stdout-only run;20GiB free guard;estimated peak below64MiB.
"""

from fractions import Fraction
from functools import lru_cache
import hashlib
import heapq
import importlib.util
import json
import math
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode=True
PIN="bb1700278361343325007cf558751af2db6df19dbf6c2090502941f1a95397be"
LN2=math.log(2)
DEADLINE=0.0


def guard():
    if time.monotonic()>DEADLINE:
        raise RuntimeError("bounded48-second invocation exhausted")


def load_counts():
    path=Path(__file__).with_name("hd_shifted_pairs_count.py")
    if hashlib.sha256(path.read_bytes()).hexdigest()!=PIN:
        raise RuntimeError("frozen profile-count helper hash mismatch")
    spec=importlib.util.spec_from_file_location("frozen_shifted_pair_counts",path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    # These caches change only this process, not the frozen source.
    module.one_pairing_count=lru_cache(None)(module.one_pairing_count)
    module.circular_run_count=lru_cache(None)(module.circular_run_count)
    return module


def terms(y):
    result=[]
    for f in range(17):
        for s in range(17-f):
            moment=sum(math.comb(s,j)*y**(min(4,f+j)+min(4,f+s-j))
                       for j in range(s+1))/2**s
            coefficient=math.comb(16,f)*math.comb(16-f,s)*moment
            result.append((f,s,math.log(coefficient)))
    return result


def polynomial(rows,logs,active):
    values=[(f,s,c+sum(logs[k]*((f,s)[axis]) for k,axis in enumerate(active)))
            for f,s,c in rows if all((f,s)[axis]==0 for axis in (0,1) if axis not in active)]
    high=max(c for _,_,c in values)
    weights=[math.exp(c-high) for _,_,c in values]
    total=sum(weights)
    means=[sum(w*(f,s)[axis] for w,(f,s,_) in zip(weights,values))/total for axis in active]
    covariance=[[sum(w*(f,s)[a]*(f,s)[b] for w,(f,s,_) in zip(weights,values))/total
                 -means[i]*means[j] for j,b in enumerate(active)] for i,a in enumerate(active)]
    return high+math.log(total),means,covariance


def objective(rows,logs,active,target):
    lp,means,covariance=polynomial(rows,logs,active)
    odds=[math.exp(x) for x in logs]
    total=1+sum(odds); probabilities=[x/total for x in odds]
    value=166*lp+1440*math.log(total)-sum(t*x for t,x in zip(target,logs))
    gradient=[166*m+1440*p-t for m,p,t in zip(means,probabilities,target)]
    hessian=[[166*covariance[i][j]+1440*((probabilities[i] if i==j else 0)
               -probabilities[i]*probabilities[j]) for j in range(len(active))] for i in range(len(active))]
    return value,gradient,hessian


def optimize(rows,full):
    singles=256-2*full; empty=4096-full-singles
    active=[axis for axis,t in enumerate((full,singles)) if t]
    target=[(full,singles)[axis] for axis in active]
    logs=[math.log(t/empty) for t in target]
    converged=False; evaluations=0
    for iteration in range(30):
        guard()
        value,g,h=objective(rows,logs,active,target); evaluations+=1
        if max(abs(x) for x in g)<1e-9:
            converged=True
            break
        if len(active)==1:
            step=[g[0]/h[0][0]]
        else:
            determinant=h[0][0]*h[1][1]-h[0][1]*h[1][0]
            if determinant<=0:
                break
            step=[(g[0]*h[1][1]-g[1]*h[0][1])/determinant,
                  (g[1]*h[0][0]-g[0]*h[1][0])/determinant]
        accepted=False
        for backoff in range(16):
            scale=2**(-backoff)
            proposed=[max(-60,min(20,x-scale*d)) for x,d in zip(logs,step)]
            candidate,_,_=objective(rows,proposed,active,target); evaluations+=1
            if candidate<=value-1e-4*scale*sum(a*b for a,b in zip(g,step)):
                logs=proposed; accepted=True; break
        if not accepted:
            break
    value,gradient,_=objective(rows,logs,active,target)
    normalizer=math.lgamma(4097)-math.lgamma(full+1)-math.lgamma(singles+1)-math.lgamma(empty+1)
    odds=[0.,0.]
    for axis,x in zip(active,logs):
        odds[axis]=math.exp(x)
    return {"full":full,"singles":singles,"u":odds[0],"v":odds[1],
            "regional_plain_log2":(value-normalizer)/LN2,"converged":converged,
            "maximum_gradient":max(abs(x) for x in gradient),"evaluations":evaluations}


def binomial_mode_log(n,p):
    if p==0 or p==1 or n==0:
        return 0.
    center=min(n,int(math.floor((n+1)*p)))
    return max(math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1)
               +k*math.log(p)+(n-k)*math.log1p(-p)
               for k in {center,max(0,center-1),min(n,center+1)})


def multinomial_mode_log(u,v,n=1440):
    if u==0:
        return binomial_mode_log(n,v/(1+v))
    if v==0:
        return binomial_mode_log(n,u/(1+u))
    pf=u/(1+u+v); conditional=v/(1+v)
    logpf=math.log(pf); logrest=math.log1p(-pf)
    logfull=n*logrest; best=-math.inf
    for k in range(n+1):
        best=max(best,logfull+binomial_mode_log(n-k,conditional))
        if k<n:
            logfull+=math.log(n-k)-math.log(k+1)+logpf-logrest
    return min(0.,best)


def logadd(a,b):
    if a==-math.inf:
        return b
    hi=max(a,b)
    return hi+math.log2(1+2**(-abs(a-b)))


def checks(rows):
    count=0
    neutral=terms(1.)
    for u,v in ((0.,.1),(.05,0.),(.04,.07)):
        active=[axis for axis,x in enumerate((u,v)) if x]
        value,_,_=polynomial(neutral,[math.log((u,v)[axis]) for axis in active],active)
        assert abs(value-16*math.log1p(u+v))<1e-12
        count+=1
    for f,s,logc in rows:
        if s==0:
            assert abs(logc-(math.log(math.comb(16,f))+2*min(4,f)*math.log(.75)))<1e-12
            count+=1
    for n in (1,2,4,6):
        for u,v in ((Fraction(0),Fraction(1,3)),(Fraction(1,4),Fraction(0)),
                    (Fraction(1,4),Fraction(1,3))):
            exact=max(Fraction(math.comb(n,f)*math.comb(n-f,s))*u**f*v**s/(1+u+v)**n
                      for f in range(n+1) for s in range(n-f+1))
            assert abs(multinomial_mode_log(float(u),float(v),n)-math.log(float(exact)))<1e-12
            count+=1
    return count


def main():
    global DEADLINE
    if len(sys.argv)!=1:
        raise ValueError("fixed profile screen accepts no parameters")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on the actual source volume")
    started=time.monotonic(); DEADLINE=started+48
    counts=load_counts(); rows=terms(.75); test_count=checks(rows)
    regional=[]
    for full in range(129):
        row=optimize(rows,full)
        row["multinomial_mode_log2"]=multinomial_mode_log(row["u"],row["v"])/LN2
        row["regional_smoothed_log2"]=row["regional_plain_log2"]+row["multinomial_mode_log2"]
        regional.append(row); guard()
    common=1+32768*math.log2(4/3)+512*math.log2(256/255)+256*math.log2(85)
    modes=("plain","smoothed")
    total={mode:-math.inf for mode in modes}
    worst={mode:[] for mode in modes}
    signs={mode:{"positive":0,"negative":0,"zero":0} for mode in modes}
    selected_runs={c:{"profiles":0,**{mode:{"minimum":math.inf,"maximum":-math.inf}
                                    for mode in modes}} for c in (1,4,16,64,128,256)}
    profile_count=0
    for f0 in range(129):
        for f1 in range(129):
            if f0+f1>255:
                continue
            multiplicity=counts.profile_count_upper(4096,256,f0,f1)
            if multiplicity==0:
                continue
            c=256-f0-f1; profile_count+=1
            for mode in modes:
                score=common+math.log2(multiplicity)+128*(regional[f0]["regional_"+mode+"_log2"]
                                                        +regional[f1]["regional_"+mode+"_log2"])
                total[mode]=logadd(total[mode],score)
                signs[mode]["positive" if score>0 else "negative" if score<0 else "zero"]+=1
                item=(score,f0,f1,c)
                if len(worst[mode])<6:
                    heapq.heappush(worst[mode],item)
                elif item>worst[mode][0]:
                    heapq.heapreplace(worst[mode],item)
                if c in selected_runs:
                    selected_runs[c][mode]["minimum"]=min(selected_runs[c][mode]["minimum"],score)
                    selected_runs[c][mode]["maximum"]=max(selected_runs[c][mode]["maximum"],score)
            if c in selected_runs:
                selected_runs[c]["profiles"]+=1
        guard()
    indices={0,64,127,128}
    for heap in worst.values():
        for _,f0,f1,_ in heap:
            indices.update((f0,f1))
    print(json.dumps({"scope":"FLOAT actual size256 profile quiet-window diagnostic",
        "whole_code_certificate":False,"directed":False,"source_size":256,
        "ordered_profile_count":profile_count,"rank_mark":"3/4",
        "log2_sum_upper_screen":total,"profile_sign_counts":signs,
        "worst_profiles":{mode:[{"log2_upper_screen":x,"full0":a,"full1":b,"runs":c}
                                for x,a,b,c in sorted(heap,reverse=True)] for mode,heap in worst.items()},
        "selected_run_counts":selected_runs,
        "selected_regional_markers":[regional[i] for i in sorted(indices)],
        "optimizer_converged_profiles":sum(r["converged"] for r in regional),
        "maximum_final_gradient":max(r["maximum_gradient"] for r in regional),
        "tests":test_count,"count_helper_sha256":PIN,
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "seconds":time.monotonic()-started,"generated_bytes":0,
        "estimated_peak_bytes":64*(1<<20)},indent=2))


if __name__=="__main__":
    main()
