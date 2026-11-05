"""Extra independent-slot smoothing for the fixed actual-profile screen.

At markers u,v and rank mark y, write W=1+u+v and B=1+y**2*u+y*v.
The quiet-half polynomial P dominates B**16 coefficientwise, because
sum_b min(4,ell_b)<=2f+s. Its normalized law is therefore a mixture:
with probability theta=B**16/P it has16 independent base slots; otherwise
it has a residual law. Each base slot is itself a mixture of the free-slot
law (1,u,v)/W, with probability p=y**2*W/B, and a residual law, since
 B=y**2*W+(1-y**2)+v*(y-y**2).

Consequently the166 quiet halves contain M extra independent free-slot
draws, with generating function
 E[t**M]=[1-theta+theta*(1-p+p*t)**16]**166.
The16 slots in one half share its base-component coin; M is not an
ordinary binomial variable. Conditional on all mixture choices, the extra
draws and1440 original free slots are independent of the residual law.
Convolution cannot increase maximum point mass.

Let a_m be the maximum atom of Multinomial(1440+m;(u,v,1)/W). For a cutoff
m>=1 and delta>=Pr[M<m], monotonicity gives the positive bound
 a_m+delta*(a_0-a_m).
The scalar Chernoff bound uses t**(-(m-1))*E[t**M], with0<t<1. This helper
tests only cutoffs0,800,1000,1200 and optimizes the scalar mark by bounded
bisection. Multinomial modes satisfy all exchange inequalities, verified
with exact integer marker weights; their values are still FLOAT here.

The screen retains the three fixed rank marks, actual source-profile
counts, common mark across both matchings, and complete-profile minima.
It is not a certificate. One48-second stdout-only job,20GiB disk guard,
estimated peak below128MiB. No frozen source edits or generated data.
"""

from fractions import Fraction
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode=True
PIN="ec52c2706499e7722482283fc194534362f70ebc058a749e9f87c21e433178eb"
CUTOFFS=(0,800,1000,1200)


def load_counts(deadline):
    path=Path(__file__).with_name("cheap_lane_quiet_profile_counts.py")
    if hashlib.sha256(path.read_bytes()).hexdigest()!=PIN:
        raise RuntimeError("frozen profile-count wrapper hash mismatch")
    spec=importlib.util.spec_from_file_location("frozen_quiet_profile_counts",path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    marks,screen=module.load_marks(deadline)
    return module,marks,screen


def exact_mode(u,v,n):
    """Return a mode; every single-coordinate exchange is checked exactly."""
    ratios=[float(x).as_integer_ratio() for x in (u,v,1.)]
    denominator=max(d for _,d in ratios)
    weights=[a*(denominator//b) for a,b in ratios]
    total=sum(weights)
    counts=[n*w//total for w in weights]
    def increment():
        best=0
        for i in (1,2):
            if weights[i]*(counts[best]+1)>weights[best]*(counts[i]+1):
                best=i
        counts[best]+=1
    for _ in range(n-sum(counts)):
        increment()
    for _ in range(32):
        change=next(((i,j) for i in range(3) for j in range(3)
                     if i!=j and counts[i]>0 and counts[i]*weights[j]>(counts[j]+1)*weights[i]),None)
        if change is None:
            break
        i,j=change; counts[i]-=1; counts[j]+=1
    else:
        counts=[0,0,0]
        for _ in range(n):
            increment()
    assert sum(counts)==n
    assert all(counts[i]*weights[j]<=(counts[j]+1)*weights[i]
               for i in range(3) for j in range(3) if i!=j and counts[i]>0)
    value=math.lgamma(n+1)-sum(math.lgamma(k+1) for k in counts)
    value+=sum(k*math.log(w/total) for k,w in zip(counts,weights) if k)
    return value,counts


def tail(theta,p,cutoff):
    """Bound Pr[M<cutoff]; the returned t is only a FLOAT witness."""
    target=cutoff-1
    if target>=166*16*theta*p:
        return 1.,1.
    low,high=0.,1.
    for _ in range(56):
        t=(low+high)/2
        base=1-p+p*t
        denominator=1-theta+theta*base**16
        mean=166*16*theta*p*t*base**15/denominator
        if mean<target:
            low=t
        else:
            high=t
    t=(low+high)/2
    logbound=166*math.log1p(theta*((1-p+p*t)**16-1))-target*math.log(t)
    return min(1.,math.exp(logbound)),t


def strengthen(row,rows,y,screen):
    u,v=row["u"],row["v"]
    active=[i for i,x in enumerate((u,v)) if x]
    logp,_,_=screen.polynomial(rows,[math.log((u,v)[i]) for i in active],active)
    b=1+y*y*u+y*v
    theta=min(1.,math.exp(16*math.log(b)-logp))
    p=y*y*(1+u+v)/b
    loga0,mode0=exact_mode(u,v,1440); a0=math.exp(loga0)
    best={"cutoff":0,"delta":0.,"t":None,"mode_counts":mode0,"log2_cap":loga0/screen.LN2}
    for cutoff in CUTOFFS[1:]:
        delta,t=tail(theta,p,cutoff)
        logam,mode=exact_mode(u,v,1440+cutoff)
        am=math.exp(logam)
        assert am<=a0*(1+1e-10)
        cap=am+delta*max(0.,a0-am)
        value=math.log2(cap)
        if value<best["log2_cap"]:
            best={"cutoff":cutoff,"delta":delta,"t":t,"mode_counts":mode,"log2_cap":value}
    return {**row,"theta":theta,"free_slot_probability":p,"expected_extra_slots":166*16*theta*p,
            "old_mode_log2":loga0/screen.LN2,"smoothing":best,
            "old_regional_log2":row["regional_plain_log2"]+loga0/screen.LN2,
            "new_regional_log2":row["regional_plain_log2"]+best["log2_cap"]}


def tests():
    coefficient_checks=mixture_checks=mode_checks=tail_checks=0
    for y in (Fraction(3,4),Fraction(2,3),Fraction(3,5)):
        for f in range(17):
            for s in range(17-f):
                moment=sum(Fraction(math.comb(s,j),2**s)*y**(min(4,f+j)+min(4,f+s-j))
                           for j in range(s+1))
                assert moment>=y**(2*f+s)
                coefficient_checks+=1
        u,v=Fraction(1,7),Fraction(1,5); w=1+u+v; b=1+y*y*u+y*v
        p=y*y*w/b; residual=(1-y*y)+v*(y-y*y)
        assert 0<p<1 and residual>0
        for a,bcoord,rr in zip((1,u,v),(1,y*y*u,y*v),(1-y*y,0,v*(y-y*y))):
            assert p*a/w+(1-p)*rr/residual==bcoord/b
            mixture_checks+=1
        # Two toy quiet halves, each with three shared-coin slots.
        theta=Fraction(7,8)
        half=[theta*math.comb(3,k)*p**k*(1-p)**(3-k) for k in range(4)]
        half[0]+=1-theta
        law=[sum(half[i]*half[k-i] for i in range(4) if 0<=k-i<4) for k in range(7)]
        assert sum(law)==1
        for t in (Fraction(1,2),Fraction(3,4)):
            pgf=(1-theta+theta*(1-p+p*t)**3)**2
            assert pgf==sum(a*t**k for k,a in enumerate(law))
            for cutoff in range(1,7):
                assert sum(law[:cutoff])<=t**(-(cutoff-1))*pgf
                tail_checks+=1
    for n in range(1,9):
        for u,v in ((Fraction(0),Fraction(1,4)),(Fraction(1,3),Fraction(0)),
                    (Fraction(1,3),Fraction(1,4)),(Fraction(2),Fraction(3))):
            got,_=exact_mode(float(u),float(v),n)
            exact=max(Fraction(math.comb(n,f)*math.comb(n-f,s))*u**f*v**s/(1+u+v)**n
                      for f in range(n+1) for s in range(n-f+1))
            assert abs(got-math.log(float(exact)))<1e-12
            mode_checks+=1
    return {"coefficientwise_base_checks":coefficient_checks,"slot_mixture_checks":mixture_checks,
            "exact_compound_tail_checks":tail_checks,"exact_small_mode_checks":mode_checks}


def main():
    if len(sys.argv)!=1:
        raise ValueError("fixed smoothing screen accepts no parameters")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on the actual source volume")
    started=time.monotonic(); counts,marks,screen=load_counts(started+48)
    checked=tests(); counter=counts.DeletionCounts(screen.load_counts()); data=[]
    for name,y in marks.MARKS:
        rows=screen.terms(y); regional=[]
        for full in range(129):
            regional.append(strengthen(screen.optimize(rows,full),rows,y,screen))
            screen.guard()
        alpha=-math.log2(y)/8
        common=-32768*math.log2(y)+alpha*(8192*math.log2(256/255)-math.log2(255))+256*math.log2(85)
        data.append((name,common,regional))
    totals={"old":-math.inf,"new":-math.inf}; worst=[]; positive=negative=profiles=0
    usage={name:0 for name,_ in marks.MARKS}
    for first in range(129):
        for second in range(129):
            if first+second>255:
                continue
            count=counter.profile(4096,256,first,second)
            if not count:
                continue
            profiles+=1
            for kind in ("old","new"):
                score,name=min((common+math.log2(count)+128*(regional[first][kind+"_regional_log2"]
                                +regional[second][kind+"_regional_log2"]),name)
                               for name,common,regional in data)
                totals[kind]=screen.logadd(totals[kind],score)
                if kind=="new":
                    positive+=score>=0; negative+=score<0; usage[name]+=1
                    marks.keep_worst(worst,(score,first,second,256-first-second,name))
        screen.guard()
    indices={name:set() for name,_ in marks.MARKS}
    for _,f0,f1,_,name in worst:
        indices[name].update((f0,f1))
    print(json.dumps({"scope":"FLOAT actual size256 quiet-window union with compound free-slot smoothing",
        "whole_code_certificate":False,"directed":False,"profile_count":profiles,
        "log2_union_bounds":totals,"positive_profiles":positive,"negative_profiles":negative,
        "mark_usage":usage,"cutoff_candidates":CUTOFFS,
        "worst_profiles":[{"log2_bound":value,"full0":f0,"full1":f1,"runs":c,"mark":name}
                           for value,f0,f1,c,name in sorted(worst,reverse=True)],
        "regional_summary":[{"mark":name,
            "cutoff_usage":{str(c):sum(r["smoothing"]["cutoff"]==c for r in regional) for c in CUTOFFS},
            "minimum_extra_smoothing_bits":min(r["old_mode_log2"]-r["smoothing"]["log2_cap"] for r in regional),
            "maximum_extra_smoothing_bits":max(r["old_mode_log2"]-r["smoothing"]["log2_cap"] for r in regional),
            "worst_regional_witnesses":[regional[i] for i in sorted(indices[name])]}
            for name,_,regional in data],
        "tests":checked,"count_wrapper_sha256":PIN,
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "seconds":time.monotonic()-started,"generated_bytes":0,"estimated_peak_bytes":128*(1<<20)
    },indent=2))


if __name__=="__main__":
    main()
