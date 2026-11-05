"""FLOAT quiet-rank bounds for all binary state paths at actual support256.

The branch has at most K births (0->1 regular steps) and at least Q quiet
(0->0 regular) steps. For a fixed binary path, quiet equations are E-defined.
Their scalar outer count is bounded by beta**256*q**(-D)/255. Apply the
fractional power alpha before averaging routes, then union over binary paths.
No field-state mask multiplicities are introduced.

For one matching's full/single counts f,s=256-2f, let U=1+u+w, P be the
quiet16-slot half polynomial, and y=q**(-alpha). At positive markers,
 S=[[v*P**2/U**32,b],[1,1]],
 H=4096*log(U)-f*log(u)-s*log(w)-log multinomial(4096;f,s).
The regular regional coefficient is entrywise at most exp(H)*S**128.
The terminal region has exp(H)*S**127*1: the last packet has no feedback,
birth, or quiet mark. Actual even/odd matching order gives
 (S0**128*S1**128)**127*S0**128*S1**127*1.
The scalar exp(128*(H0+H1)) and source factor are paid once. Marks satisfy
0<b<=1,v>=1, so the branch correction is b**(-K)*v**(-Q), outside alpha.
Zero route markers are used only for coefficients whose degree is zero.

Uniform optional smoothing fixes a path first. A quiet half extracts IID
free slots with compound PGF J(t)=1-theta+theta*(1-p+p*t)**16; a free half
has PGF t**16<=J(t). Across256 halves, Pr[M<m]<=t**(-(m-1))*J(t)**256.
Conditional on M, the extracted categories are IID with law(1,u,w)/U and
independent of residuals. Thus every regional entry, including terminal,
can be multiplied by rho=a_m+delta*(1-a_m), where a_m is the maximum atom
of Multinomial(m;(u,w,1)/U). There is no guaranteed1440-free-slot baseline.

Four actual profiles and nine fixed K/Q branches only; no all-profile or
whole-code certificate. Bounded single-start optimization selects witnesses;
reported bounds use the actual finite positive transfer, not a Perron proxy.
FLOAT only,48-second guard,20GiB free, estimated peak below128MiB, no files.
"""

from fractions import Fraction
import hashlib
import importlib.util
import itertools
import json
import math
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode=True
PIN="473428234351fd0813a4d5cf3e9bcf4db96179ef7e1cd0b884c232cd0e8a2a73"
PROFILES=((128,127),(115,55),(64,64),(0,0))
KS=(256,512,1024)
QS=(16384,20480,24576)
LN2=math.log(2)
LNQ=math.log(256)
LOG_SOURCE=256*(160*LNQ-32*math.log(255))-math.log(255)
DEADLINE=0.


def guard():
    if time.monotonic()>DEADLINE:
        raise RuntimeError("48-second pilot deadline exhausted; no retry")


def load():
    path=Path(__file__).with_name("cheap_lane_quiet_profile_smoothing.py")
    if hashlib.sha256(path.read_bytes()).hexdigest()!=PIN:
        raise RuntimeError("frozen smoothing helper hash mismatch")
    spec=importlib.util.spec_from_file_location("frozen_path_smoothing",path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    counts,_,screen=module.load_counts(DEADLINE)
    return module,counts.DeletionCounts(screen.load_counts())


def add(a,b):
    if a<b:
        a,b=b,a
    return a if b==-math.inf else a+math.log1p(math.exp(b-a))


def mul(a,b):
    return (add(a[0]+b[0],a[1]+b[2]),add(a[0]+b[1],a[1]+b[3]),
            add(a[2]+b[0],a[3]+b[2]),add(a[2]+b[1],a[3]+b[3]))


def power(a,n):
    out=(0.,-math.inf,-math.inf,0.)
    while n:
        if n&1:
            out=mul(out,a)
        n//=2
        if n:
            a=mul(a,a)
    return out


def chain(a,b,packets=128,pairs=128):
    left=power(a,packets); right=power(b,packets)
    total=mul(mul(power(mul(left,right),pairs-1),left),power(b,packets-1))
    return add(total[0],total[1])


def tests():
    checked=0
    a=((Fraction(3,2),Fraction(1,5)),(Fraction(1),Fraction(1)))
    b=((Fraction(7,4),Fraction(1,5)),(Fraction(1),Fraction(1)))
    for packets,pairs in ((1,1),(2,1),(1,2),(2,2)):
        steps=[a if region%2==0 else b for region in range(2*pairs)
               for _ in range(packets)]
        steps=steps[:-1]
        exact=Fraction(0)
        for tail in itertools.product((0,1),repeat=len(steps)):
            states=(0,)+tail; value=Fraction(1)
            for i,matrix in enumerate(steps):
                value*=matrix[states[i]][states[i+1]]
            exact+=value; checked+=1
        got=chain(tuple(math.log(float(x)) for row in a for x in row),
                  tuple(math.log(float(x)) for row in b for x in row),packets,pairs)
        assert abs(got-math.log(float(exact)))<1e-12
    # Post-fractional marks dominate the selected branch path by path.
    birthmark=Fraction(1,5); quietmark=Fraction(7,4)
    markchecks=0
    for tail in itertools.product((0,1),repeat=7):
        states=(0,)+tail
        births=sum(a==0 and b==1 for a,b in zip(states,states[1:]))
        quiet=sum(a==0 and b==0 for a,b in zip(states,states[1:]))
        for limit in (1,2,3):
            for minimum in (1,3,5):
                if births<=limit and quiet>=minimum:
                    assert birthmark**(births-limit)*quietmark**(quiet-minimum)>=1
                    markchecks+=1
    # Every path has256 half slots, even when its last packet is free.
    theta=Fraction(3,4); p=Fraction(2,3)
    for t in (Fraction(1,3),Fraction(3,4)):
        assert t**16<=1-theta+theta*(1-p+p*t)**16
    return {"exact_path_terms":checked,"branch_mark_checks":markchecks,
            "charged_regular_steps":2*128*128-1,"free_half_pgf_checks":2}


def make_terms(np):
    degrees=[]; coeff=[]
    for f in range(17):
        for s in range(17-f):
            row=[0.]*9
            scale=math.comb(16,f)*math.comb(16-f,s)/2**s
            for j in range(s+1):
                row[min(4,f+j)+min(4,f+s-j)]+=scale*math.comb(s,j)
            degrees.append((f,s)); coeff.append(row)
    return np.asarray(degrees,dtype=float),np.asarray(coeff,dtype=float)


def smooth_row(row,x,helper):
    u,w=row["u"],row["w"]; y=math.exp(-x)
    base=1+y*y*u+y*w
    theta=min(1.,math.exp(16*math.log(base)-row["logP"]))
    p=min(1.,y*y*(1+u+w)/base)
    mean=4096*theta*p
    best={"log2_cap":0.,"cutoff":0,"t":1.,"delta":1.,"mode_counts":[0,0,0]}
    for fraction in (.4,.6,.8):
        cutoff=max(1,int(mean*fraction)); target=cutoff-1
        lo,hi=0.,1.
        for _ in range(40):
            t=(lo+hi)/2; inner=1-p+p*t
            tilted=4096*theta*p*t*inner**15/(1-theta+theta*inner**16)
            if tilted<target:
                lo=t
            else:
                hi=t
        t=(lo+hi)/2
        delta=min(1.,math.exp(256*math.log1p(theta*((1-p+p*t)**16-1))-target*math.log(t)))
        logatom,mode=helper.exact_mode(u,w,cutoff)
        atom=math.exp(logatom); cap=atom+delta*(1-atom)
        value=math.log2(cap)
        if value<best["log2_cap"]:
            best={"log2_cap":value,"cutoff":cutoff,"t":t,"delta":delta,"mode_counts":mode}
    return {**best,"theta":theta,"slot_p":p,"mean_extracted":mean}


def optimize(profile,k,q,count,helper,np,minimize,degrees,coeff):
    targets=[(f,256-2*f) for f in profile]
    masks=[np.logical_and(degrees[:,0]==0 if f==0 else True,
                         degrees[:,1]==0 if s==0 else True) for f,s in targets]
    metadata=[]; variables=[.4054651081081644,-8.,1.]
    bounds=[(.0001,LNQ),(-40.,0.),(0.,20.)]
    for (f,s),mask in zip(targets,masks):
        active=[axis for axis,t in enumerate((f,s)) if t]
        indices=[]; empty=4096-f-s
        for axis in active:
            indices.append(len(variables)); variables.append(math.log((f,s)[axis]/empty)); bounds.append((-30.,12.))
        normalizer=math.log(math.comb(4096,f)*math.comb(4096-f,s))
        metadata.append((active,indices,normalizer,degrees[mask],coeff[mask]))
    evaluations=0
    def evaluate(values,detail=False):
        nonlocal evaluations
        guard(); evaluations+=1
        x,logb,logv=map(float,values[:3]); monomials=np.exp(-x*np.arange(9))
        rows=[]; matrices=[]
        for (f,s),(active,indices,norm,deg,cfs) in zip(targets,metadata):
            logs=[0.,0.]; odds=[0.,0.]
            for axis,index in zip(active,indices):
                logs[axis]=float(values[index]); odds[axis]=math.exp(logs[axis])
            logU=math.log1p(sum(odds))
            terms=np.log(cfs@monomials)+deg@np.asarray(logs)
            high=float(terms.max()); logP=high+math.log(float(np.exp(terms-high).sum()))
            h=4096*logU-f*logs[0]-s*logs[1]-norm
            matrix=(logv+2*logP-32*logU,logb,0.,0.)
            matrices.append(matrix); rows.append({"full":f,"singles":s,"u":odds[0],"w":odds[1],
                                                  "logP":logP,"H_bits":h/LN2})
        source=x/LNQ*LOG_SOURCE
        routing=128*sum(row["H_bits"] for row in rows)*LN2
        paths=chain(*matrices); birth=-k*logb; quiet=-q*logv; source_count=math.log(count)
        score=source+routing+paths+birth+quiet+source_count
        if not detail:
            return score/LN2
        smoothing=[smooth_row(row,x,helper) for row in rows]
        smoothing_bits=128*sum(row["log2_cap"] for row in smoothing)
        return {"profile":profile,"K":k,"Q":q,"bare_log2_bound":score/LN2,
            "smoothed_log2_bound":score/LN2+smoothing_bits,"alpha":x/LNQ,"y":math.exp(-x),
            "birth_mark":math.exp(logb),"quiet_mark":math.exp(logv),
            "components_bits":{"source":source/LN2,"routing_normalizer":routing/LN2,
              "finite_binary_paths":paths/LN2,"birth_correction":birth/LN2,"quiet_correction":quiet/LN2,
              "source_count":source_count/LN2,"uniform_smoothing":smoothing_bits},
            "regional_witnesses":[{**row,"smoothing":sm} for row,sm in zip(rows,smoothing)]}
    # Fixed deterministic starting point matching the intended transition counts.
    ratio=q/32768; births=k/32768
    lam=(1-ratio-births)/(1-ratio-2*births)
    variables[1]=math.log(lam*lam*births*births/((ratio+births)*(1-ratio-births)))
    initial=evaluate(variables,True)
    logg=sum(2*r["logP"]-32*math.log1p(r["u"]+r["w"]) for r in initial["regional_witnesses"])/2
    variables[2]=max(0.,min(20.,math.log(lam*ratio/(ratio+births))-logg))
    result=minimize(evaluate,np.asarray(variables),method="L-BFGS-B",bounds=bounds,
                    options={"maxiter":28,"maxfun":260,"ftol":1e-11,"gtol":1e-6,"maxls":12})
    final=evaluate(result.x,True)
    final.update({"optimizer_success":bool(result.success),"optimizer_message":str(result.message),
                  "evaluations":evaluations,"iterations":int(result.nit)})
    return final


def main():
    global DEADLINE
    if len(sys.argv)!=1:
        raise ValueError("this fixed pilot accepts no sweep arguments")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); DEADLINE=started+48
    import numpy as np
    from scipy.optimize import minimize
    helper,counter=load(); checked=tests(); degrees,coeff=make_terms(np); results=[]
    for profile in PROFILES:
        count=counter.profile(4096,256,*profile)
        assert count>0 and sum(profile)<=255
        for k in KS:
            for q in QS:
                results.append(optimize(profile,k,q,count,helper,np,minimize,degrees,coeff))
    best_grid=[]
    for k in KS:
        for q in QS:
            subset=[r for r in results if r["K"]==k and r["Q"]==q]
            best_grid.append({"K":k,"Q":q,
               "selected_profile_log2_union":math.log2(sum(2**(r["smoothed_log2_bound"]-max(x["smoothed_log2_bound"] for x in subset)) for r in subset))+max(x["smoothed_log2_bound"] for x in subset),
               "minimum_profile_bound":min(r["smoothed_log2_bound"] for r in subset),
               "maximum_profile_bound":max(r["smoothed_log2_bound"] for r in subset)})
    best=[min((r for r in results if tuple(r["profile"])==profile),key=lambda r:r["smoothed_log2_bound"]) for profile in PROFILES]
    print(json.dumps({"scope":"FLOAT selected actual support256 profiles; all binary paths with births<=K and quiet regular steps>=Q",
        "whole_code_certificate":False,"directed":False,"selected_profiles":PROFILES,"branch_grid":best_grid,
        "best_per_profile":best,"all_scores":[{key:r[key] for key in ("profile","K","Q","bare_log2_bound","smoothed_log2_bound","alpha","evaluations","optimizer_success")} for r in results],
        "mark_witnesses":len(results),"total_evaluations":sum(r["evaluations"] for r in results),"tests":checked,
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"smoothing_pin":PIN,
        "seconds":time.monotonic()-started,"generated_bytes":0,"estimated_peak_bytes":128*(1<<20)},indent=2))


if __name__=="__main__":
    main()
