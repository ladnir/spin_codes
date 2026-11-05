"""Finite-population bound for closed-set quartet rows.

Retain pending input/selection categories; maximize only internal run ages.
Exact small tests and FLOAT full-width local diagnostics. Serial Peach180s,
<128MiB arrays,no files. Not a regional or whole-code certificate.
"""
from fractions import Fraction as F
import itertools
import math
import os
from pathlib import Path
import shutil
import threading
import time
import numpy as np
from quartet_closed_transfer import TYPES,local
from quartet_two_cohort import compositions
from quartet_history_screen import logchoose


PROPOSAL=np.array([0.,-.039394433547311924,-1.9293817778852493,3.5448453647630407,
                   -.03625712677466449,3.5448453647630904,-.04243334838197499,
                   -1.5382677895381995,-1.1471538011850797,3.544845364762959,
                   -.049691309338134336,-1.1481576278443342])


def quartet_caps(logh=PROPOSAL,logx=-.03796846262217181,logb=2.5,logR=0.,logS=0.):
    """Return log2 caps indexed by counts of choices0,1,2 in a quartet."""
    h={t:logh[i] for i,t in enumerate(TYPES)};result={}
    for hist in itertools.combinations_with_replacement(TYPES,4):
        values=[]
        for (dest,unpaid,covered,activity,deaths),ways in local(hist):
            common=math.log2(ways)+sum(h[t] for t in dest)-sum(h[t] for t in hist)
            common+=logb*(707*activity-1905*deaths)/2560
            common+=logR*sum(t%3!=0 for t in dest)+logS*activity
            for c in range(covered+1):
                values.append(common+math.log2(math.comb(covered,c))+(unpaid+c)*logx-logb*6*c/5)
        if not values:continue
        peak=max(values);value=peak+math.log2(sum(2**(v-peak) for v in values))
        category=tuple(sum(t%3==i for t in hist) for i in range(3))
        result[category]=max(result.get(category,-math.inf),value)
    assert len(result)==15
    return result


def multinomial(ns):
    answer=1;remaining=sum(ns)
    for n in ns:answer*=math.comb(remaining,n);remaining-=n
    return answer


def tilted_terms(caps,counts):
    """Positive saddle-point proposal; exact coefficient identity needs no optimum."""
    width=sum(counts);groups=width//4;active=[i for i,n in enumerate(counts) if n]
    # Omit the largest category from the coefficient array.
    omitted=max(active,key=lambda i:counts[i]);kept=[i for i in active if i!=omitted]
    categories=[m for m in caps if all(m[i]==0 for i in range(3) if i not in active)]
    powers=np.array([[m[i] for i in kept] for m in categories],dtype=float)
    base=np.array([caps[m]+math.log2(multinomial(m)) for m in categories])
    x=np.zeros(len(kept));target=np.array([counts[i]/groups for i in kept])
    def evaluate(x):
        terms=base+powers@x;peak=max(terms);raw=np.exp2(terms-peak);p=raw/sum(raw)
        total=peak+math.log2(sum(raw));mean=p@powers
        centered=powers-mean;cov=(centered.T*p)@centered
        return total,p,mean,cov
    for _ in range(40):
        total,p,mean,cov=evaluate(x)
        if not len(x) or max(abs(mean-target))<1e-11:break
        direction=np.linalg.solve(math.log(2)*cov+1e-12*np.eye(len(x)),target-mean)
        direction/=max(1.,max(abs(direction))/4)
        old=total-target@x;rate=1.
        for _ in range(30):
            proposal=x+rate*direction
            if evaluate(proposal)[0]-target@proposal<=old+1e-13:break
            rate*=.5
        x=proposal
    return categories,kept,x,evaluate


def population_bound(caps,counts,offset=None):
    width=sum(counts);assert width%4==0
    categories,kept,x,evaluate=tilted_terms(caps,counts)
    if offset is not None:x=x+np.asarray(offset)[:len(x)]
    logQ,p,mean,cov=evaluate(x);groups=width//4
    shape=tuple(counts[i]+1 for i in kept) or (1,)
    assert math.prod(shape)*16<64*(1<<20)
    current=np.zeros(shape);current.flat[0]=1.
    # Terms are positive. Coefficients above the requested degrees cannot
    # return to the target; their omission is exact, not a tail truncation.
    for step in range(groups):
        after=np.zeros_like(current)
        for category,prob in zip(categories,p):
            degree=tuple(category[i] for i in kept) or (0,)
            if any(d>=n for d,n in zip(degree,shape)):continue
            source=tuple(slice(0,n-d) for d,n in zip(degree,shape))
            destination=tuple(slice(d,n) for d,n in zip(degree,shape))
            after[destination]+=prob*current[source]
        current=after
    mass=float(current[tuple(counts[i] for i in kept) or (0,)])
    assert mass>0 and math.isfinite(mass)
    cauchy=groups*logQ-sum(counts[i]*xi for i,xi in zip(kept,x))-math.log2(multinomial(counts))
    exact_float=cauchy+math.log2(mass)
    return exact_float,cauchy,math.log2(mass)


def tests():
    coefficients=0
    integer_caps={m:1+3*m[1]+7*m[2]+2*m[1]*m[2] for m in compositions(4,3)}
    logs={m:math.log2(v) for m,v in integer_caps.items()}
    for width in (4,8):
        for counts in compositions(width,3):
            if width==4:exact=F(integer_caps[counts])
            else:
                total=0
                for m,v in integer_caps.items():
                    other=tuple(n-a for n,a in zip(counts,m))
                    if min(other)<0:continue
                    total+=multinomial(m)*multinomial(other)*v*integer_caps[other]
                exact=F(total,multinomial(counts))
            value,upper,mass=population_bound(logs,counts)
            assert abs(value-math.log2(exact))<1e-9 and math.log2(exact)<=upper+1e-9
            shifted=population_bound(logs,counts,offset=(.25,-.375))[0]
            assert abs(value-shifted)<1e-9;coefficients+=1
    print(dict(exact_population_coefficients=coefficients,independent_tilt_replays=coefficients),flush=True)


def diagnostic():
    start=time.monotonic();caps=quartet_caps();independent=128*max(caps.values())
    for allowed,selected in ((384,0),(384,32),(384,128),(384,256),(408,32),(408,128),(512,256),(512,512)):
        counts=(512-allowed,allowed-selected,selected)
        value,upper,mass=population_bound(caps,counts)
        replay=population_bound(caps,counts,offset=(.125,-.125))[0]
        assert abs(value-replay)<1e-8 and value<=independent+1e-8
        print(dict(allowed=allowed,selected=selected,independent_quartet_log2=independent,
                   population_log2=value,cauchy_log2=upper,coefficient_mass_log2=mass,
                   local_only=True,directed=False),flush=True)
    print(dict(seconds=time.monotonic()-start,files_written=0),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(180,lambda:os._exit(124));timer.daemon=True;timer.start()
    try:tests();diagnostic()
    finally:timer.cancel()


if __name__=='__main__':main()
