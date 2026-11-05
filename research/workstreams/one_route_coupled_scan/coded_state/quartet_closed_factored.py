"""Exact factorization of the seven-type kernel, summing all destinations.

Full outgoing coverage at each tested incoming population. FLOAT full-width
bounds and exact local checks, not a global certificate. Serial Peach180s,
<=256MiB, no generated files.
"""
from fractions import Fraction as F
from collections import defaultdict
import itertools
import math
import os
from pathlib import Path
import shutil
import threading
import time
import numpy as np
from quartet_closed_quotient import kernel
from quartet_closed_population import multinomial
from quartet_two_cohort import compositions


def nonfull(h,activation,x,y):
    base=h[1]+activation*(1+x*y)*h[4]
    return sum(math.comb(4,j)*(x*h[0])**j*base**(4-j) for j in range(1,5))


def quartet_value(counts,h,activation,death,x,y):
    """Weighted destination sum for one labelled incoming quartet."""
    selected=sum(counts[4:]);blocked=counts[3];t=1+x*y
    if selected and blocked:return 0
    if selected:
        values=(h[1]+activation*h[4],h[2],h[3],0,
                activation*h[5],activation*h[6],death*h[3]+activation*h[6])
        return t**4*math.prod(v**n for v,n in zip(values,counts) if n)
    answer=nonfull(h,activation,x,y)
    if not blocked:
        values=(h[1]+activation*t*h[4],h[2],h[3])
        answer+=math.prod(v**n for v,n in zip(values,counts[:3]) if n)
    return answer


def source_fit(pop,h,activation,death,x,y):
    width=sum(pop);groups=width//4
    active=[i for i,n in enumerate(pop) if n];omitted=max(active,key=lambda i:pop[i]);kept=[i for i in active if i!=omitted]
    powers=[];base=[]
    for m in compositions(4,7):
        if any(m[i] for i in range(7) if i not in active):continue
        value=quartet_value(m,h,activation,death,x,y)
        if value<=0:continue
        base.append(math.log2(value*multinomial(m))-sum(n*math.log2(v) for n,v in zip(m,h)))
        powers.append([m[i] for i in kept])
    base=np.array(base);powers=np.array(powers,dtype=float);target=np.array([pop[i]/groups for i in kept]);tilts=np.zeros(len(kept))
    def evaluate(tilts):
        values=base+powers@tilts;peak=max(values);raw=np.exp2(values-peak);p=raw/sum(raw)
        total=peak+math.log2(sum(raw));mean=p@powers;centered=powers-mean
        return total,mean,(centered.T*p)@centered
    for _ in range(100):
        total,mean,cov=evaluate(tilts)
        if not len(tilts) or max(abs(mean-target))<1e-10:break
        direction=np.linalg.solve(math.log(2)*cov+1e-10*np.eye(len(tilts)),target-mean)
        direction/=max(1.,max(abs(direction))/4);rate=1.;old=total-target@tilts
        for _ in range(25):
            proposed=np.clip(tilts+rate*direction,-100,100)
            if evaluate(proposed)[0]-target@proposed<=old+1e-11:break
            rate*=.5
        tilts=proposed
    total,mean,cov=evaluate(tilts);z=np.zeros(7);z[active]=1.
    for i,t in zip(kept,tilts):z[i]=2**t
    cauchy=groups*(total-target@tilts)-math.log2(multinomial(pop))
    return z,total,cauchy,float(max(abs(mean-target),default=0.))


def row_bound(pop,h,activation,death,x,y):
    """Extract selected and blocked counts; retain exact selected-age split."""
    z,logQ,cauchy,residual=source_fit(pop,h,activation,death,x,y)
    q=z/np.asarray(h);t=1+x*y;n=nonfull(h,activation,x,y)
    shared=sum(q[:3]);ordinary=(h[1]+activation*h[4])*q[0]+h[2]*q[1]+h[3]*q[2]
    paid=(h[1]+activation*t*h[4])*q[0]+h[2]*q[1]+h[3]*q[2]
    live_terms=np.array([activation*h[5]*q[4],activation*h[6]*q[5],(death*h[3]+activation*h[6])*q[6]])
    live=sum(live_terms)
    weights={(0,0):n*shared**4+paid**4}
    for j in range(1,5):
        weights[0,j]=n*math.comb(4,j)*shared**(4-j)*q[3]**j
        weights[j,0]=t**4*math.comb(4,j)*ordinary**(4-j)*live**j
    total=sum(weights.values());assert abs(math.log2(total)-logQ)<1e-7
    weights={k:v/total for k,v in weights.items() if v}
    selected=sum(pop[4:]);blocked=pop[3];target=(selected,blocked)
    current=np.zeros((selected+1,blocked+1));current[0,0]=1.
    for _ in range(sum(pop)//4):
        after=np.zeros_like(current)
        for (a,b),value in weights.items():
            if a<=selected and b<=blocked:after[a:,b:]+=value*current[:selected+1-a,:blocked+1-b]
        current=after
    mass=float(current[target])
    if mass<=0:return -math.inf,cauchy,residual
    split=math.log2(multinomial(pop[4:]))
    if selected:
        split+=sum(count*math.log2(value/live) for count,value in zip(pop[4:],live_terms) if count)
    return cauchy+math.log2(mass)+split,cauchy,residual


def tests():
    checked=0
    for h,a,d,x,y in ((tuple(map(F,(1,2,3,4,5,6,7))),F(3,2),F(1,4),F(1,3),F(1,8)),
                       ((F(1),)*7,F(2),F(1),F(1),F(1))):
        for hist in itertools.combinations_with_replacement(range(7),4):
            counts=tuple(hist.count(t) for t in range(7));exact=F(0)
            for (dest,u,v,activity,deaths),ways in kernel(hist):
                exact+=ways*math.prod(h[t] for t in dest)*a**activity*d**deaths*x**u*(1+x*y)**v
            assert quartet_value(counts,h,a,d,x,y)==exact;checked+=1
    # Exact width8 source coefficients after summing ALL destinations.
    cats=list(compositions(4,7));table={m:quartet_value(m,(1,)*7,2,1,1,1)*multinomial(m) for m in cats}
    profiles=((8,0,0,0,0,0,0),(0,0,0,8,0,0,0),(0,0,0,0,0,0,8),
              (2,2,0,0,2,0,2),(2,0,2,2,0,0,2),(0,0,0,4,0,0,4),
              (1,2,1,0,1,1,2),(2,1,1,1,1,1,1))
    for pop in profiles:
        value=0
        for m,a in table.items():
            other=tuple(n-v for n,v in zip(pop,m))
            if other in table:value+=a*table[other]
        exact=F(value,multinomial(pop));upper,_,_=row_bound(pop,(1.,)*7,2.,1.,1.,1.)
        assert exact==0 or math.log2(exact)<=upper+1e-8
    print(dict(exact_quartet_factorizations=checked,all_destination_width8_checks=len(profiles)),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(180,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:
        tests();logb=2.5;x=2**(-.03796846262217181)
        h=np.exp2([0.,-.039394433547311924,-.03625712677466449,-.04243334838197499,
                   -1.9293817778852493,-1.5382677895381995,-1.1471538011850797])
        a=2**(707*logb/2560);d=2**(-1905*logb/2560);y=2**(-6*logb/5)
        profiles=((512,0,0,0,0,0,0),(104,176,72,32,16,16,96),
                  (104,120,64,32,16,16,160),(104,56,32,64,32,32,192),
                  (0,0,0,0,0,0,512))
        for pop in profiles:
            value,cauchy,residual=row_bound(pop,h,a,d,x,y)
            regional=16*value-2100*math.log2(x)-math.log2(math.comb(8192,2100))
            print(dict(population=pop,all_destinations=True,log2_row_upper=value,
                       source_cauchy=cauchy,regional_row_envelope=regional,
                       residual=residual,seconds=time.monotonic()-start),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
