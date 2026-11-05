"""Population Cauchy bound with an exact original-inactive coefficient.

Local exact toy checks and FLOAT width512 diagnostic, not a certificate.
Serial Peach180s,<128MiB arrays,no generated files.
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
from quartet_closed_population import PROPOSAL,multinomial
from quartet_two_cohort import compositions


def coefficient_caps(logh=PROPOSAL,logb=2.5,logR=0.,logS=0.):
    h={t:logh[i] for i,t in enumerate(TYPES)};caps={}
    for hist in itertools.combinations_with_replacement(TYPES,4):
        row=np.full(5,-np.inf)
        for (dest,unpaid,covered,activity,deaths),ways in local(hist):
            value=math.log2(ways)+sum(h[t] for t in dest)-sum(h[t] for t in hist)
            value+=logb*(707*activity-1905*deaths)/2560
            value+=logR*sum(t%3!=0 for t in dest)+logS*activity
            for paid in range(covered+1):
                term=value+math.log2(math.comb(covered,paid))-logb*6*paid/5
                row[unpaid+paid]=np.logaddexp2(row[unpaid+paid],term)
        category=tuple(sum(t%3==i for t in hist) for i in range(3))
        caps[category]=np.maximum(caps.get(category,np.full(5,-np.inf)),row)
    assert len(caps)==15 and all(np.all(np.isfinite(v)) for v in caps.values())
    return caps


def bound(caps,counts,inactive,replay_shift=0.):
    width=sum(counts);groups=width//4;assert width%4==0
    active=[i for i,n in enumerate(counts) if n];omitted=max(active,key=lambda i:counts[i]);kept=[i for i in active if i!=omitted]
    terms=[];powers=[];js=[]
    for m,values in caps.items():
        if any(m[i] for i in range(3) if i not in active):continue
        for j,value in enumerate(values):
            if (inactive==0 and j!=0) or (inactive==width and j!=4):continue
            terms.append(value+math.log2(multinomial(m)))
            powers.append([m[i] for i in kept]+([j] if 0<inactive<width else []));js.append(j)
    base=np.array(terms);powers=np.array(powers,dtype=float);js=np.array(js,dtype=int)
    target=np.array([counts[i]/groups for i in kept]+([inactive/groups] if 0<inactive<width else []))
    x=np.zeros(len(target))
    def evaluate(x):
        values=base+powers@x;peak=max(values);raw=np.exp2(values-peak);p=raw/sum(raw)
        total=peak+math.log2(sum(raw));mean=p@powers;centered=powers-mean
        return total,p,mean,(centered.T*p)@centered
    for _ in range(40):
        total,p,mean,cov=evaluate(x)
        if not len(x) or max(abs(mean-target))<1e-11:break
        direction=np.linalg.solve(math.log(2)*cov+1e-12*np.eye(len(x)),target-mean)
        direction/=max(1.,max(abs(direction))/4);rate=1.;old=total-target@x
        for _ in range(30):
            proposal=x+rate*direction
            if evaluate(proposal)[0]-target@proposal<=old+1e-13:break
            rate*=.5
        x=proposal
    if 0<inactive<width:x[-1]+=replay_shift
    total,p,mean,cov=evaluate(x);q=np.bincount(js,weights=p,minlength=5)
    probabilities=np.zeros(inactive+1);probabilities[0]=1
    for _ in range(groups):
        after=np.zeros_like(probabilities)
        for j,weight in enumerate(q):
            if j<=inactive:after[j:]+=weight*probabilities[:inactive+1-j]
        probabilities=after
    mass=float(probabilities[inactive]);assert mass>0
    cauchy=groups*(total-target@x)-math.log2(multinomial(counts))
    return cauchy+math.log2(mass),cauchy,math.log2(mass)


def tests():
    integer={m:[1+3*m[1]+7*m[2]+2*m[1]*m[2]+j*(m[0]+1) for j in range(5)] for m in compositions(4,3)}
    caps={m:np.log2(values) for m,values in integer.items()};checked=0
    for width in (4,8):
        for counts in compositions(width,3):
            for inactive in range(width+1):
                if width==4:exact=F(integer[counts][inactive])
                else:
                    numerator=0
                    for m,values in integer.items():
                        other=tuple(n-a for n,a in zip(counts,m))
                        if min(other)<0:continue
                        for j,value in enumerate(values):
                            if 0<=inactive-j<=4:numerator+=multinomial(m)*multinomial(other)*value*integer[other][inactive-j]
                    exact=F(numerator,multinomial(counts))
                value,cauchy,mass=bound(caps,counts,inactive)
                assert math.log2(exact)<=value+1e-9 and value<=cauchy+1e-9
                assert abs(value-bound(caps,counts,inactive,.3)[0])<1e-9;checked+=1
    print(dict(exact_inactive_coefficient_envelopes=checked,inactive_tilt_replays=checked),flush=True)


def diagnostic():
    start=time.monotonic();caps=coefficient_caps()
    for allowed,selected in ((384,0),(384,32),(384,128),(408,32),(408,128),(512,256),(512,512)):
        for inactive in (128,131,160):
            counts=(512-allowed,allowed-selected,selected)
            value,cauchy,mass=bound(caps,counts,inactive)
            replay=bound(caps,counts,inactive,.25)[0];assert abs(value-replay)<1e-8
            print(dict(allowed=allowed,selected=selected,inactive=inactive,
                       conditional_coefficient_log2=value,inactive_mass_log2=mass,
                       local_only=True,directed=False),flush=True)
    print(dict(seconds=time.monotonic()-start,files_written=0),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(180,lambda:os._exit(124));timer.daemon=True;timer.start()
    try:tests();diagnostic()
    finally:timer.cancel()


if __name__=='__main__':main()
