"""All-destination two-update upper bound at a fixed original-inactive count.

Keeps both partitions' B/S exclusion. The first local transition is enlarged
coefficientwise to independent coordinate transitions; this loss is explicit.
Source types are bounded by positive markers; the pair's inactivity count is
extracted. FLOAT proposals and separate directed replay, never a full cert.
Serial180s, <128MiB, source-only, no files.
"""
from collections import defaultdict
from fractions import Fraction as F
import itertools
from functools import lru_cache
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
from quartet_closed_directed import coefficient_upper,power_upper
from quartet_two_update import log_upper


def add(*polys):
    out=[0]*max(map(len,polys))
    for p in polys:
        for i,v in enumerate(p):out[i]+=v
    return out


def mul(a,b):
    out=[0]*(len(a)+len(b)-1)
    for i,u in enumerate(a):
        for j,v in enumerate(b):out[i+j]+=u*v
    return out


def scale(p,c):return [v*c for v in p]


def badd(*polys):
    out=defaultdict(lambda:0)
    for p in polys:
        for k,v in p.items():out[k]+=v
    return dict(out)


def bscale(p,c):return {k:v*c for k,v in p.items() if v*c}


def bmul(p,q):
    out=defaultdict(lambda:0)
    for (f,i),v in p.items():
        for (g,j),w in q.items():out[f+g,i+j]+=v*w
    return dict(out)


def bpower(p,n):
    out={(0,0):1}
    for _ in range(n):out=bmul(out,p)
    return out


def flagged_first(z,a,d,y):
    ordinary=sum(z[:4]);t={(0,0):1,(0,1):y};full={(1,0):1,(1,1):y}
    return ({(0,1):ordinary},badd({(0,0):ordinary},bscale(full,z[0])),
            bscale(full,z[1]),bscale(full,z[2]+d*z[6]),
            badd(bscale(t,a*ordinary),bscale(full,a*z[0])),
            bscale(full,a*z[4]),bscale(full,a*(z[5]+z[6])))


def flagged_polynomial(z,h,a,d,y):
    t={(0,0):1,(0,1):y};base=badd({(0,0):h[1]},bscale(t,a*h[4]));nonfull={}
    for j in range(1,5):nonfull=badd(nonfull,bscale(bmul({(0,j):h[0]**j},bpower(base,4-j)),math.comb(4,j)))
    v=flagged_first(z,a,d,y);shared=badd(*v[:4])
    ordinary=badd(bscale(v[0],h[1]+a*h[4]),bscale(v[1],h[2]),bscale(v[2],h[3]))
    live=badd(bscale(v[4],a*h[5]),bscale(v[5],a*h[6]),bscale(v[6],d*h[3]+a*h[6]))
    paid=badd(bmul(base,v[0]),bscale(v[1],h[2]),bscale(v[2],h[3]));selected={}
    for j in range(1,5):selected=badd(selected,bscale(bmul(bpower(ordinary,4-j),bpower(live,j)),math.comb(4,j)))
    result=badd(bmul(nonfull,bpower(shared,4)),bmul(bpower(t,4),selected),bpower(paid,4))
    assert all(f<=4 and i<=8 and v>=0 for (f,i),v in result.items())
    return result


@lru_cache(maxsize=1)
def nonfull_occupancy():
    table=[[1]]
    for g in range(1,129):
        after=[0]*(4*g+1)
        for i,v in enumerate(table[-1]):
            for j in range(1,5):after[i+j]+=v*math.comb(4,j)
        table.append(after)
    return table


def flag_factor(f,inactive,width=512):
    if f%4:return F(0)
    groups=width//4;full=f//4;rest=width-f;zeros=min(inactive,rest)
    if zeros<groups-full:return F(0)
    return F(math.comb(groups,full)*nonfull_occupancy()[groups-full][zeros],
             math.comb(width,f)*math.comb(rest,zeros))


def power(p,n):
    out=[1]
    for _ in range(n):out=mul(out,p)
    return out


def first(z,a,d,y):
    """Relaxed first-coordinate type generating polynomials, variable x."""
    t=[1,y]
    return ([0,sum(z[:4])],
            add(scale(t,z[0]),[z[1]+z[2]+z[3]]),
            scale(t,z[1]),scale(t,z[2]+d*z[6]),
            scale(t,a*sum(z[:4])),scale(t,a*z[4]),scale(t,a*z[5]+a*z[6]))


def pair_polynomial(z,h,a,d,y):
    """Q(z,x), after exact second-quartet summation. Degree in x<=8."""
    t=[1,y];base=add([h[1]],scale(t,a*h[4]))
    nonfull=[0]
    for j in range(1,5):nonfull=add(nonfull,scale(mul([0]*j+[h[0]**j],power(base,4-j)),math.comb(4,j)))
    v=first(z,a,d,y)
    shared=add(*v[:4]);ordinary=add(scale(v[0],h[1]+a*h[4]),scale(v[1],h[2]),scale(v[2],h[3]))
    live=add(scale(v[4],a*h[5]),scale(v[5],a*h[6]),scale(v[6],d*h[3]+a*h[6]))
    paid=add(mul(base,v[0]),scale(v[1],h[2]),scale(v[2],h[3]))
    selected=[0]
    for j in range(1,5):selected=add(selected,scale(mul(power(ordinary,4-j),power(live,j)),math.comb(4,j)))
    result=add(mul(nonfull,power(shared,4)),mul(power(t,4),selected),power(paid,4))
    assert len(result)<=9 and min(result)>=0
    return result


def point(z,h,a,d,x,y):
    t=1+x*y;base=h[1]+a*t*h[4]
    nonfull=sum(math.comb(4,j)*(x*h[0])**j*base**(4-j) for j in range(1,5))
    v=(x*sum(z[:4]),t*z[0]+sum(z[1:4]),t*z[1],t*(z[2]+d*z[6]),
       a*t*sum(z[:4]),a*t*z[4],a*t*(z[5]+z[6]))
    shared=sum(v[:4]);ordinary=(h[1]+a*h[4])*v[0]+h[2]*v[1]+h[3]*v[2]
    live=a*h[5]*v[4]+a*h[6]*v[5]+(d*h[3]+a*h[6])*v[6]
    paid=base*v[0]+h[2]*v[1]+h[3]*v[2]
    return nonfull*shared**4+t**4*sum(math.comb(4,j)*ordinary**(4-j)*live**j for j in range(1,5))+paid**4


def fit(pop,h,inactive):
    active=[i for i,n in enumerate(pop) if n];omitted=max(active,key=lambda i:pop[i]);kept=[i for i in active if i!=omitted]
    target=np.array([pop[i]/128 for i in kept]+[inactive/128]);v=np.zeros(len(target))
    a=2**(707/1024);d=2**(-1905/1024);y=1/8
    def evaluate(v):
        z=np.zeros(7);z[active]=1
        for i,u in zip(kept,v):z[i]=2**u
        return math.log2(point(z,h,a,d,2**v[-1],y))-target@v
    def gradient(v):
        g=np.zeros(len(v));eps=2e-5
        for i in range(len(v)):
            p=v.copy();m=v.copy();p[i]+=eps;m[i]-=eps
            g[i]=(evaluate(p)-evaluate(m))/(2*eps)
        return evaluate(v),g
    inverse=np.eye(len(v));value,grad=gradient(v)
    for _ in range(120):
        direction=-inverse@grad
        if grad@direction>=0:direction=-grad;inverse=np.eye(len(v))
        direction/=max(1.,max(abs(direction))/3);rate=1.
        for _ in range(25):
            nxt=np.clip(v+rate*direction,-40,40);step=nxt-v
            if evaluate(nxt)<=value+1e-4*grad@step:break
            rate*=.5
        if max(abs(step))<1e-9:break
        nv,ng=gradient(nxt);delta=ng-grad;curve=delta@step
        if curve>1e-12:
            correction=np.eye(len(v))-np.outer(step,delta)/curve
            inverse=correction@inverse@correction.T+np.outer(step,step)/curve
        else:inverse=np.eye(len(v))
        v,value,grad=nxt,nv,ng
        if max(abs(grad))<1e-8:break
    z=np.zeros(7);z[active]=1
    for i,u in zip(kept,v):z[i]=2**u
    return z,2**v[-1],float(max(abs(grad)))


def avoidance_proposal(blocked,selected,width=512):
    counts=(width-blocked-selected,blocked,selected);active=[i for i,n in enumerate(counts) if n]
    if not blocked or not selected:return 0.,[F(1)]*3
    terms=[]
    for b in range(5):
        for s in range(5-b):
            if b and s:continue
            m=(4-b-s,b,s)
            if any(m[i] and not counts[i] for i in range(3)):continue
            terms.append((m,multinomial(m)))
    omit=max(active,key=lambda i:counts[i]);kept=[i for i in active if i!=omit];v=np.zeros(len(kept));target=np.array([counts[i]/(width/4) for i in kept])
    powers=np.array([[m[i] for i in kept] for m,w in terms]);base=np.array([math.log2(w) for m,w in terms])
    for _ in range(60):
        logs=base+powers@v;peak=max(logs);p=np.exp2(logs-peak);p/=sum(p);mean=p@powers;cent=powers-mean
        if max(abs(mean-target))<1e-10:break
        move=np.linalg.solve(math.log(2)*(cent.T*p)@cent+1e-10*np.eye(len(v)),target-mean)
        move/=max(1.,max(abs(move))/3);v=np.clip(v+move,-40,40)
    markers=[F(int(i in active)) for i in range(3)]
    for i,u in zip(kept,v):markers[i]=power_upper(int(round(1024*u)))
    weights={(m[1],m[2]):F(w)*math.prod(z**n for z,n in zip(markers,m)) for m,w in terms}
    total=sum(weights.values());probs={k:w/total for k,w in weights.items()}
    mass=coefficient_upper(probs,(blocked,selected),width//4)
    if not mass:return -math.inf,markers
    upper=log_upper(width//4,total)+log_upper(1,mass)+log_upper(-1,multinomial(counts))
    upper+=sum(log_upper(-n,z) for n,z in zip(counts,markers) if n)
    return min(F(0),upper),markers


def evaluate(pop,h,inactive,directed=False):
    z,x,residual=fit(pop,h,inactive)
    geometry,_=avoidance_proposal(pop[3],sum(pop[4:]))
    if geometry==-math.inf:return -math.inf,residual
    if directed:
        h=[F(max(1,round(v*(1<<24))),1<<24) for v in h]
        z=[power_upper(round(1024*math.log2(v))) if v else F(0) for v in z]
        x=power_upper(round(1024*math.log2(x)));a=power_upper(707);d=power_upper(-1905);y=F(1,8)
        poly=pair_polynomial(z,h,a,d,y);weighted=[v*x**i for i,v in enumerate(poly)];total=sum(weighted)
        mass=coefficient_upper({(i,):v/total for i,v in enumerate(weighted) if v},(inactive,),128)
        assert mass>0
        value=geometry+log_upper(128,total)+log_upper(1,mass)+log_upper(-inactive,x)+log_upper(-1,multinomial(pop))
        value+=sum(log_upper(-n,zz)+log_upper(-n,hh) for n,zz,hh in zip(pop,z,h) if n)
        return value,residual
    poly=pair_polynomial(z,h,2**(707/1024),2**(-1905/1024),1/8)
    weighted=np.array([v*x**i for i,v in enumerate(poly)]);total=sum(weighted);p=weighted/total;current=np.array([1.])
    for _ in range(128):current=np.convolve(current,p)[:inactive+1]
    assert current[inactive]>0
    value=float(geometry)+128*math.log2(total)+math.log2(current[inactive])-inactive*math.log2(x)-math.log2(multinomial(pop))
    value-=sum(n*(math.log2(zz)+math.log2(hh)) for n,zz,hh in zip(pop,z,h) if n)
    return value,residual


def flagged_bound(pop,h,inactive):
    z,x,residual=fit(pop,h,inactive)
    h=[F(max(1,round(v*(1<<24))),1<<24) for v in h]
    z=[power_upper(round(1024*math.log2(v))) if v else F(0) for v in z]
    x=power_upper(round(1024*math.log2(x)));a=power_upper(707);d=power_upper(-1905);y=F(1,8)
    poly=flagged_polynomial(z,h,a,d,y)
    weighted={key:value*x**key[1] for key,value in poly.items() if value};total=sum(weighted.values())
    masses=coefficient_upper({key:v/total for key,v in weighted.items()},(512-pop[3],inactive),128,return_axis=0)
    mass=sum(masses[f]*flag_factor(f,inactive) for f in range(((sum(pop[4:])+3)//4)*4,513-pop[3],4))
    assert mass>0
    value=log_upper(128,total)+log_upper(1,mass)+log_upper(-inactive,x)+log_upper(-1,multinomial(pop))
    value+=sum(log_upper(-n,zz)+log_upper(-n,hh) for n,zz,hh in zip(pop,z,h) if n)
    return value,residual


def flagged_float(pop,h,inactive,upward=False,proposal=None):
    z,x,residual=fit(pop,h,inactive) if proposal is None else proposal
    if upward:
        h=[F(max(1,round(v*(1<<24))),1<<24) for v in h]
        z=[power_upper(round(1024*math.log2(v))) if v else F(0) for v in z]
        x=power_upper(round(1024*math.log2(x)));a=power_upper(707);d=power_upper(-1905);y=F(1,8)
    else:a=2**(707/1024);d=2**(-1905/1024);y=1/8
    poly=flagged_polynomial(z,h,a,d,y);weights={key:value*x**key[1] for key,value in poly.items() if value};total=sum(weights.values())
    targets=(512-pop[3],inactive)
    current=relative_convolution({key:v/total for key,v in weights.items()},targets,128,upward)
    mass=sum((F.from_float(float(current[f,inactive])) if upward else float(current[f,inactive]))*flag_factor(f,inactive)
             for f in range(((sum(pop[4:])+3)//4)*4,513-pop[3],4))
    if not mass:return -math.inf
    if upward:
        value=log_upper(128,total)+log_upper(1,mass)+log_upper(-inactive,x)+log_upper(-1,multinomial(pop))
        value+=sum(log_upper(-n,zz)+log_upper(-n,hh) for n,zz,hh in zip(pop,z,h) if n)
    else:value=128*math.log2(total)+math.log2(mass)-inactive*math.log2(x)-math.log2(multinomial(pop))-sum(n*(math.log2(zz)+math.log2(hh)) for n,zz,hh in zip(pop,z,h) if n)
    return value


def float_upper(value):
    value=F(value);rounded=float(value)
    if F.from_float(rounded)<value:rounded=float(np.nextafter(rounded,np.inf))
    assert F.from_float(rounded)>=value
    return rounded


def relative_convolution(probabilities,targets,groups,upward=False):
    # Requires IEEE754 binary64 operations with gradual underflow.
    assert np.finfo(np.float64).nmant==52
    assert float(np.multiply(np.float64(2.**-1022),np.float64(.5)))==2.**-1023
    assert float(np.multiply(np.float64(2.**-1023),np.float64(1.)))==2.**-1023
    assert float(np.add(np.float64(2.**-1023),np.float64(2.**-1023)))==2.**-1022
    assert all(v>=0 for v in probabilities.values())
    terms=[(key,float_upper(v) if upward else float(v)) for key,v in probabilities.items() if v]
    df=max(key[0] for key,v in terms);di=max(key[1] for key,v in terms)
    current=np.ones((1,1))
    for _ in range(groups):
        rows,cols=current.shape;nr=min(targets[0]+1,rows+df);nc=min(targets[1]+1,cols+di)
        after=np.zeros((nr,nc))
        for (f,i),v in terms:
            r=min(rows,nr-f);c=min(cols,nc-i)
            if r<=0 or c<=0:continue
            product=np.multiply(current[:r,:c],v)
            if upward:np.nextafter(product,np.inf,out=product)
            dest=after[f:f+r,i:i+c];np.add(dest,product,out=dest)
            if upward:np.nextafter(dest,np.inf,out=dest)
        assert np.all(np.isfinite(after));current=after
    return current


def exit_diagnostic(pop,h,inactive,proposal):
    """Tilted relaxation's mean endpoint, not an encoder-state prediction."""
    z,x,_=proposal;a=2**(707/1024);d=2**(-1905/1024);y=1/8;fmax=512-pop[3]
    def local(h):return {key:value*x**key[1] for key,value in flagged_polynomial(z,h,a,d,y).items() if value}
    poly=local(h);total=sum(poly.values())
    current=relative_convolution({key:v/total for key,v in poly.items()},(fmax,inactive),127)
    factors=np.array([float(flag_factor(f,inactive)) if f>=sum(pop[4:]) else 0. for f in range(fmax+1)])
    adjoint={key:float(current[:fmax+1-key[0],inactive-key[1]]@factors[key[0]:]) for key in poly if key[1]<=inactive}
    mass=sum(adjoint[key]*v for key,v in poly.items() if key in adjoint);assert mass>0
    mean=np.zeros(7);eps=2e-5
    for i in range(7):
        plus=np.array(h);minus=np.array(h);plus[i]*=2**eps;minus[i]*=2**(-eps)
        lp=local(plus);lm=local(minus)
        deriv=sum(weight*(lp[key]-lm[key])/(2*eps) for key,weight in adjoint.items())
        mean[i]=128*deriv/(math.log(2)*mass)
    assert min(mean)>=-1e-6 and abs(sum(mean)-512)<1e-5
    print(dict(tilted_relaxation_endpoint=mean.tolist(),source=pop,
               original_count=inactive,mean_sum=float(sum(mean)),
               not_an_encoder_state_distribution=True),flush=True)


def exact_pair_tests():
    from quartet_closed_pair_population import exact_transition,exact_one
    h=tuple(F(8+i,8) for i in range(7));a=F(3,2);d=F(1,4);y=F(1,8)
    profiles=((8,0,0,0,0,0,0),(0,8,0,0,0,0,0),(0,0,8,0,0,0,0),
              (0,0,0,0,0,0,8),(0,0,0,4,0,0,4),(2,1,1,1,1,1,1))
    checked=0
    for pop in profiles:
        width=sum(pop);denom=math.prod(v**n for v,n in zip(h,pop));values=[]
        for node in range(2*width+1):
            values.append(sum(v*exact_one(mid,h,a,d,F(node),y) for mid,v in exact_transition(pop,a,d,F(node),y).items())/denom)
        # Newton interpolation at 0..2w, converted exactly to monomials.
        exact=[F(0)]*(2*width+1);falling=[F(1)];differences=values
        for k in range(2*width+1):
            for j,c in enumerate(falling):exact[j]+=differences[0]*c/math.factorial(k)
            differences=[b-a for a,b in zip(differences,differences[1:])]
            falling=mul(falling,[-k,1])
        assert min(exact)>=0
        z=[F(int(n>0)) for n in pop];poly=bpower(flagged_polynomial(z,h,a,d,y),width//4)
        for inactive,true in enumerate(exact):
            upper=sum(value*flag_factor(f,inactive,width) for (f,i),value in poly.items()
                      if i==inactive and sum(pop[4:])<=f<=width-pop[3])/multinomial(pop)/denom
            assert upper>=true;checked+=1
    print(dict(exact_mixed_pair_quota_coefficient_checks=checked,source_profiles=len(profiles)),flush=True)


def refit_flagged(pop,h,inactive):
    """Fit the actual joint coefficient, not the unconditioned local sum."""
    initial,x,_=fit(pop,h,inactive);active=[i for i,n in enumerate(pop) if n];omit=max(active,key=lambda i:pop[i]);kept=[i for i in active if i!=omit]
    if not kept:return initial,x,0.
    v=np.array([math.log2(initial[i]) for i in kept]);target=np.array([pop[i]/128 for i in kept]);fmax=512-pop[3]
    factors=np.array([float(flag_factor(f,inactive)) if f>=sum(pop[4:]) else 0. for f in range(fmax+1)])
    a=2**(707/1024);d=2**(-1905/1024);y=1/8;start=time.monotonic()
    def local(v):
        z=np.zeros(7);z[active]=1
        for i,u in zip(kept,v):z[i]=2**u
        return {key:value*x**key[1] for key,value in flagged_polynomial(z,h,a,d,y).items() if value}
    def evaluate(v,need_gradient=True):
        poly=local(v);total=sum(poly.values())
        current=relative_convolution({key:value/total for key,value in poly.items()},(fmax,inactive),127)
        adjoint={}
        for f,i in poly:
            if i<=inactive:adjoint[f,i]=float(current[:fmax+1-f,inactive-i]@factors[f:])
        mass=sum(adjoint[key]*value for key,value in poly.items() if key in adjoint);assert mass>0
        value=math.log2(total)+math.log2(mass/total)/128-target@v
        if not need_gradient:return value,None
        grad=np.zeros(len(v));eps=2e-5
        for i in range(len(v)):
            p=v.copy();m=v.copy();p[i]+=eps;m[i]-=eps;plus=local(p);minus=local(m)
            deriv=sum(weight*(plus[key]-minus[key])/(2*eps) for key,weight in adjoint.items())
            grad[i]=deriv/(math.log(2)*mass)-target[i]
        return value,grad
    inverse=np.eye(len(v));value,grad=evaluate(v)
    for iteration in range(35):
        direction=-inverse@grad
        if grad@direction>=0:direction=-grad;inverse=np.eye(len(v))
        direction/=max(1.,max(abs(direction))/2);rate=1.
        for _ in range(15):
            nxt=np.clip(v+rate*direction,-40,40);step=nxt-v
            if evaluate(nxt,False)[0]<=value+1e-4*grad@step:break
            rate*=.5
        if max(abs(step))<1e-8:break
        nv,ng=evaluate(nxt);delta=ng-grad;curve=delta@step
        if curve>1e-12:
            correction=np.eye(len(v))-np.outer(step,delta)/curve
            inverse=correction@inverse@correction.T+np.outer(step,step)/curve
        else:inverse=np.eye(len(v))
        v,value,grad=nxt,nv,ng
        if iteration%10==9:print(dict(refit_iteration=iteration+1,gradient_residual=float(max(abs(grad))),seconds=time.monotonic()-start),flush=True)
        if max(abs(grad))<1e-7:break
    z=np.zeros(7);z[active]=1
    for i,u in zip(kept,v):z[i]=2**u
    return z,x,float(max(abs(grad)))


def tests():
    a=F(3,2);d=F(1,4);y=F(1,8);checked=0
    for hist in itertools.combinations_with_replacement(range(7),4):
        direct=defaultdict(F)
        for (dest,u,v,act,deaths),ways in kernel(hist):
            for paid in range(v+1):direct[tuple(dest.count(i) for i in range(7)),u+paid]+=ways*math.comb(v,paid)*a**act*d**deaths*y**paid
        if 3 in hist and any(i>=4 for i in hist):assert not direct;checked+=1;continue
        relaxed={((0,)*7,0):F(1)}
        for source in hist:
            z=[F(int(i==source)) for i in range(7)];choices=first(z,a,d,y);after=defaultdict(F)
            for (pop,i),weight in relaxed.items():
                for dest,poly in enumerate(choices):
                    new=list(pop);new[dest]+=1
                    for di,v in enumerate(poly):
                        if v:after[tuple(new),i+di]+=weight*v
            relaxed=dict(after)
        assert all(value<=relaxed.get(key,0) for key,value in direct.items());checked+=1
    for seed in range(1,6):
        z=[F(seed+i,8) for i in range(7)];h=[F(8+i,8) for i in range(7)];x=F(seed,3)
        poly=pair_polynomial(z,h,a,d,y)
        assert sum(v*x**i for i,v in enumerate(poly))==point(z,h,a,d,x,y);checked+=1
    print(dict(exact_first_quartet_domination_and_substitution_checks=checked),flush=True)
    checked=0
    for hist in itertools.combinations_with_replacement(range(7),4):
        direct=defaultdict(F)
        for (dest,u,v,act,deaths),ways in kernel(hist):
            for paid in range(v+1):direct[tuple(dest.count(i) for i in range(7)),u+paid]+=ways*math.comb(v,paid)*a**act*d**deaths*y**paid
        relaxed={((0,)*7,0,0):F(1)}
        for source in hist:
            z=[F(int(i==source)) for i in range(7)];choices=flagged_first(z,a,d,y);after=defaultdict(F)
            for (pop,f,i),weight in relaxed.items():
                for dest,poly in enumerate(choices):
                    new=list(pop);new[dest]+=1
                    for (df,di),v in poly.items():
                        if v:after[tuple(new),f+df,i+di]+=weight*v
            relaxed=dict(after)
        counted=defaultdict(F)
        for (pop,f,i),v in relaxed.items():
            if f in (0,4) and (f or pop[0]):counted[pop,i]+=v
        assert all(value<=counted.get(key,0) for key,value in direct.items());checked+=1
    for width in (4,8,12):
        for f in range(0,width+1,4):
            previous=F(0)
            for i in range(width-f+1):
                value=flag_factor(f,i,width);assert previous<=value<=1;previous=value;checked+=1
    print(dict(exact_flagged_quartet_and_occupancy_checks=checked),flush=True)
    checked=0
    probs={(0,0):F(1,7),(1,0):F(2,7),(0,1):F(3,7),(2,2):F(1,7)};exact={(0,0):F(1)}
    for groups in range(1,5):
        after=defaultdict(F)
        for (f,i),v in exact.items():
            for (df,di),w in probs.items():after[f+df,i+di]+=v*w
        exact=dict(after);upper=relative_convolution(probs,(2*groups,2*groups),groups,True)
        for key,v in exact.items():assert F.from_float(float(upper[key]))>=v;checked+=1
    for tiny in (F(1,1<<1023),F(1,1<<1074),F(1,1<<1075)):
        p={(0,0):tiny,(1,1):1-tiny};upper=relative_convolution(p,(2,2),2,True)
        assert F.from_float(float(upper[0,0]))>=tiny*tiny
        assert F.from_float(float(upper[1,1]))>=2*tiny*(1-tiny);checked+=2
    print(dict(exact_relative_rounding_checks=checked),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(180,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:
        tests();h=np.exp2([0.,-.022766994810343107,-.03148991732888965,-.038897622883545944,
                         -1.9738559013803745,-1.4394120075915136,-1.4317249625558666])
        profiles=((0,0,0,64,0,0,448),(104,176,72,32,16,16,96),(104,120,64,32,16,16,160),(104,56,32,64,32,32,192))
        import sys
        if '--exact-tests' in sys.argv:exact_pair_tests();return
        if '--exit' in sys.argv:
            pop=profiles[1]
            proposal=(np.array([.5322266218627301,1.,.45883641227852656,.2540744589849891,
                                .8243003414401759,.5692101436924697,2.6201648825483423]),.7465851247925207,0.)
            exit_diagnostic(pop,h,262,proposal)
            for inactive in (0,128,512):
                value=flagged_float(pop,h,inactive,True,proposal)
                print(dict(edge_original_count=inactive,directed_pair_log2_upper=str(value),
                           display=float(value),seconds=time.monotonic()-start),flush=True)
            return
        if '--refit' in sys.argv:
            pop=profiles[1];proposal=refit_flagged(pop,h,262)
            for inactive in (262,263):
                value=flagged_float(pop,h,inactive,True,proposal)
                print(dict(source=pop,inactive=inactive,refitted_joint_coefficient=True,
                           source_markers=list(map(float,proposal[0])),inactive_marker=float(proposal[1]),
                           directed_pair_log2_upper=str(value),display=float(value),
                           seconds=time.monotonic()-start),flush=True)
            return
        if '--relative' in sys.argv:
            for pop in profiles:
                proposal=flagged_float(pop,h,262);value=flagged_float(pop,h,262,True)
                print(dict(source=pop,inactive=262,float_proposal=proposal,
                           outward_binary64_pair_log2_upper=str(value),display=float(value),
                           seconds=time.monotonic()-start),flush=True)
            return
        for pop in profiles:
            values=[]
            for inactive in (262,263):
                value,residual=evaluate(pop,h,inactive,True);values.append(value)
                print(dict(source=pop,pair_original_inactive=inactive,directed_raw_pair_log2_upper=str(value),
                           display=float(value),all_destinations=True,source_fit_residual=residual,
                           first_transition_relaxed=True,seconds=time.monotonic()-start),flush=True)
            diagnostic=4*sum(values)+log_upper(-1,math.comb(8192,2100))
            print(dict(source=pop,balanced_pair_label_diagnostic_only=float(diagnostic),
                       not_a_regional_operator=True,not_a_certificate=True),flush=True)
            for inactive in (262,263):
                value,residual=flagged_bound(pop,h,inactive)
                print(dict(source=pop,pair_original_inactive=inactive,shared_full_quartet_flag=True,
                           directed_raw_pair_log2_upper=str(value),display=float(value),
                           all_destinations=True,seconds=time.monotonic()-start),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
