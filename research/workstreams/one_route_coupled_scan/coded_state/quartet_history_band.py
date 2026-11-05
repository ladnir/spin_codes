"""Directed, varying-history band for the unchanged field-quartet code.

Imports the frozen local count; never changes that checkpoint. Positive
fixed-point coefficients, a telescoping potential, and regional quota caps
sum all time-varying histories in the stated band. Not whole-code coverage.
No output files; <16 MiB DP arrays; serial runs with a 60-second watchdog.
"""
import argparse
from fractions import Fraction as F
import math
import os
from pathlib import Path
import shutil
import threading
import time

import numpy as np
import quartet_two_update as q

BRIDGE_MARKERS=('112627/167924','2042267/1022811','97905/432928','173812/965177','13592/664659')
BRIDGE_UPPER=F(-106287443107715149,10995116277760)


def configure():
    q.THETA=F(7,40)
    q.P,q.ZERO,q.HOLD,q.NEXT,q.DEATH=map(F,q.FROZEN_MARKERS)
    q.H,q.X,q.Y=q.HOLD*q.NEXT,q.NEXT,q.DEATH/q.NEXT
    q.ROOTS.clear()


def coefficient_box(probabilities,target,precision=31,power=128):
    """All positive coefficients in a bounded box, rounded upward.

    This is the frozen DP with its final array retained rather than one cell.
    At most current, next, and one sliced-product temporary coexist.
    """
    scale=1<<precision
    integers={k:(v.numerator*scale+v.denominator-1)//v.denominator
              for k,v in probabilities.items()}
    shape=tuple(x+1 for x in target)
    assert math.prod(shape)*8*3 < 16*(1<<20)
    current=np.zeros(shape,dtype=np.uint64);current[0,0,0]=scale
    nxt=np.zeros_like(current);maximum=scale;total=sum(integers.values())
    for step in range(power):
        assert maximum*total+scale-1 < 1<<64
        nxt.fill(0)
        for (a,b,d),value in integers.items():
            if any(k>n for k,n in zip((a,b,d),target)):continue
            extent=tuple(min(n-k,4*step)+1 for k,n in zip((a,b,d),target))
            source=tuple(slice(0,e) for e in extent)
            dest=tuple(slice(k,k+e) for k,e in zip((a,b,d),extent))
            nxt[dest]+=current[source]*np.uint64(value)
        nxt+=np.uint64(scale-1);nxt>>=np.uint64(precision)
        maximum=(maximum*total+scale-1)//scale
        assert int(nxt.max())<=maximum
        current,nxt=nxt,current
    return current


def cheap_log_upper(n,precision=31):
    """log2(n/2^precision) upper on a 1/16 grid, exact integer powers."""
    assert n>0
    power=n**16
    ceil=power.bit_length()-(power&(power-1)==0)
    assert power<=1<<ceil
    return F(ceil,16)-precision


def row_constants(r,threshold,zero_free=False):
    # The binomial factor is deliberately kept: all exact inactive row sets.
    pi=math.comb(q.M,r)*q.P**r*(1-q.P)**(q.M-r)
    base=q.log_upper(1,math.comb(q.M,r))+q.THETA+32*q.THETA*(-q.K+128*r)
    if not zero_free:
        base+=32*q.THETA*threshold+q.log_upper(q.THETA*threshold,q.ZERO)
    return base,q.log_interval(pi)[0]


def band_bound(box,total,weights,deaths,rows,distance,quota=True,starts=(1,1),zero_free=False):
    lo,hi=weights;dl,dh=deaths;rl,rh=rows
    assert 0<lo<=hi<=512 and 0<=dl<=dh<=min(lo,512-hi)
    threshold=distance*q.N//100
    assert 1<=starts[0]<=starts[1]<=2048
    # If A+D-W<0, the same formula overbounds the zero-selection-free
    # count because q*ell>=1 and every feasible local polynomial has P>=1.
    assert zero_free or (1<<32)*q.ZERO>=1
    assert hi<box.shape[0] and hi-lo+dh<box.shape[1] and dh<box.shape[2]
    # f(w)=X^w ell^(theta*w). Multiplying T(w,v) by f(v)/f(w)
    # removes v from every marker except the extracted coefficient.
    lp=q.log_interval(q.X)
    assert lp[1]-32*q.THETA<0 if zero_free else lp[1]+q.THETA*q.log_interval(q.ZERO)[1]<0
    initial=max(q.log_upper(1,math.comb(512,w))
                +q.log_upper(w,(1-q.P)*q.X) for w in range(lo,hi+1))
    initial+=q.log_upper(1,hi-lo+1)+q.log_upper(-hi,q.X)
    initial+=32*q.THETA*hi if zero_free else q.log_upper(-q.THETA*hi,q.ZERO)
    common=q.log_upper(128,total)
    worst=None;where=None
    for w in range(lo,hi+1):
        base=common+q.log_upper(-1,math.comb(512,w))+q.log_upper(-w,q.H)
        base+=32*q.THETA*w if zero_free else q.log_upper(-q.THETA*w,q.ZERO)
        for d in range(dl,dh+1):
            # Sum every possible next weight before upper-bounding deaths.
            mass=sum(int(box[w,v-w+d,d]) for v in range(lo,hi+1)
                     if 0<=v-w+d<=512-w)
            if not mass:continue
            term=base+cheap_log_upper(mass)+q.log_upper(-d,q.X*q.Y)
            if not zero_free:term+=-32*q.THETA*d+q.log_upper(-q.THETA*d,q.ZERO)
            if worst is None or term>worst:worst,where=term,(w,d)
    step=worst+q.log_upper(1,dh-dl+1)
    # In r, -(255)*log binom(M,r) plus affine terms is convex.
    # Therefore endpoints bound the entire integer interval.
    service=None
    if quota:
        dual=list(map(F,q.FROZEN_DUAL))
        service=min(sum(c*x for c,x in zip(dual,[128,w,v-w+d,d]))
                    for w in (lo,hi) for v in (lo,hi) for d in (dl,dh))
        assert service>0
    endpoint_rows=[row_constants(r,threshold,zero_free) for r in (rl,rh)]
    cap={n:q.log_upper(1,q.atom_cap(n*service)) if n and quota else F(0)
         for n in range(17)}
    upper=None;worst_start=None
    for start in range(starts[0],starts[1]+1):
        steps=q.L-start+1;regions=(steps+15)//16
        # Every complete zero prefix region has numerator quota mass pi,
        # which cancels its conditioning cost exactly. No route is changed.
        fixed=max(base-regions*pi_lower for base,pi_lower in endpoint_rows)
        first_transitions=(steps-1)%16
        value=(fixed+initial+(steps-1)*step+q.log_upper(1,rh-rl+1)
               +cap[first_transitions]+(regions-1)*cap[16])
        if upper is None or value>upper:upper,worst_start=value,start
    upper+=q.log_upper(1,starts[1]-starts[0]+1)
    return dict(weights=weights,deaths=deaths,inactive_rows=rows,distance_percent=distance,
                all_varying_histories=True,whole_code=False,upper=str(upper),
                upper_float=float(upper),worst_step=where,quota_service=str(service),
                activation_steps=starts,worst_start=worst_start,zero_free=zero_free)


def marker_proposal(zero_free=True):
    """FLOAT proposal only; the caller replays its rational markers."""
    from scipy.optimize import minimize
    from scipy.special import logsumexp
    from diagnose_pairing_sufficiency import local_records,choose_bits
    records=[]
    for a,rows in enumerate(local_records(4,True)):
        for row in rows:
            if not row[4]:continue
            t,c,v,d=row[:4];b=v-a+d
            records.append((a,b,d,t,c,math.log2(math.comb(4,a)),*row[4:]))
    data=np.array(records)
    lo,hi,dl,dh=216,232,24,40;m=float(32*q.THETA)
    def objective(values):
        odds,h,x,y=values[:4]
        ell=0 if zero_free else 2**values[4]
        logell=0 if zero_free else values[4]
        lp=-float(np.logaddexp2(0,-odds));lc=-float(np.logaddexp2(0,odds))
        terms=data[:,0]*h+data[:,1]*x+data[:,2]*y+data[:,3]*lp+(4-data[:,3])*lc-m*data[:,4]+data[:,5]
        terms+=float(q.THETA)*np.log2(data[:,6:] @ np.array([ell**i for i in range(5)]))
        total=float(logsumexp(terms*math.log(2))/math.log(2))
        step=max(128*total-choose_bits(512,w)-w*h-d*(x+y)
                 +(m*w if zero_free else -float(q.THETA)*w*logell-m*d-float(q.THETA)*d*logell)
                 for w in (lo,hi) for d in (dl,dh))+math.log2((hi-lo+1)*(dh-dl+1))
        initial=max(choose_bits(512,w)+w*(lc+x) for w in range(lo,hi+1))
        initial+=math.log2(hi-lo+1)+max(-w*(x-m if zero_free else x+float(q.THETA)*logell) for w in (lo,hi))
        starts=(2048,) if zero_free else (129,640)
        threshold=45*q.N//100
        result=max(choose_bits(8192,r)-((q.L-start+16)//16)*(choose_bits(8192,r)+r*lp+(8192-r)*lc)
                   +float(q.THETA)+m*(-q.K+128*r)+(q.L-start)*step
                   +(0 if zero_free else m*threshold+float(q.THETA)*threshold*logell)
                   for r in (1850,2100) for start in starts)
        return (result+initial+math.log2(251))/q.K
    initial=[math.log2(float(q.P/(1-q.P))),math.log2(float(q.H)),math.log2(float(q.X)),math.log2(float(q.Y))]
    if not zero_free:initial.append(math.log2(float(q.ZERO)))
    fit=minimize(objective,initial,method='Nelder-Mead',options={'maxiter':360,'xatol':1e-7,'fatol':1e-9})
    q.P=F(1/(1+2**-float(fit.x[0]))).limit_denominator(1<<20)
    q.H,q.X,q.Y=[F(2**float(v)).limit_denominator(1<<20) for v in fit.x[1:4]]
    q.ZERO=F(0) if zero_free else F(2**float(fit.x[4])).limit_denominator(1<<20)
    q.ROOTS.clear()
    print(dict(proposal_upper=float(fit.fun*q.K),zero_free=zero_free,p=str(q.P),h=str(q.H),x=str(q.X),y=str(q.Y),
               ell=str(q.ZERO),optimizer_success=bool(fit.success),float_only=True),flush=True)


def tests():
    import itertools
    from collections import defaultdict
    probabilities={(0,0,0):F(1,3),(1,0,0):F(1,3),(0,1,0):F(1,3)}
    box=coefficient_box(probabilities,(3,3,0),precision=12,power=3)
    for a in range(4):
        for b in range(4):
            expected=F(math.factorial(3),math.factorial(a)*math.factorial(b)*math.factorial(3-a-b)*27) if a+b<=3 else F(0)
            assert F(int(box[a,b,0]),1<<12)>=expected
    for n in (1,2,3,7,100,287823,1<<31):
        assert cheap_log_upper(n)>=q.log_interval(F(n,1<<31))[0]
    # Exact potential telescoping and row-bound composition, including
    # all nonconstant three-step paths and a free terminal state.
    matrix=((F(1,4),F(1,3)),(F(1,5),F(1,6)));potential=(F(2),F(3))
    beta=max(sum(matrix[i][j]*potential[j] for j in range(2))/potential[i] for i in range(2))
    for i in range(2):
        actual=sum(matrix[i][j]*matrix[j][k]*matrix[k][l]
                   for j in range(2) for k in range(2) for l in range(2))
        assert actual<=beta**3*potential[i]/min(potential)
    for n in range(2,30):
        for r in range(1,n):
            assert math.comb(n,r)**2>=math.comb(n,r-1)*math.comb(n,r+1)
    # Independent partition enumeration with varying next weight, retaining
    # target labels and the actual local compatibility polynomial.
    def weight(before,after,target):
        coeff=q.choices(2,before,after,target,True)
        if not any(coeff):return F(0)
        c=target.bit_count() if before else (target&after).bit_count()
        t=target.bit_count()
        return F(1,3)**t*F(2,3)**(2-t)*F(1,4)**c*sum(n*F(1,3)**z for z,n in enumerate(coeff))
    poly=defaultdict(F)
    for a in range(3):
        before=(1<<a)-1
        for after,target in itertools.product(range(4),repeat=2):
            b=(after&~before).bit_count();d=(before&~after).bit_count()
            poly[a,b,d]+=math.comb(2,a)*weight(before,after,target)
    square=defaultdict(F)
    for x,v in poly.items():
        for y,u in poly.items():square[tuple(a+b for a,b in zip(x,y))]+=v*u
    direct=defaultdict(F);incoming={0,1}
    for group in itertools.combinations(range(4),2):
        groups=(group,tuple(x for x in range(4) if x not in group))
        for after,target in itertools.product(range(16),repeat=2):
            value=F(1)
            for block in groups:
                local=lambda mask:sum(((mask>>x)&1)<<j for j,x in enumerate(block))
                value*=weight(local(3),local(after),local(target))
            b=(after&~3).bit_count();d=(3&~after).bit_count()
            direct[b,d]+=value/F(6)
    for b,d in itertools.product(range(3),repeat=2):
        assert direct[b,d]==square[2,b,d]/math.comb(4,2)
    # All-zero prefix regions have no local restriction on target labels.
    for r in range(5):
        mass=sum(F(1,3)**r*F(2,3)**(4-r) for mask in range(16) if mask.bit_count()==r)
        assert mass==math.comb(4,r)*F(1,3)**r*F(2,3)**(4-r)
    for z in range(-8,1):assert (4*F(1,3))**(-z)>=1
    print(dict(box_checks=16,log_checks=7,potential_checks=2,binomial_concavity_checks=406,
               labeled_partition_cases=1536,varying_transition_identities=9,
               zero_prefix_quota_checks=5,negative_zero_count_checks=9),flush=True)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--tests-only',action='store_true')
    parser.add_argument('--late-proposal',action='store_true')
    parser.add_argument('--bridge-proposal',action='store_true')
    parser.add_argument('--frozen',action='store_true')
    args=parser.parse_args()
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(60,lambda:os._exit(124));timer.daemon=True;timer.start()
    start=time.monotonic()
    try:
        configure();tests()
        if args.tests_only:return
        if args.frozen:
            q.P,q.H,q.X,q.Y,q.ZERO=map(F,BRIDGE_MARKERS);q.ROOTS.clear()
            probabilities,total,_=q.local()
            box=coefficient_box(probabilities,(232,56,40))
            result=band_bound(box,total,(216,232),(24,40),(1850,2100),45,
                              quota=False,starts=(1,2048))
            assert F(result['upper'])==BRIDGE_UPPER and BRIDGE_UPPER < -9666
            print(result,flush=True)
            print(dict(seconds=time.monotonic()-start,files_written=0),flush=True)
            return
        if args.late_proposal or args.bridge_proposal:
            marker_proposal(zero_free=args.late_proposal)
            probabilities,total,_=q.local()
            box=coefficient_box(probabilities,(232,56,40))
            for starts in (((2048,2048),(640,2048),(1,2048)) if args.late_proposal
                           else ((129,640),(129,1023),(1,2048))):
                print(band_bound(box,total,(216,232),(24,40),(1850,2100),45,
                                 quota=False,starts=starts,zero_free=args.late_proposal),flush=True)
            return
        probabilities,total,records=q.local()
        # Verifies all 35 exact dual inequalities before reusing the dual.
        q.quota_bound(records,frozen=True)
        box=coefficient_box(probabilities,(232,56,40))
        assert F(int(box[224,32,32]),1<<31)==F(287823,1<<31)
        for weights,deaths,rows in (((220,228),(28,36),(1900,2048)),
                                    ((216,232),(24,40),(1850,2100))):
            for distance in (45,46):
                print(band_bound(box,total,weights,deaths,rows,distance),flush=True)
        for starts in ((1,16),(1,128),(1,2048)):
            print(band_bound(box,total,(216,232),(24,40),(1850,2100),45,starts=starts),flush=True)
        print(dict(seconds=time.monotonic()-start,files_written=0),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
