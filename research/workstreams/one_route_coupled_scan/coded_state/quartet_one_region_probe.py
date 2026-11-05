"""Full-width regional transfer diagnostic; NOT a directed certificate.

Sequential importance sampling of the exact seven-type positive kernel.
Carries every sampled population, the original-inactive quota, and pending
types across region boundaries. No population truncation or early terminal.
The normalizing-constant estimator is unbiased in exact arithmetic; its
logarithm and finite sample output are not upper bounds. Serial, file-free.
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
from quartet_closed_quotient import kernel


H=np.exp2([0.,-.022766994810343107,-.03148991732888965,
           -.038897622883545944,-1.9738559013803745,
           -1.4394120075915136,-1.4317249625558666])
PROFILES=((104,176,72,32,16,16,96),(104,120,64,32,16,16,160),
          (512,0,0,0,0,0,0))
A=2**(707/1024)
D=2**(-1905/1024)


def tables():
    histories=list(itertools.combinations_with_replacement(range(7),4))
    hist=np.array([[row.count(i) for i in range(7)] for row in histories],dtype=np.int16)
    rows={};maxlen=0;total=np.zeros((len(hist),5))
    for i,row in enumerate(histories):
        choices=[[] for _ in range(5)]
        for (dest,u,v,a,d),ways in kernel(row):
            new=np.array([dest.count(t) for t in range(7)],dtype=np.int16)
            potential=float(np.prod(H**(new-hist[i])))
            for c in range(v+1):
                value=ways*math.comb(v,c)*A**a*D**d*(1/8)**c*potential
                choices[u+c].append((new,value,c,a,d))
        for j,entries in enumerate(choices):
            rows[i,j]=entries;maxlen=max(maxlen,len(entries))
            total[i,j]=sum(v for new,v,c,a,d in entries)
    dest=np.zeros((len(hist),5,maxlen,7),dtype=np.int16)
    stats=np.zeros((len(hist),5,maxlen,3),dtype=np.int16)
    cdf=np.ones((len(hist),5,maxlen))
    for (i,j),entries in rows.items():
        cumulative=0.
        for k,(new,value,c,a,d) in enumerate(entries):
            cumulative+=value/total[i,j];cdf[i,j,k]=cumulative
            dest[i,j,k]=new;stats[i,j,k]=(c,a,d)
        if entries:cdf[i,j,len(entries)-1:]=1.
    choose=np.zeros((8193,5));choose[:,0]=1
    for k in range(1,5):choose[k:,k]=choose[k:,k-1]*(np.arange(k,8193)-k+1)/k
    return hist,total,dest,stats,cdf,choose


def quota_probabilities(remaining,positions,choose):
    # Probability of each particular inactive subset with cardinality j,
    # not the hypergeometric probability of cardinality j itself.
    return np.array([choose[remaining,j]*choose[positions-remaining,4-j]
                     /(choose[positions,4]*math.comb(4,j)) for j in range(5)]).T


def resample(weights,rng):
    # Independent multinomial resampling, after normalizing all live mass.
    cumulative=np.cumsum(weights);cumulative/=cumulative[-1];cumulative[-1]=1.
    return np.searchsorted(cumulative,rng.random(len(weights)),side='right')


def probe(pop,particles,seed,regions=1):
    hist,total,dest,stats,cdf,choose=tables();rng=np.random.default_rng(seed)
    population=np.tile(np.array(pop,dtype=np.int16),(particles,1))
    accumulated=np.zeros((particles,3),dtype=np.int32)
    logz=0.;start=time.monotonic();minimum_ess=particles
    width=sum(pop);quota=2100 if width==512 else 16
    assert width in (4,512) and particles<=1024 and regions<=2
    for region in range(regions):
        remaining=np.full(particles,quota,dtype=np.int32);positions=16*width
        before=logz
        for step in range(16):
            available=population.copy();outgoing=np.zeros_like(population)
            for group in range(width//4):
                slots=width-4*group
                probability=np.ones((particles,len(hist)))
                for t in range(7):probability*=choose[available[:,t,None],hist[None,:,t]]
                probability/=choose[slots,4]
                qp=quota_probabilities(remaining,positions,choose)
                local=qp@total.T
                group_weights=probability*local;weights=group_weights.sum(axis=1)
                mean=float(weights.mean())
                if mean==0:
                    print(dict(source=pop,particles=particles,seed=seed,region=region,
                               step=step,group=group,particle_extinction=True,
                               not_a_zero_operator=True),flush=True)
                    return
                assert math.isfinite(mean) and np.all(weights>=0)
                logz+=math.log2(mean)
                ess=float(weights.sum()**2/(weights@weights));minimum_ess=min(minimum_ess,ess)
                parents=resample(weights,rng)
                available=available[parents];outgoing=outgoing[parents]
                remaining=remaining[parents];accumulated=accumulated[parents]
                group_weights=group_weights[parents];local=local[parents];qp=qp[parents]
                cumulative=np.cumsum(group_weights,axis=1)
                pick=np.sum(cumulative<rng.random(particles)[:,None]*cumulative[:,-1,None],axis=1)
                assert np.all(pick<len(hist))
                js=qp*total[pick];jc=np.cumsum(js,axis=1)
                j=np.sum(jc<rng.random(particles)[:,None]*jc[:,-1,None],axis=1)
                outcome=np.sum(cdf[pick,j]<rng.random(particles)[:,None],axis=1)
                available-=hist[pick];outgoing+=dest[pick,j,outcome]
                remaining-=j;accumulated+=stats[pick,j,outcome]
                positions-=4
                assert np.all(available>=0) and np.all(remaining>=0) and np.all(remaining<=positions)
            assert not available.any() and np.all(outgoing.sum(axis=1)==width)
            population=outgoing
        assert positions==0 and not remaining.any()
        print(dict(source=pop,particles=particles,seed=seed,region=region+1,
                   float_region_log2_estimate=logz-before,float_chain_log2_estimate=logz,
                   endpoint_mean=population.mean(axis=0).tolist(),
                   selected_range=[int(population[:,4:].sum(axis=1).min()),int(population[:,4:].sum(axis=1).max())],
                   tilted_paid_activity_deaths_mean=accumulated.mean(axis=0).tolist(),
                   minimum_quartet_ess=minimum_ess,seconds=time.monotonic()-start,
                   quota_exact=True,population_truncation=False,terminal_applied=False,
                   importance_estimate_not_upper_bound=True),flush=True)


def exact_toy(regions=2):
    """Complete width4 operator, sixteen updates per region, rational weights.

    No destination population is dropped. This tests composition and the
    estimator, not the requested width512 bound. Raw integer scaling avoids
    rounding through the entire calculation.
    """
    from collections import defaultdict
    scale=1<<48;rows={}
    for hist in itertools.combinations_with_replacement(range(7),4):
        row=defaultdict(int)
        for (dest,u,v,a,d),ways in kernel(hist):
            for c in range(v+1):
                weight=scale*ways*math.comb(v,c)*F(13,8)**a*F(9,32)**d*F(1,8)**c
                assert weight.denominator==1
                row[dest,u+c]+=weight.numerator
        rows[hist]=tuple(row.items())
    state={((0,0,0,0),0):1};normal=1;peak=1;start=time.monotonic()
    for region in range(regions):
        for step in range(16):
            after=defaultdict(int)
            for (src,used),value in state.items():
                for (dest,j),weight in rows[src]:
                    quota=used+j
                    if quota>16 or quota+4*(15-step)<16:continue
                    after[dest,0 if step==15 else quota]+=value*weight
            state=dict(after);normal*=scale;peak=max(peak,len(state))
        normal*=math.comb(64,16)
        value=F(sum(state.values()),normal)
        exponent=value.numerator.bit_length()-value.denominator.bit_length()
        if value<F(2)**exponent:exponent-=1
        assert F(2)**exponent<=value<F(2)**(exponent+1)
        print(dict(exact_width4_region=region+1,steps_per_region=16,quota=16,
                   exact_chain_log2_integer_bracket=[exponent,exponent+1],
                   display=math.log2(value.numerator)-math.log2(value.denominator),
                   endpoint_populations=len(state),peak_states=peak,
                   all_populations_and_quota_allocations=True,
                   terminal_applied=False,seconds=time.monotonic()-start),flush=True)


def tests():
    hist,total,dest,stats,cdf,choose=tables();checked=0
    for n in ((4,2,1,1,0,0,0),(0,0,0,4,0,0,4),(2,1,1,1,1,1,1)):
        exact=sum(F(math.prod(math.comb(x,int(c)) for x,c in zip(n,h)),math.comb(8,4))
                  for h in hist if all(x>=c for x,c in zip(n,h)))
        assert exact==1;checked+=1
    for positions in range(4,17):
        for r in range(positions+1):
            prob=quota_probabilities(np.array([r]),positions,choose)[0]
            for j,p in enumerate(prob):
                exact=F(math.comb(positions-4,r-j),math.comb(positions,r)) if 0<=r-j<=positions-4 else F(0)
                assert abs(p-float(exact))<1e-14;checked+=1
            assert abs(sum(math.comb(4,j)*p for j,p in enumerate(prob))-1)<1e-14
    # Retain outgoing types at an internal boundary. This is an exact width4,
    # four-update rational composition with quota2 every two updates.
    def exact(regional=True,terminal_inside=False):
        from collections import defaultdict
        from quartet_closed_quotient import accepted
        state={((0,0,0,0),0):F(1)}
        for step in range(4):
            after=defaultdict(F)
            for (src,quota),value in state.items():
                for (out,u,v,a,d),ways in kernel(src):
                    for c in range(v+1):
                        j=u+c;running=quota+j
                        if regional and running>2:continue
                        if regional and step%2 and running!=2:continue
                        if not regional and running>4:continue
                        if terminal_inside and step==1 and not accepted(out):continue
                        weight=F(ways*math.comb(v,c))*F(3,2)**a*F(1,4)**d*F(1,8)**c
                        after[out,0 if regional and step%2 else running]+=value*weight
            state=dict(after)
        return {key:value for key,value in state.items() if regional or key[1]==4}
    regional=exact();global_only=exact(False);extra_terminal=exact(True,True)
    assert sum(regional.values())!=sum(global_only.values())
    assert sum(regional.values())!=sum(extra_terminal.values())
    print(dict(exact_normalization_checks=checked,exact_small_boundary_outputs=len(regional),
               global_quota_negative_control=True,early_terminal_negative_control=True),flush=True)


def empty_payment_floor():
    """Exact lower bound for the all-Z row, no terminal imposed on operator.

    Its subset of terminal-accepting, zero-payment canonical paths suffices.
    Every nonterminal selected run has length>=3 and nonnegative score.
    The terminal run contributes at least lambda per last selected bit.
    Every allowed final input is selected at the free terminal. Therefore
    a^A*d^D*Phi_h(end) >= (2/5)^(512-J_last), using a>=8/5 and h>=1/4.
    Average the true hypergeometric J_last; the bound is independent of
    all quartet partitions. This lower bound is not a failure lower bound.
    """
    assert 3*707>1905 and (F(8,5)**1024)<=2**707 and min(H)>=.25 and H[0]==1.
    value=sum(F(math.comb(512,j)*math.comb(7680,2100-j),math.comb(8192,2100))
              *F(2,5)**(512-j) for j in range(513))
    if np.all(H==1):value=F(1)
    exponent=value.numerator.bit_length()-value.denominator.bit_length()
    if value<F(2)**exponent:exponent-=1
    assert F(2)**exponent<=value<F(2)**(exponent+1)
    print(dict(all_Z_row_certified_log2_lower_integer=exponent,
               zero_payment_canonical_subset=True,exact_region_quota=True,
               not_a_code_failure_lower_bound=True),flush=True)
    return exponent


def suffix_floors():
    """Directed full-width regional lower witnesses, exact quota and endpoint.

    Pay precisely the original-inactive positions in a fixed full suffix.
    At least three suffix updates force the final population to all L2.
    These witnesses are closed; no coordinate partition restriction is used.
    Different suffix lengths are NOT summed (their paid sets can overlap).
    """
    a=F(math.floor(A*(1<<24)),1<<24)
    assert a**1024<=2**707
    h=F.from_float(float(H[6]));normal=math.comb(8192,2100);best=None
    for updates in range(3,17):
        length=512*updates;outside=8192-length
        low=max(0,2100-outside);high=min(length,2100)
        left=math.comb(length,low);right=math.comb(outside,2100-low);numerator=0
        for c in range(low,high+1):
            numerator+=(left*right)<<(3*(2100-c))
            if c<high:
                left=left*(length-c)//(c+1)
                right=right*(2100-c)//(outside-(2100-c)+1)
        value=h**512*a**length*F(numerator,normal*(1<<(3*2100)))
        exponent=value.numerator.bit_length()-value.denominator.bit_length()
        if value<F(2)**exponent:exponent-=1
        assert F(2)**exponent<=value<F(2)**(exponent+1)
        if best is None or exponent>best[0]:best=(exponent,updates)
    print(dict(full_width_suffix_certified_log2_lower_integer=best[0],suffix_updates=best[1],
               start_population=[512,0,0,0,0,0,0],end_population=[0,0,0,0,0,0,512],
               quota=2100,steps=16,not_a_code_failure_lower_bound=True),flush=True)


def suffix_tests():
    from quartet_closed_transfer import advance
    from quartet_closed_quotient import PROJECT
    rng=np.random.default_rng(192);checks=0
    for width in (4,8):
        permutations=[rng.permutation(width).reshape(-1,4) for _ in range(4)]
        for suffix in (3,4):
            for indices in itertools.combinations(range(4*width),2):
                inactive=np.zeros((4,width),dtype=bool);inactive.flat[list(indices)]=True
                paid=inactive.copy();paid[:4-suffix]=False
                allowed=~inactive|paid;full=[]
                for t,p in enumerate(permutations):
                    flags=np.zeros(width,dtype=bool);flags[p]=np.all(allowed[t,p],axis=1)[:,None]
                    full.append(flags)
                z=allowed.copy()
                for t in range(3):z[t]&=full[t+1]
                selected=np.zeros_like(z)
                for i in range(width):
                    t=0
                    while t<4:
                        if not z[t,i]:t+=1;continue
                        end=t
                        while end<4 and z[end,i]:end+=1
                        if end==4 or end-t>=3:selected[t:end,i]=True
                        t=end
                covered=selected.copy()
                for t in range(1,4):
                    p=permutations[t];covered[t,p]|=np.any(selected[t-1,p],axis=1)[:,None]
                assert np.array_equal(inactive&covered,paid)
                deaths=np.count_nonzero(selected[:-1]&~selected[1:])
                assert 707*np.count_nonzero(selected)-1905*deaths>=707*width*suffix
                for i in range(width):
                    age=0
                    for t in range(3):age=advance(age,bool(z[t,i]),bool(selected[t,i]))
                    assert PROJECT[3*age+2]==6
                checks+=1
    print(dict(exact_suffix_closure_score_and_endpoint_checks=checks),flush=True)


def main():
    global A,D
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profile',type=int,default=0,choices=range(3))
    parser.add_argument('--particles',type=int,default=128,choices=(64,128,256,512,1024))
    parser.add_argument('--seed',type=int,default=1)
    parser.add_argument('--regions',type=int,default=1,choices=(1,2))
    parser.add_argument('--flat',action='store_true')
    parser.add_argument('--toy',action='store_true')
    parser.add_argument('--floors-only',action='store_true')
    args=parser.parse_args()
    if args.flat:H[:]=1.
    if args.toy:H[:]=1.;A=13/8;D=9/32
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(180,lambda:os._exit(124));timer.daemon=True;timer.start()
    try:
        tests()
        if args.toy:
            exact_toy(args.regions)
            probe((4,0,0,0,0,0,0),args.particles,args.seed,args.regions)
        else:
            empty_payment_floor()
            if args.floors_only:suffix_tests();suffix_floors()
            else:probe(PROFILES[args.profile],args.particles,args.seed,args.regions)
    finally:timer.cancel()


if __name__=='__main__':main()
