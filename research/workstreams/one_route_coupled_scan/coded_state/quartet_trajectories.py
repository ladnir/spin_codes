"""Bounded varying-activity quartet diagnostic, plus directed family replays.

Extends the joint count in RANDOM_PAIRING_SUFFICIENCY.md. Each tested weight
and death schedule is fixed, not selected after seeing the code. The bound
unions all labeled supports/messages consistent with that schedule. It does
not cover arbitrary schedules or give a whole-code certificate. No files.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
import itertools
import math
import os
from pathlib import Path
import shutil
import threading
import time

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp
from diagnose_pairing_sufficiency import choices as baseline_choices, local_records, choose_bits, LN2

N, K, S, L, M, G = 1 << 21, 1 << 20, 512, 4096, 8192, 4
RECORDS = local_records(G, True)
choices=baseline_choices
# Compact selected certificates, not raw optimizer/timing dumps. Transition
# marker triples follow the sorted (incoming weight, next weight, deaths) list.
FROZEN = {
    'alternate-240-272': ('708/971','122381/399341',
        ('3665126/538619','169780/396939','23127/718958','18611349/946322','128591/942344','28109/913427'),
        '-7510553/64'),
    'alternate-224-288': ('668/887','411007/835602',
        ('3583751/958530','830877/995885','5269/76017','17120015/983674','144000/946867','37685/597186'),
        '-80139747/256'),
    '224-at-45': ('387/577','153435/995816',
        ('6101401/666363','75953/338894','10966/285397'),'-1830657/32'),
}


def mixed_choices(g,before,after,target,strict):
    """NEW setup condition: every A+delta D, binary diagonal D, is superregular.

    Not a consequence of accepting only A and A+delta I. A field-multiplier
    witness and rejection bound must be checked separately for this variant.
    """
    original=baseline_choices(g,before,after,target,strict)
    if not any(original) or not before:return original
    result=[0]*(g+1)
    for zeros in range(1<<g):
        if zeros & ~(before&after):continue
        rows=target & (~after | zeros)
        if rows.bit_count()>=before.bit_count():continue
        result[zeros.bit_count()]+=1
    assert all(a<=b for a,b in zip(result,original))
    return tuple(result)


def mixed_records():
    groups=[]
    for a in range(G+1):
        before=(1<<a)-1;records=[]
        for after,target in itertools.product(range(1<<G),repeat=2):
            coeff=mixed_choices(G,before,after,target,True)
            if not any(coeff):continue
            covered=target.bit_count() if a else (target&after).bit_count()
            records.append((target.bit_count(),covered,after.bit_count(),(before&~after).bit_count(),*coeff))
        groups.append(np.array(records,dtype=float))
    return groups


def profile(name, cycle, deaths, rows=1974, percent=46):
    weights = [cycle[t % len(cycle)] for t in range(L)]
    transitions = Counter((weights[t-1], weights[t], deaths[(t-1) % len(deaths)])
                          for t in range(1, L))
    for (u,v,d), count in transitions.items():
        assert count > 0 and 0 <= d <= u and 0 <= u-d <= v
        assert v-u+d <= S-u
    A = sum(weights)
    D = sum(count*d for (u,v,d), count in transitions.items())
    z = max(0, A+D-percent*N//100)
    possible = sum(count*(u-d) for (u,v,d), count in transitions.items())
    return dict(name=name, first=weights[0], transitions=sorted(transitions.items()),
                A=A, D=D, z=z, rows=rows, impossible=z>possible)


def objective(p, free_theta=False, active_floor=False):
    # Keep the already replayable theta=3/16. Optimize the positive markers,
    # not the construction. Incoming-coefficient Chernoff bound avoids arrays.
    r, z = p['rows'], p['z']
    row_bits = choose_bits(M,r)
    def evaluate(values):
        m = float(values[-1]) if free_theta else 6
        theta = m/32
        active = float(values[2+3*len(p['transitions'])]) if active_floor else 0
        odds, zero = values[:2]
        lp = -float(np.logaddexp2(0,-odds))
        lc = -float(np.logaddexp2(0,odds))
        pi = row_bits+r*lp+(M-r)*lc
        total = (row_bits+choose_bits(S,p['first'])+p['first']*lc-256*pi+theta
                 +m*(p['A']-K+128*r-z)-theta*z*zero
                 +active*(p['first']-129*(M-r)))
        for i,((u,v,d), count) in enumerate(p['transitions']):
            incoming,nxt,death = values[2+3*i:5+3*i]
            terms=[]
            for a,data in enumerate(RECORDS):
                poly=data[:,4:] @ np.exp2(zero*np.arange(G+1))
                logs=(data[:,0]*lp+(G-data[:,0])*lc-m*data[:,1]
                      +nxt*data[:,2]+death*data[:,3]+theta*np.log2(poly)
                      +active*(G-data[:,0] if a else data[:,2]))
                terms.append(float(logsumexp(logs*LN2)/LN2)+math.log2(math.comb(G,a))+a*incoming)
            kernel=(S//G)*float(logsumexp(np.array(terms)*LN2)/LN2)
            kernel-=choose_bits(S,u)+u*incoming+v*nxt+d*death
            total+=count*kernel
        return total
    return evaluate


def replay(p, values, active_floor=False, frozen=None):
    """All rounding is outward, using rational markers and integer powers."""
    theta=F(3,16);m=6
    length=3*len(p['transitions'])
    if frozen is None:
        probability=F(1/(1+2**-float(values[0]))).limit_denominator(1000)
        zero=F(2**float(values[1])).limit_denominator(1000000)
        markers=[F(2**float(x)).limit_denominator(1000000) for x in values[2:2+length]]
    else:
        assert not active_floor
        probability,zero=map(F,frozen[:2]);markers=list(map(F,frozen[2]))
        assert len(markers)==length
    active=F(2**float(values[2+length])).limit_denominator(1000000) if active_floor else F(1)
    assert active>=1
    assert 0<probability<1 and zero>0 and all(x>0 for x in markers)
    roots={}
    def root(value):
        if value not in roots:
            cube=value**3;num=cube.numerator << 640;den=cube.denominator
            lo,hi=0,1<<42
            while hi**16*den<num:hi*=2
            while hi-lo>1:
                mid=(lo+hi)//2
                if mid**16*den>=num:hi=mid
                else:lo=mid
            roots[value]=F(hi,1<<40)
            assert roots[value]**16>=cube
        return roots[value]
    def term(c,value):
        value=F(value);assert value>0
        if value.numerator&(value.numerator-1)==0 and value.denominator&(value.denominator-1)==0:
            return c*(value.numerator.bit_length()-value.denominator.bit_length())
        bits=(value.numerator**16).bit_length()-(value.denominator**16).bit_length()
        return c*F(bits+(1 if c>=0 else -1),16)
    assert term(-1000000,F(1))==0 and term(-3,F(1,8))==9
    r,z=p['rows'],p['z']
    pi=math.comb(M,r)*probability**r*(1-probability)**(M-r)
    total=(term(1,math.comb(M,r))+term(1,math.comb(S,p['first']))
           +term(p['first'],1-probability)+term(-256,pi)
           +theta+m*(p['A']-K+128*r-z)+term(-theta*z,zero)
           +term(p['first']-129*(M-r),active))
    for i,((u,v,d), count) in enumerate(p['transitions']):
        h,nxt,death=markers[3*i:3*i+3];poly=F(0)
        for a in range(G+1):
            before=(1<<a)-1;local=F(0)
            for after,target in itertools.product(range(1<<G),repeat=2):
                coefficients=choices(G,before,after,target,True)
                if not any(coefficients):continue
                pol=sum(c*zero**j for j,c in enumerate(coefficients))
                tc=target.bit_count();covered=tc if a else (target&after).bit_count()
                local+=(probability**tc*(1-probability)**(G-tc)*F(1,1<<(m*covered))
                        *nxt**after.bit_count()*death**(before&~after).bit_count()*root(pol)
                        *active**(G-tc if a else after.bit_count()))
            poly+=math.comb(G,a)*local*h**a
        kernel=poly**(S//G)/(math.comb(S,u)*h**u*nxt**v*death**d)
        total+=term(count,kernel)
    return dict(log2_upper_rational=str(total),log2_upper=float(total),directed=True,
                p=str(probability),zero=str(zero),markers=list(map(str,markers)),active_marker=str(active),
                whole_code=False)


def tests():
    from diagnose_pairing_sufficiency import objective as old_objective
    p=profile('test',[256],[32])
    assert p['A']==K and p['D']==4095*32 and p['z']==214927
    varying=profile('test',[240,272],[16,48])
    assert varying['A']==K and varying['D']==131024
    assert profile('reset',[0,256],[0,256])['impossible']
    # Replacing an exact incoming coefficient by evaluation at positive h
    # can only increase the original constant-profile bound.
    for h in (2,12,64):
        v=[math.log2(8/3),math.log2(5/16),math.log2(h),math.log2(5/18),math.log2(1/24)]
        assert objective(p)(v)>=old_objective(4,True)([6,v[0],v[1],v[3],v[4]])-1e-6
        assert abs(objective(p,active_floor=True)([*v,0])-objective(p)(v))<1e-6
    print(dict(profile_checks='pass',original_coefficient_comparisons=3),flush=True)


def ideal_bound(p):
    """Integer upper exponent for fully uniform 512x512 raw inner matrices.

    Only for nonempty states throughout. Uses deliberately loose support and
    zero-set counts <=2^N and <=2^(A-w1), respectively. Not our quartet code.
    """
    assert p['first']>0 and all(u>0 and v>0 for (u,v,d),count in p['transitions'])
    exponent=K-128*p['rows']+p['A']-p['first']-p['z']-S*(L-1)
    return M+N+p['A']-p['first']+32*exponent


def main():
    global choices,RECORDS
    parser=argparse.ArgumentParser();parser.add_argument('--replay',action='store_true')
    parser.add_argument('--mixed-acceptance',action='store_true')
    parser.add_argument('--hard-ladder',action='store_true')
    parser.add_argument('--frozen',action='store_true')
    args=parser.parse_args()
    if args.mixed_acceptance:choices=mixed_choices;RECORDS=mixed_records()
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(60,lambda:os._exit(124));timer.daemon=True;timer.start()
    start=time.monotonic()
    try:
        if not args.mixed_acceptance:tests()
        if args.frozen:
            assert not args.mixed_acceptance
            selected=[profile('alternate-240-272',[240,272],[16,48]),
                      profile('alternate-224-288',[224,288],[16,80]),
                      profile('224-at-45',[224],[32],percent=45)]
            for p in selected:
                result=replay(p,(),frozen=FROZEN[p['name']])
                assert F(result['log2_upper_rational'])==F(FROZEN[p['name']][3])
                assert result['log2_upper'] < -40
                print(dict(name=p['name'],**result),flush=True)
            return
        profiles=[profile('constant-224',[224],[32]),profile('constant-256',[256],[32]),
                  profile('constant-288',[288],[32]),profile('alternate-240-272',[240,272],[16,48]),
                  profile('alternate-224-288',[224,288],[16,80]),
                  profile('constant-256-sticky',[256],[0]),
                  profile('constant-256-replace',[256],[128]),
                  profile('reset-0-256',[0,256],[0,256]),
                  profile('constant-128',[128],[16]),profile('constant-384',[384],[32])]
        if args.hard_ladder:
            profiles=[profile('224-at-'+str(d),[224],[32],percent=d) for d in (45,46,47,49)]
        for p in profiles:
            if p['impossible']:
                print(dict(name=p['name'],structurally_impossible=True),flush=True);continue
            f=objective(p)
            initial=[math.log2(8/3),math.log2(5/16)]
            for transition,count in p['transitions']:
                initial += [math.log2(12),math.log2(5/18),math.log2(1/24)]
            bounds=[(-10,10),(-12,12)]+[(-16,16)]*(len(initial)-2)
            result=minimize(lambda x:f(x)/K,initial,method='L-BFGS-B',bounds=bounds,
                            options={'maxiter':100,'ftol':1e-12,'gtol':1e-7,'maxls':25})
            record=dict(name=p['name'],A=p['A'],D=p['D'],z=p['z'],float_log2_upper=f(result.x),
                        converged=bool(result.success),evaluations=result.nfev,whole_code=False,
                        mixed_acceptance=args.mixed_acceptance)
            if args.replay:record.update(replay(p,result.x))
            if args.mixed_acceptance:record['requires_stronger_setup_acceptance']=True
            if record['float_log2_upper']>=0:
                record['ideal_dense_inner_log2_upper']=ideal_bound(p)
                flexible=objective(p,True)
                tuned=minimize(lambda x:flexible(x)/K,[*result.x,6],method='L-BFGS-B',
                               bounds=[*bounds,(0.001,32)],options={'maxiter':160,'ftol':1e-12,'gtol':1e-7})
                record.update(float_free_theta_upper=flexible(tuned.x),theta=float(tuned.x[-1])/32)
                floor=objective(p,active_floor=True)
                floor_result=minimize(lambda x:floor(x)/K,[*result.x,0],method='L-BFGS-B',
                                      bounds=[*bounds,(0,12)],options={'maxiter':160,'ftol':1e-12,'gtol':1e-7})
                record.update(float_active_floor_upper=floor(floor_result.x),active_marker_log2=float(floor_result.x[-1]))
                if args.replay:record['active_floor_replay']=replay(p,floor_result.x,active_floor=True)
            print(record,flush=True)
        print(dict(seconds=time.monotonic()-start,schedules=len(profiles),scope='fixed schedules, one inactive-row count'),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
