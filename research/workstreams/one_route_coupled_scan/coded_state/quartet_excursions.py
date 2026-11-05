"""Bounded exceptional transitions for the unchanged field-quartet code.

Exact positive local envelopes cover every weight, including zero. A good
transition uses the frozen band coefficient DP. No files, no optimizer,
<16 MiB coefficient arrays, 60-second watchdog, serial numerical runs.
"""
from fractions import Fraction as F
import itertools
import math
import os
from pathlib import Path
import shutil
import threading
import time

import quartet_two_update as q
import quartet_history_band as band

FROZEN_EXCEPTION_UPPER=F(-3901196831485731,43980465111040)
FROZEN_RESET_UPPER=F(-13146227725709835,2199023255552)


def good_step(box,total):
    """Directed log2 beta: G f <= beta f on all 513 state weights."""
    common=q.log_upper(128,total);worst=None;where=None
    for w in range(216,233):
        base=(common+q.log_upper(-1,math.comb(512,w))
              +q.log_upper(-w,q.H)+q.log_upper(-q.THETA*w,q.ZERO))
        for d in range(24,41):
            mass=sum(int(box[w,v-w+d,d]) for v in range(216,233))
            term=(base+band.cheap_log_upper(mass)+q.log_upper(-d,q.X*q.Y)
                  -32*q.THETA*d+q.log_upper(-q.THETA*d,q.ZERO))
            if worst is None or term>worst:worst,where=term,(w,d)
    return worst+q.log_upper(1,17),where


def full_step():
    """All-weight envelope without enumerating a 513^3 coefficient cube.

    Local incoming a has tilted row sum L_a * ell^(-theta*a).
    Their maximum to power 128 bounds the actual random-partition row sum.
    Integer powers compare all five rational/algebraic expressions exactly.
    """
    assert q.THETA==F(7,40)
    death=q.root_upper(F(1,(1<<32)*q.ZERO))/q.X
    locals=[]
    for a in range(5):
        before=(1<<a)-1;value=F(0)
        for after,target in itertools.product(range(16),repeat=2):
            coefficients=q.choices(4,before,after,target,True)
            if not any(coefficients):continue
            b=(after&~before).bit_count();d=(before&~after).bit_count()
            t=target.bit_count();c=t if a else (target&after).bit_count()
            poly=sum(n*q.ZERO**j for j,n in enumerate(coefficients))
            value+=(q.P**t*(1-q.P)**(4-t)*q.root_upper(poly/F(1<<(32*c)))
                    *q.X**b*death**d)
        locals.append(value)
    powers=[value**40/q.ZERO**(7*a) for a,value in enumerate(locals)]
    index=max(range(5),key=lambda a:powers[a])
    # The maximizer is checked, never assumed to be the empty quartet.
    upper=q.log_upper(F(128,40),powers[index])
    assert locals[0]==(1+(1-q.P)*q.X)**4
    return upper,index,[q.log_upper(F(1,40),v) for v in powers]


def certificate(beta,eta,budget):
    assert 0<=budget<=32 and eta>=beta
    # Initial and final weights may be ANY nonzero / arbitrary weights.
    # Omitting the -1 in the binomial activation sum only increases it.
    assert q.log_interval(q.X)[1]+q.THETA*q.log_interval(q.ZERO)[1]<0
    initial=(q.log_upper(512,1+(1-q.P)*q.X)
             +q.log_upper(-512,q.X)+q.log_upper(-512*q.THETA,q.ZERO))
    rows=[band.row_constants(r,45*q.N//100) for r in (1850,2100)]
    # For n>=2048 and budget<=32, binom(n,j)*(eta/beta)^j increases
    # through the retained budget. Use n<=4095 in this count only.
    count=q.log_upper(1,budget+1)+q.log_upper(1,math.comb(4095,budget))
    upper=None;where=None
    for start in range(1,2049):
        n=q.L-start;regions=(n+16)//16
        fixed=max(base-regions*pi for base,pi in rows)
        value=fixed+initial+(n-budget)*beta+budget*eta+count
        if upper is None or value>upper:upper,where=value,start
    upper+=q.log_upper(1,251)+11  # All row counts and activation times.
    return upper,where


def bernstein_power(values,power,precision=512):
    """Normalized coefficients of (sum binom(g,a) values[a] h^a)^power.

    All values lie in [0,1]. Hypergeometric weights form a convex combination,
    so normalized coefficients stay <=1. Arbitrary-precision integers and
    exact binomial denominators avoid underflow at the extreme degrees.
    """
    assert all(0<=v<=1 for v in values)
    g=len(values)-1;scale=1<<precision
    coeff=[(v.numerator*scale+v.denominator-1)//v.denominator for v in values]
    previous=[scale];old_choose=[1]
    for step in range(power):
        degree=g*(step+1)
        new_choose=[1]
        for u in range(degree):new_choose.append(new_choose[-1]*(degree-u)//(u+1))
        current=[]
        for u in range(degree+1):
            numerator=sum(math.comb(g,a)*old_choose[u-a]*coeff[a]*previous[u-a]
                          for a in range(g+1) if 0<=u-a<len(previous))
            denominator=new_choose[u]*scale
            value=(numerator+denominator-1)//denominator
            assert 0<=value<=scale
            current.append(value)
        previous,old_choose=current,new_choose
    return previous


def short_resets(beta):
    """Unlimited short zero runs, with exact reset and reactivation prices."""
    reset=[F(1)]
    for a in range(1,5):
        # For after=empty, the actual local compatibility condition is |T|<a.
        reset.append(sum(math.comb(4,t)*q.P**t*(1-q.P)**(4-t)
                         *q.root_upper(F(1,1<<(32*t))) for t in range(a)))
    coefficients=bernstein_power(reset,128)
    b=max(band.cheap_log_upper(coefficients[w],512)-32*q.THETA*w
          +q.log_upper(-w,q.X)+q.log_upper(-2*q.THETA*w,q.ZERO)
          for w in range(216,233))
    activation=sum(math.comb(512,w)*((1-q.P)*q.X)**w for w in range(216,233))
    a=q.log_upper(1,activation)
    ceil=lambda v:-((-v.numerator)//v.denominator)
    ae,be,ge=map(ceil,(a,b,beta));le=-270
    assert ge<le<0
    # Scalar matrix entries are these explicit powers of two, not estimates.
    allowed=-1
    for loops in range(9):
        test=F(2)**(ge-le)+sum(F(2)**(ae+be-(k+2)*le) for k in range(loops+1))
        if test>1:break
        allowed=loops
    assert allowed>=0
    assert allowed==2  # Three consecutive zero states, two zero-to-zero edges.
    # Every zero-age potential is >=2^(ae-le). The minimum product potential
    # is therefore still f(232), not an unjustified terminal reset condition.
    assert ae-le>=232*(q.log_interval(q.X)[1]+q.THETA*q.log_interval(q.ZERO)[1])
    initial=ae+q.log_upper(-232,q.X)+q.log_upper(-232*q.THETA,q.ZERO)
    rows=[band.row_constants(r,45*q.N//100) for r in (1850,2100)]
    upper=None;where=None
    for start in range(1,2049):
        n=q.L-start;regions=(n+16)//16
        fixed=max(base-regions*pi for base,pi in rows)
        value=fixed+initial+n*le
        if upper is None or value>upper:upper,where=value,start
    upper+=q.log_upper(1,251)+11
    assert upper==FROZEN_RESET_UPPER and upper < -5978
    return dict(reset_log2_upper=str(b),reset_float=float(b),activation_log2_upper=str(a),
                activation_float=float(a),dyadic_exponents=(ge,ae,be,le),
                maximum_zero_run=allowed+1,unlimited_episodes=True,
                log2_upper=str(upper),upper_float=float(upper),worst_activation=where,
                whole_code=False,directed=True)


def tests():
    # Noncommuting matrices: exact sum of every path with <=h bad edges.
    # The full envelope may also contain good edges; this only overcounts.
    good=((F(1,4),F(1,8)),(F(0),F(1,6)))
    bad=((F(1,9),F(1,7)),(F(1,5),F(1,8)))
    f=(F(2),F(3));full=tuple(tuple(good[i][j]+bad[i][j] for j in range(2)) for i in range(2))
    beta=max(sum(good[i][j]*f[j] for j in range(2))/f[i] for i in range(2))
    eta=max(sum(full[i][j]*f[j] for j in range(2))/f[i] for i in range(2))
    checks=0
    for n in range(1,6):
        for h in range(n+1):
            upper=sum(math.comb(n,j)*beta**(n-j)*eta**j for j in range(h+1))
            for initial in range(2):
                exact=F(0)
                for path in itertools.product(range(2),repeat=n):
                    for types in itertools.product(range(2),repeat=n):
                        if sum(types)>h:continue
                        value=F(1);i=initial
                        for j,t in zip(path,types):
                            value*=(bad if t else good)[i][j];i=j
                        exact+=value
                assert exact<=upper*f[initial]/min(f);checks+=1
    # Degree-four local envelope, including zero occupancy and full resets.
    for occupancies in itertools.product(range(5),repeat=3):
        local=[F(5-a,7+a) for a in range(5)]
        assert math.prod(local[a] for a in occupancies)<=max(local)**3
    # Independent exact polynomial powers validate Bernstein normalization.
    from collections import defaultdict
    bchecks=0
    for values in ((F(1),F(1,3),F(1,7)),(F(0),F(1,2),F(1)),(F(1),F(1),F(1))):
        base=[math.comb(2,a)*v for a,v in enumerate(values)];poly=[F(1)]
        for power in range(1,5):
            out=defaultdict(F)
            for i,x in enumerate(poly):
                for j,y in enumerate(base):out[i+j]+=x*y
            poly=[out[i] for i in range(2*power+1)]
            upper=bernstein_power(values,power,precision=32)
            for u,value in enumerate(poly):
                assert F(upper[u],1<<32)>=value/math.comb(2*power,u);bchecks+=1
    # Every local reset mask is equivalent to the retained target-count test.
    for a in range(1,5):
        for target in range(16):
            assert bool(any(q.choices(4,(1<<a)-1,0,target,True)))==(target.bit_count()<a)
    # Short-reset automaton: exact matrix powers include every return and
    # termination at ANY zero age, not only trajectories ending in the band.
    lam=F(1,2);a=F(1,8);b=F(1,16);g=F(1,8);h=2
    psi=[F(1)]+[a*sum(lam**(-k) for k in range(1,h-i+2)) for i in range(h+1)]
    matrix=[[F(0)]*4 for _ in range(4)]
    matrix[0][0]=g;matrix[0][1]=b
    for i in range(1,4):
        matrix[i][0]=a
        if i<3:matrix[i][i+1]=1
    assert all(sum(matrix[i][j]*psi[j] for j in range(4))<=lam*psi[i] for i in range(4))
    vector=[F(1)]*4
    for length in range(1,9):
        vector=[sum(matrix[i][j]*vector[j] for j in range(4)) for i in range(4)]
        assert all(vector[i]<=lam**length*psi[i]/min(psi) for i in range(4))
    print(dict(noncommuting_path_checks=checks,local_occupancy_envelopes=125,
               bernstein_cells=bchecks,reset_masks=64,zero_age_row_checks=4,
               zero_age_terminal_checks=32),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(60,lambda:os._exit(124));timer.daemon=True;timer.start()
    start=time.monotonic()
    try:
        tests();band.configure()
        q.P,q.H,q.X,q.Y,q.ZERO=map(F,band.BRIDGE_MARKERS);q.ROOTS.clear()
        probabilities,total,_=q.local()
        box=band.coefficient_box(probabilities,(232,56,40))
        old=band.band_bound(box,total,(216,232),(24,40),(1850,2100),45,
                            quota=False,starts=(1,2048))
        assert F(old['upper'])==band.BRIDGE_UPPER
        beta,where=good_step(box,total);eta,group,local=full_step()
        print(dict(log2_beta=str(beta),beta_float=float(beta),worst_good=where,
                   log2_eta=str(eta),eta_float=float(eta),worst_incoming_group=group,
                   local_log2_uppers=list(map(float,local))),flush=True)
        last=None
        for budget in range(33):
            upper,activation=certificate(beta,eta,budget)
            if budget in (0,8,16,24,32) or upper>=-40:
                print(dict(exception_budget=budget,log2_upper=str(upper),upper_float=float(upper),
                           worst_activation=activation,distance_percent=45,
                           whole_code=False,directed=True),flush=True)
            if upper>=-40:break
            last=(budget,upper)
        assert last==(25,FROZEN_EXCEPTION_UPPER) and last[1] < -88
        print(dict(last_40bit_budget=last[0] if last else None,
                   last_log2_upper=str(last[1]) if last else None,
                   last_upper_float=float(last[1]) if last else None,
                   seconds=time.monotonic()-start,files_written=0),flush=True)
        print(short_resets(beta),flush=True)
        print(dict(seconds=time.monotonic()-start,files_written=0),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
