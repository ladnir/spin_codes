"""Exact geometry tail, not a whole-code distance certificate.

Uniform over all inactive sets of one requested size. Serial,30s watchdog,
small scalar integers/rationals (<4 MiB numeric storage), no files.
"""
import argparse
from fractions import Fraction as F
from itertools import product
import math
import os
from pathlib import Path
import shutil
import threading

from quartet_two_update import log_upper,log_interval


def factors(p,v):
    a=1-p;g=p+a*v
    return 1-a**4+a**4*g**4,g


def tests():
    # Three width-four updates, fixed quartets, and independent labels.
    # v=(3/2)^2 makes both sides rational before squaring Cauchy--Schwarz.
    p=F(1,3);v=F(9,4);total=F(0);conditioned=F(0);mass=F(0)
    for bits in product((0,1),repeat=12):
        inactive=sum(bits);prob=p**inactive*(1-p)**(12-inactive)
        active=[4-sum(bits[4*t:4*t+4]) for t in range(3)]
        z=active[2]+sum(active[t] for t in (0,1) if active[t+1]==4)
        total+=prob*F(3,2)**z
        if inactive==4:conditioned+=prob*F(3,2)**z;mass+=prob
    local,last=factors(p,v)
    assert total**2<=local**2*last**4
    assert mass==math.comb(12,4)*p**4*(1-p)**8
    assert (conditioned/mass)**2<=local**2*last**4/mass**2
    checks=0
    for z in range(9):
        exponents=[a-max(0,(a-z+4)//5) for a in range(32)]
        assert all(x<=y for x,y in zip(exponents,exponents[1:]));checks+=31
    concavity=0
    for m in range(2,33):
        for r in range(1,m):
            assert math.comb(m,r)**2>=math.comb(m,r-1)*math.comb(m,r+1)
            concavity+=1
    tangent=0
    for size in range(2,13):
        for anchor in range(1,size):
            x=F(anchor,size)
            for k in range(size+1):
                assert math.comb(size,k)*x**k*(1-x)**(size-k)<=1;tangent+=1
    print(dict(exact_label_patterns=4096,parity_moment_checks=2,
               coverage_exponent_monotonicity_checks=checks,
               exact_binomial_log_concavity_checks=concavity,
               exact_entropy_tangent_checks=tangent),flush=True)


def geometry_exponent(r,p,v):
    choose=math.comb(8192,r);n=8192-r;local,last=factors(p,v)
    return (log_upper(2-512,choose)-512*r*log_interval(p)[0]
            -512*n*log_interval(1-p)[0]+log_upper(128*4095,local)
            +log_upper(512,last))


def certificate(r,p,v,bits=60):
    assert 0<r<8192 and 0<p<1 and v>1
    # Square the complete row-set union bound before directed logarithms.
    exponent=geometry_exponent(r,p,v)
    lower=log_interval(v)[0];ratio=(exponent+2*bits)/lower
    cutoff=-((-ratio.numerator)//ratio.denominator)
    upper=exponent-cutoff*lower
    assert upper<=-2*bits
    print(dict(inactive_rows=r,p=str(p),v=str(v),Z_at_least=cutoff,
               all_row_sets_failure_at_most=f'2^-{bits}',squared_log2_upper=str(upper),
               consequence=f'A <= {cutoff-1} + 5*C_T outside that event',
               arbitrary_histories=True,free_terminal=True,
               distance_certificate=False,arithmetic='directed rational',files_written=0),flush=True)
    return cutoff


def activity_cutoff(r,z,bits=60):
    """Uniform low-activity class outside the separately charged Z event."""
    N=1<<21;n=8192-r;row_term=log_upper(1,math.comb(8192,r))
    def upper(a):
        inactive_floor=max(0,(a-z+4)//5)
        entropy=log_upper(N,N)+log_upper(-a,a)+log_upper(-(N-a),N-a)
        return 1+row_term+log_upper(1,a)+entropy+32*(a-128*n-inactive_floor)
    lo,hi=1,N//2
    assert upper(lo)<=-bits
    while lo<hi:
        mid=(lo+hi+1)//2
        if upper(mid)<=-bits:lo=mid
        else:hi=mid-1
    bound=upper(lo);assert bound<=-bits
    print(dict(inactive_rows=r,activity_at_most=lo,good_geometry_field_log2_upper=str(bound),
               combined_failure_at_most=f'2^-{bits-1}',
               includes_geometry_and_field_unions=True,all_row_sets_of_this_size=True,
               all_row_counts=False,arbitrary_histories=True,free_terminal=True,
               whole_code_distance=False,arithmetic='directed rational',files_written=0),flush=True)
    return lo


def interval_proposal(left,right,bits):
    """Bounded FLOAT proposal, never used as a certificate."""
    endpoints=[]
    for r in (left,right):
        n=8192-r
        choose=(math.lgamma(8193)-math.lgamma(r+1)-math.lgamma(n+1))/math.log(2)
        endpoints.append((r,n,choose))
    center=F(left+right,2*8192);best=None
    for shift in range(0,61,5):
        p=center+F(shift,1000)
        if not 0<p<1:continue
        for extra in range(2,25,2):
            v=F(64+extra,64);local,last=factors(p,v)
            scores=[((2-512)*choose-512*r*math.log2(float(p))-512*n*math.log2(float(1-p))
                     +128*4095*math.log2(float(local))+512*math.log2(float(last))+2*bits)
                    /math.log2(float(v)) for r,n,choose in endpoints]
            score=max(scores)
            if best is None or score<best[0]:best=(score,p,v)
    return best[1:]


def interval_certificate(left,right,p,v,bits=74,verbose=True):
    """Exact endpoint geometry and tangent field bounds for every r in range."""
    assert 1<=left<=right<8192
    exponent=max(geometry_exponent(r,p,v) for r in (left,right))
    lower=log_interval(v)[0];ratio=(exponent+2*bits)/lower
    cutoff=-((-ratio.numerator)//ratio.denominator);z=cutoff-1
    geometry_bound=exponent-cutoff*lower
    assert geometry_bound<=-2*bits
    N=1<<21;y=F(left+right,2*8192);mid_n=8192-F(left+right,2)
    def field_upper(b):
        x=b*mid_n/N;max_a=b*(8192-left)
        assert 0<x<1 and max_a<=N//2
        points=[F(left),F(right)]
        crossover=8192-F(z,b)
        if left<crossover<right:points.append(crossover)
        bounds=[]
        for r in points:
            n=8192-r;a=b*n
            row=log_upper(-r,y)+log_upper(-n,1-y)
            entropy=log_upper(-a,x)+log_upper(-(N-a),1-x)
            e=min(a-128*n,F(4,5)*a-128*n+F(z,5))
            bounds.append(1+log_upper(1,max_a)+row+entropy+32*e)
        return max(bounds)
    lo,hi=1,min(159,N//2//(8192-left))
    assert field_upper(lo)<=-bits
    # Search proposes b; the final exact check, not optimality, is used.
    while lo<hi:
        mid=(lo+hi+1)//2
        if field_upper(mid)<=-bits:lo=mid
        else:hi=mid-1
    bound=field_upper(lo);assert bound<=-bits
    record=dict(inactive_interval=[left,right],p=str(p),v=str(v),Z_upper=z,
                activity_per_active_row=lo,per_row_count_geometry_bits=bits,
                per_row_count_field_bits=bits,field_log2_upper=str(bound),
                geometry_squared_log2_upper=str(geometry_bound),
                includes_every_intermediate_row_count=True,distance_certificate=False)
    if verbose:print(record,flush=True)
    return record


def interval_family(left,right,width,search=False):
    assert 1<=left<=right<8192 and 1<=width<=128
    records=[]
    for start in range(left,right+1,width):
        end=min(right,start+width-1)
        center=F(start+end,2*8192)
        p,v=(interval_proposal(start,end,74) if search else
             (min(center+F(3,100),(1+center)/2),F(39,32)))
        verbose=(right-left<512 or not records or start<=2100<=end or end==right)
        records.append(interval_certificate(start,end,p,v,verbose=verbose))
        if not verbose and len(records)%16==0:
            print(dict(phase='directed row intervals',completed=len(records),through=end),flush=True)
    union_terms=2*(right-left+1)
    loss=(union_terms-1).bit_length()
    from quartet_low_activity import remaining_activity_floor
    improved=sum(record['activity_per_active_row']*(8192-r)>=remaining_activity_floor(8192-r)
                 for record in records for r in range(record['inactive_interval'][0],record['inactive_interval'][1]+1))
    print(dict(inactive_rows=[left,right],intervals=len(records),
               combined_failure_at_most=f'2^-{74-loss}',union_terms=union_terms,
               row_counts_with_stronger_activity_cutoff=improved,
               frozen_parameter_rule=not search,
               all_row_sets=True,all_histories=True,free_terminal=True,
               all_row_counts=False,whole_code_distance=False,
               final_arithmetic='directed rational',files_written=0),flush=True)
    return records


def certified_row_bounds(r):
    """One slice of the replayed 105..8191 family; shared allowance 2^-60."""
    assert 105<=r<=8191
    left=105+64*((r-105)//64);right=min(8191,left+63)
    center=F(left+right,2*8192)
    p=min(center+F(3,100),(1+center)/2)
    return interval_certificate(left,right,p,F(39,32),verbose=False)


def screen(r):
    # Proposal only. The selected rational parameters must pass certificate().
    n=8192-r;choose=(math.lgamma(8193)-math.lgamma(r+1)-math.lgamma(n+1))/math.log(2)
    best=None
    for numerator in range(250,311,5):
        p=F(numerator,1000)
        for extra in range(2,25,2):
            v=F(64+extra,64);local,last=factors(p,v)
            score=((2-512)*choose-512*r*math.log2(float(p))-512*n*math.log2(float(1-p))
                   +128*4095*math.log2(float(local))+512*math.log2(float(last))+120)/math.log2(float(v))
            if best is None or score<best[0]:best=(score,p,v)
    print(dict(proposed_float_cutoff=best[0],p=str(best[1]),v=str(best[2]),not_a_certificate=True),flush=True)
    return best[1:]


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(30,lambda:os._exit(124));timer.daemon=True;timer.start()
    try:
        parser=argparse.ArgumentParser()
        parser.add_argument('--tests-only',action='store_true')
        parser.add_argument('--screen',action='store_true')
        parser.add_argument('--rows',type=int,default=2100)
        parser.add_argument('--p',default='29/100')
        parser.add_argument('--v',default='39/32')
        parser.add_argument('--intervals',nargs=2,type=int)
        parser.add_argument('--width',type=int,default=32)
        args=parser.parse_args();tests()
        if args.tests_only:return
        if args.intervals:
            interval_family(*args.intervals,args.width,search=args.screen)
            return
        p,v=screen(args.rows) if args.screen else (F(args.p),F(args.v))
        cutoff=certificate(args.rows,p,v)
        activity_cutoff(args.rows,cutoff-1)
    finally:timer.cancel()


if __name__=='__main__':main()
