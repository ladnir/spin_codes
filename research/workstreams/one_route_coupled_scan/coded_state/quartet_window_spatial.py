"""Exact quartet occupation formula, FLOAT short-window moment diagnostic.

No temporal factorization is asserted. Serial,30s watchdog,<32MiB,no files.
"""
from fractions import Fraction as F
import itertools
import math
import os
from pathlib import Path
import shutil
import threading
import time


def occupied_moments(width,p):
    assert width%4==0
    a=p.numerator**4;d=p.denominator**4
    local=[d]+[math.comb(4,j)*a for j in range(1,5)]
    polynomial=[1]
    for _ in range(width//4):
        after=[0]*(len(polynomial)+4)
        for j,x in enumerate(polynomial):
            for k,y in enumerate(local):after[j+k]+=x*y
        polynomial=after
    denominator=d**(width//4)
    return [F(value,denominator*math.comb(width,k)) for k,value in enumerate(polynomial)]


def exact_tests():
    p=F(2,3);moments=occupied_moments(8,p);checks=0
    partitions=[]
    for rest in itertools.combinations(range(1,8),3):
        group={0,*rest};partitions.append((group,set(range(8))-group))
    for k in range(9):
        selected=set(range(k))
        direct=sum(p**(4*sum(bool(part&selected) for part in parts)) for parts in partitions)/len(partitions)
        assert direct==moments[k];checks+=len(partitions)
    # Independent check of the marked-subset expansion on a whole quartet.
    width=4;p=F(2,3);alpha_factor=F(3,2);full_factor=F(2)
    moments=occupied_moments(width,p)
    formula=F(0)
    for k in range(width+1):
        inner=sum(math.comb(k,h)*(alpha_factor-1)**(k-h)*(full_factor-alpha_factor)**h*moments[h] for h in range(k+1))
        formula+=math.comb(width,k)*p**k*moments[k]**3*inner
    # First input has a 4-bit active mask. Next four quartet-full bits are
    # independent with probability p^4; other inputs do not affect the event.
    direct=F(0);paths=0
    for initial in range(16):
        count=initial.bit_count();probability=p**count*(1-p)**(4-count)
        for full in itertools.product((0,1),repeat=4):
            prob=probability*math.prod(p**4 if bit else 1-p**4 for bit in full)
            factor=full_factor if all(full) else alpha_factor if all(full[:3]) else F(1)
            direct+=prob*factor**count;paths+=1
    assert direct==formula
    print(dict(exact_partition_subset_checks=checks,exact_spatial_category_paths=paths),flush=True)


def logadd(values):
    largest=max(values)
    return largest+math.log2(sum(2**(x-largest) for x in values))


def screen(p,slope,tau,logb):
    width=512;alpha=3*slope-tau;assert 0<=alpha<=slope
    G=occupied_moments(width,p)
    # Only this final evaluation is floating point. G itself is exact.
    logs=[math.log2(x.numerator)-math.log2(x.denominator) for x in G]
    a=math.expm1(float(logb*alpha)*math.log(2))
    c=2**float(logb*slope)-2**float(logb*alpha)
    assert a>0 and c>0
    outer=[]
    for k in range(width+1):
        inner=logadd([math.log2(math.comb(k,h))+(k-h)*math.log2(a)+h*math.log2(c)+logs[h] for h in range(k+1)])
        outer.append(math.log2(math.comb(width,k))+k*math.log2(float(p))+3*logs[k]+inner)
    actual=logadd(outer)
    independent=width*math.log2(1+a*float(p)**13+c*float(p)**17)
    four_a=math.expm1(float(4*logb*alpha)*math.log(2))
    four_c=2**float(4*logb*slope)-2**float(4*logb*alpha)
    holder=width/4*math.log2(1+four_a*float(p)**13+four_c*float(p)**17)
    print(dict(active_probability=str(p),slope=str(slope),death_weight=str(tau),log2_base=str(logb),
               FLOAT_log2_one_start_actual_moment=actual,FLOAT_independent_column_moment=independent,
               FLOAT_read4_upper=holder,FLOAT_log2_coupling_excess=actual-independent,
               exact_partition_formula=True,terminal_terms_not_included=True,
               overlapping_start_times_not_independent=True,not_a_directed_certificate=True,
               whole_code_distance=False),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(30,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:
        exact_tests()
        screen(1-F(2100,8192),F(1,4),F(37,50),F(1))
        screen(1-F(2100,8192)*F(6,7),F(707,2560),F(381,512),F(3))
        print(dict(seconds=time.monotonic()-start,files_written=0),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
