"""A directed regional transfer allowing moving seven-type populations.

Restricts populations to a stated finite set, but sums every path within
that set and every original-inactive allocation with the regional quota.
Not full coverage. Serial Peach180s, <=256MiB, no generated files.
"""
from collections import defaultdict
from fractions import Fraction as F
import itertools
import math
import os
from pathlib import Path
import shutil
import threading
import time
import numpy as np
from quartet_closed_joint import bound
from quartet_closed_quotient import coefficient_table
from quartet_closed_directed import atoms,tilted_atoms,coefficient_upper,power_upper
from quartet_closed_population import multinomial
from quartet_two_update import log_upper


def populations(center,radius):
    z,v0,v1,b,l0,l1,l2=center;result=[]
    for a,c in itertools.product(range(-radius,radius+1),repeat=2):
        p=(z,v0-a-c,v1,b,l0+a,l1+c,l2)
        assert min(p)>=0
        result.append(p)
    return result


def matrix_coefficient(integer_edges,count,target,steps,precision=31):
    """Outward positive recurrence for a substochastic polynomial matrix."""
    scale=1<<precision;total=np.zeros(count,dtype=object)
    for left,right,poly in integer_edges:total[left]+=sum(int(x) for x in poly)
    rowmax=int(max(total));maximum=scale
    current=np.zeros((target+1,count),dtype=np.uint64);current[0,:]=scale
    for _ in range(steps):
        assert maximum*rowmax+scale-1<1<<64
        after=np.zeros_like(current)
        for left,right,poly in integer_edges:
            for degree,coefficient in enumerate(poly[:target+1]):
                if coefficient:after[degree:,left]+=coefficient*current[:target+1-degree,right]
        after+=np.uint64(scale-1);after>>=np.uint64(precision)
        maximum=(maximum*rowmax+scale-1)//scale
        assert int(after.max())<=maximum
        current=after
    return [F(int(x),scale) for x in current[target]]


def regional(center,radius=1,balanced=False,age_modulus=None):
    start=time.monotonic();pops=populations(center,radius);size=len(pops)
    fit=bound(coefficient_table(),center,center,131,return_fit=True)
    weighted,total,shift,markers=tilted_atoms(atoms(),fit)
    probabilities=defaultdict(F)
    for (old,new,j),value in weighted:
        degrees=(j-new[0],old[4],new[4]) if age_modulus is None else (j-new[0],old[4],old[5],new[4])
        probabilities[degrees]+=value/total
    max_young=max(p[4] for p in pops)
    target=(512-center[0],max_young,max_young) if age_modulus is None else (512-center[0],max_young,age_modulus-1,max_young)
    grid,scale=coefficient_upper(probabilities,target,128,return_array=True,modular_axes=() if age_modulus is None else (2,))
    def population_power(old,new):
        degrees=[old[i] for i in fit['oldkeep']]+[new[i] for i in fit['newkeep']]
        return sum(d*x for d,x in zip(degrees,markers))
    center_power=population_power(center,center);normal=multinomial(center)
    exact_edges=[]
    for left,old in enumerate(pops):
        for right,new in enumerate(pops):
            # Exact seven-type conservation, not a fitted constraint.
            if new[5]!=old[4]:continue
            deaths=old[5]+old[6]-new[6]
            if not 0<=deaths<=min(old[6],new[3]):continue
            ratio=power_upper(center_power-population_power(old,new))*F(normal,multinomial(old))
            values=grid[:,old[4],new[4]] if age_modulus is None else grid[:,old[4],old[5]%age_modulus,new[4]]
            exact_edges.append((left,right,ratio,values))
    potential=[F(1)]*size
    if balanced:
        aggregate=np.zeros((size,size))
        for left,right,ratio,values in exact_edges:
            aggregate[left,right]=float(ratio*F(sum(int(v) for v in values),scale))
        proposed=np.ones(size)
        for _ in range(256):
            proposed=aggregate@proposed
            proposed/=max(proposed)
        assert min(proposed)>0
        potential=[F(max(1,int(round(float(v)*(1<<40)))),1<<40) for v in proposed]
    # Any positive potential is valid; convergence of the proposal is not
    # assumed. This potential must be carried through regional boundaries.
    exact_edges=[(left,right,ratio*potential[right]/potential[left],values)
                 for left,right,ratio,values in exact_edges]
    row_sums=[F(0)]*size
    for left,right,ratio,values in exact_edges:
        row_sums[left]+=ratio*F(sum(int(v) for v in values),scale)
    normalizer=max(row_sums);assert normalizer>0
    integer_edges=[]
    for left,right,ratio,values in exact_edges:
        factor=ratio/normalizer
        # values already use the same fixed-point scale. This is the
        # exact ceiling of scale*(ratio/normalizer)*(values/scale).
        coefficients=np.array([(int(v)*factor.numerator+factor.denominator-1)//factor.denominator for v in values],dtype=np.uint64)
        integer_edges.append((left,right,coefficients))
    mass=matrix_coefficient(integer_edges,size,2100-16*center[0],16)
    scalar=F(128*shift)-F(center_power,1024)+log_upper(128,total)+log_upper(-1,normal)
    common=16*scalar-F(2100*markers[-1],1024)+log_upper(16,normalizer)+log_upper(-1,math.comb(8192,2100))
    results=[common+log_upper(1,value) for value in mass]
    worst=max(range(size),key=lambda i:results[i])
    print(dict(center=center,radius=radius,populations=size,edges=len(integer_edges),balanced=balanced,age_modulus=age_modulus,
               all_internal_paths=True,all_original_quota_allocations=True,
               worst_start=pops[worst],directed_region_log2_upper=str(results[worst]),
               display=float(results[worst]),potential_range_log2=float(log_upper(1,max(potential)/min(potential))),
               minimum_scaled_quota_mass=float(min(mass)),seconds=time.monotonic()-start,
               full_certificate=False),flush=True)
    return results


def tests():
    precision=20;scale=1<<precision
    edges=[(0,0,np.array([scale//4,scale//4],dtype=np.uint64)),
           (0,1,np.array([scale//8,scale//8],dtype=np.uint64)),
           (1,0,np.array([scale//2],dtype=np.uint64)),
           (1,1,np.array([0,scale//4],dtype=np.uint64))]
    exact={(0,i):F(1) for i in range(2)};checked=0
    for steps in range(1,5):
        after=defaultdict(F)
        for left,right,poly in edges:
            for degree,weight in enumerate(poly):
                for (old_degree,state),value in exact.items():
                    if state==right:after[old_degree+degree,left]+=F(int(weight),scale)*value
        exact=dict(after)
        for degree in range(steps+1):
            upper=matrix_coefficient(edges,2,degree,steps,precision)
            assert all(exact.get((degree,i),F(0))<=v for i,v in enumerate(upper));checked+=2
    print(dict(exact_polynomial_matrix_comparisons=checked),flush=True)
    probs={(0,0):F(1,7),(1,0):F(2,7),(0,1):F(3,7),(2,4):F(1,7)}
    exact={(0,0):F(1)};checked=0
    for steps in range(1,5):
        after=defaultdict(F)
        for old,a in exact.items():
            for new,b in probs.items():after[old[0]+new[0],(old[1]+new[1])%4]+=a*b
        exact=dict(after)
        values,scale=coefficient_upper(probs,(2*steps,3),steps,return_array=True,modular_axes=(1,))
        for key,value in exact.items():
            assert value<=F(int(values[key]),scale);checked+=1
    print(dict(exact_modular_coefficient_comparisons=checked),flush=True)


def main():
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',choices=('baseline','first','second','wide','score'),default='score')
    args=parser.parse_args()
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(180,lambda:os._exit(124));timer.daemon=True;timer.start()
    try:
        tests()
        first=(104,156,60,64,32,32,64);second=(104,176,72,32,16,16,96)
        cases={'first':(first,1,4),'second':(second,1,4),'wide':(second,2,8),
               'score':((104,120,64,32,16,16,160),2,None)}
        if args.case=='baseline':
            for center in (first,second):regional(center,1,True)
        else:
            center,radius,modulus=cases[args.case]
            regional(center,radius,True,modulus)
    finally:timer.cancel()


if __name__=='__main__':main()
