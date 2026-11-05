"""Directed seven-type conditional coefficients at b=2^(5/2).

Floating optimization selects arbitrary dyadic log markers only. All
trusted coefficients, outward DP, and logarithm bounds use integers and
rationals. Serial Peach180s, <=256MiB, no generated files. No full cert.
"""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
import itertools
import math
import os
from pathlib import Path
import shutil
import threading
import time
import numpy as np
from quartet_closed_joint import bound
from quartet_closed_quotient import kernel,coefficient_table
from quartet_closed_population import multinomial
from quartet_two_update import log_upper


@lru_cache(maxsize=1)
def dyadic_root():
    """Upper endpoint for 2^(1/1024) on a 2^-80 grid, using integer roots."""
    target=1<<(80*1024+1);value=target
    for _ in range(10):value=math.isqrt(value-1)+1
    assert (value-1)**1024<target<=value**1024
    return value


@lru_cache(maxsize=1024)
def fractional_power(exponent):
    assert 0<=exponent<1024
    scale=1<<80;result=scale;base=dyadic_root();n=exponent
    while n:
        if n&1:result=(result*base+scale-1)//scale
        n>>=1
        if n:base=(base*base+scale-1)//scale
    # Round once more to a compact dyadic rational, always upward.
    return F((result+(1<<32)-1)>>32,1<<48)


def power_upper(exponent):
    whole,fraction=divmod(exponent,1024)
    return fractional_power(fraction)*F(2)**whole


def atoms():
    """Each coefficient is an integer times2^(exponent/1024)."""
    out={}
    for hist in itertools.combinations_with_replacement(range(7),4):
        old=tuple(hist.count(t) for t in range(7))
        for (dest,u,v,a,d),ways in kernel(hist):
            new=tuple(dest.count(t) for t in range(7))
            for paid in range(v+1):
                key=old,new,u+paid
                exponent=707*a-1905*d-3072*paid
                weight=multinomial(old)*ways*math.comb(v,paid)
                if key in out:
                    previous,previous_exp=out[key];assert previous_exp==exponent
                    weight+=previous
                out[key]=weight,exponent
    return out


def coefficient_upper(probabilities,targets,groups,precision=31,return_axis=None,return_array=False,modular_axes=()):
    """Positive outward fixed-point recurrence; explicit no-overflow bound."""
    scale=1<<precision
    integers={k:(p.numerator*scale+p.denominator-1)//p.denominator for k,p in probabilities.items()}
    shape=tuple(t+1 for t in targets);assert math.prod(shape)*24<256*(1<<20)
    current=np.zeros(shape,dtype=np.uint64);current[(0,)*len(shape)]=scale
    maximum=scale;total=sum(integers.values())
    terms=[]
    for degrees,value in integers.items():
        if any(d>t and i not in modular_axes for i,(d,t) in enumerate(zip(degrees,targets))):continue
        choices=[]
        for i,(d,t) in enumerate(zip(degrees,targets)):
            if i not in modular_axes:choices.append(((slice(0,t+1-d),slice(d,t+1)),));continue
            d%=t+1
            choices.append(((slice(None),slice(None)),) if d==0 else
                           ((slice(0,t+1-d),slice(d,t+1)),(slice(t+1-d,t+1),slice(0,d))))
        for pairs in itertools.product(*choices):
            source,dest=zip(*pairs)
            terms.append((source,dest,np.uint64(value)))
    for _ in range(groups):
        assert maximum*total+scale-1<1<<64
        after=np.zeros_like(current)
        for source,dest,value in terms:after[dest]+=value*current[source]
        after+=np.uint64(scale-1);after>>=np.uint64(precision)
        maximum=(maximum*total+scale-1)//scale
        assert int(after.max())<=maximum
        current=after
    if return_array:return current,scale
    if return_axis is not None:
        selection=list(targets);selection[return_axis]=slice(None)
        return [F(int(x),scale) for x in current[tuple(selection)]]
    return F(int(current[targets]),scale)


def tilted_atoms(local_atoms,fit):
    """Exact upper local weights at arbitrary rounded proposal markers."""
    assert np.all(np.isfinite(fit['tilts']))
    markers=[int(round(1024*float(x))) for x in fit['tilts']]
    weighted=[]
    for key,powers in zip(fit['keys'],fit['powers']):
        integer,field_power=local_atoms[key]
        exponent=field_power+sum(int(p)*x for p,x in zip(powers,markers))
        weighted.append((key,integer,exponent))
    shift=max(e for k,c,e in weighted)//1024
    weighted=[(key,integer*power_upper(exponent-1024*shift)) for key,integer,exponent in weighted]
    return weighted,sum(value for key,value in weighted),shift,markers


def directed(local_atoms,fit,incoming,outgoing,inactive,coordinates,polynomial=False):
    def degree(key,coordinate):
        old,new,j=key;kind,index=coordinate
        if kind=='old':return old[index]
        if kind=='new':return new[index]
        if kind=='paid':return j-new[0]
        assert kind=='inactive';return j
    weighted,total,shift,markers=tilted_atoms(local_atoms,fit)
    target_values=[incoming[i] for i in fit['oldkeep']]+[outgoing[i] for i in fit['newkeep']]
    if fit['interior']:target_values.append(inactive)
    target_power=sum(t*x for t,x in zip(target_values,markers))
    terms=defaultdict(F)
    for key,value in weighted:
        degrees=tuple(degree(key,c) for c in coordinates)
        assert min(degrees)>=0
        terms[degrees]+=value
    probabilities={k:v/total for k,v in terms.items()}
    targets=tuple(degree((incoming,outgoing,inactive),c) for c in coordinates)
    if polynomial:
        assert coordinates[0]==('paid',0) and fit['interior']
        targets=(sum(outgoing)-outgoing[0],)+targets[1:]
        masses=coefficient_upper(probabilities,targets,fit['groups'],return_axis=0)
        # The original-inactive marker is removed only after the exact
        # regional quota is selected, not separately at each mean count.
        population_power=target_power-inactive*markers[-1]
        scalar=F(fit['groups']*shift)-F(population_power,1024)
        scalar+=log_upper(fit['groups'],total)+log_upper(-1,multinomial(incoming))
        return scalar,F(markers[-1],1024),masses
    mass=coefficient_upper(probabilities,targets,fit['groups'])
    result=F(fit['groups']*shift)-F(target_power,1024)
    result+=log_upper(fit['groups'],total)+log_upper(1,mass)+log_upper(-1,multinomial(incoming))
    return result,mass,len(terms)


def tests():
    for f in (0,1,2,255,511,512,767,1023):
        upper=fractional_power(f);assert upper**1024>=F(2)**f
    probs={(0,0):F(1,7),(1,0):F(2,7),(0,1):F(3,7),(1,1):F(1,7)}
    exact={(0,0):F(1)};checked=0
    for groups in range(1,5):
        after=defaultdict(F)
        for left,a in exact.items():
            for right,b in probs.items():after[tuple(x+y for x,y in zip(left,right))]+=a*b
        exact=dict(after)
        for target,value in exact.items():
            upper=coefficient_upper(probs,target,groups)
            assert value<=upper and upper-value<F(1,1<<24);checked+=1
        vector=coefficient_upper(probs,(groups,groups//2),groups,return_axis=0)
        assert all(exact.get((i,groups//2),F(0))<=value for i,value in enumerate(vector))
    print(dict(directed_power_checks=8,outward_coefficient_checks=checked),flush=True)


def stationary_regions():
    """All original-inactive allocations, but a fixed population every step."""
    start=time.monotonic();table=coefficient_table();local_atoms=atoms()
    coordinates=(('paid',0),('old',4),('new',4))
    for pop in ((104,156,60,64,32,32,64),(104,176,72,32,16,16,96)):
        fit=bound(table,pop,pop,131,return_fit=True)
        scalar,marker,masses=directed(local_atoms,fit,pop,pop,131,coordinates,polynomial=True)
        total=sum(masses);paid=2100-16*pop[0]
        normalized={(i,):p/total for i,p in enumerate(masses)}
        mass=coefficient_upper(normalized,(paid,),16)
        result=16*scalar-2100*marker+log_upper(16,total)+log_upper(1,mass)+log_upper(-1,math.comb(8192,2100))
        print(dict(population=pop,all_quota_allocations=True,original_quota=2100,
                   directed_region_log2_upper=str(result),display=float(result),
                   seconds=time.monotonic()-start,full_certificate=False),flush=True)


def main():
    import sys
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(180,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:
        tests()
        if '--regions' in sys.argv:
            stationary_regions();return
        table=coefficient_table();local_atoms=atoms()
        for pop in ((104,156,60,64,32,32,64),(104,176,72,32,16,16,96)):
            fits={i:bound(table,pop,pop,i,return_fit=True) for i in (131,132)}
            for coordinates in ((('inactive',0),('new',4)),(('paid',0),('old',4),('new',4))):
                values={i:directed(local_atoms,fits[i],pop,pop,i,coordinates) for i in (131,132)}
                regional=12*values[131][0]+4*values[132][0]+log_upper(1,math.comb(16,4))+log_upper(-1,math.comb(8192,2100))
                print(dict(population=pop,coordinates=coordinates,
                           directed_region_log2_upper=str(regional),display=float(regional),
                           terms=values[131][2],seconds=time.monotonic()-start,
                           full_certificate=False),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
