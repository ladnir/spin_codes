"""Joint positive coefficient extraction for the exact seven-type kernel.

Selected finite-width FLOAT bounds, not a whole-code certificate. Serial
Peach180s, <=512MiB, no generated files. The inactive count is original.
"""
from collections import defaultdict
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


def positive_coefficient(terms,targets,groups):
    """Retain exactly the selected nonnegative degrees; no tail estimate."""
    shape=tuple(t+1 for t in targets)
    assert math.prod(shape)*16<256*(1<<20)
    current=np.zeros(shape);current[(0,)*len(shape)]=1.
    slices=[]
    for degrees,value in terms.items():
        if any(d>t for d,t in zip(degrees,targets)):continue
        source=tuple(slice(0,t+1-d) for d,t in zip(degrees,targets))
        dest=tuple(slice(d,t+1) for d,t in zip(degrees,targets))
        slices.append((source,dest,value))
    for _ in range(groups):
        after=np.zeros_like(current)
        for source,dest,value in slices:after[dest]+=value*current[source]
        current=after
    return float(current[targets])


def extract(fit,incoming,outgoing,inactive,coordinates):
    """Coordinates are ('old',i), ('new',i), ('paid',0), ('inactive',0).

    paid = original inactivity - outgoing Z count. Since that Z count is
    fixed in the requested coefficient, extracting paid is valid even
    when other population coefficients are bounded by positive evaluation.
    """
    def degree(key,coordinate):
        old,new,j=key;kind,index=coordinate
        if kind=='old':return old[index]
        if kind=='new':return new[index]
        if kind=='paid':return j-new[0]
        assert kind=='inactive';return j
    targetkey=(incoming,outgoing,inactive)
    targets=tuple(degree(targetkey,c) for c in coordinates)
    assert min(targets)>=0
    terms=defaultdict(float)
    for key,p in zip(fit['keys'],fit['probabilities']):
        degrees=tuple(degree(key,c) for c in coordinates)
        assert min(degrees)>=0
        terms[degrees]+=p
    mass=positive_coefficient(terms,targets,fit['groups'])
    # Do not interpret numerical underflow as a zero coefficient.
    if mass<=0:raise ArithmeticError('Coefficient mass underflowed')
    return fit['cauchy']+math.log2(mass),math.log2(mass),len(terms)


def tests():
    # Exact integer multivariate recurrence, checked against enumeration.
    terms={(0,0):2,(1,0):3,(0,1):5,(1,1):7,(2,0):11}
    count=0
    for groups in (1,2,3):
        exact=defaultdict(int)
        for sequence in itertools.product(terms,repeat=groups):
            powers=tuple(sum(p[j] for p in sequence) for j in range(2))
            exact[powers]+=math.prod(terms[p] for p in sequence)
        for targets,value in exact.items():
            assert positive_coefficient(terms,targets,groups)==value;count+=1
    # With b=1, quartet weights are integers. Full extraction at width4
    # must reproduce the exact requested kernel coefficient.
    table=coefficient_table(0.);checked=0
    for key,logvalue in itertools.islice(table.items(),0,None,37):
        old,new,j=key
        fit=bound(table,old,new,j,return_fit=True)
        coords=tuple(('old',k) for k in range(6))+tuple(('new',k) for k in range(6))+(('inactive',0),)
        value,_,_=extract(fit,old,new,j,coords)
        # The table excludes the incoming multinomial; it cancels here.
        assert abs(value-logvalue)<1e-8;checked+=1
    print(dict(exact_positive_recurrences=count,full_quartet_extractions=checked),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(180,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:
        tests();table=coefficient_table()
        for pop in ((104,156,60,64,32,32,64),(104,176,72,32,16,16,96)):
            fits={i:bound(table,pop,pop,i,return_fit=True) for i in (131,132)}
            for coordinates in ((('paid',0),),(('inactive',0),('new',4)),
                                (('paid',0),('old',4),('new',4)),
                                (('paid',0),('old',4),('old',5),('new',4))):
                values=[extract(fits[i],pop,pop,i,coordinates) for i in (131,132)]
                regional=12*values[0][0]+4*values[1][0]+math.log2(math.comb(16,4))-math.log2(math.comb(8192,2100))
                print(dict(population=pop,coordinates=coordinates,conditional_region_log2=float(regional),
                           log2_masses=[v[1] for v in values],terms=values[0][2],
                           seconds=time.monotonic()-start,full_certificate=False),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
