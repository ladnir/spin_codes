"""Restricted certificate for dense XOR-zero groups with mostly weight-2 columns.

Outer counts use exact shortening bounds. Independent lane shuffles give
a direct output moment, without assumptions about the state distribution.
This does not cover general four-word groups or the complete SPIN code.
"""
import argparse
from fractions import Fraction as Q
from itertools import product
from math import comb,sqrt
from flint import arb,ctx
from joint_support import span


def count_cap(spectrum,dimensions,k,exceptions):
    n=len(spectrum)-1
    assert 0<=exceptions<n
    total=sum(count*(1<<dimensions[v])*sum(comb(n-v,j) for j in range(min(exceptions,n-v)+1))
              for v,count in enumerate(spectrum))
    return min(1<<(3*k),(1<<k)*total)


def lane_moment(z):return (1+4*z*z+z**4)/6


def contribution(q,exceptions,count,threshold=209715):
    windows=(256-exceptions)*q
    if threshold>=2*windows:z_num=1;z_den=1
    else:
        # This floating calculation proposes a rational witness only.
        # The upper bound is valid at every rational 0<z<=1.
        a=threshold/windows
        x=2*a/(sqrt((8-4*a)**2+4*(4-a)*a)+(8-4*a))
        z_den=1<<24;z_num=max(1,min(z_den,round(sqrt(x)*z_den)))
    assert 0<z_num<=z_den
    z=arb(z_num)/z_den
    value=arb(comb(2048,q))*arb(count)**q*z**(-threshold)*lane_moment(z)**windows
    return value,(z_num,z_den)


def self_test():
    checks=0
    for z in (Q(0),Q(1,8),Q(1,4),Q(1,2),Q(3,4),Q(1)):
        masks=[x for x in range(16) if x.bit_count()==2]
        for state in range(16):
            actual=sum(z**(state^x).bit_count() for x in masks)/6
            assert actual<=lane_moment(z);checks+=1
    # Identity and counting tests include codes not containing all ones.
    for basis,n in (([1,2,4],4),([15,51,85],7),([0x97,0x4b,0x2d,0x1e],8)):
        code=span(basis);k=len(basis);ones=(1<<n)-1
        spectrum=[sum(w.bit_count()==v for w in code) for v in range(n+1)]
        dimensions=[0]*(n+1)
        for support in range(1<<n):
            d=sum(w&~support==0 for w in code).bit_length()-1
            dimensions[support.bit_count()]=max(dimensions[support.bit_count()],d)
        shells=[0]*(n+1)
        for a,b,c in product(code,repeat=3):
            fourth=a^b^c;d=a^b;h=ones^a^c
            exceptions=sum(sum((x>>j)&1 for x in (a,b,c,fourth))!=2 for j in range(n))
            assert exceptions==(h&~d&ones).bit_count()
            shells[exceptions]+=1
        for e in range(n):
            actual=sum(shells[:e+1])
            assert actual<=count_cap(spectrum,dimensions,k,e)
            inflated=[v*(2 if w else 1) for w,v in enumerate(spectrum)]
            assert actual<=count_cap(inflated,dimensions,k,e)
            checks+=2
    print('Even-column exact lane-moment/identity/count checks:',checks,'passed',flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exceptions',type=int,nargs='+',default=[0,4,8,16])
    parser.add_argument('--precision',type=int,default=192)
    parser.add_argument('--budget',type=int,default=60)
    parser.add_argument('--positive-shortening',action='store_true')
    parser.add_argument('--dual-shortening',action='store_true')
    args=parser.parse_args();ctx.prec=args.precision;self_test()
    from bch_joint_support import authenticated_caps
    from shortened_bound import dimension_caps
    spectrum=authenticated_caps();dimensions=dimension_caps()
    if args.positive_shortening:
        from shortening_polynomial import improve_dimensions
        dimensions=improve_dimensions(dimensions)
    if args.dual_shortening:
        from dual_shortening import improve_dimensions, self_test as dual_test
        dual_test();dimensions=improve_dimensions(dimensions)
    # Imported legacy proof modules can change flint's process-global context.
    ctx.prec=args.precision
    assert ctx.prec == args.precision
    for exceptions in args.exceptions:
        count=count_cap(spectrum,dimensions,128,exceptions)
        print('Exceptions per group',exceptions,'log2 count upper',arb(count).log()/arb(2).log(),flush=True)
        total=arb(0);first=None;accepted=None;witness=None
        for q in range(2048,0,-1):
            value,z=contribution(q,exceptions,count)
            total+=value
            if total<arb(2)**(-args.budget):first=q;accepted=total;witness=z
        if first is None:
            value,z=contribution(2048,exceptions,count)
            print('NOT CLOSED even at q=2048: log2 bound',value.log()/arb(2).log(),flush=True)
        else:
            print('VERIFIED restricted even-column family q',first,'through 2048; precision',args.precision,
                  'upper',accepted,'margin',-accepted.log()/arb(2).log(),'first witness',witness,flush=True)
    print('Only XOR-zero groups with the stated per-group exception cap; no full-code certificate.',flush=True)


if __name__=='__main__':main()
