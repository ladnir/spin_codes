"""Exact population identity when allowed support equals selected support.

Integer quartet-occupancy tests and FLOAT512 examples. No full certificate.
"""
import itertools
import math
from fractions import Fraction
from quartet_closed_transfer import TYPES,local
from quartet_closed_population import multinomial


def empty_counts(width,zeros):
    """Counts of zero-position subsets, indexed by all-zero quartets."""
    state={(0,0):1}
    for _ in range(width//4):
        after={}
        for (used,empty),count in state.items():
            for j in range(min(4,zeros-used)+1):
                key=(used+j,empty+(j==4))
                after[key]=after.get(key,0)+count*math.comb(4,j)
        state=after
    return {e:count for (n,e),count in state.items() if n==zeros}


def routing(width,oldzeros,newzeros):
    """Expected number of permitted outgoing zero subsets."""
    if newzeros>oldzeros:return Fraction(0)
    count=sum(v*math.comb(4*e,newzeros) for e,v in empty_counts(width,oldzeros).items() if 4*e>=newzeros)
    return Fraction(count,math.comb(width,oldzeros))


def diagonal_log(selected,inactive,logb=2.5,width=512):
    zeros=width-selected
    if zeros%4 or inactive<zeros or inactive>width:return -math.inf
    return (math.log2(math.comb(width//4,zeros//4))-math.log2(math.comb(width,zeros))
            +math.log2(math.comb(selected,inactive-zeros))
            +logb*(707*selected/2560-6*(inactive-zeros)/5))


def tests():
    checked=0
    for width in (4,8,12):
        for oldzeros in range(width+1):
            exact={}
            for positions in itertools.combinations(range(width),oldzeros):
                old=set(positions)
                free=set().union(*(set(range(j,j+4)) for j in range(0,width,4) if all(i in old for i in range(j,j+4))))
                for newzeros in range(len(free)+1):
                    exact[newzeros]=exact.get(newzeros,0)+math.comb(len(free),newzeros)
            for newzeros in range(width+1):
                assert routing(width,oldzeros,newzeros)==Fraction(exact.get(newzeros,0),math.comb(width,oldzeros))
                if newzeros==oldzeros:
                    expected=Fraction(math.comb(width//4,oldzeros//4),math.comb(width,oldzeros)) if oldzeros%4==0 else 0
                    assert routing(width,oldzeros,newzeros)==expected
                checked+=1
    forbidden=deaths=0
    # Type7 (unselected age2) and type16 (a just-ended selected run)
    # require a non-full next input quartet. A selected pending coordinate
    # requires that quartet to be full. They cannot share a quartet.
    for hist in itertools.combinations_with_replacement(TYPES,4):
        if any(t in (7,16) for t in hist) and any(t%3==2 for t in hist):
            assert not local(hist);forbidden+=1
        for (dest,u,v,a,d),ways in local(hist):
            assert sum(t==16 for t in dest)==d;deaths+=1
    print(dict(exact_occupancy_coefficients=checked,
               forbidden_mixed_quartets=forbidden,death_label_identities=deaths),flush=True)


def main():
    tests()
    for selected in (381,384,408,448,480,512):
        vals=[diagonal_log(selected,i) for i in (131,132)]
        print(dict(selected=selected,exact_formula_float_logs=vals),flush=True)


if __name__=='__main__':main()
