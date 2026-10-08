"""Exact dimension obstruction for an unpunctured regional MDS construction.

A [b,k] MDS outer repeats over R independent message rows. Each of b
coordinate regions has an independent uniform permutation of the R rows.
A causal linear inner consumes width-symbol packets and carries s symbols
of state. There is no additional route. For a fixed prefix of W complete
packet positions per region, shortening the outer into these windows leaves
S=sum_i max(0,t_i-(b-k)) dimensions. Enforcing zero state after every
window costs at most b*s dimensions. Thus S>b*s exhibits a nonzero message
whose output is supported in the windows, for EVERY choice of local maps
and diagonals. Initial state must be zero, and zero input/state emits zero.

Each t_i is Binomial(b,W/R). No independence between rows is needed.
Since sum_i t_i=bW and (t-(b-k))_+ <= kt/b, we have S<=kW. Therefore
Pr[S>b*s] >= (E[S]-b*s)/(kW-b*s), when positive.
All quantities in this lower bound are exact Fractions. This is a failure
probability lower bound, not an unsuccessful upper-bound certificate.
"""
import argparse
from fractions import Fraction as F
from itertools import combinations, product
import json
import math


def mean_row_dimension(length,dimension,p):
    if not 1<=dimension<length or not 0<=p<=1:
        raise ValueError('invalid MDS parameters or window fraction')
    p=F(p)
    return sum((t-(length-dimension))*math.comb(length,t)*p**t*(1-p)**(length-t)
               for t in range(length-dimension+1,length+1))


def lower_bound(length,dimension,rows,window,state):
    if rows<1 or not 0<=window<=rows or state<0:
        raise ValueError('invalid row/window/state counts')
    mean=rows*mean_row_dimension(length,dimension,F(window,rows))
    constraints=length*state
    maximum=window*dimension
    probability=max(F(0),(mean-constraints)/(maximum-constraints)) if maximum>constraints else F(0)
    return mean,probability


def screen(length=64,dimension=32,state=4,width=4,distance=F(2,5),k_log=20):
    distance=F(distance)
    if not 0<distance<1 or width<1:
        raise ValueError('invalid distance or packet width')
    unit=math.lcm(width,64)
    rows=((2**k_log+dimension*unit-1)//(dimension*unit))*unit
    strict_count=(distance.numerator*rows-1)//distance.denominator
    window=(strict_count//width)*width
    mean,p=lower_bound(length,dimension,rows,window,state)
    return dict(outer=[length,dimension],state=state,width=width,requested_k=2**k_log,
        padded_k=rows*dimension,rows=rows,output_n=rows*length,rate=str(F(dimension,length)),
        target_distance=str(distance),window=window,window_fraction=str(F(window,rows)),
        bad_weight_at_most=window*length,zero_state_constraints=length*state,
        maximum_surviving_dimension=dimension*window,
        expected_surviving_dimension=float(mean),failure_probability_lower_bound=str(p),
        failure_probability_lower_bound_float=float(p),
        margin_upper_bits=(-math.log2(p.numerator)+math.log2(p.denominator)) if p else None,
        exact=True,distribution='independent uniform regional permutations',
        scope='lower bound for this window event, independent of field size and inner coefficients')


def enumerate_small(length,dimension,rows,window,state):
    choices=list(combinations(range(rows),window))
    mean=0;failures=0;total=len(choices)**length
    if total>1000000:
        raise ValueError('enumeration is for small regression tests only')
    for regions in product(choices,repeat=length):
        count=[0]*rows
        for region in regions:
            for row in region:count[row]+=1
        size=sum(max(0,t-(length-dimension)) for t in count)
        mean+=size
        failures+=size>length*state
    return F(mean,total),F(failures,total)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--outer',type=int,nargs=2,default=[64,32])
    p.add_argument('--state',type=int,default=4)
    p.add_argument('--width',type=int,default=4)
    p.add_argument('--distance',type=F,nargs='+',default=[F(7,20),F(2,5)])
    a=p.parse_args()
    for d in a.distance:print(json.dumps(screen(*a.outer,a.state,a.width,d)),flush=True)
