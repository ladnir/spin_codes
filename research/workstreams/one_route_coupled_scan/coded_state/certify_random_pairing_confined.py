"""Directed restricted certificate: two constant-support suffix families.

GF(2^32), delta .46, unchanged RS/route/SL2 construction except random
matchings. Not a whole-code certificate. Serial <=30s, source-only.
"""
from fractions import Fraction as F
import math
import time
from diagnose_random_pairing import crossing_law

SCALE=128


def logs(value):
    """Directed log2 endpoints in units 1/SCALE, by integer powers only."""
    assert value>0
    a=value.numerator**SCALE;b=value.denominator**SCALE
    difference=a.bit_length()-b.bit_length()
    return difference-1,difference+1


def ceildiv(a,b):return -(-a//b)


def certify(R,w,cutoff,tilt):
    n=512;M=8192;q=1<<32;steps=16*R
    weight=460*(1<<21)//1000
    masks=(math.comb(n,w)-1).bit_length()
    budget=80+17+masks
    # Exact Bin(R,w/n) marginal. Regional selected sets are independent
    # uniform subsets; row functions are negatively associated, not independent.
    moment=sum((F(math.comb(R,b)*w**b*(n-w)**(R-b),n**R)
                *tilt**max(0,b-cutoff) for b in range(R+1)),F(0))
    lo_tilt,_=logs(tilt);_,hi_moment=logs(moment)
    assert lo_tilt>0
    cap=ceildiv(M*hi_moment+budget*SCALE,lo_tilt)
    assert M*hi_moment-cap*lo_tilt<=-budget*SCALE
    # U(b)=1+2^(R-4096+31b); log U is 31-Lipschitz and increasing.
    exponent=R-4096+31*cutoff
    _,epsilon=logs(1+F(2)**exponent)
    row_enumerator_scaled=M*epsilon+31*cap*SCALE
    maximum=max(prob/F(math.comb(n-w,k)) for k,prob in crossing_law(n,w).items())
    _,matching_upper=logs(maximum)
    zeros=max(0,steps*w-weight)
    assert zeros<=(steps-1)*w
    zero_choices=(math.comb((steps-1)*w,zeros)-1).bit_length()
    good_scaled=((1+masks+zero_choices-32*zeros)*SCALE
                 +(steps-1)*matching_upper+row_enumerator_scaled)
    good_bits=ceildiv(good_scaled,SCALE)
    assert good_bits<=-97
    # For this R,w, route failure after all masks is <=2^-97; the
    # good-route bad-message contribution is <=2^good_bits<=2^-97.
    return dict(R=R,w=w,cutoff=cutoff,tilt=str(tilt),excess_cap=cap,
                good_route_upper_bits=good_bits,route_union_upper_bits=-97,
                restricted_failure_upper_bits=-96,output_threshold=weight,
                masks_and_messages='all',whole_code_certificate=False)


def main():
    began=time.monotonic()
    for args in ((192,256,120,F(9,8)),(256,256,124,F(33,32))):
        print(certify(*args),flush=True)
        if time.monotonic()-began>25:raise TimeoutError('restricted replay budget')
    print(dict(seconds=time.monotonic()-began,joint_two_family_upper_bits=-95))


if __name__=='__main__':main()
