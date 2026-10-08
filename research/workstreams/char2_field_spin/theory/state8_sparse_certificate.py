"""Rigorous Arb replay of the frozen GRS256/176 active-row cover.

This covers exactly 1<=A<=600, not the entire code on its own. It retains
the coarse width16/state8 return envelope from the floating checkpoint.
The stronger joint-return construction bound only decreases that envelope.
Rational marks are selected from saved floating witnesses; no floating
quantity participates in the accepted probability arithmetic.
"""

import argparse
from fractions import Fraction
import json
import math

from flint import arb, arb_poly, ctx, fmpq

from cauchy_firstmoment_screen import geometry
from sparse_active_interval_screen import ANCHORS, FROZEN_INTERVALS, interpolated_anchors
from superregular_return_screen import local_matrix


def ball(value):
    value = Fraction(value)
    return arb(fmpq(value.numerator, value.denominator))


def poly_product(left, right, degree):
    a,b,c,d = left
    e,f,g,h = right
    return ((a*e+b*g).truncate(degree+1), (a*f+b*h).truncate(degree+1),
            (c*e+d*g).truncate(degree+1), (c*f+d*h).truncate(degree+1))


def polynomial_power(local, packets, degree):
    """2x2 matrix polynomial power, truncated only above the needed degree."""
    out = (arb_poly([1]),arb_poly([]),arb_poly([]),arb_poly([1]))
    factor = tuple(arb_poly([ball(local[w][i][j])*math.comb(16,w)
                            for w in range(min(16,degree)+1)])
                   for i,j in ((0,0),(0,1),(1,0),(1,1)))
    while packets:
        if packets&1:
            out = poly_product(out,factor,degree)
        packets >>= 1
        if packets:
            factor = poly_product(factor,factor,degree)
    return out


def matrix_product(left,right):
    a,b,c,d = left
    e,f,g,h = right
    return a*e+b*g,a*f+b*h,c*e+d*g,c*f+d*h


def terminal_power(matrix, exponent):
    result = tuple(map(arb,(1,0,0,1)))
    while exponent:
        if exponent&1:
            result = matrix_product(result,matrix)
        exponent >>= 1
        if exponent:
            matrix = matrix_product(matrix,matrix)
    return result[0]+result[1]


def regional(coefficients, active, rows, tau):
    weights = [ball(Fraction(math.comb(active,h),math.comb(rows,h))*tau**h)
               for h in range(active+1)]
    return tuple(sum((polynomial[h]*weights[h] for h in range(active+1)),arb(0))
                 for polynomial in coefficients)


def source_maximum(lower, upper, rows, parity, tau):
    # Consecutive factor ratio is (L-A)/(A+1)/tau^p. Verify this exact
    # sufficient condition before using the left endpoint as maximizer.
    if tau**parity >= Fraction(rows-lower,lower+1):
        return ball(math.comb(rows,lower))/ball(tau)**(parity*lower+1)
    return max((ball(math.comb(rows,a))/ball(tau)**(parity*a+1)).upper()
               for a in range(lower,upper+1))


def rational_marks(anchor, point, n):
    y = Fraction(math.exp(-math.exp(point[0])*anchor/n)).limit_denominator(10**10)
    tau = min(Fraction(255),Fraction(math.exp(point[1])).limit_denominator(10**8))
    if not 0<y<=1 or not 1<=tau<=255:
        raise ValueError('invalid rational marks')
    return y,tau


def integer_bits(upper):
    if not upper.is_finite() or not upper>0:
        raise ArithmeticError('a finite positive upper bound is required')
    bits = math.floor(float(-upper.log()/arb(2).log()))
    while not upper<=arb(2)**(-bits):
        bits -= 1
    while upper<=arb(2)**(-bits-1):
        bits += 1
    return bits


def certify_anchor(anchor, bins, point, precision=512):
    g = geometry(20,256,176)
    y,tau = rational_marks(anchor,point,g['n'])
    local = [local_matrix(w,y,tau) for w in range(17)]
    degree = max(hi for _,hi in bins)
    old_precision = ctx.prec
    try:
        ctx.prec = precision
        coefficients = polynomial_power(local,g['rows']//16,degree)
        output_mark = ball(y)**(-g['threshold'])
        entries = []
        total = arb(0)
        for lo,hi in bins:
            matrix = regional(coefficients,hi,g['rows'],tau)
            value = (source_maximum(lo,hi,g['rows'],80,tau)*output_mark
                     *terminal_power(matrix,256)*(hi-lo+1))
            upper = value.upper()
            if not upper.is_finite() or not 0<upper<arb(2)**(-80):
                raise ArithmeticError(f'insufficient directed replay for {lo}..{hi} at {precision} bits')
            entries.append(dict(lower=lo,upper=hi,anchor=anchor,
                                certified_integer_bits=integer_bits(upper)))
            total += upper
        return total.upper(), dict(anchor=anchor,precision_bits=precision,
                                   y=str(y),tau=str(tau),entries=entries,
                                   upper_ball=str(total.upper()),
                                   certified_integer_bits=integer_bits(total.upper()))
    finally:
        ctx.prec = old_precision


def certify(progress=None, only_anchor=None):
    if (FROZEN_INTERVALS[0][0]!=1 or FROZEN_INTERVALS[-1][1]!=600
        or any(left[1]+1!=right[0] for left,right in zip(FROZEN_INTERVALS,FROZEN_INTERVALS[1:]))):
        raise ArithmeticError('frozen cover has gaps, overlaps, or wrong endpoints')
    anchors = interpolated_anchors(ANCHORS)
    grouped = {}
    for lo,hi,anchor in FROZEN_INTERVALS:
        grouped.setdefault(anchor,[]).append((lo,hi))
    if only_anchor is not None:
        grouped = {only_anchor:grouped[only_anchor]}
    old_precision = ctx.prec
    try:
        ctx.prec = 256
        total = arb(0)
        receipts = []
        for anchor,bins in grouped.items():
            for precision in (512,1024,2048,4096,8192,16384):
                try:
                    upper,receipt = certify_anchor(anchor,bins,anchors[anchor],precision)
                    break
                except ArithmeticError:
                    if precision==16384:
                        raise
            total += upper
            receipts.append(receipt)
            if progress:
                progress(receipt)
        upper = total.upper()
        bits = integer_bits(upper)
        if bits<40:
            raise ArithmeticError('union fails forty-bit target')
        exact_upper = upper.fmpq()
        g = geometry(20,256,176)
        return dict(outer=[256,176],q=256,width=16,state=8,k_log=20,
                    rows=g['rows'],k=g['k'],n=g['n'],output_threshold=g['threshold'],
                    active_range=[1,600] if only_anchor is None else None,
                    complete_sparse_range=only_anchor is None,
                    interval_count=sum(len(receipt['entries']) for receipt in receipts),
                    failure_upper_ball=str(upper),certified_integer_bits=bits,
                    failure_upper_dyadic=[str(exact_upper.p),str(exact_upper.q)],
                    rigorous=True,whole_code=False,
                    scope='specified active-row intervals only',receipts=receipts)
    finally:
        ctx.prec = old_precision


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--anchor',type=int)
    parser.add_argument('--quiet',action='store_true')
    args = parser.parse_args()
    result = certify(None if args.quiet else lambda value:print(json.dumps(value),flush=True),args.anchor)
    del result['receipts']
    print(json.dumps(result),flush=True)
