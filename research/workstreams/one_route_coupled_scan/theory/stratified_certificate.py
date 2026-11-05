"""Exact rational route moments and outward-rounded complete finite cover.

The floating module proposes only rectangles and rational moment markers.
Every accepted rectangle is rechecked independently with exact route moments
and outward-rounded Decimal logarithms. No artifacts or caches are written.
"""

import argparse
from decimal import Decimal, Context, ROUND_CEILING, ROUND_FLOOR, localcontext
from fractions import Fraction
from functools import lru_cache
import math
import time

from stratified_fractional import propose_cover


def ceil_fraction(x):
    return -((-x.numerator)//x.denominator)


@lru_cache(maxsize=4096)
def exact_excess(B, c, selected, N, tau_num, tau_den):
    """Exact maximum of the padded-binomial E[tau^shortening-1]."""
    r = min(2*c, B)
    mu = Fraction(selected*B, N)
    best = Fraction(0)
    for f in range(max(0, ceil_fraction(mu-r)), min(B-r, mu.numerator//mu.denominator)+1):
        largest = f+r-B//2
        if largest <= 0:
            continue
        p = (mu-f)/r
        if p == 0 or p == 1:
            d = max(0, f+(r if p else 0)-B//2)
            candidate = Fraction(tau_num, tau_den)**d-1
        else:
            a, den = p.numerator, p.denominator
            b = den-a
            tau_base = tau_den**largest
            total = 0
            for j in range(max(0, B//2-f+1), r+1):
                d = f+j-B//2
                probability_numerator = math.comb(r,j)*a**j*b**(r-j)
                excess_numerator = tau_num**d*tau_den**(largest-d)-tau_base
                total += probability_numerator*excess_numerator
            candidate = Fraction(total, den**r*tau_base)
        best = max(best, candidate)
    return best


class Upward:
    def __init__(self, precision):
        self.up = Context(prec=precision, rounding=ROUND_CEILING)
        self.down = Context(prec=precision, rounding=ROUND_FLOOR)
        self._logs = {}

    def log_int(self, n):
        if n not in self._logs:
            # Decimal.ln is correctly rounded to nearest; adjacent numbers
            # enclose its exact value irrespective of the context rounding.
            nearest = self.up.ln(Decimal(n))
            self._logs[n] = (self.down.next_minus(nearest), self.up.next_plus(nearest))
        return self._logs[n]

    def log_fraction(self, x):
        nl, nu = self.log_int(x.numerator)
        dl, du = self.log_int(x.denominator)
        return self.down.subtract(nl, du), self.up.subtract(nu, dl)

    def entropy_choose(self, n, k):
        if k == 0 or k == n:
            return Decimal(0)
        if not 0 < k < n:
            raise ValueError((n,k))
        result = self.up.multiply(Decimal(n), self.log_int(n)[1])
        result = self.up.add(result, self.up.multiply(Decimal(-k), self.log_int(k)[0]))
        return self.up.add(result, self.up.multiply(Decimal(k-n), self.log_int(n-k)[0]))

    def times_log(self, coefficient, base):
        lo, hi = self.log_fraction(base)
        value = self.up.multiply(Decimal(coefficient.numerator), lo if coefficient < 0 else hi)
        return self.up.divide(value, Decimal(coefficient.denominator))

    def routed(self, M, excess, tau, require_positive):
        if excess == 0:
            return Decimal('-Infinity') if require_positive else Decimal(0)
        full = self.up.multiply(Decimal(M), self.log_fraction(1+excess)[1])
        if not require_positive:
            return full
        # exp(M log(1+x))-1 <= M*x*exp(M*x). This avoids loss of
        # significance when x is tiny. Take the smaller of two valid bounds.
        value = self.up.add(self.log_int(M)[1], self.log_fraction(excess)[1])
        xupper = self.up.divide(Decimal(excess.numerator), Decimal(excess.denominator))
        value = self.up.add(value, self.up.multiply(Decimal(M), xupper))
        value = self.up.subtract(value, self.log_fraction(1-1/tau)[0])
        return min(full, value)

    def rectangle(self, N, B, L, e, item):
        ul, uh, cl, ch, marker, _ = item
        if ul > uh or cl > ch:
            raise ValueError('Invalid proposed rectangle')
        if uh*B+2*ch*N <= (B//2)*N:
            return Decimal('-Infinity')
        free = uh-L*(cl-1)
        if free < cl:
            raise ValueError('Proposal retained an infeasible rectangle')
        n1 = free-1
        n2 = N-ul+ch
        count = self.up.add(self.entropy_choose(n1,min(ch-1,n1//2)),
                            self.entropy_choose(n2,min(ch,n2//2)))
        w = min(Fraction(free), Fraction(e)+Fraction(free,128))
        coefficient = ch-min(ul,N//2)+w
        closest = min(uh,N//2) if ul <= N//2 else N-ul
        tau = Fraction(*marker)
        if not 1 < tau <= 2**128:
            raise ValueError('Invalid fractional moment')
        excess = exact_excess(B,ch,closest,N,*marker)
        route = self.routed(N//B,excess,tau,uh <= N//2)
        if route.is_infinite():
            return route
        result = self.up.add(count,self.times_log(coefficient,tau))
        result = self.up.add(result,route)
        ln2 = self.log_int(2)[1 if result < 0 else 0]
        return self.up.divide(result,ln2)


def certify(B=256, L=256, distance=Fraction(45,100), seconds=180):
    N = 2**21
    e = (distance*N).numerator//(distance*N).denominator
    # e is checked again here; floating distance only chooses a partition.
    boxes = propose_cover(B,L,float(distance),threshold=-100)
    for precision in (80,120):
        bound = Upward(precision)
        worst = Decimal('-Infinity')
        started = time.monotonic()
        with localcontext(bound.up):
            for index,item in enumerate(boxes):
                if time.monotonic()-started > seconds:
                    raise RuntimeError('Directed replay exceeded its time cap')
                value = bound.rectangle(N,B,L,e,item)
                if value > -100:
                    raise ArithmeticError(f'Rectangle fails: {item}, upper={value}')
                worst = max(worst,value)
                if (index+1)%400 == 0:
                    print(dict(precision=precision,verified=index+1,total=len(boxes)),flush=True)
            states = N*((N+L)//(L+1))
            # There are fewer than 2^ceil(log2(states)) candidate (u,c) pairs.
            union_bits = 100-(states-1).bit_length()
            print(dict(kind='directed_complete_cover', B=B,L=L,N=N,
                       threshold=e,distance=str(distance),precision=precision,
                       boxes=len(boxes),worst_log2_upper=str(worst),
                       certified_margin_at_least_bits=union_bits,
                       seconds=round(time.monotonic()-started,3)),flush=True)
    return boxes


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--B',type=int,default=256,choices=(64,128,256))
    ap.add_argument('--L',type=int,default=256)
    ap.add_argument('--distance',type=Fraction,default=Fraction(45,100))
    args = ap.parse_args()
    certify(args.B,args.L,args.distance)
