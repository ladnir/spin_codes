"""Directed full cover of ONE trajectory class. No whole-code claim or files.

Decimal ln is correctly rounded (half-even); adjacent representable numbers
enclose it. All arithmetic operations use directed endpoints.
"""
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction
import math
import time
from feedback_support_probe import N, K, M, U, W, MASK_N, HMAX, point
from feedback_support import binomial_moment


class Interval:
    def __init__(self, lo, hi=None):
        self.lo = Decimal(lo)
        self.hi = self.lo if hi is None else Decimal(hi)

    @staticmethod
    def cast(x):
        return x if isinstance(x,Interval) else Interval(x)

    def __add__(self, other):
        other = self.cast(other)
        with localcontext() as c:
            c.rounding = ROUND_FLOOR
            lo = self.lo+other.lo
            c.rounding = ROUND_CEILING
            hi = self.hi+other.hi
        return Interval(lo,hi)
    __radd__ = __add__

    def __neg__(self):
        # copy_negate is exact, unlike unary minus under the active context.
        return Interval(self.hi.copy_negate(),self.lo.copy_negate())

    def __sub__(self,other):
        return self+-self.cast(other)

    def __mul__(self,other):
        other = self.cast(other)
        with localcontext() as c:
            c.rounding = ROUND_FLOOR
            lo = min(a*b for a in (self.lo,self.hi) for b in (other.lo,other.hi))
            c.rounding = ROUND_CEILING
            hi = max(a*b for a in (self.lo,self.hi) for b in (other.lo,other.hi))
        return Interval(lo,hi)
    __rmul__ = __mul__

    def __truediv__(self,other):
        other = self.cast(other)
        assert other.lo>0
        with localcontext() as c:
            c.rounding = ROUND_FLOOR
            lo = 1/other.hi
            c.rounding = ROUND_CEILING
            hi = 1/other.lo
        return self*Interval(lo,hi)

    def log(self):
        assert self.lo>0
        return Interval(self.lo.ln().next_minus(),self.hi.ln().next_plus())


def rational(x):
    x = Fraction(x)
    return Interval(x.numerator)/x.denominator


def entropy_choose(n,k):
    """Natural logarithm of the standard entropy upper bound for binom(n,k)."""
    if k==0 or k==n:
        return Interval(0)
    return n*Interval(n).log()-k*Interval(k).log()-(n-k)*Interval(n-k).log()


def poly_value(n,x):
    result = Interval(0)
    for j in range(63,-1,-1):
        result = result*x+math.comb(n,j)
    return result


def mask_cap(a,b):
    choices = [entropy_choose(MASK_N,b).hi]
    for tilt in (Fraction(1,256),Fraction(1,64),Fraction(1,16),Fraction(1,8),
                 Fraction(1,4),Fraction(1,2),Fraction(1),Fraction(2),
                 Fraction(4),Fraction(16),Fraction(256)):
        x = rational(tilt)
        h = b if tilt<1 else a
        value = 8191*poly_value(256,x).log()+poly_value(128,x).log()-h*x.log()
        choices.append(value.hi)
    if a==b==HMAX:
        choices.append((8191*Interval(math.comb(256,63)).log()
                        +Interval(math.comb(128,63)).log()).hi)
    return min(choices)


def route_cap(h,tau):
    """Natural log upper bound for E[tau^D(S^c)], |S^c|=128+h."""
    return M*rational(binomial_moment(256,128,128+h,N,tau)).log()


def box(a,b,tau):
    assert 0<=a<=b<=HMAX and isinstance(tau,int) and 1<tau<=2**128
    reduction = (64*(U-W)+126)//127
    va,vb = min(W,U-a-reduction),min(W,U-b-reduction)
    # binom(U,v) increases toward U/2; the dimension term uses va.
    support = Interval(U+1).log()+entropy_choose(U,max(U//2,vb))
    alpha = Interval(tau).log()/(128*Interval(2).log())
    value = Interval(mask_cap(a,b))+Interval(tau).log()*(va-K)+alpha*support+route_cap(b,tau)
    return (value/Interval(2).log()).hi


def certify(precision=80, seconds=55):
    start = time.monotonic()
    pending, accepted = [(0,HMAX)], []
    with localcontext() as ctx:
        ctx.prec = precision
        while pending:
            if time.monotonic()-start>seconds:
                raise RuntimeError('Certificate time cap: no result')
            a,b = pending.pop()
            tau = point((a+b)//2)[1]  # FLOAT proposal only; box() is decisive.
            upper = box(a,b,tau)
            if upper<=-80:
                accepted.append((a,b,tau,upper))
            else:
                if a==b:
                    raise RuntimeError(f'Unclosed hole count {a}: {upper}')
                mid = (a+b)//2
                pending.extend(((mid+1,b),(a,mid)))
        accepted.sort()
        assert accepted[0][0]==0 and accepted[-1][1]==HMAX
        assert all(x[1]+1==y[0] for x,y in zip(accepted,accepted[1:]))
        assert sum(b-a+1 for a,b,_,_ in accepted)==HMAX+1
        assert HMAX+1 < 2**19
        print(f'{precision} digits: {len(accepted)} boxes cover all {HMAX+1} hole counts; '
              f'every count <=2^-80, class union <2^-61. '
              f'Largest box upper log2={max(x[3] for x in accepted)}',flush=True)
        worst = max(accepted,key=lambda entry:entry[3])
        print(f'Weakest box: H={worst[0]}..{worst[1]}, tau={worst[2]}.',flush=True)
    print(f'Elapsed {time.monotonic()-start:.3f}s. Restricted class only; no files.')
    return accepted


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--precision',type=int,default=80)
    parser.add_argument('--seconds',type=float,default=55)
    args = parser.parse_args()
    certify(args.precision,args.seconds)
