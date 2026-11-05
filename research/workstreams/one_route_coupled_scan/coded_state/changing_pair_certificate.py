"""Directed bound for at most d changing pairs, not a whole-code certificate."""
from decimal import localcontext
from fractions import Fraction
import math
import time
from changing_pair_probe import P, point
from feedback_support_probe import N, K, W, MASK_N
from feedback_support_certificate import Interval, rational, entropy_choose, poly_value, route_cap


def mask_cap(d,a,b):
    # In this verifier's domain b < MASK_N/2, so entropy_choose increases.
    choices = [entropy_choose(MASK_N,b).hi]
    for tilt in (Fraction(1,256),Fraction(1,64),Fraction(1,16),Fraction(1,8),
                 Fraction(1,4),Fraction(1,2),Fraction(1),Fraction(2),
                 Fraction(4),Fraction(16),Fraction(256)):
        x = rational(tilt)
        h = b if tilt<1 else a
        choices.append(((P-d)*poly_value(256,x).log()
                        +(2*d+1)*poly_value(128,x).log()-h*x.log()).hi)
    if a==b==63*(8192+d):
        choices.append(((P-d)*Interval(math.comb(256,63)).log()
                        +(2*d+1)*Interval(math.comb(128,63)).log()).hi)
    return min(choices)


def box(d,a,b,tau):
    assert 0<=d<=256 and 0<=a<=b<=63*(8192+d)<MASK_N//2
    assert isinstance(tau,int) and 1<tau<=2**128
    U = 128*(8193+d)
    reduction = (64*(U-W)+126)//127
    va,vb = min(W,U-a-reduction),min(W,U-b-reduction)
    support = Interval(U+1).log()+entropy_choose(U,max(U//2,vb))
    alpha = Interval(tau).log()/(128*Interval(2).log())
    # The choice of exceptional positions, like the mask, is outside alpha.
    value = (Interval(math.comb(P,d)).log()+Interval(mask_cap(d,a,b))
             +Interval(tau).log()*(va-K)+alpha*support+route_cap(b,tau))
    return (value/Interval(2).log()).hi


def certify(d=32,precision=80,seconds=45):
    if not 0<=d<=256:
        raise ValueError('Bounded verifier supports 0<=d<=256')
    started = time.monotonic()
    hmax = 63*(8192+d)
    pending,accepted = [(0,hmax)],[]
    with localcontext() as ctx:
        ctx.prec = precision
        while pending:
            if time.monotonic()-started>seconds:
                raise RuntimeError('Certificate time cap: no result')
            a,b = pending.pop()
            tau = point(d,(a+b)//2)[1]  # FLOAT proposes, interval decides.
            upper = box(d,a,b,tau)
            if upper<=-80:
                accepted.append((a,b,tau,upper))
            else:
                if a==b:
                    raise RuntimeError(f'Unclosed d={d}, H={a}: log2 upper {upper}')
                mid = (a+b)//2
                pending.extend(((mid+1,b),(a,mid)))
        accepted.sort()
        assert accepted[0][0]==0 and accepted[-1][1]==hmax
        assert all(x[1]+1==y[0] for x,y in zip(accepted,accepted[1:]))
        assert sum(b-a+1 for a,b,_,_ in accepted)==hmax+1
        union_bits = 80-(hmax+1).bit_length()
        worst = max(accepted,key=lambda x:x[3])
        print(f'd<={d}, {precision} digits: {len(accepted)} boxes cover all '
              f'{hmax+1} hole totals; class failure <2^-{union_bits}.',flush=True)
        print(f'Weakest box H={worst[0]}..{worst[1]}, tau={worst[2]}, '
              f'log2 upper={worst[3]}.',flush=True)
    print(f'Elapsed {time.monotonic()-started:.3f}s; restricted class, no files.')
    return accepted


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--changes',type=int,default=32)
    parser.add_argument('--precision',type=int,default=80)
    parser.add_argument('--seconds',type=float,default=45)
    args = parser.parse_args()
    certify(args.changes,args.precision,args.seconds)
