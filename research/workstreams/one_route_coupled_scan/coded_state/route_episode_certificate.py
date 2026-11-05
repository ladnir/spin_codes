"""Complete p=32..64 episode cover; floating proposal plus directed replay.

Scope and inequalities: ROUTE_EPISODES.md. No generated files or benchmarks.
"""
import argparse
from decimal import Decimal
from fractions import Fraction
from functools import lru_cache
import math
from pathlib import Path
import sys
import time

theory = Path(__file__).resolve().parent.parent / 'theory'
if theory.is_dir():
    sys.path.insert(0, str(theory))
from stratified_certificate import Upward
from stratified_fractional import entropy_logchoose, LN2, NEG_INF
from route_episode_probe import stable_excess

N, M, n = 1 << 21, 8192, 16384
MARKERS = ((17,16),(9,8),(5,4),(3,2)) + tuple((1 << k,1) for k in (1,2,4,8,16,32,64,128))


@lru_cache(maxsize=4096)
def float_log_factor(v, partials, num, den):
    theta = math.log(num/den)
    excess = stable_excess(partials, v/64, theta)
    if excess == NEG_INF:
        return NEG_INF
    full = excess+math.log1p(math.exp(-excess)) if excess>0 else math.log1p(math.exp(excess))
    return min(full, excess-math.log(-math.expm1(-theta)))


def shape(box, terminal):
    pl,ph,vl,vh,cl,ch = box
    delta = int(terminal)
    vh = min(vh,n-(1-delta))
    vl = max(vl,2*cl-delta)
    ch = min(ch,(vh+delta)//2)
    return (pl,ph,vl,vh,cl,ch) if vl<=vh and cl<=ch else None


def parameters(box, terminal, weight):
    pl,ph,vl,vh,cl,ch = box
    delta = int(terminal)
    partials = min(2*ch-delta,256)
    zero = (vh+63)//64 + partials-1 <=128
    comp_n = vh-(cl-delta)-1
    gap_n = n-vl+ch-1
    comp_k = min(ch-1,comp_n//2)
    gap_k = min(ch-delta,gap_n//2)
    witness = (max(0,128*vl-319*ch+128*delta-weight)+1)//2
    coefficient = Fraction(1,64)-1-64*(cl-delta)-Fraction(107*witness,128)
    return zero,partials,comp_n,comp_k,gap_n,gap_k,coefficient


def float_box(box, terminal, weight):
    pl,ph,vl,vh,cl,ch = box
    zero,r,cn,ck,gn,gk,coefficient = parameters(box,terminal,weight)
    if zero:
        return NEG_INF, MARKERS[0]
    count = entropy_logchoose(M,ph)+entropy_logchoose(cn,ck)+entropy_logchoose(gn,gk)
    best = (float('inf'),None)
    for marker in MARKERS:
        num,den = marker
        log_f = float_log_factor(vh,r,num,den)
        selected_p = ph if log_f>=0 else pl
        value = (count+selected_p*log_f+float(coefficient)*math.log(num/den))/LN2
        if value < best[0]:
            best = value,marker
    return best


def feasible_count(box, terminal):
    pl,ph,vl,vh,cl,ch = box
    delta = int(terminal)
    stop = min(ch,(vh+delta)//2)
    flat = min(stop,(vl+delta)//2)
    count = max(0,flat-cl+1)*(vh-vl+1)
    start = max(cl,flat+1)
    if start<=stop:
        amount = stop-start+1
        count += amount*(vh+delta+1)-amount*(start+stop)
    return (ph-pl+1)*count


def cover(weight, seconds=60, term_bits=80):
    started = time.monotonic()
    boxes = []
    checked = 0
    expected = delivered = 0
    for terminal in (False,True):
        root = shape((32,64,1,n,1,(n+int(terminal))//2),terminal)
        expected += feasible_count(root,terminal)
        stack = [root]
        while stack:
            if time.monotonic()-started>seconds or checked>=60000:
                raise RuntimeError(f'Cover cap: checked={checked}, boxes={len(boxes)}, pending={len(stack)}')
            box = shape(stack.pop(),terminal)
            if box is None:
                continue
            checked += 1
            value,marker = float_box(box,terminal,weight)
            if value <=-term_bits:
                boxes.append((terminal,box,marker,value))
                delivered += feasible_count(box,terminal)
                continue
            pl,ph,vl,vh,cl,ch = box
            if pl==ph and vl==vh and cl==ch:
                raise RuntimeError(f'Unclosed exact cell: {terminal,box,value,marker}')
            spans = (16*(ph-pl),vh-vl,2*(ch-cl))
            axis = max(range(3),key=lambda i:spans[i])
            index = (0,2,4)[axis]
            lo,hi = box[index:index+2]
            mid = (lo+hi)//2
            left,right = list(box),list(box)
            left[index+1],right[index] = mid,mid+1
            stack.extend((tuple(right),tuple(left)))
    assert expected == delivered
    worst = max(boxes,key=lambda row:row[3])
    print(f'Floating complete cover: boxes={len(boxes)}, checks={checked}, cells={delivered}, '
          f'worst={worst}, seconds={time.monotonic()-started:.3f}',flush=True)
    return boxes,expected


@lru_cache(maxsize=256)
def exact_factor(v, r, num, den):
    mu = Fraction(v,64)
    tau = Fraction(num,den)
    best = Fraction(0)
    ceil_mu_r = -((-(mu-r).numerator)//(mu-r).denominator)
    for full in range(max(0,ceil_mu_r),min(256-r,mu.numerator//mu.denominator)+1):
        largest = full+r-128
        if largest<=0:
            continue
        probability = (mu-full)/r
        if probability in (0,1):
            d = max(0,full+(r if probability else 0)-128)
            candidate = tau**d-1
        else:
            a,denominator = probability.numerator,probability.denominator
            b = denominator-a
            tau_base = den**largest
            total = 0
            for j in range(max(0,129-full),r+1):
                d = full+j-128
                prob_num = math.comb(r,j)*a**j*b**(r-j)
                excess_num = num**d*den**(largest-d)-tau_base
                total += prob_num*excess_num
            candidate = Fraction(total,denominator**r*tau_base)
        best = max(best,candidate)
    return min(1+best,best*tau/(tau-1))


def directed(boxes, weight, precision, seconds=90, term_bits=80):
    started = time.monotonic()
    arithmetic = Upward(precision)
    worst = Decimal('-Infinity')
    for terminal,box,marker,_ in boxes:
        if time.monotonic()-started>seconds:
            raise RuntimeError('Directed replay cap')
        pl,ph,vl,vh,cl,ch = box
        zero,r,cn,ck,gn,gk,coefficient = parameters(box,terminal,weight)
        if zero:
            continue
        factor = exact_factor(vh,r,*marker)
        if not factor:
            continue
        count = arithmetic.up.add(arithmetic.entropy_choose(M,ph),arithmetic.entropy_choose(cn,ck))
        count = arithmetic.up.add(count,arithmetic.entropy_choose(gn,gk))
        log_factor_upper = arithmetic.log_fraction(factor)[1]
        selected_p = ph if factor>=1 else pl
        value = arithmetic.up.add(count,arithmetic.up.multiply(Decimal(selected_p),log_factor_upper))
        value = arithmetic.up.add(value,arithmetic.times_log(coefficient,Fraction(*marker)))
        log2 = arithmetic.up.divide(value,arithmetic.log_int(2)[1 if value<0 else 0])
        if log2 > -term_bits:
            raise ArithmeticError(f'Directed failure: {terminal,box,marker,log2}')
        worst = max(worst,log2)
    print(f'Directed {precision} digits: all boxes PASS, worst upper log2={worst}, '
          f'seconds={time.monotonic()-started:.3f}',flush=True)


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--float-only',action='store_true')
    parser.add_argument('--seconds',type=int,default=60)
    args = parser.parse_args()
    weight = 491*N//1000
    boxes,cells = cover(weight,args.seconds)
    if not args.float_only:
        for precision in (80,120):
            directed(boxes,weight,precision)
        margin = 80-(cells-1).bit_length()
        assert margin>=40
        print(f'Complete p=32..64, >491/1000 source-class certificate: failure <=2^-{margin}')
