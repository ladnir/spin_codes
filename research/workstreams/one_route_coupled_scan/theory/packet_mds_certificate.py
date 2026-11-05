"""Bounded packet-route extension of the frozen persistent MDS certificate.

The float cover is only a proposal. Optional directed replays use the frozen
exact regional moments. No data files are written. See PACKET_ROUTE.md. Only the route moment changes.
"""

import argparse
from decimal import Decimal, localcontext
from fractions import Fraction
import math
import time

import numpy as np

from stratified_fractional import excess_log, entropy_logchoose, LN2, NEG_INF
from stratified_certificate import Upward, exact_excess

PACKET = 4  # CLI-fixed ensemble: independent outer rows, not payload lanes.

MARKERS = ((17,16), (9,8), (5,4), (3,2), (2,1), (4,1), (8,1))
FINE_MARKERS = ((17,16), (9,8), (19,16), (5,4), (21,16), (11,8),
                (23,16), (3,2), (2,1), (4,1), (8,1))


def c_limit(n, b, e):
    return min((n+1)//2, (e-1)//(b+1)+1)


def float_box(N, B, b, bits, e, vl, vh, cl, ch):
    if b*vh*B+2*ch*N <= (B//2)*N:
        return NEG_INF, 'inverse', (5,4)
    n = N//b
    n1, n2 = vh-cl, n-vl+ch
    count = (entropy_logchoose(n1, min(ch-1, n1//2))
             + entropy_logchoose(n2, min(ch, n2//2)))
    free = b*vh
    volume = min(float(free), e+free/bits)
    if e < free:
        volume = min(volume, e+entropy_logchoose(free, e)/(bits*LN2))
    volume -= b*(cl-1)
    al, ah = max(0,b*(vl-ch)), b*(vh-cl)
    ul, uh = b*vl, b*vh
    closest_a = min(ah,N//2) if al <= N//2 else N-al
    closest_u = min(uh,N//2) if ul <= N//2 else N-ul
    best = (float('inf'), 'inverse', (5,4))
    for marker in MARKERS:
        theta = math.log(marker[0]/marker[1])
        for kind, selected, coefficient in (
            ('inverse', closest_a, volume-min(al,N//2)),
            ('forward', closest_u, max(0,uh-N//2)-b*(cl-1)),
        ):
            excess = excess_log(B,B//2,ch,selected*B/N,theta*PACKET)
            routed = (N//(B*PACKET))*float(np.logaddexp(0.0,excess))
            if kind == 'forward' and uh <= N//2:
                corrected = NEG_INF if excess == NEG_INF else float('inf') if excess>700 else (
                    math.log(N//(B*PACKET))+excess+(N//(B*PACKET))*math.exp(excess)
                    -math.log(math.expm1(theta)))
                routed = min(routed,corrected)
            value = (count+theta*coefficient+routed)/LN2
            if value < best[0]:
                best = value,kind,marker
    return best


def cover(N,B,b,bits,distance,term_bits,seconds):
    if b%PACKET or (N//B)%PACKET:
        raise ValueError('Packets must divide batch size and regional length')
    n=N//b
    e=int(distance*N)
    queue=[(1,n,1,c_limit(n,b,e))]
    boxes=[]
    started=time.monotonic()
    checked=0
    worst=(NEG_INF,None)
    while queue:
        if time.monotonic()-started > seconds or checked >= 50000:
            raise RuntimeError(f'Cover cap: {checked} checks; {len(queue)} pending')
        vl,vh,cl,ch=queue.pop()
        vl=max(vl,2*cl-1)
        ch=min(ch,(vh+1)//2,c_limit(n,b,e))
        if vl>vh or cl>ch:
            continue
        checked+=1
        value,kind,marker=float_box(N,B,b,bits,e,vl,vh,cl,ch)
        if value<=-term_bits:
            boxes.append((vl,vh,cl,ch,kind,marker,value))
            if value>worst[0]:
                worst=value,(vl,vh,cl,ch,kind,marker)
            continue
        if vl==vh and cl==ch:
            raise RuntimeError(f'Point does not close: {(vl,cl,value,kind,marker)}')
        if vh-vl >= 2*(ch-cl) and vl<vh:
            mid=(vl+vh)//2
            queue.extend(((vl,mid,cl,ch),(mid+1,vh,cl,ch)))
        else:
            mid=(cl+ch)//2
            queue.extend(((vl,vh,cl,mid),(vl,vh,mid+1,ch)))
    print(dict(kind='floating_packet_cover_not_certificate',N=N,B=B,b=b,
               field_bits=bits,packet=PACKET,distance=str(distance),term_bits=term_bits,
               boxes=len(boxes),checks=checked,worst=worst,
               forward_boxes=sum(x[4]=='forward' for x in boxes),
               seconds=round(time.monotonic()-started,3)),flush=True)
    return boxes


class PersistentUpward(Upward):
    def projective_route(self,M,excess,tau):
        if excess == 0:
            return Decimal('-Infinity')
        ordinary=self.routed(M,excess,tau,False)
        # (exp(M log(1+x))-1)/(tau-1) <= M*x*exp(M*x)/(tau-1).
        corrected=self.up.add(self.log_int(M)[1],self.log_fraction(excess)[1])
        xupper=self.up.divide(Decimal(excess.numerator),Decimal(excess.denominator))
        corrected=self.up.add(corrected,self.up.multiply(Decimal(M),xupper))
        corrected=self.up.subtract(corrected,self.log_fraction(tau-1)[0])
        # Ordinary (nonprojective) counting is also valid, independently.
        return min(ordinary,corrected)

    def rectangle(self,N,B,b,bits,e,item):
        vl,vh,cl,ch,kind,marker,_=item
        if b*vh*B+2*ch*N <= (B//2)*N:
            return Decimal('-Infinity')
        n=N//b
        n1,n2=vh-cl,n-vl+ch
        count=self.up.add(self.entropy_choose(n1,min(ch-1,n1//2)),
                          self.entropy_choose(n2,min(ch,n2//2)))
        tau=Fraction(*marker)
        if kind=='inverse':
            free=b*vh
            volume=min(Fraction(free),Fraction(e)+Fraction(free,bits))
            if e<free:
                den=self.down.multiply(Decimal(bits),self.log_int(2)[0])
                entropy=self.up.divide(self.entropy_choose(free,e),den)
                volume=min(volume,Fraction(self.up.add(Decimal(e),entropy)))
            al,ah=max(0,b*(vl-ch)),b*(vh-cl)
            selected=min(ah,N//2) if al<=N//2 else N-al
            coefficient=volume-b*(cl-1)-min(al,N//2)
            excess=exact_excess(B,ch,selected,N,marker[0]**PACKET,marker[1]**PACKET)
            routed=self.routed(N//(B*PACKET),excess,tau,False)
        elif kind=='forward':
            ul,uh=b*vl,b*vh
            selected=min(uh,N//2) if ul<=N//2 else N-ul
            coefficient=Fraction(max(0,uh-N//2)-b*(cl-1))
            excess=exact_excess(B,ch,selected,N,marker[0]**PACKET,marker[1]**PACKET)
            routed=(self.projective_route(N//(B*PACKET),excess,tau) if uh<=N//2
                    else self.routed(N//(B*PACKET),excess,tau,False))
        else:
            raise ValueError(kind)
        if routed.is_infinite():
            return routed
        result=self.up.add(count,self.times_log(coefficient,tau))
        result=self.up.add(result,routed)
        return self.up.divide(result,self.log_int(2)[1 if result<0 else 0])


def certify(N,B,b,bits,distance,term_bits,directed,seconds,sum_boxes=False):
    boxes=cover(N,B,b,bits,distance,term_bits,seconds)
    if not directed:
        return boxes
    e=int(distance*N)
    for precision in (80,120):
        bound=PersistentUpward(precision)
        worst=Decimal('-Infinity')
        # Exact dyadic upper sum: each rectangle's per-class log2 bound
        # is rounded up to an integer exponent, then multiplied by its
        # integer area. Capping the denominator at 2^256 only enlarges it.
        mass_numerator=0
        started=time.monotonic()
        with localcontext(bound.up):
            for item in boxes:
                if time.monotonic()-started>seconds:
                    raise RuntimeError('Directed replay time cap')
                value=bound.rectangle(N,B,b,bits,e,item)
                if value>-term_bits:
                    raise ArithmeticError(f'Rectangle fails: {item}, {value}')
                worst=max(worst,value)
                if sum_boxes and value.is_finite():
                    area=(item[1]-item[0]+1)*(item[3]-item[2]+1)
                    exponent=min(256,int((-value).to_integral_value(rounding='ROUND_FLOOR')))
                    mass_numerator+=area << (256-exponent)
            states=(N//b)*c_limit(N//b,b,e)
            print(dict(kind='directed_packet_complete_cover',N=N,B=B,b=b,
                       field_bits=bits,packet=PACKET,threshold=e,precision=precision,
                       boxes=len(boxes),worst_log2_upper=str(worst),
                       certified_margin_at_least_bits=term_bits-(states-1).bit_length(),
                       seconds=round(time.monotonic()-started,3)),flush=True)
            if sum_boxes:
                print(dict(kind='directed_dyadic_rectangle_sum',precision=precision,
                           packet=PACKET,B=B,b=b,threshold=e,boxes=len(boxes),
                           numerator=mass_numerator,denominator_power=256,
                           certified_margin_at_least_bits=(256-(mass_numerator-1).bit_length()
                                                           if mass_numerator else 'infinite'),
                           closes_40_bits=mass_numerator <= (1<<216)),flush=True)
    return boxes


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--packet',type=int,default=4,choices=(1,2,4,8))
    ap.add_argument('--B',type=int,default=64,choices=(64,128,256))
    ap.add_argument('--batch',type=int,default=16,choices=(8,16,32,64,128))
    ap.add_argument('--field-bits',type=int,default=128,choices=(32,64,128))
    ap.add_argument('--logK',type=int,default=20,choices=range(12,21))
    ap.add_argument('--distance',type=Fraction,default=Fraction(35,100))
    ap.add_argument('--term-bits',type=int,default=100)
    ap.add_argument('--seconds',type=int,default=180)
    ap.add_argument('--float-only',action='store_true')
    ap.add_argument('--fine-markers',action='store_true',
                    help='Refine the bounded rational marker grid near one')
    ap.add_argument('--sum-boxes',action='store_true',
                    help='Also sum exact dyadic rectangle upper bounds, using their own areas')
    args=ap.parse_args()
    PACKET=args.packet
    if args.fine_markers:
        MARKERS=FINE_MARKERS
    if not 0<args.distance<Fraction(1,2) or not 1<=args.term_bits<=1000:
        ap.error('Require 0<distance<1/2 and 1<=term-bits<=1000')
    if not 1<=args.seconds<=180:
        ap.error('Time cap must be 1..180 seconds')
    certify(2**(args.logK+1),args.B,args.batch,args.field_bits,args.distance,
            args.term_bits,not args.float_only,args.seconds,args.sum_boxes)
