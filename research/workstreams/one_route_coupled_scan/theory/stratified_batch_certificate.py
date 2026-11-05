"""Finite Toeplitz-batch bound using the exact inverse random-coordinate set.

Float only proposes a complete rectangular cover. The optional directed
replay uses exact rational moments and outward arithmetic from the scalar
certificate. The uniform full-field Toeplitz diagonals are not cheap masks.
"""

import argparse
from decimal import Decimal, localcontext
from fractions import Fraction
import math
import time

import numpy as np

from stratified_fractional import excess_log, entropy_logchoose, LN2, NEG_INF
from stratified_certificate import Upward, exact_excess

MARKERS = ((17,16),(9,8),(5,4),(3,2),(2,1),(4,1),(8,1))


def float_box(N,B,L,b,e,vl,vh,cl,ch,entropy_volume=False):
    n,m = N//b,L//b
    if b*vh*B+2*ch*N <= (B//2)*N:
        return NEG_INF,(5,4)
    free = b*vh-L*(cl-1)
    n1 = vh-m*(cl-1)-1
    n2 = n-vl+ch
    count = entropy_logchoose(n1,min(ch-1,n1//2))+entropy_logchoose(n2,min(ch,n2//2))
    al,ah = max(0,b*(vl-ch)),b*(vh-cl)
    closest = min(ah,N//2) if al <= N//2 else N-al
    volume = min(float(free),e+free/128)
    if entropy_volume and e<free:
        volume = min(volume,e+entropy_logchoose(free,e)/(128*LN2))
    coefficient = volume-min(al,N//2)
    best = (float('inf'),(5,4))
    for marker in MARKERS:
        theta = math.log(marker[0]/marker[1])
        excess = excess_log(B,B//2,ch,closest*B/N,theta)
        routed = (N//B)*float(np.logaddexp(0.0,excess))
        value = (count+theta*coefficient+routed)/LN2
        if value < best[0]:
            best = value,marker
    return best


def cover(B,L,b,distance,term_bits=100,entropy_volume=False,seconds=120):
    N=2**21
    if L%b or N%b:
        raise ValueError('Batch must divide memory and output length')
    n,m=N//b,L//b
    e=(distance*N).numerator//(distance*N).denominator
    queue=[(1,n,1,(n+m)//(m+1))]
    boxes=[]
    started=time.monotonic()
    checked=0
    worst=(NEG_INF,None)
    while queue:
        if time.monotonic()-started > seconds or checked > 50000:
            raise RuntimeError(f'Cover cap reached: {checked} checks; {len(queue)} pending')
        vl,vh,cl,ch=queue.pop()
        vl=max(vl,m*(cl-1)+cl)
        ch=min(ch,(vh+m)//(m+1))
        if vl>vh or cl>ch:
            continue
        checked+=1
        value,marker=float_box(N,B,L,b,e,vl,vh,cl,ch,entropy_volume)
        if value<=-term_bits:
            boxes.append((vl,vh,cl,ch,marker,value))
            if value>worst[0]:
                worst=value,(vl,vh,cl,ch,marker)
            continue
        if vl==vh and cl==ch:
            raise RuntimeError(f'Point does not close: {(vl,cl,value,marker)}')
        if vh-vl >= m*(ch-cl) and vl<vh:
            mid=(vl+vh)//2
            queue.extend(((vl,mid,cl,ch),(mid+1,vh,cl,ch)))
        else:
            mid=(cl+ch)//2
            queue.extend(((vl,vh,cl,mid),(vl,vh,mid+1,ch)))
    print(dict(kind='floating_batch_complete_cover_not_certificate',B=B,L=L,b=b,
               distance=str(distance),term_bits=term_bits,entropy_volume=entropy_volume,
               boxes=len(boxes),checks=checked,worst=worst,
               seconds=round(time.monotonic()-started,3)),flush=True)
    return boxes


class BatchUpward(Upward):
    def rectangle(self,N,B,L,b,e,item,entropy_volume=False):
        vl,vh,cl,ch,marker,_=item
        if b*vh*B+2*ch*N <= (B//2)*N:
            return Decimal('-Infinity')
        n,m=N//b,L//b
        free=b*vh-L*(cl-1)
        n1=vh-m*(cl-1)-1
        n2=n-vl+ch
        count=self.up.add(self.entropy_choose(n1,min(ch-1,n1//2)),
                          self.entropy_choose(n2,min(ch,n2//2)))
        al,ah=max(0,b*(vl-ch)),b*(vh-cl)
        closest=min(ah,N//2) if al<=N//2 else N-al
        volume=min(Fraction(free),Fraction(e)+Fraction(free,128))
        if entropy_volume and e<free:
            denominator=self.down.multiply(Decimal(128),self.log_int(2)[0])
            entropy=self.up.divide(self.entropy_choose(free,e),denominator)
            volume=min(volume,Fraction(self.up.add(Decimal(e),entropy)))
        coefficient=volume-min(al,N//2)
        tau=Fraction(*marker)
        excess=exact_excess(B,ch,closest,N,*marker)
        routed=self.routed(N//B,excess,tau,False)
        result=self.up.add(count,self.times_log(coefficient,tau))
        result=self.up.add(result,routed)
        return self.up.divide(result,self.log_int(2)[1 if result<0 else 0])


def certify(B,L,b,distance,do_directed=True,term_bits=100,entropy_volume=False,seconds=180):
    boxes=cover(B,L,b,distance,term_bits,entropy_volume)
    if not do_directed:
        return
    N=2**21
    n,m=N//b,L//b
    e=(distance*N).numerator//(distance*N).denominator
    for precision in (80,120):
        bound=BatchUpward(precision)
        worst=Decimal('-Infinity')
        started=time.monotonic()
        with localcontext(bound.up):
            for index,item in enumerate(boxes):
                if time.monotonic()-started>seconds:
                    raise RuntimeError('Directed batch replay exceeded time cap')
                value=bound.rectangle(N,B,L,b,e,item,entropy_volume)
                if value>-term_bits:
                    raise ArithmeticError(f'Rectangle fails: {item}, {value}')
                worst=max(worst,value)
            states=n*((n+m)//(m+1))
            print(dict(kind='directed_batch_complete_cover',B=B,L=L,b=b,N=N,
                       threshold=e,distance=str(distance),precision=precision,
                       boxes=len(boxes),worst_log2_upper=str(worst),
                       term_bits=term_bits,entropy_volume=entropy_volume,
                       certified_margin_at_least_bits=term_bits-(states-1).bit_length(),
                       seconds=round(time.monotonic()-started,3)),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--B',type=int,default=256,choices=(64,128,256))
    ap.add_argument('--L',type=int,default=256)
    ap.add_argument('--batch',type=int,default=64)
    ap.add_argument('--distance',type=Fraction,default=Fraction(45,100))
    ap.add_argument('--float-only',action='store_true')
    ap.add_argument('--term-bits',type=int,default=100)
    ap.add_argument('--entropy-volume',action='store_true')
    args=ap.parse_args()
    if args.L<1 or args.batch<1 or args.term_bits<1 or not 0<args.distance<Fraction(1,2):
        ap.error('Require positive memory, batch, term bits, and 0 < distance < 1/2')
    certify(args.B,args.L,args.batch,args.distance,not args.float_only,args.term_bits,args.entropy_volume)
