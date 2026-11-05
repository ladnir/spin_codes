"""Directed grouped-reset certificates; bounded source-only serial runs.

The proposal is FLOAT. Replay uses exact rational route moments and
outward-rounded logarithms. No generated files or frozen-helper edits.
"""
import argparse
from decimal import Decimal, localcontext
from fractions import Fraction
import time

from striped_reset_probe import cover
from ideal_reset_certificate import partition, controls, retained


class GroupedUpward(retained.PersistentUpward):
    def rectangle_grouped(self, N, B, b, d, e, item, density_bits=0, field_bits=128, sum_rank=False):
        vl, vh, cl, ch, kind, marker, _ = item
        if b*vh*B+2*ch*N <= (B//2)*N:
            return Decimal('-Infinity')
        n1, n2 = vh-cl, N//b-vl+ch
        count = self.up.add(self.entropy_choose(n1, min(ch-1, n1//2)),
                            self.entropy_choose(n2, min(ch, n2//2)))
        tau = Fraction(*marker)
        if field_bits < 16 or not 1 < tau <= 2**field_bits:
            raise ValueError('Field or fractional marker outside valid range')
        if kind == 'inverse':
            free = b*vh
            volume = min(Fraction(free), Fraction(e)+Fraction(free,field_bits))
            if e < free:
                den = self.down.multiply(Decimal(field_bits), self.log_int(2)[0])
                entropy = self.up.divide(self.entropy_choose(free,e), den)
                volume = min(volume, Fraction(self.up.add(Decimal(e),entropy)))
            al, ah = max(0,b*(vl-ch)), b*(vh-cl)
            selected = min(ah,N//2) if al <= N//2 else N-al
            volume -= d*(cl-1)
            if sum_rank:
                if d != 1: raise ValueError('Sum-rank refinement requires running sum')
                r = cl-1
                volume = min(volume,Fraction(free-b*r),
                             Fraction(e+b*(vh-2*r),2)+Fraction(free,field_bits))
            coefficient = volume-min(al,N//2)
            # Point-density inflation for every noninitial live batch.
            coefficient += Fraction(density_bits*(vh-cl),field_bits)
            excess = retained.exact_excess(B,ch,selected,N,*marker)
            routed = self.routed(N//B,excess,tau,False)
        elif kind == 'forward':
            ul, uh = b*vl, b*vh
            selected = min(uh,N//2) if ul <= N//2 else N-ul
            coefficient = Fraction(max(0,uh-N//2)-b*(cl-1))
            coefficient += Fraction(density_bits*(ch-1),field_bits)
            excess = retained.exact_excess(B,ch,selected,N,*marker)
            routed = (self.projective_route(N//B,excess,tau) if uh <= N//2
                      else self.routed(N//B,excess,tau,False))
        else:
            raise ValueError(kind)
        if routed.is_infinite():
            return routed
        result = self.up.add(count,self.times_log(coefficient,tau))
        result = self.up.add(result,routed)
        return self.up.divide(result,self.log_int(2)[1 if result < 0 else 0])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--outer',type=int,choices=(256,512),required=True)
    ap.add_argument('--group',type=int,choices=(1,2,4,8,16,32,64),default=32)
    ap.add_argument('--distance-permille',type=int,default=450)
    ap.add_argument('--precision',type=int,choices=(80,120),required=True)
    ap.add_argument('--seconds',type=int,choices=range(1,121),default=60)
    ap.add_argument('--density-bits',type=int,default=0)
    ap.add_argument('--field-bits',type=int,choices=(16,32,64,128),default=128,
                    help='Base field for counting; extension to GF128 preserves distance.')
    ap.add_argument('--sum-rank',action='store_true',help='Retain active-coordinate reset rank; group1 only')
    args = ap.parse_args()
    if not 0 < args.distance_permille < 500 or args.density_bits < 0:
        ap.error('Invalid distance/density')
    N, B, d = 2**21,args.outer,args.group
    b = 512 if B == 256 else 256
    start = time.monotonic()
    if args.field_bits < 32 and d != 1:
        ap.error('Only running-sum reset is supported below GF32')
    if args.sum_rank and d != 1: ap.error('--sum-rank requires --group 1')
    result = cover(B,b,d,args.seconds,args.distance_permille,args.density_bits,args.field_bits,args.sum_rank)
    if result['status'] != 'FLOAT_COMPLETE_COVER_NOT_CERTIFICATE':
        raise ArithmeticError(result)
    # The floating deterministic-exclusion label is checked by integer
    # arithmetic during replay; canonicalize only for the coverage reader.
    boxes = [(vl,vh,cl,ch,'inverse' if kind=='outer_excludes' else kind,marker,val)
             for vl,vh,cl,ch,kind,marker,val in result['boxes']]
    deadline = time.monotonic()+args.seconds
    pairs,digest = partition(boxes,N//b,result['cmax'],deadline)
    controls(boxes,N//b,result['cmax'],deadline)
    print(dict(kind='exact_partition',pairs=pairs,sha256=digest,
               corrupted_controls=True),flush=True)
    bound = GroupedUpward(args.precision)
    worst = Decimal('-Infinity')
    e = args.distance_permille*N//1000
    with localcontext(bound.up):
        for index,item in enumerate(boxes):
            if time.monotonic()>deadline:
                raise RuntimeError('Replay cap; NO certificate')
            value = bound.rectangle_grouped(N,B,b,d,e,item,args.density_bits,args.field_bits,args.sum_rank)
            if value > -result['term_bits']:
                raise ArithmeticError((item,value))
            worst = max(worst,value)
            if (index+1)%100 == 0:
                print(dict(replayed=index+1,total=len(boxes)),flush=True)
    if time.monotonic()>deadline:
        raise RuntimeError('Final cap; NO certificate')
    print(dict(kind='directed_grouped_reset_certificate',outer=B,state=b,group=d,
               N=N,threshold=e,distance_permille=args.distance_permille,
               precision=args.precision,boxes=len(boxes),pairs=pairs,
               partition_sha256=digest,worst_log2_upper=str(worst),
               term_bits=result['term_bits'],margin_bits=40,
               conservative_integer_margin=int(worst.copy_negate())-(N//b*result['cmax']-1).bit_length(),
               density_bits=args.density_bits,field_bits=args.field_bits,sum_rank=args.sum_rank,
               seconds=round(time.monotonic()-start,3)),flush=True)


if __name__ == '__main__':
    main()
