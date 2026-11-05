"""Bounded regression checks for the directed grouped wrapper, not a certificate."""
from decimal import localcontext
from fractions import Fraction
import time
from grouped_reset_certificate import GroupedUpward,retained
from ideal_reset_certificate import partition,controls


def main():
    start=time.monotonic()
    new=GroupedUpward(80)
    old=retained.PersistentUpward(80)
    checked=0
    with localcontext(new.up):
        for v,c in ((2200,64),(3200,128),(4000,1000)):
            for kind in ('forward','inverse'):
                item=(v,v,c,c,kind,(5,4),0)
                exact=new.rectangle_grouped(2**21,256,512,512,943718,item)
                assert exact == old.rectangle(2**21,256,512,128,943718,item)
                weak=new.rectangle_grouped(2**21,256,512,32,943718,item)
                assert weak >= exact
                inflated=new.rectangle_grouped(2**21,256,512,32,943718,item,1)
                assert inflated >= weak
                lower_distance=new.rectangle_grouped(2**21,256,512,32,900000,item)
                assert lower_distance <= weak
                for bits in (32,64,128):
                    generic=new.rectangle_grouped(2**21,256,512,512,943718,item,field_bits=bits)
                    assert generic == old.rectangle(2**21,256,512,bits,943718,item)
                    baseline=new.rectangle_grouped(2**21,256,512,1,943718,item,field_bits=bits)
                    stronger=new.rectangle_grouped(2**21,256,512,1,943718,item,field_bits=bits,sum_rank=True)
                    assert stronger <= baseline
                checked+=1
    for n in range(1,17):
        cmax=(n+1)//2
        boxes=[(2*c-1,n,c,c,'inverse',(5,4),0) for c in range(1,cmax+1)]
        total,_=partition(boxes,n,cmax,start+15)
        assert total == sum(n-2*c+2 for c in range(1,cmax+1))
        controls(boxes,n,cmax,start+15)
    # The generic group1 reset is exactly a running sum. Every nonempty
    # closed support has at least two elements and at least one equation.
    from itertools import product
    patterns=0
    for word in product(range(4),repeat=4):
        total=0
        for value in word: total^=value
        if total==0 and any(word):
            assert sum(value!=0 for value in word)>=2
            patterns+=1
    assert time.monotonic()-start<15
    print(dict(kind='directed_wrapper_regression',rectangle_cases=checked,
               partitions=16,corruption_controls=48,nonzero_closed_patterns=patterns))


if __name__ == '__main__':
    main()
