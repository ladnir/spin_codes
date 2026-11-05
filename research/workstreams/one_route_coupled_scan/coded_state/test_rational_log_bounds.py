"""Independent exact-series, known-power, and high-precision smoke checks."""
from fractions import Fraction as F
from decimal import Decimal,localcontext
from random import Random
from rational_log_bounds import log2_bounds,ln_normalized,Q,TERMS


def run():
    rng=Random(905);powers=0;series=0;checks=0
    for exponent in range(-256,257):
        v=F(2)**exponent;lo,hi=log2_bounds(v)
        assert lo==hi==exponent*(1<<24);powers+=1
    for n,d in ((2,1),(3,2),(7,4),(17,16),(101,100)):
        z=F(n-d,n+d)
        exact=2*sum((z**(2*j+1)/(2*j+1) for j in range(TERMS)),F(0))
        tail=2*z**(2*TERMS+1)/((2*TERMS+1)*(1-z*z))
        lo,hi=ln_normalized(n,d)
        assert F(lo,Q)<=exact and F(hi,Q)>=exact+tail;series+=1
    with localcontext() as ctx:
        ctx.prec=75;ln2=Decimal(2).ln()
        for _ in range(256):
            n=rng.getrandbits(rng.randrange(1,1025))+1;d=rng.getrandbits(rng.randrange(1,1025))+1
            v=F(n,d);lo,hi=log2_bounds(v)
            reference=(Decimal(n).ln()-Decimal(d).ln())/ln2*(1<<24)
            assert Decimal(lo)<=reference<=Decimal(hi)
            assert hi-lo<=2;checks+=1
            il,iu=log2_bounds(1/v)
            assert lo+il<=0<=hi+iu
    print(dict(exact_powers=powers,exact_series_enclosures=series,decimal_smoke_checks=checks))


if __name__=='__main__':run()
