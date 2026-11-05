"""Exact rational first moments of the balanced auxiliary hole process.

P=[x^256]K_x / binom(512,256). Normalize each row by P1 before
reporting its next hole population and persistent overlap. This is a
moment-weighted auxiliary process, not an encoder state trajectory.
Actual routing fixes a quota only over sixteen updates: it does NOT
force each update to be balanced. Therefore these means do not upper
bound the regional operator or prove a distance certificate.

Fixed ten-point diagnostic, no search, no files, parent-only execution.
The exact mass and means use Fraction; decimals are display values.
"""

from fractions import Fraction
from math import comb
import os
import shutil
import threading


def quota_drift_main():
    from functools import lru_cache
    import resource,math
    if shutil.disk_usage('/tmp').free <20*(1<<30): raise RuntimeError('20GiB guard')
    resource.setrlimit(resource.RLIMIT_AS,(1<<30,1<<30))
    timer=threading.Timer(55,lambda:os._exit(124));timer.daemon=True;timer.start()
    n=512; quota=256;m=128
    @lru_cache(None)
    def Q(w,n,h):
        return sum(comb(w,r)*comb(n-w,h-r)*2**(w-r)
                   for r in range(max(0,h-(n-w)),min(w,h)+1))
    try:
        for a,w in ((64,8),(96,24),(128,32),(128,48),(160,40),(160,64),(192,64),(192,96),(256,128),(512,512)):
            tot=Q(w,n,quota); mass=Fraction(0);mg=Fraction(0);mug=Fraction(0)
            for o in range(5):
              for e in range(5-o):
                z=4-o-e;u=o+e
                if o>w or e>a-w or z>n-a:continue
                typep=Fraction(comb(w,o)*comb(a-w,e)*comb(n-a,z),comb(n,4))
                for r in range(o+1):
                  for s in range(e+1):
                    for t in range(z+1):
                        j=r+s+t
                        lp=Fraction(comb(o,r)*comb(e,s)*comb(z,t)*2**(o-r)*Q(w-o,n-4,quota-j),tot)
                        B=2*j-r;loss=max(0,4+u-B)
                        g=Fraction(1) if j==0 else Fraction(2**loss-1,2**(4+u)-1)
                        p=typep*lp
                        mass+=p;mg+=p*g;mug+=p*g*u
            assert mass==1
            print(dict(state=(a,w),log2_exact_quota_row_mass=math.log2(tot)-math.log2(comb(n,quota)),mean_next_a=float(4*m*mg),mean_next_w=float(m*mug),mean_exact_rational=True,full_operator_bound=False),flush=True)
    finally:timer.cancel()


if __name__ == "__main__":
    quota_drift_main()

