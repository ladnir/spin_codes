"""Exact toy checks for an idealized pair-randomness comparison, not a cert.

Uniform SL(2,q) maps each nonzero pair uniformly onto nonzero pairs. Rejection
of non-superregular A and A+delta I preserves this up to a point-density bound.
No encoder changes. All setup is fixed for all messages.
"""
from collections import Counter
from fractions import Fraction as F
import time


def field(q,poly):
    def mul(a,b):
        out=0
        while b:
            if b&1:out^=a
            b>>=1;a<<=1
            if a&q:a^=poly
        return out
    table=[[mul(a,b) for b in range(q)] for a in range(q)]
    return table


def run():
    start=time.monotonic();receipts=[]
    for q,poly in ((4,7),(8,11),(16,19)):
        t=field(q,poly);delta=2
        mats=[(a,b,c,d) for a in range(q) for b in range(q)
              for c in range(q) for d in range(q) if t[a][d]^t[b][c]==1]
        good=[A for A in mats if all(A) and A[0]!=delta and A[3]!=delta
              and (t[A[0]^delta][A[3]^delta]^t[A[1]][A[2]])!=0]
        assert len(mats)==q*(q*q-1)
        assert len(good)>=q*(q*q-7*q+2)
        tested=0;max_fiber=0
        for x in range(q):
            for y in range(q):
                if not (x or y):continue
                def image(A):
                    a,b,c,d=A
                    return (t[a][x]^t[b][y],t[c][x]^t[d][y])
                census=Counter(map(image,mats));assert len(census)==q*q-1 and set(census.values())=={q}
                restricted=Counter(map(image,good));assert max(restricted.values(),default=0)<=q
                max_fiber=max(max_fiber,max(restricted.values(),default=0));tested+=1
                if time.monotonic()-start>30:raise TimeoutError('30s toy-check cap')
        # Every retained matrix has a nonzero upper-right pivot. Three shears
        # L(v) U(b) L(w) reconstruct it; all divisions are setup work.
        for a,b,c,d in good:
            inv=next(i for i in range(1,q) if t[b][i]==1)
            v=t[d^1][inv];w=t[a^1][inv]
            assert (1^t[b][w],b,v^t[1^t[v][b]][w],1^t[v][b])==(a,b,c,d)
        receipts.append(dict(q=q,matrices=len(mats),retained=len(good),inputs=tested,
                             restricted_max_fiber=max_fiber,three_shear_checks=len(good)))
    q=1<<32;N=1<<21
    # Product of point-density losses is bounded by Bernoulli's inequality.
    inflation=1/(1-F(7*N,2*q));assert inflation<2
    print(dict(exact_local_checks=receipts,inflation_upper=str(inflation),
               seconds=time.monotonic()-start,whole_code_certificate=False))


if __name__=='__main__':run()
