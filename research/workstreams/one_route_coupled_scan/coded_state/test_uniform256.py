"""Exact finite tests of new proof interfaces. Not a K20 certificate."""
from collections import Counter
from fractions import Fraction
import itertools
import random
import time
from check_maps import Field
from ideal_reset_certificate import partition,controls


def main():
    started=time.monotonic()
    # Running-sum closed-segment enumerator, including the zero output.
    checks=0
    for q in (2,4):
        for length in range(2,7):
            counts=Counter()
            for word in itertools.product(range(q),repeat=length):
                total=0
                for x in word: total^=x
                if total==0: counts[sum(x!=0 for x in word)]+=1
            for z in (Fraction(1,8),Fraction(1,2),Fraction(3,4),Fraction(1)):
                actual=sum(number*z**w for w,number in counts.items())
                S=1+(q-1)*z
                exact=(S**length+(q-1)*(1-z)**length)/q
                upper=(1+(q-1)*z*z)*S**(length-2)
                assert actual==exact and actual<=upper
                checks+=1
    print(f'{checks} exact parity enumerators/weighted bounds PASS.',flush=True)

    # Distinct-time Vandermonde reset columns, also across batch boundaries.
    field=Field(4,0x13)
    points=field.points[:12]
    ranks=0
    for state in (2,4):
        cols=[]
        for x in points:
            column=[1]
            for _ in range(1,state): column.append(field.table[column[-1]][x])
            cols.append(column)
        for size in range(1,state+2):
            for chosen in itertools.combinations(range(12),size):
                matrix=[[cols[j][i] for j in chosen] for i in range(state)]
                assert field.rank(matrix)==min(size,state)
                ranks+=1
        # Identical-lane outputs at different times cancel a running sum,
        # but not this globally distinct-column syndrome.
        assert [a^b for a,b in zip(cols[0],cols[state])]!=[0]*state
    print(f'{ranks} exact reset-column ranks; repeated-lane negative control PASS.',flush=True)

    # Directly count the smallest MDS closed segment and compare support bound.
    mul=field.table
    counts=Counter()
    for word in itertools.product(range(16),repeat=4):
        s0=field.xor(word)
        s1=field.xor(mul[x][y] for x,y in zip(points[:4],word))
        if s0==s1==0: counts[sum(x!=0 for x in word)]+=1
    assert counts[1]==counts[2]==0 and counts[3]>0
    for weight in (3,4):
        from math import comb
        assert counts[weight]<=comb(4,weight)*16**(weight-2)
    print('Exhaustive 65536-output MDS reset count PASS.',flush=True)

    # Actual b=s interface: full feedback, invertible injection, inverse stream.
    gf4=Field(2,0x7)
    for state in itertools.product(range(4),repeat=2):
        if state==(0,0): continue
        counts=Counter()
        for tape in itertools.product(range(4),repeat=4):
            feedback=gf4.apply([list(tape[:2]),list(tape[2:])],state)
            counts[tuple(feedback)]+=1
        assert len(counts)==16 and set(counts.values())=={16}
    rng=random.Random(20261022)
    for _ in range(64):
        X=[[rng.randrange(16) for _ in range(2)] for _ in range(6)]
        T=[[[rng.randrange(16) for _ in range(2)] for _ in range(2)] for _ in range(6)]
        H=[[[1,1],points[2*t:2*t+2]] for t in range(6)]
        state=[0,0]; outputs=[]
        for t in range(6):
            y=[x^a for x,a in zip(X[t],field.apply(T[t],state))]
            state=[a^b for a,b in zip(state,field.apply(H[t],y))]
            outputs.append(y)
        state=[0,0]; recovered=[]
        for t,y in enumerate(outputs):
            recovered.append([a^b for a,b in zip(y,field.apply(T[t],state))])
            state=[a^b for a,b in zip(state,field.apply(H[t],y))]
        assert recovered==X
    print('15 exhaustive full-feedback laws and 64 complete inverse fixtures PASS.',flush=True)

    # Check every possible small c slice independently of the cover generator.
    for n in range(1,17):
        cmax=(n+1)//2
        boxes=[(2*c-1,n,c,c,'inverse',(5,4),0) for c in range(1,cmax+1)]
        assert partition(boxes,n,cmax,time.monotonic()+2)[0]==sum(n-2*c+2 for c in range(1,cmax+1))
        controls(boxes,n,cmax,time.monotonic()+2)
    assert time.monotonic()-started<15
    print('16 exact partitions and 48 corrupted-cover controls PASS.',flush=True)


if __name__=='__main__': main()
