"""Exact small-mask checks for mixed-history inequalities; no output files."""
from itertools import product
from math import comb
import random
import time
from layered_mixed_gate import choose_log2,occupation_dual


def pairs(width,bit):
    return [(x,x^(1<<bit)) for x in range(width) if not x&(1<<bit)]


def types(S,ps):
    return [(int(a in S),int(b in S)) for a,b in ps]


def tests():
    started=time.monotonic();width=4;defect_cases=orientation_cases=local_cases=0
    supports=[{i for i in range(width) if mask&(1<<i)} for mask in range(1<<width)]
    for masks in product(range(16),repeat=3):
        states=[supports[m] for m in masks]
        ty=[types(S,pairs(width,t%2)) for t,S in enumerate(states)]
        h=sum(x==y for row in ty for x,y in row)
        r=sum((x,y)==(0,1) for row in ty for x,y in row)
        V=sum(map(len,states));P=6;u=16
        U=len(states[0])+sum(len(states[t]|states[t+1]) for t in range(2))
        assert U>=3*u//4-width-min(r,P-r)-2*h
        assert U>=V>=P-h
        defect_cases+=1
        if any(pair==(0,0) for row in ty for pair in row):continue
        B=sum(len(states[t+1]-states[t]) for t in range(2))
        l=sum((x,y)==(1,0) for row in ty for x,y in row)
        holes=l+r;remaining=holes;weighted=0
        for weight in range(3):
            take=min(remaining,3*comb(2,weight));remaining-=take;weighted+=weight*take
        assert 2*min(l,r)>=weighted-(B+width//2)
        orientation_cases+=1

    for S in supports:
        for T in supports:
            old=pairs(width,0);new=pairs(width,1)
            h_old=sum(x==y for x,y in types(S,old))
            h_new=sum(x==y for x,y in types(T,new))
            full=[j for j,(a,b) in enumerate(old) if a in S and b in S]
            for selected in range(1<<len(full)):
                exceptional={full[j] for j in range(len(full)) if selected&(1<<j)}
                F=set()
                for j,(a,b) in enumerate(old):
                    if j in exceptional:continue
                    if b in S:F.update((a,b))
                    elif a in S:F.add(b)
                heads={b for a,b in old}
                charge=len((F&heads)-T)+len(exceptional)
                assert charge>=width//4-h_old-h_new
                local_cases+=1

    discrepancy=0
    for v in range(512):
        zeros=9-v.bit_count()
        for phase in range(9):
            for length in range(19):
                hits=sum(not (v&(1<<((phase+t)%9))) for t in range(length))
                assert abs(9*hits-zeros*length)<=20
                discrepancy+=1

    rng=random.Random(802);joint_cases=0
    for _ in range(4000):
        states=[];keys=[]
        for t in range(3):
            ps=pairs(4,t%2)
            choices=[rng.randrange(1,4) for p in ps]
            S={p[k] for p,c in zip(ps,choices) for k in range(2) if c&(1<<k)}
            ex={j for j,c in enumerate(choices) if c==3 and rng.randrange(2)}
            states.append(S);keys.append(ex)
        ty=[types(S,pairs(4,t%2)) for t,S in enumerate(states)]
        l=sum(pair==(1,0) for row in ty for pair in row)
        r=sum(pair==(0,1) for row in ty for pair in row)
        f=sum(pair==(1,1) for row in ty for pair in row);e=sum(map(len,keys))
        V=l+r+2*f;D=l+2*r+2*f-2*e
        assert V-D-e==e-r
        B=sum(len(states[t+1]-states[t]) for t in range(2))
        W=len(states[0]);C=0;charge=e
        for t in range(2):
            intersection=states[t]&states[t+1]
            zeros={x for x in intersection if rng.randrange(2)}
            Y=(states[t]|states[t+1])-zeros;C+=len(zeros);W+=len(Y)
            F=set()
            for j,(a,b) in enumerate(pairs(4,t%2)):
                if j in keys[t]:continue
                if b in states[t]:F.update((a,b))
                elif a in states[t]:F.add(b)
            charge+=len(F-(Y&states[t+1]))
        assert C>=V+B-2-W
        assert 2*B<=W
        assert charge>=l+r+e-B-2
        assert charge>=D+e+B-W-4
        joint_cases+=1
    binomial_cases=0
    for n in range(1,100):
        for k in range(n+1):
            assert comb(n,k)<=1<<choose_log2(n,k)
            binomial_cases+=1
    dual_cases=0
    for T in range(1,7):
        for a,h,m in ((4,1,2),(4,3,5),(8,4,1)):
            exact=0
            for v in range(8):
                exact+=min(sum((a*int(not(v&(1<<(t%3))))-h)*z[t] for t in range(T))
                           +m*sum(z[t] and not z[t+1] for t in range(T-1))
                           for z in product((0,1),repeat=T))
            assert occupation_dual(T,a,h,m,bits=3)==exact
            dual_cases+=1
    ordered_cases=0
    from itertools import combinations
    halfsets=list(combinations(range(4),2))
    for sets in product(halfsets,repeat=3):
        degrees=[sum(v in S for S in sets) for v in range(4)]
        for s in range(5):
            actual=sum(sorted(degrees,reverse=True)[:s])
            assert actual==min(k*s+sum(max(0,d-k) for d in degrees) for k in range(4))
            ordered_cases+=1
    swap_cases=0
    for x,y,k in product(range(8),range(1,9),range(9)):
        change=max(0,x+1-k)+max(0,y-1-k)-max(0,x-k)-max(0,y-k)
        assert abs(change)<=1;swap_cases+=1
    assert time.monotonic()-started<15
    return dict(defect_support_checks=defect_cases,small_cycle_orientation_checks=orientation_cases,
                local_defect_charge_checks=local_cases,nine_bit_interval_checks=discrepancy,
                joint_birth_rank_root_checks=joint_cases,binomial_upper_checks=binomial_cases,
                occupation_dual_checks=dual_cases,ordered_head_identity_checks=ordered_cases,
                head_swap_range_checks=swap_cases)


if __name__=='__main__':print(tests())
