"""Bounded proof diagnostics for the actual layered inner's one-sided states.

All algebra and directed integer bounds are separate from native performance.
"""
from math import comb, isqrt
from fractions import Fraction
import argparse
import random
from itertools import product
import time
import numpy as np
from layered_route_gate import graph_envelope
from test_layered_inner import field


def head_envelope(M,R):
    choose = 1
    bits = []
    for s in range(M+1):
        if s:
            choose = choose*(M-s+1)//s
        bits.append(choose.bit_length())
    result = []
    for s in range(M+1):
        n = R*s
        # Per-subset tail <2^(-bitlength(C(M,s))-77), then size/suffix union.
        numerator = 7*n*(bits[s]+77)
        root = isqrt(numerator//20)
        root += 20*root*root < numerator
        result.append(min(n, (n+1)//2+root, R*M//2))
    return result


def dimension(M,R,P,r,heads,edges):
    s = np.arange(M+1,dtype=np.int64)
    # F consists of all heads in the live batches, plus r selected tails.
    values = np.minimum.reduce((np.full(M+1,P+r,dtype=np.int64),
                                P+edges, heads+r, R*s))-128*s
    return min((P+r)//2,max(0,R*M-2**20),max(0,int(values.max())))


def point(R,r,permille=400):
    M,b,N = 8192,512,2**21
    u=R*M;P=(u-b)//2;W=permille*N//1000
    heads=np.asarray(head_envelope(M,R),dtype=np.int64)
    edges=np.asarray(graph_envelope(M,R),dtype=np.int64)
    d=dimension(M,R,P,r,heads,edges)
    rank=max(0,P-W,3*u//4-b-min(r,P-r)-W)
    inverse=2*P+32*(b-r-rank+d)+1
    forward=P+32*(max(0,u-N//2)-u//4+b//2)+1
    return dict(R=R,r=r,dimension=d,zero_output_rank=rank,
                inverse_log2_upper=inverse,forward_log2_upper=forward)


def mask_tests():
    cases=0
    width=8
    for i in range(3):
        j=(i+1)%3
        pairs_i=[(x,x^(1<<i)) for x in range(width) if not x&(1<<i)]
        pairs_j=[(x,x^(1<<j)) for x in range(width) if not x&(1<<j)]
        for a in range(1<<(width//2)):
            S={p[(a>>k)&1] for k,p in enumerate(pairs_i)}
            r=a.bit_count()
            for c in range(1<<(width//2)):
                T={p[(c>>k)&1] for k,p in enumerate(pairs_j)}
                rnext=c.bit_count()
                assert len({x for x in range(width) if x&(1<<i)}-T)==width//4
                assert len(S & T)<=width//4+min(r,rnext,width//2-r,width//2-rnext)
                cases+=1
    return cases


def exact_tests():
    moments=0
    for M in range(2,15,2):
        for s in range(M+1):
            for j in range(s+1):
                actual=Fraction(comb(M-j,M//2-j),comb(M,M//2)) if j<=M//2 else 0
                assert actual<=Fraction(1,2**j)
                moments+=1
    rng=random.Random(72)
    graphs=0
    # Exhaust all tail selections for several small, oriented route fixtures.
    for _ in range(12):
        M,R,k=6,3,2
        edges=[]
        for region in range(R):
            order=list(range(M));rng.shuffle(order)
            edges+=list(zip(order[::2],order[1::2]))
        H=[0]*(M+1);E=[0]*(M+1)
        for mask in range(1<<M):
            s=mask.bit_count()
            H[s]=max(H[s],sum(bool(mask&(1<<b)) for a,b in edges))
            E[s]=max(E[s],sum(bool(mask&(1<<a)) and bool(mask&(1<<b)) for a,b in edges))
        P=len(edges)
        for mask in range(1<<P):
            r=mask.bit_count();degrees=[0]*M
            for j,(a,b) in enumerate(edges):
                degrees[b]+=1
                degrees[a]+=bool(mask&(1<<j))
            actual=sum(max(0,d-k) for d in degrees)
            cap=max(0,max(min(P+r,P+E[s],H[s]+r,R*s)-k*s for s in range(M+1)))
            assert actual<=cap
            graphs+=1
    return dict(head_factorial_checks=moments,oriented_graph_checks=graphs,
                consecutive_wiring_checks=mask_tests(),
                temporal_rank_checks=temporal_tests())


def temporal_tests():
    table=field(2,7)
    def rank(rows):
        rows=[row[:] for row in rows];pivot=0
        for col in range(6):
            found=next((j for j in range(pivot,len(rows)) if rows[j][col]),None)
            if found is None:continue
            rows[pivot],rows[found]=rows[found],rows[pivot]
            inv=next(x for x in range(1,4) if table[x][rows[pivot][col]]==1)
            rows[pivot]=[table[inv][v] for v in rows[pivot]]
            for j in range(pivot+1,len(rows)):
                scalar=rows[j][col]
                rows[j]=[v^table[scalar][w] for v,w in zip(rows[j],rows[pivot])]
            pivot+=1
        return pivot
    cases=0
    for orientations in product(range(4),repeat=3):
        supports=[]
        for t,mask in enumerate(orientations):
            bit=t%2
            pairs=[(x,x^(1<<bit)) for x in range(4) if not x&(1<<bit)]
            supports.append({p[(mask>>j)&1] for j,p in enumerate(pairs)})
        r=sum(mask.bit_count() for mask in orientations)
        edges=[((t,x),(t+1,x)) for t in range(2) for x in supports[t]&supports[t+1]]
        assert len(edges)<=2+min(r,6-r)
        potential=len(supports[0])+sum(len(supports[t]|supports[t+1]) for t in range(2))
        assert potential>=8-min(r,6-r)
        indices={v:j for j,v in enumerate((t,x) for t in range(3) for x in sorted(supports[t]))}
        for selected in range(1<<len(edges)):
            rows=[]
            for j,(a,b) in enumerate(edges):
                if selected&(1<<j):
                    row=[0]*6;row[indices[a]]=2;row[indices[b]]=1
                    rows.append(row)
            assert rank(rows)==selected.bit_count()
            cases+=1
    return cases


def cover(permille=400,seconds=50,defects=0):
    """Cover every integer (R,r), using monotonicity, not sampled points."""
    started=time.monotonic()
    M,b,N=8192,512,2**21
    W=permille*N//1000
    classes=intervals=forward_regions=0
    worst=None
    for R in range(129,257):
        u=R*M;P=(u-b)//2
        assert 0<=defects<=P
        # At most H defective (00 or 11) pairs: (P+1)^H location encodings,
        # two support choices per pair, and at most 2^H exceptional branches.
        extra=defects*P.bit_length()
        classes+=P+1
        forward=P+extra+defects+32*(max(0,u-N//2)-u//4+b//2+2*defects)+1
        if forward<=-80:
            forward_regions+=1
            continue
        heads=np.asarray(head_envelope(M,R),dtype=np.int64)
        edges=np.asarray(graph_envelope(M,R),dtype=np.int64)
        cached={}
        def dim(r):
            if r not in cached:
                cached[r]=dimension(M,R,P,r,heads,edges)
            return cached[r]
        stack=[(0,P)];accepted=[]
        while stack:
            if time.monotonic()-started>seconds:
                raise RuntimeError('Coverage time cap; no partial coverage claim')
            lo,hi=stack.pop()
            # d(r) is nondecreasing. min(r,P-r) attains its maximum nearest P/2.
            center=max(lo,min(hi,P//2))
            rank_min=max(0,P-defects-W,
                         3*u//4-b-min(center,P-center)-2*defects-W)
            upper=2*P+extra+2*defects+32*(b+defects-lo-rank_min+dim(min(P,hi+defects)))+1
            if upper<=-80:
                accepted.append((lo,hi));intervals+=1
                if worst is None or upper>worst[0]:worst=(upper,R,lo,hi)
            elif lo==hi:
                return dict(complete=False,permille=permille,defect_budget=defects,failed_class=(R,lo),
                            inverse_log2_upper=upper,forward_log2_upper=forward)
            else:
                mid=(lo+hi)//2
                stack.extend(((mid+1,hi),(lo,mid)))
        accepted.sort()
        cursor=0
        for lo,hi in accepted:
            assert lo==cursor and hi>=lo
            cursor=hi+1
        assert cursor==P+1
    assert classes<2**27
    return dict(complete=True,permille=permille,defect_budget=defects,weight_threshold=W,
                classes=classes,forward_regions=forward_regions,
                inverse_intervals=intervals,worst_interval=worst,
                combined_failure_below_power_of_two=-52,
                whole_code_certificate=False)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--permille',type=int,default=400)
    parser.add_argument('--defects',type=int,default=0)
    args=parser.parse_args()
    print(exact_tests(),flush=True)
    print(cover(args.permille,defects=args.defects),flush=True)


if __name__=='__main__':main()
