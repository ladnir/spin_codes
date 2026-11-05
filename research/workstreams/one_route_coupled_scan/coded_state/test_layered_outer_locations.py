"""Exact location/outer compatibility checks, including negative controls."""
from fractions import Fraction as F
from itertools import combinations,product
from math import comb
import random
import time
from layered_fixed_rank import matrix_rank,path_nullity
from test_layered_inner import field,setup,pair
from layered_outer_locations import positive_rows,rows,start,TERMS


def dot(a,c,mul):
    out=0
    for x,y in zip(a,c):out^=mul[x][y]
    return out


def nullspace(A,n,mul):
    A=[r[:] for r in A];pivots=[];k=0
    for j in range(n):
        p=next((i for i in range(k,len(A)) if A[i][j]),None)
        if p is None:continue
        A[k],A[p]=A[p],A[k]
        inv=next(x for x in range(1,len(mul)) if mul[x][A[k][j]]==1)
        A[k]=[mul[inv][x] for x in A[k]]
        for i in range(len(A)):
            if i!=k:
                c=A[i][j];A[i]=[x^mul[c][y] for x,y in zip(A[i],A[k])]
        pivots.append(j);k+=1
    out=[]
    for j in range(n):
        if j in pivots:continue
        v=[0]*n;v[j]=1
        for i,p in enumerate(pivots):v[p]=A[i][j]
        assert all(dot(r,v,mul)==0 for r in A)
        out.append(v)
    return out


def fixtures():
    mul=field(2,7);keys,_=setup(mul);rng=random.Random(803)
    width,L,M,K=4,4,4,8;n=width*L;gains=[];zeros_counts=[];split_edges=0
    for trial in range(128):
        G=[[0]*K for _ in range(n)]
        for region in range(4):
            order=list(range(M));rng.shuffle(order)
            for j,block in enumerate(order):
                G[region*M+j][2*block]=1;G[region*M+j][2*block+1]=region
        Y=[0]*n
        for j in rng.sample(range(n),1+trial%9):Y[j]=rng.randrange(1,4)
        Q=[[0]*width]
        for t in range(L):Q.append([mul[2][Q[-1][j]]^Y[t*width+j] for j in range(width)])
        inactive={t*width+j for t in range(L-1) for j in range(width) if not Q[t+1][j]}
        ze={t*width+j for t in range(1,L-1) for j in range(width)
            if not Y[t*width+j] and Q[t][j] and Q[t+1][j]}
        fixed={j:0 for j in range(width)}
        for t in range(1,L):
            for j in range(width):
                if j>>(t%2)&1:continue
                k=j^(1<<(t%2));x,y=Q[t][j],Q[t][k]
                p,_=pair(mul,rng.choice(keys),x,y)
                if not x and not y:fixed[t*width+j]=fixed[t*width+k]=0
                elif x and not y:fixed[t*width+j]=1
                elif not p:fixed[t*width+j]=0;fixed[t*width+k]=1
        T=[]
        for j in inactive:
            row=[0]*n;row[j]=1;T.append(row)
        for j in ze:
            row=[0]*n;row[j]=1;row[j-width]=2;T.append(row)
        A=[];Gfix=[]
        for j,g in fixed.items():
            row=[0]*n;row[j]=1
            if j>=width:row[j-width]=g
            A.append(row);Gfix.append(G[j])
        extra=[]
        for (pos,g),row in zip(fixed.items(),A):
            t,j=divmod(pos,width)
            if 1<=t<L-1 and g==1 and Q[t][j] and Q[t+1][j] and Y[pos] and Q[t][j]==Q[t+1][j]:
                extra.append(row)
        assert matrix_rank(T+extra,n,mul)-matrix_rank(T,n,mul)==len(extra)
        assert matrix_rank(T+extra+A,n,mul)==matrix_rank(T+A,n,mul)
        split_edges+=len(extra);T+=extra
        basis=nullspace(T,n,mul)
        image=[[dot(row,v,mul) for v in basis] for row in A]
        relations=nullspace([list(row) for row in zip(*image)],len(A),mul)
        quotient=[[dot(v,col,mul) for col in zip(*Gfix)] for v in relations]
        compatible=K-matrix_rank(quotient,K,mul)
        kappa=path_nullity(width,L,inactive,ze,fixed)
        joint=[r+[0]*K for r in T]+[r+g for r,g in zip(A,Gfix)]
        exact=n+K-matrix_rank(joint,n+K,mul)
        assert exact==kappa+compatible
        zero_rows=[g for r,g in zip(image,Gfix) if not any(r)]
        zero_rank=matrix_rank(zero_rows,K,mul)
        assert compatible<=K-zero_rank
        gains.append(K-zero_rank-compatible);zeros_counts.append(len(zero_rows))
    assert max(gains)>0
    # Same number of zeros, different locations: RS[4,2] in four blocks.
    # Four zeros in one row remove 2 dimensions; one per row removes 4.
    assert sum(min(j,2) for j in (4,0,0,0))==2
    assert sum(min(j,2) for j in (1,1,1,1))==4
    return dict(fixtures=len(gains),joint_dimension_identity=True,
                extra_relation_rank_cases=sum(g>0 for g in gains),
                max_rank_beyond_forced_zeros=max(gains),max_forced_zeros=max(zeros_counts),
                independently_checked_split_edges=split_edges)


def colour_checks():
    checks=0;M=4;zero_factor=F(1,4)
    for s in range(1,M+1):
        sets=[set(v) for v in combinations(range(M),s)]
        for pattern in product(range(3),repeat=M):
            J={i for i,v in enumerate(pattern) if v==1};U={i for i,v in enumerate(pattern) if v==2}
            actual=sum((zero_factor**len(J&S) for S in sets if U<=S),F(0))/len(sets)
            for p in (F(1,4),F(1,2),F(3,4)):
                conditioning=comb(M,s)*p**s*(1-p)**(M-s)
                bound=p**len(U)*(1-p+p*zero_factor)**len(J)/conditioning
                assert actual<=bound
                checks+=1
    return checks


def split_checks():
    # A gamma=1 edge at nonzero output has X=0 iff Qnext=Q.
    # This relation lowers actual state dimension by one, but is already
    # present in the auxiliary kernel, so it does not lower kappa again.
    mul=field(2,7);count=0
    for a,c in product(range(1,4),repeat=2):
        y=c^mul[2][a];x=c^a
        if not y:assert x
        else:assert (not x)==(a==c)
        count+=1
    extra=[[1,1]];fixed=[[1,1]]
    assert matrix_rank(extra,2,mul)==1
    assert matrix_rank(fixed+extra,2,mul)==matrix_rank(fixed,2,mul)
    return count


def actual_local_checks():
    mul=field(2,7);keys,_=setup(mul);checked=0
    for x,y,xx,yy in product(range(4),repeat=4):
        for key in keys:
            p,q=pair(mul,key,x,y)
            gammas=((0,0) if not (x or y) else ((1,-1) if x and not y else
                    ((0,1) if x and y and not p else (-1,-1))))
            for a,c,g,v in zip((x,y),(xx,yy),gammas,(p^xx,q^yy)):
                if g<0:continue
                assert v==(c^mul[g][a])
                output=c^mul[2][a]
                if g==0:assert bool(v)==bool(c)
                elif bool(a)!=bool(c) or (a and c and not output):assert v
                elif a and c:assert (not v)==(a==c)
                else:assert not v
                checked+=1
    return checked


def positive_checks():
    import numpy as np
    # Local float path is only a proposer. Check it against exact fractions
    # with nontrivial small weights, both matching depths and the zero split.
    markers=[F(2,3),F(3,4),F(2),F(1,3),F(2,3),F(1,2)]+list(map(F,(2,3,3,4,5,6,7,8)))
    checks=0
    for s in (1024,4096,8192):
        for two in (False,True):
            _,rr,ee,_=positive_rows(markers,29,F(s,8192),two=two)
            ff,gg=rows(start(markers),s,29,two=two)
            assert np.max(np.abs(ff-np.log2(np.array(rr,dtype=float))))<1e-11
            assert np.max(np.abs(gg-np.log2(np.array(ee,dtype=float))))<1e-11
            checks+=ff.size+gg.size
    assert all(len(t)==12 for t in TERMS)
    return checks


if __name__=='__main__':
    began=time.monotonic()
    print(fixtures())
    print(dict(exact_colour_bounds=colour_checks(),split_cases=split_checks(),
               actual_local_classifications=actual_local_checks(),
               exact_positive_rows=positive_checks(),whole_code_certificate=False))
    assert time.monotonic()-began<30
