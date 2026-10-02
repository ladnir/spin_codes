"""Exact vertex clipping in fixed dimension; no float feasibility decisions."""
from fractions import Fraction as Q
from itertools import product


def dot(a,b):return sum(x*y for x,y in zip(a,b))


def rank(rows):
    if not rows:return 0
    rows=[list(map(Q,row)) for row in rows];r=0
    for j in range(len(rows[0])):
        pivot=next((i for i in range(r,len(rows)) if rows[i][j]),None)
        if pivot is None:continue
        rows[r],rows[pivot]=rows[pivot],rows[r]
        divisor=rows[r][j];rows[r]=[x/divisor for x in rows[r]]
        for i in range(r+1,len(rows)):
            factor=rows[i][j]
            if factor:rows[i]=[x-factor*y for x,y in zip(rows[i],rows[r])]
        r+=1
        if r==len(rows):break
    return r


def vertices(cell,planes):
    """Clip all box edges by each halfspace, retaining exact face incidences.

    Two distinct vertices are adjacent exactly when their common active
    faces have rank dimension-1. Degenerate vertices retain every active
    face, including constraints that did not remove any previous vertex.
    """
    d=len(cell)//2
    if len(cell)!=2*d or not d or any(lo>=hi for lo,hi in zip(cell[::2],cell[1::2])):
        raise ValueError('positive-width box required')
    normals=[]
    for j in range(d):
        for sign in (-1,1):normals.append(tuple(Q(sign*int(j==k)) for k in range(d)))
    current={tuple(p):frozenset(2*j+int(x==cell[2*j+1]) for j,x in enumerate(p))
             for p in product(*zip(cell[::2],cell[1::2]))}
    for normal,bound in planes:
        normal=tuple(map(Q,normal));bound=Q(bound)
        if len(normal)!=d:raise ValueError('halfspace dimension mismatch')
        if sum(max(a*lo,a*hi) for a,lo,hi in zip(normal,cell[::2],cell[1::2]))<=bound:continue
        index=len(normals);normals.append(normal)
        values={p:dot(normal,p)-bound for p in current}
        inside=[p for p,v in values.items() if v<=0];outside=[p for p,v in values.items() if v>0]
        if not inside:return ()
        next_={p:current[p]|({index} if values[p]==0 else set()) for p in inside}
        edge_cache={}
        for p in inside:
            if values[p]==0:continue
            for q in outside:
                common=current[p]&current[q]
                if len(common)<d-1:continue
                if common not in edge_cache:edge_cache[common]=rank([normals[j] for j in common])>=d-1
                if not edge_cache[common]:continue
                point=tuple((values[q]*x-values[p]*y)/(values[q]-values[p]) for x,y in zip(p,q))
                faces=common|{index}
                next_[point]=next_.get(point,frozenset())|faces
        current={p:frozenset(f) for p,f in next_.items()}
    return tuple(sorted(current))
