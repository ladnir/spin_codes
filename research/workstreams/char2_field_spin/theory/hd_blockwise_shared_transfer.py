"""Whole-joint-image rank law for a fresh invertible rank-one mixer.

All geometry is conditioned before the current mixer. The operator counts a
return event once per joint subspace, not once per field-valued source word.
No whole-code distance certificate or native timing is claimed.
"""

from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations,product

from hd_blockwise_maps import mul as mul8


def multiply(a,b,q):
    if q==256:
        return mul8(a,b)
    if q==65536:
        # AES GF256[u]/(u^2+u+0x20), matching WideDiagonal16.h.
        a0,a1=a&255,a>>8
        b0,b1=b&255,b>>8
        high=mul8(a1,b1)
        return ((mul8(a0,b0)^mul8(0x20,high))
                |((mul8(a0,b1)^mul8(a1,b0)^high)<<8))
    modulus={2:3,4:7,8:11,16:19}[q]
    result=0
    while b:
        if b&1:
            result^=a
        b>>=1
        a<<=1
        if a&q:
            a^=modulus
    return result


def power(a,n,q):
    result=1
    while n:
        if n&1:
            result=multiply(result,a,q)
        a=multiply(a,a,q)
        n>>=1
    return result


def dot(a,b,q):
    result=0
    for x,y in zip(a,b):
        result^=multiply(x,y,q)
    return result


def rref(matrix,q):
    if not matrix:
        return (),()
    work=[list(row) for row in matrix]
    if len({len(row) for row in work})!=1:
        raise ValueError('ragged matrix')
    pivots=[]
    for column in range(len(work[0])):
        pivot=next((i for i in range(len(pivots),len(work)) if work[i][column]),None)
        if pivot is None:
            continue
        row=len(pivots)
        work[row],work[pivot]=work[pivot],work[row]
        inverse=power(work[row][column],q-2,q)
        work[row]=[multiply(inverse,x,q) for x in work[row]]
        for i in range(len(work)):
            scalar=work[i][column]
            if i!=row and scalar:
                work[i]=[x^multiply(scalar,y,q) for x,y in zip(work[i],work[row])]
        pivots.append(column)
        if len(pivots)==len(work):
            break
    return tuple(tuple(row) for row in work),tuple(pivots)


def rank(matrix,q):
    return len(rref(matrix,q)[1])


def kernel(matrix,columns,q):
    if not matrix:
        return tuple(tuple(int(i==j) for j in range(columns)) for i in range(columns))
    reduced,pivots=rref(matrix,q)
    free=[j for j in range(columns) if j not in pivots]
    answer=[]
    for j in free:
        vector=[0]*columns
        vector[j]=1
        for i,pivot in enumerate(pivots):
            vector[pivot]=reduced[i][j]
        answer.append(vector)
    return tuple(tuple(v[i] for v in answer) for i in range(columns))


def mm(a,b,q):
    if not a:
        return ()
    if not b or not b[0]:
        return tuple(() for _ in a)
    return tuple(tuple(dot(row,column,q) for column in zip(*b)) for row in a)


def columns_basis(matrix,q):
    pivots=rref(matrix,q)[1]
    return tuple(tuple(row[j] for j in pivots) for row in matrix)


def geometry(joint,q=65536):
    """joint has rows (Q,C), with arbitrary source columns; quotient duplicates."""
    joint=tuple(tuple(row) for row in joint)
    if not joint or len(joint)%2 or len({len(row) for row in joint})!=1:
        raise ValueError('joint state/feedback matrix must have 2s rows')
    if any(not 0<=x<q for row in joint for x in row):
        raise ValueError('coefficient outside the field')
    joint=columns_basis(joint,q)
    s=len(joint)//2
    d=len(joint[0])
    state,feedback=joint[:s],joint[s:]
    difference=tuple(tuple(a^b for a,b in zip(x,y)) for x,y in zip(state,feedback))
    r=rank(difference,q)
    diagonal=columns_basis(mm(state,kernel(difference,d,q),q),q)
    e=len(diagonal[0])
    assert d==r+e
    k=e+r-rank(tuple(x+y for x,y in zip(diagonal,difference)),q)
    a0=rank(state,q)-e
    a1=rank(feedback,q)-e
    # Quotient both state and feedback by E, without choosing a quotient basis.
    augmented=(tuple(joint[i]+diagonal[i]+(0,)*e for i in range(s))
               +tuple(joint[s+i]+(0,)*e+diagonal[i] for i in range(s)))
    a01=rank(augmented,q)-2*e
    assert 0<=k<=min(e,r) and max(a0,a1)<=a01<=r
    return dict(field=q,state=s,joint_dimension=d,difference_rank=r,
                diagonal_dimension=e,diagonal_difference_intersection=k,
                quotient_state_rank=a0,quotient_feedback_rank=a1,
                quotient_joint_rank=a01)


def transition_bounds(profile,tau=1):
    """Positive upper transition and a capped live-return event for one U."""
    q,s=profile['field'],profile['state']
    r,e=profile['difference_rank'],profile['diagonal_dimension']
    k=profile['diagonal_difference_intersection']
    a0,a1,a01=(profile[key] for key in
              ('quotient_state_rank','quotient_feedback_rank','quotient_joint_rank'))
    tau=F(tau)
    if tau<1:
        raise ValueError('rank-moment mark must be at least one')
    nonzero=q**s-1
    up=(F((q**s-q**r)*(q**e-1),q**e*nonzero)
        -F(q**e-q**k,q**e*(q-1)*nonzero))
    good=q**r-q**(r-a0)-q**(r-a1)+q**(r-a01)
    down_lower=F(good,q**(e+1)*nonzero)
    down_upper=min(F(good,q**e*(q-1)*nonzero),1-up)
    assert 0<=up<=1 and 0<=down_lower<=down_upper<=1-up
    positive={}
    if down_upper:
        positive[r-1]=down_upper
    positive[r]=1-up-down_lower
    if up:
        positive[r+1]=up
    live_return=(down_upper if e==0 else 1-up if e==1 else F(1))
    # A normalized stochastic domination suffices for nonincreasing rank payoffs.
    moment=tau**(-r)*(1+up*(1/tau-1)+down_upper*(tau-1))
    return dict(profile=profile,up=up,down_lower=down_lower,down_upper=down_upper,
                good_directions=good,positive_rank_weights=positive,
                positive_mass=sum(positive.values(),F(0)),
                capped_live_return_upper=live_return,rank_moment_upper=moment,
                tau=tau,output_constraints_before_current_mixer=True,
                source_word_union=False,whole_certificate=False)


def successors_exact(joint,q):
    """Small-field oracle for the actual conditioned (u,v) sampler."""
    joint=columns_basis(tuple(tuple(row) for row in joint),q)
    s=len(joint)//2
    if q**s>64:
        raise ValueError('bounded exhaustive sampler only')
    state,feedback=joint[:s],joint[s:]
    difference=tuple(tuple(a^b for a,b in zip(x,y)) for x,y in zip(state,feedback))
    counts=Counter()
    vectors=tuple(product(range(q),repeat=s))
    for u in vectors:
        if not any(u):
            continue
        for v in vectors:
            if dot(v,u,q)==1:
                continue
            observed=tuple(dot(v,column,q) for column in zip(*state))
            successor=tuple(tuple(x^multiply(u[i],y,q) for x,y in zip(row,observed))
                            for i,row in enumerate(difference))
            counts[rank(successor,q)]+=1
    total=sum(counts.values())
    return {r:F(n,total) for r,n in counts.items()}


def all_subspaces(ambient,q):
    """Canonical column bases via all reduced row-echelon bases; tiny only."""
    if ambient>4 or q>4:
        raise ValueError('bounded subspace oracle only')
    for dimension in range(ambient+1):
        for pivots in combinations(range(ambient),dimension):
            slots=[(i,j) for i,pivot in enumerate(pivots)
                   for j in range(pivot+1,ambient) if j not in pivots]
            for entries in product(range(q),repeat=len(slots)):
                rows=[[int(j==pivot) for j in range(ambient)] for pivot in pivots]
                for (i,j),x in zip(slots,entries):
                    rows[i][j]=x
                yield tuple(tuple(rows[i][j] for i in range(dimension)) for j in range(ambient))


def beta64_line_receipt():
    state=(34,232,1,0)
    feedback=(72,27,211,244)
    profile=geometry(tuple((x,) for x in state+feedback))
    result=transition_bounds(profile,tau=65536)
    # These two vectors are noncollinear over GF256 and its extension.
    assert rank((state,feedback),65536)==2
    result['exact_return_probability']=F(65535,65536*(65536**4-1))
    result['state_vector']=state
    result['feedback_vector']=feedback
    return result


if __name__=='__main__':
    import json
    print(json.dumps(beta64_line_receipt(),default=str,indent=2))
