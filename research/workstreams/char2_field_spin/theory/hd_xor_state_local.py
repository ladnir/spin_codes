"""Exact local support/weight kernel before diagonal-XOR-diagonal refresh.

One actual RS[16,4] block emits Y=X+A Q and feeds back W=Q+A^T X.
Q is independently uniform nonzero on its exact coordinate support S;
X is uniform-field on U, including zeros. Only |U|<=2 is implemented.
No selected coefficient-column subset is assumed MDS. The entire A is
checked MDS before using its distance to limit output-zero subsets.

Put V=Y_U. Translation makes V uniform and independent of Q, while
Y_out=A_out,S Q and W=(I+C_U A_U)[:,S]Q+C_U V. For nonempty S, Y_out
has at most h-1 zeros. Enumerate those zero constraints, remove forced
zero V variables, and enforce Q nonzero by integer inclusion-exclusion.
Exact W support follows by integer subset inversion. The final positive
zero-set coefficients and ordinary weight enumerators have no rounding.
This is a local law, not a distance certificate or regional operator.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations,product
from math import comb

import numpy as np

from hd_blockwise_shared_transfer import multiply,rank,rref,mm
from hd_bucket_feedback import maps


@lru_cache(5)
def field_tables(q):
    if type(q) is not int or q not in (2,4,8,16,256):
        raise ValueError('Supported fields are2,4,8,16,256')
    table=np.array([[multiply(i,j,q) for j in range(q)] for i in range(q)],dtype=np.uint8)
    inverse=np.zeros(q,dtype=np.uint8)
    for i in range(1,q):inverse[i]=int(np.flatnonzero(table[i]==1)[0])
    return table,inverse


def batch_ranks(matrices,q):
    """Exact elimination for a bounded batch of small field matrices."""
    work=np.asarray(matrices)
    if (work.ndim!=3 or len(work)>8192 or work.shape[1]>8 or work.shape[2]>6
            or not np.issubdtype(work.dtype,np.integer) or np.any(work<0)
            or np.any(work>=q)):
        raise ValueError('Invalid or over-budget rank batch')
    work=work.astype(np.uint8,copy=True)
    table,inverse=field_tables(q)
    count,rows,columns=work.shape
    ranks=np.zeros(count,dtype=np.int64)
    for column in range(columns):
        eligible=(np.arange(rows)[None,:]>=ranks[:,None])&(work[:,:,column]!=0)
        chosen=np.flatnonzero(eligible.any(axis=1))
        if not len(chosen):continue
        pivot=np.argmax(eligible[chosen],axis=1)
        target=ranks[chosen]
        old=work[chosen,target,:].copy()
        work[chosen,target,:]=work[chosen,pivot,:]
        work[chosen,pivot,:]=old
        normalized=table[inverse[work[chosen,target,column]][:,None],work[chosen,target,:]]
        factors=work[chosen,:,column].copy()
        factors[np.arange(len(chosen)),target]=0
        work[chosen]^=table[factors[:,:,None],normalized[:,None,:]]
        work[chosen,target,:]=normalized
        ranks[chosen]+=1
    return ranks


@lru_cache(1)
def current_emission():
    _,feedback=maps(local_rank=4)
    return tuple(tuple(feedback[j][i] for j in range(4)) for i in range(16))


@lru_cache(8)
def checked_emission(emission,q):
    field_tables(q)
    n=len(emission)
    h=len(emission[0]) if n else 0
    if (not 1<=h<=4 or not h<=n<=16 or any(len(row)!=h for row in emission)
            or any(type(x) is not int or not 0<=x<q for row in emission for x in row)):
        raise ValueError('A local emission of at most16 rows and4 columns is required')
    restrictions=np.array([tuple(emission[i] for i in indices)
                           for indices in combinations(range(n),h)],dtype=np.uint8)
    if np.any(batch_ranks(restrictions,q)!=h):
        raise ValueError('Full emission is not MDS; zero-set cutoff would be unsound')
    return n,h


def _nonzero_kernel_count(basis,s,q):
    total=0
    for keep in range(1<<s):
        columns=[i for i in range(s) if keep>>i&1]
        restricted=tuple(tuple(row[i] for i in columns) for row in basis)
        total+=(-1)**(s-len(columns))*q**(len(columns)-rank(restricted,q))
    return total


def _zero_spaces(emission,outside,S,q):
    """Group actual zero constraints by rowspace, with size multiplicities."""
    h=len(emission[0]);s=len(S)
    spaces={}
    examined=0
    for size in range(h):
        for positions in combinations(outside,size):
            examined+=1
            rows=tuple(tuple(emission[i][j] for j in S) for i in positions)
            reduced,pivots=rref(rows,q)
            basis=tuple(reduced[:len(pivots)])
            if basis not in spaces:
                spaces[basis]=[0]*h if _nonzero_kernel_count(basis,s,q)>0 else None
            if spaces[basis] is not None:spaces[basis][size]+=1
    kept=tuple((basis,tuple(counts)) for basis,counts in spaces.items() if counts is not None)
    return kept,examined


def _empty_state_counts(emission,U,q):
    """At most q**2 explicit inputs; Y=X when Q=0."""
    n,h=len(emission),len(emission[0])
    values=np.array(tuple(product(range(q),repeat=len(U))),dtype=np.uint8)
    table,_=field_tables(q)
    feedback=np.zeros((len(values),h),dtype=np.uint8)
    for j,position in enumerate(U):
        feedback^=table[values[:,j,None],np.array(emission[position])[None,:]]
    support=((feedback!=0)*(1<<np.arange(h))).sum(axis=1)
    weight=(values!=0).sum(axis=1)
    counts=np.zeros((1<<h,n+1),dtype=np.int64)
    np.add.at(counts,(support,weight),1)
    positive=np.array([[sum(int(row[w])*comb(n-w,k) for w in range(n+1) if n-w>=k)
                        for k in range(n+1)] for row in counts],dtype=np.int64)
    return counts,positive


def local_counts(input_support,state_mask,*,emission=None,field_order=256):
    """Exact counts indexed by (support(W),weight(Y)), with denominator.

    Also returns positive zero-set coefficients b[T,k]: the exact kernel
    entry is sum_k b[T,k] z**(n-k)(1-z)**k / denominator. Final ordinary
    weight counts are recovered by exact integer polynomial conversion.
    """
    q=field_order
    a=current_emission() if emission is None else tuple(tuple(row) for row in emission)
    n,h=checked_emission(a,q)
    U=tuple(input_support)
    if (len(U)>2 or len(set(U))!=len(U) or any(type(i) is not int or not 0<=i<n for i in U)
            or type(state_mask) is not int or not 0<=state_mask<(1<<h)):
        raise ValueError('Exact support required; input cutoff is two')
    U=tuple(sorted(U));S=tuple(i for i in range(h) if state_mask>>i&1)
    denominator=q**len(U)*(q-1)**len(S)
    if not S:
        counts,positive=_empty_state_counts(a,U,q)
        rank_systems=zero_subsets=0
    else:
        outside=tuple(i for i in range(n) if i not in U)
        spaces,zero_subsets=_zero_spaces(a,outside,S,q)
        multiplicities=np.array([row[1] for row in spaces],dtype=np.int64)
        basis=np.zeros((len(spaces),h-1,len(S)),dtype=np.uint8)
        for i,(rows,_) in enumerate(spaces):
            if rows:basis[i,:len(rows),:]=rows
        # Actual coefficient columns of I+C_U A_U, never an MDS surrogate.
        c_u=tuple(tuple(a[position][j] for position in U) for j in range(h))
        a_u=tuple(a[position] for position in U)
        correction=mm(c_u,a_u,q) if U else tuple((0,)*h for _ in range(h))
        transformed=tuple(tuple(int(i==j)^correction[i][j] for j in range(h)) for i in range(h))
        W=np.array([tuple(transformed[i][j] for j in S)+tuple(a[j][i] for j in U)
                    for i in range(h)],dtype=np.uint8)
        zero_flags=np.array([[(mask>>i)&1 for i in range(h)] for mask in range(1<<h)],dtype=np.uint8)
        sums=np.zeros((h+len(U),1<<h),dtype=np.int64)
        powers=np.array([q**i for i in range(7)],dtype=np.int64)
        rank_systems=0
        for input_zero in range(1<<len(U)):
            Vcols=[len(S)+i for i in range(len(U)) if not input_zero>>i&1]
            offset=input_zero.bit_count()
            for keep in range(1<<len(S)):
                Qcols=[i for i in range(len(S)) if keep>>i&1]
                columns=Qcols+Vcols
                sign=(-1)**(len(S)-len(Qcols))
                for start in range(0,len(spaces),256):
                    stop=min(start+256,len(spaces));batch=stop-start
                    constraints=np.zeros((batch,1<<h,2*h-1,len(columns)),dtype=np.uint8)
                    constraints[:,:,:h,:]=zero_flags[None,:,:,None]*W[:,columns][None,None,:,:]
                    constraints[:,:,h:,:len(Qcols)]=basis[start:stop][:,:,Qcols][:,None,:,:]
                    ranks=batch_ranks(constraints.reshape((batch*(1<<h),2*h-1,len(columns))),q).reshape((batch,1<<h))
                    solutions=powers[len(columns)-ranks]
                    # Absolute sum over all Q/V deletion terms is at most
                    # 697*(q+1)**(s+u). Even a16x subset-inversion bound
                    # is below2**63 for q<=256 and s+u<=6.
                    sums[offset:offset+h]+=sign*(solutions.T@multiplicities[start:stop]).T
                    rank_systems+=solutions.size
        exact=sums.copy()
        for bit in range(h):
            for mask in range(1<<h):
                if not mask>>bit&1:exact[:,mask]-=exact[:,mask|(1<<bit)]
        positive=exact[:n+1,np.arange(1<<h)^((1<<h)-1)].T
        if np.any(positive<0):raise ArithmeticError('Negative exact zero-set count')
        counts=np.zeros((1<<h,n+1),dtype=np.int64)
        for k in range(positive.shape[1]):
            for j in range(k+1):
                counts[:,n-k+j]+=((-1)**j)*comb(k,j)*positive[:,k]
    if np.any(counts<0) or int(counts.sum())!=denominator:
        raise ArithmeticError('Incomplete exact local enumerator')
    return dict(field_order=q,width=n,state=h,input_support=U,state_mask=state_mask,
        denominator=denominator,counts=tuple(tuple(int(x) for x in row) for row in counts),
        positive_zero_coefficients=tuple(tuple(int(x) for x in row) for row in positive),
        rank_systems=rank_systems,outside_zero_subsets=zero_subsets,
        exact_local=True,whole_code_certificate=False,input_law='uniform-field')


def marked_kernel(data,z):
    """Evaluate the exact nonnegative weight enumerators at a rational mark."""
    z=F(z)
    if not 0<=z<=1:raise ValueError('Weight mark must be in [0,1]')
    return tuple(sum((F(count,data['denominator'])*z**weight
                      for weight,count in enumerate(row)),F(0)) for row in data['counts'])


if __name__=='__main__':
    import argparse,json
    from math import log2
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--examples',action='store_true',help='Run seven bounded native-map local cases')
    args=parser.parse_args()
    if args.examples:
        for U,S in (((),15),((0,),1),((5,),3),((0,),15),((0,1),3),((0,5),15),((0,5),0)):
            data=local_counts(U,S)
            marked=marked_kernel(data,F(1,2))
            weights=[sum(row[w] for row in data['counts']) for w in range(17)]
            print(json.dumps(dict(input_support=U,state_mask=S,denominator=data['denominator'],
                rank_systems=data['rank_systems'],weight_counts=weights,
                possible_W=[i for i,row in enumerate(data['counts']) if sum(row)],
                zero_W_count=sum(data['counts'][0]),log2_marked_total=log2(sum(marked)),
                log2_marked_zero=None if not marked[0] else log2(marked[0]),
                exact_local=True,whole_code_certificate=False)),flush=True)
    else:parser.print_help()
