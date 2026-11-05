"""Exact toy audit of the fixed-coordinate consistency correction.

GF4 tests and a linear-time path-component kernel; no K20 certificate.
"""
from itertools import product
import random
import time
from test_layered_inner import field,setup,pair


def matrix_rank(rows,width,mul):
    rows=[row[:] for row in rows];r=0;q=len(mul)
    for col in range(width):
        pivot=next((j for j in range(r,len(rows)) if rows[j][col]),None)
        if pivot is None:continue
        rows[r],rows[pivot]=rows[pivot],rows[r]
        inv=next(a for a in range(1,q) if mul[a][rows[r][col]]==1)
        rows[r]=[mul[inv][x] for x in rows[r]]
        for j in range(r+1,len(rows)):
            c=rows[j][col]
            if c:rows[j]=[x^mul[c][y] for x,y in zip(rows[j],rows[r])]
        r+=1
    return r


def path_nullity(width,batches,inactive,zero_edges,fixed):
    """Kernel dimension of state zeros, output-zero edges, and fixed X=0.

    fixed maps physical positions to gamma in {0,1}: X=Qnext+gamma*Q.
    zero_edges are Y=0 positions with two active adjacent state variables.
    """
    n=width*batches;parent=list(range(n));pinned=[False]*n
    def root(i):
        while parent[i]!=i:
            parent[i]=parent[parent[i]];i=parent[i]
        return i
    def pin(i):pinned[root(i)]=True
    def join(i,j):
        a,c=root(i),root(j)
        if a!=c:parent[c]=a;pinned[a]|=pinned[c]
    for i in inactive:pin(i)
    # Only adjacent-time edges occur, so the only cycles are parallel edges.
    for position in zero_edges:
        assert position>=width
        join(position-width,position)
    for position,gamma in fixed.items():
        if position<width or gamma==0:
            pin(position)
        else:
            join(position-width,position)
            if position in zero_edges:pin(position)  # delta != 1
    return sum(parent[i]==i and not pinned[i] for i in range(n))


def sweep_nullity(width,batches,inactive,zero_edges,fixed):
    """Same kernel using one live-component bit per coordinate path.

    This is a counting interface, not an encoder or a certificate. A component
    is counted only after it ends: a later fixed equation can still pin it.
    """
    total=0
    for j in range(width):
        live=False
        for t in range(batches):
            pos=t*width+j;gamma=fixed.get(pos)
            zero=pos in zero_edges
            linked=t>0 and (zero or gamma==1)
            pinned=pos in inactive or gamma==0 or (t==0 and gamma is not None)
            pinned|=zero and gamma==1
            if linked:live=live and not pinned
            else:
                total+=live
                live=not pinned
        total+=live
    return total


def hard_class_target():
    """Required kernel cap, NOT a proved kernel cap."""
    a,D,e,K_R,dbar=790492,638768,368520,524288,163923
    old_bound=734997
    pattern_bits=old_bound-1-32*(a+dbar-D-e)
    cap=D+e-K_R+(-100-pattern_bits-1)//32
    bound=pattern_bits+1+32*(cap+K_R-D-e)
    assert cap==407155 and bound<=-100 and bound+32>-100
    return dict(pattern_bits=pattern_bits,required_kernel_upper=cap,
                conditional_log2_upper=bound,kernel_cap_proved=False)


def test():
    started=time.monotonic();mul=field(2,7);keys,_=setup(mul);rng=random.Random(602)
    b,L,M=4,4,4;n=b*L;K=8;cases=0;positive=0;best=None
    for trial in range(256):
        # RS[4,2] per block, one independent regional block permutation.
        G=[[0]*K for _ in range(n)]
        for region in range(4):
            order=list(range(M));rng.shuffle(order)
            for j,block in enumerate(order):
                G[region*M+j][2*block]=1;G[region*M+j][2*block+1]=region
        Y=[0]*n
        for j in rng.sample(range(n),1+trial%8):Y[j]=rng.randrange(1,4)
        states=[[0]*b]
        for t in range(L):states.append([mul[2][states[-1][j]]^Y[t*b+j] for j in range(b)])
        inactive={t*b+j for t in range(L-1) for j in range(b) if states[t+1][j]==0}
        zero_edges={t*b+j for t in range(1,L-1) for j in range(b)
                    if Y[t*b+j]==0 and states[t][j] and states[t+1][j]}
        fixed={j:0 for j in range(b)};exceptional=0;X=[0]*n;F=[]
        X[:b]=states[1]
        for t in range(1,L):
            bit=t%2
            for j in range(b):
                if j&(1<<bit):continue
                k=j^(1<<bit);x,y=states[t][j],states[t][k]
                key=rng.choice(keys);p,z=pair(mul,key,x,y)
                X[t*b+j]=p^states[t+1][j];X[t*b+k]=z^states[t+1][k]
                if not x and not y:
                    fixed[t*b+j]=fixed[t*b+k]=0
                elif x and not y:
                    fixed[t*b+j]=1;F.append(t*b+k)
                elif not p:
                    assert x and y
                    fixed[t*b+j]=0;fixed[t*b+k]=1;exceptional+=1
                else:F.extend((t*b+j,t*b+k))
        D=len(F);assert D+len(fixed)==n
        flat=[v for state in states[1:] for v in state]
        T=[];A=[]
        for j in inactive:
            row=[0]*n;row[j]=1;T.append(row)
        for j in zero_edges:
            row=[0]*n;row[j]=1;row[j-b]=2;T.append(row)
        for j,gamma in fixed.items():
            row=[0]*n;row[j]=1
            if j>=b:row[j-b]=gamma
            value=0
            for c,v in zip(row,flat):value^=mul[c][v]
            assert value==X[j]
            A.append(row)
        a=n-matrix_rank(T,n,mul)
        kappa=path_nullity(b,L,inactive,zero_edges,fixed)
        assert kappa==sweep_nullity(b,L,inactive,zero_edges,fixed)
        assert kappa==n-matrix_rank(T+A,n,mul)
        assert kappa<=D
        projection=[G[j] for j in fixed]
        d=K-matrix_rank(projection,K,mul)
        joint=[row+[0]*K for row in T]
        joint += [row+g for row,g in zip(A,projection)]
        exact_dimension=n+K-matrix_rank(joint,n+K,mul)-d
        upper=min(a,kappa+K-d)
        assert exact_dimension<=upper<=a
        rho=max(0,a-kappa-K+d)
        positive+=rho>0;cases+=1
        receipt=dict(output_weight=sum(bool(x) for x in Y),state_dimension=a,
                     free_coordinates=D,exceptions=exceptional,outer_shortening=d,
                     fixed_zero_kernel=kappa,extra_rank=rho,
                     exact_admissible_state_dimension=exact_dimension,
                     admissible_state_dimension_upper=upper)
        if best is None or rho>best['extra_rank']:best=receipt
    # Exhaustively compare graph and matrix kernels on all width-one length3
    # patterns, including duplicate pins, zeros, and parallel delta/1 edges.
    kernels=0
    for gammas in product((-1,0,1),repeat=3):
        fixed={j:g for j,g in enumerate(gammas) if g>=0}
        for zeros in range(4):
            ze={j for j in (1,2) if zeros&(1<<(j-1))}
            for mask in range(8):
                inactive={j for j in range(3) if mask&(1<<j)};rows=[]
                for j in inactive:
                    row=[0]*3;row[j]=1;rows.append(row)
                for j in ze:
                    row=[0]*3;row[j]=1;row[j-1]=2;rows.append(row)
                for j,g in fixed.items():
                    row=[0]*3;row[j]=1
                    if j:row[j-1]=g
                    rows.append(row)
                assert path_nullity(1,3,inactive,ze,fixed)==3-matrix_rank(rows,3,mul)
                assert sweep_nullity(1,3,inactive,ze,fixed)==path_nullity(1,3,inactive,ze,fixed)
                kernels+=1
    # Negative control: two fixed equations along a zero-output run impose
    # only one extra constraint after the temporal equations, not two.
    temporal=[[2,1,0],[0,2,1]]
    fixed_rows=[[1,0,0],[0,1,0]]
    assert matrix_rank(temporal+fixed_rows,3,mul)-matrix_rank(temporal,3,mul)==1
    assert time.monotonic()-started<15
    return dict(rs_affine_fixtures=cases,positive_correction_fixtures=positive,
                exhaustive_kernel_checks=kernels,largest_toy_correction=best,
                hard_class_target=hard_class_target(),
                whole_code_certificate=False)


if __name__=='__main__':print(test())
