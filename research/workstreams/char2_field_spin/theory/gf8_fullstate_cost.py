"""Exact algebra/cost reference for the full-state8 GF256 candidate.

No native kernel or timing. The original packet has A=G_128,C=G_192,
where G_a[i,j]=1/(a xor i xor j), state8,width8,eta2. A scalar state
change normalizes both XOR-circulants and changes eta to13, exactly.
"""
from __future__ import annotations

from itertools import combinations
import json
from random import Random

from gf8_packet_maps import (MUL,INV,determinant,matmul,matvec,rank,
                             transpose_four_records,xor_sum)


A_KERNEL=tuple(INV[128^i] for i in range(8))
C_KERNEL=tuple(INV[192^i] for i in range(8))
A_SUM=xor_sum(A_KERNEL)
C_SUM=xor_sum(C_KERNEL)
NORMAL_A=tuple(MUL[INV[A_SUM]][v] for v in A_KERNEL)
NORMAL_C=tuple(MUL[INV[C_SUM]][v] for v in C_KERNEL)
NORMAL_ETA=MUL[2][MUL[A_SUM][C_SUM]]


def circulant(kernel):
    return tuple(tuple(kernel[i^j] for j in range(8)) for i in range(8))


def zeta8(values):
    result=list(values)
    for bit in (1,2,4):
        for index in range(8):
            if not index&bit:result[index]^=result[index|bit]
    return tuple(result)


def convolution19(kernel,values):
    """Nineteen explicit product sites; both zeta transforms use12 XORs."""
    g=zeta8(kernel)
    if g[0]!=1:raise ValueError('normalized coefficient sum must be one')
    z=list(zeta8(values))
    # Descending masks make this destructive schedule safe: later outputs
    # use only lower-index original z registers, never an overwritten value.
    z[7]^=(MUL[g[1]][z[6]]^MUL[g[2]][z[5]]^MUL[g[3]][z[4]]^
           MUL[g[4]][z[3]]^MUL[g[5]][z[2]]^MUL[g[6]][z[1]]^MUL[g[7]][z[0]])
    z[6]^=MUL[g[2]][z[4]]^MUL[g[4]][z[2]]^MUL[g[6]][z[0]]
    z[5]^=MUL[g[1]][z[4]]^MUL[g[4]][z[1]]^MUL[g[5]][z[0]]
    z[4]^=MUL[g[4]][z[0]]
    z[3]^=MUL[g[1]][z[2]]^MUL[g[2]][z[1]]^MUL[g[3]][z[0]]
    z[2]^=MUL[g[2]][z[0]]
    z[1]^=MUL[g[1]][z[0]]
    return zeta8(z)


def four_convolutions(kernel,packets):
    """Two explicit4x4 record transposes in each direction.

    Each of the eight intermediate ZMMs holds the same coordinate of four
    independent packets. Thus all four scalar circuits share19GFNI sites.
    This interpreter checks record ordering, not a compiler allocation.
    """
    if len(packets)!=4 or any(len(v)!=8 for v in packets):raise ValueError('four8-record packets required')
    packed=(transpose_four_records(tuple(tuple(v[:4]) for v in packets))+
            transpose_four_records(tuple(tuple(v[4:]) for v in packets)))
    transformed=tuple(convolution19(kernel,tuple(packed[i][lane] for i in range(8)))
                      for lane in range(4))
    packed_output=tuple(zip(*transformed))
    low=transpose_four_records(packed_output[:4])
    high=transpose_four_records(packed_output[4:])
    return tuple(a+b for a,b in zip(low,high))


def packet(values,state,dx,ds,normalized=True,adjoint=False):
    if normalized:
        a,c,eta=NORMAL_A,NORMAL_C,NORMAL_ETA
        apply=lambda kernel,x:convolution19(kernel,x)
    else:
        a,c,eta=A_KERNEL,C_KERNEL,2
        apply=lambda kernel,x:matvec(circulant(kernel),x)
    if adjoint:
        emitted=apply(c,state)
        injected=apply(a,values)
        return (tuple(MUL[d][v^e] for d,v,e in zip(dx,values,emitted)),
                tuple(MUL[d][h^MUL[eta][m]] for d,h,m in zip(ds,state,injected)))
    x=tuple(MUL[d][v] for d,v in zip(dx,values))
    r=tuple(MUL[d][v] for d,v in zip(ds,state))
    emitted=apply(a,r)
    injected=apply(c,x)
    return (tuple(v^MUL[eta][e] for v,e in zip(x,emitted)),
            tuple(v^e for v,e in zip(r,injected)))


def four_packet_adjoint(packets,state,input_diagonals,state_diagonals):
    """Packets are listed in reverse traversal order, newest to oldest."""
    moments=four_convolutions(NORMAL_A,packets)
    history=[]
    for moment,diagonal in zip(moments,state_diagonals):
        history.append(state)
        state=tuple(MUL[d][h^MUL[NORMAL_ETA][m]] for d,h,m in zip(diagonal,state,moment))
    emissions=four_convolutions(NORMAL_C,history)
    outputs=tuple(tuple(MUL[d][v^e] for d,v,e in zip(diagonal,values,emitted))
                  for diagonal,values,emitted in zip(input_diagonals,packets,emissions))
    return outputs,state


def cost():
    return dict(width=8,state=8,batch_packets=4,outputs_per_batch=32,
                each_four_convolutions_gfni512=19,both_maps_gfni512=38,
                input_diagonal_gfni512=8,state_diagonal_gfni512=8,eta_gfni512=8,
                batch_gfni512=62,batch_boolean=90,batch_record_shuffles=64,
                per16_gfni512=31,per16_boolean=45,per16_record_shuffles=32,
                fixed_nontrivial_coefficients=14,
                proposed_peak_live_zmm_at_most=28,
                register_bound_scope='destructive hand schedule, not compiler-validated',
                native_implemented=False,whole_code_certified=False)


def check():
    assert (A_SUM,C_SUM,NORMAL_ETA)==(15,74,13)
    assert xor_sum(NORMAL_A)==xor_sum(NORMAL_C)==1
    assert zeta8(NORMAL_A)==(1,179,167,39,148,101,182,93)
    assert zeta8(NORMAL_C)==(1,61,163,207,40,50,208,83)
    a,c=circulant(A_KERNEL),circulant(C_KERNEL)
    ac=matmul(a,c)
    # Partial fractions reduce the product to one Cauchy XOR-circulant.
    assert ac==circulant(tuple(MUL[A_SUM^C_SUM][INV[64^i]] for i in range(8)))
    returned=tuple(tuple(int(i==j)^MUL[2][ac[i][j]] for j in range(8)) for i in range(8))
    minors=0
    for size in range(1,9):
        for rows in combinations(range(8),size):
            for columns in combinations(range(8),size):
                assert determinant(tuple(tuple(returned[i][j] for j in columns) for i in rows))
                minors+=1
    assert minors==12869
    extended=tuple(tuple(int(i==j) for j in range(8))+tuple(MUL[2][v] for v in a[i]) for i in range(8))
    extended+=tuple(c[i]+tuple(int(i==j) for j in range(8)) for i in range(8))
    assert rank(extended)==16
    normal_a,normal_c=circulant(NORMAL_A),circulant(NORMAL_C)
    normal_extended=tuple(tuple(int(i==j) for j in range(8))+
                          tuple(MUL[NORMAL_ETA][v] for v in normal_a[i]) for i in range(8))
    normal_extended+=tuple(normal_c[i]+tuple(int(i==j) for j in range(8)) for i in range(8))
    circuit_cases=0
    for kernel in (NORMAL_A,NORMAL_C):
        dense=circulant(kernel)
        for coordinate in range(8):
            for scalar in range(256):
                values=tuple(scalar if i==coordinate else 0 for i in range(8))
                assert convolution19(kernel,values)==matvec(dense,values)
                circuit_cases+=1
    rng=Random(81913)
    basis_cases=0
    for trial in range(16*256+32):
        x=[0]*8;h=[0]*8
        if trial<16*256:
            coordinate,scalar=divmod(trial,256)
            (x if coordinate<8 else h)[coordinate%8]=scalar
        else:x=[rng.randrange(256) for _ in range(8)];h=[rng.randrange(256) for _ in range(8)]
        dx=tuple(rng.randrange(1,256) for _ in range(8));ds=tuple(rng.randrange(1,256) for _ in range(8))
        normalized=packet(x,h,dx,ds)
        # Independent dense16x16 reference, including the literal matrix
        # transpose. This does not derive its adjoint from packet().
        dense=tuple(tuple(MUL[v][d] for v,d in zip(row,dx+ds)) for row in normal_extended)
        flat=matvec(dense,tuple(x)+tuple(h))
        assert normalized==(flat[:8],flat[8:])
        original=packet(x,tuple(MUL[C_SUM][v] for v in h),dx,ds,normalized=False)
        assert normalized[0]==original[0]
        assert tuple(MUL[C_SUM][v] for v in normalized[1])==original[1]
        normalized=packet(x,h,dx,ds,adjoint=True)
        flat=matvec(tuple(zip(*dense)),tuple(x)+tuple(h))
        assert normalized==(flat[:8],flat[8:])
        original=packet(x,tuple(MUL[INV[C_SUM]][v] for v in h),dx,ds,normalized=False,adjoint=True)
        assert normalized[0]==original[0]
        assert normalized[1]==tuple(MUL[C_SUM][v] for v in original[1])
        basis_cases+=2
    batches=0
    for _ in range(128):
        packets=tuple(tuple(rng.randrange(256) for _ in range(8)) for _ in range(4))
        h=tuple(rng.randrange(256) for _ in range(8))
        dx=tuple(tuple(rng.randrange(1,256) for _ in range(8)) for _ in range(4))
        ds=tuple(tuple(rng.randrange(1,256) for _ in range(8)) for _ in range(4))
        got,final=four_packet_adjoint(packets,h,dx,ds)
        expected=[];state=h
        for i in range(4):
            output,state=packet(packets[i],state,dx[i],ds[i],adjoint=True)
            expected.append(output)
        assert got==tuple(expected) and final==state
        batches+=1
    return dict(return_minors=minors,circuit_basis_scalars=circuit_cases,
                fullstep_basis_scalar_random=basis_cases,dense_fullstep_and_transpose=basis_cases,
                four_packet_batches=batches,
                exact_checks=minors+circuit_cases+basis_cases+batches)


if __name__=='__main__':
    print(json.dumps(dict(checks=check(),cost=cost(),row_sums=[A_SUM,C_SUM],normalized_eta=NORMAL_ETA)))
