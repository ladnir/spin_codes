"""Exact width8/state4 full-Cauchy packet candidate; no native kernel.

The12x12 Cauchy matrix is partitioned [[F,A0],[C0,E]]. Replacing its
bottom rows by E^-1 times those rows gives [[F,A0],[C,I]], C=E^-1C0.
The public interface returns F,A=A0/2,C with eta2. Both the birth map F
and the return map F+2AC are fully superregular.
"""
from __future__ import annotations

from itertools import combinations
import json
from random import Random

from gf8_packet_maps import (MUL,INV,determinant,matmul,matvec,rank,
                             transpose_four_records,xor_sum,zeta4,five_product_circulant)


WIDTH=8
STATE=4
ETA=2
ROWS=tuple(range(12))
COLUMNS=tuple(range(128,140))


def identity(size):
    return tuple(tuple(int(i==j) for j in range(size)) for i in range(size))


def add(left,right):
    return tuple(tuple(a^b for a,b in zip(lrow,rrow)) for lrow,rrow in zip(left,right))


def scale(matrix,scalar):
    return tuple(tuple(MUL[scalar][v] for v in row) for row in matrix)


T=tuple(tuple(INV[u^v] for v in COLUMNS) for u in ROWS)
F_ORIGINAL=tuple(row[:8] for row in T[:8])
EFFECTIVE_A_ORIGINAL=tuple(row[8:] for row in T[:8])
C0=tuple(row[:8] for row in T[8:])
E=tuple(row[8:] for row in T[8:])
E_SUM=xor_sum(E[0])
E_INVERSE=scale(E,INV[MUL[E_SUM][E_SUM]])
C_ORIGINAL=matmul(E_INVERSE,C0)

# Constant row/column scales preserve every rank and Hamming support bound.
# They give unit-sum kernels for both state maps, enabling five products.
OUTPUT_SCALES=tuple(INV[xor_sum(EFFECTIVE_A_ORIGINAL[i])] for i in (0,4))
INPUT_SCALES=tuple(INV[xor_sum(C_ORIGINAL[0][j:j+4])] for j in (0,4))
F_NORMAL=tuple(tuple(MUL[OUTPUT_SCALES[i//4]][MUL[v][INPUT_SCALES[j//4]]]
                     for j,v in enumerate(row)) for i,row in enumerate(F_ORIGINAL))
EFFECTIVE_A_NORMAL=tuple(tuple(MUL[OUTPUT_SCALES[i//4]][v] for v in row)
                             for i,row in enumerate(EFFECTIVE_A_ORIGINAL))
C_NORMAL=tuple(tuple(MUL[v][INPUT_SCALES[j//4]] for j,v in enumerate(row))
                   for row in C_ORIGINAL)


def matrices(normalized=True):
    """Return F(8x8),A(8x4),C(4x8),eta2; packet y=Fx+eta*A*r."""
    f,a,c=((F_NORMAL,EFFECTIVE_A_NORMAL,C_NORMAL) if normalized else
           (F_ORIGINAL,EFFECTIVE_A_ORIGINAL,C_ORIGINAL))
    return f,scale(a,INV[ETA]),c,ETA


def extended(normalized=True):
    f,a,c,eta=matrices(normalized)
    effective=scale(a,eta)
    return tuple(f[i]+effective[i] for i in range(8))+tuple(c[i]+identity(4)[i] for i in range(4))


def return_matrix(normalized=True):
    f,a,c,eta=matrices(normalized)
    return add(f,scale(matmul(a,c),eta))


def selected_rank(output_rows,input_columns,state_columns,include_return=False,normalized=True):
    """Exact numerical rank for a fixed support pattern.

    The proven formula is min(len(output_rows)+4*include_return,
                              len(input_columns)+len(state_columns)).
    All four outgoing-state rows are included on return, not a subset.
    """
    rows=tuple(output_rows)+(tuple(range(8,12)) if include_return else ())
    columns=tuple(input_columns)+tuple(8+j for j in state_columns)
    matrix=extended(normalized)
    return rank(tuple(tuple(matrix[i][j] for j in columns) for i in rows))


def convolution9(kernel,values):
    """Four-point XOR-circulant with nine explicit multiplication sites."""
    g0,g1,g2,g3=zeta4(kernel)
    z0,z1,z2,z3=zeta4(values)
    return zeta4((MUL[g0][z0],
                  MUL[g0][z1]^MUL[g1][z0],
                  MUL[g0][z2]^MUL[g2][z0],
                  MUL[g0][z3]^MUL[g1][z2]^MUL[g2][z1]^MUL[g3][z0]))


def four_maps(kernels,values,unit_sum):
    """Model four scalar circuits packed into four128-bit ZMM lanes."""
    packed=transpose_four_records(values)
    scalar=five_product_circulant if unit_sum else convolution9
    result=tuple(scalar(k,tuple(packed[i][lane] for i in range(4)))
                 for lane,k in enumerate(kernels))
    packed_result=tuple(zip(*result))
    return transpose_four_records(packed_result)


def block_kernel(matrix,row,column):
    kernel=tuple(matrix[row][column:column+4])
    assert all(matrix[row+i][column+j]==kernel[i^j] for i in range(4) for j in range(4))
    return kernel


def packed_maps(values,state,adjoint=False):
    """Fixed14-GFNI map schedule, separate from the three diagonals."""
    f,a,c=F_NORMAL,EFFECTIVE_A_NORMAL,C_NORMAL
    if adjoint:
        f=tuple(zip(*f))
        a,c=tuple(zip(*c)),tuple(zip(*a))
    mix_kernels=tuple(block_kernel(f,i,j) for i,j in ((0,0),(0,4),(4,0),(4,4)))
    mixed=four_maps(mix_kernels,(tuple(values[:4]),tuple(values[4:]))*2,False)
    pair_kernels=(block_kernel(c,0,0),block_kernel(c,0,4),
                  block_kernel(a,0,0),block_kernel(a,4,0))
    paired=four_maps(pair_kernels,(tuple(values[:4]),tuple(values[4:]),tuple(state),tuple(state)),True)
    output=(tuple(x^y^z for x,y,z in zip(mixed[0],mixed[1],paired[2]))+
            tuple(x^y^z for x,y,z in zip(mixed[2],mixed[3],paired[3])))
    next_state=tuple(x^y^z for x,y,z in zip(state,paired[0],paired[1]))
    return output,next_state


def packet(values,state,dx,ds,adjoint=False):
    if adjoint:
        output,next_state=packed_maps(values,state,True)
        return (tuple(MUL[d][v] for d,v in zip(dx,output)),
                tuple(MUL[d][v] for d,v in zip(ds,next_state)))
    x=tuple(MUL[d][v] for d,v in zip(dx,values))
    r=tuple(MUL[d][v] for d,v in zip(ds,state))
    return packed_maps(x,r)


def all_minors(matrix):
    count=0
    for size in range(1,min(len(matrix),len(matrix[0]))+1):
        for rows in combinations(range(len(matrix)),size):
            for columns in combinations(range(len(matrix[0])),size):
                if not determinant(tuple(tuple(matrix[i][j] for j in columns) for i in rows)):
                    raise AssertionError((rows,columns))
                count+=1
    return count


def cost():
    return dict(width=8,state=4,eta=2,per8_pair_maps_gfni512=5,
                per8_input_mix_gfni512=9,per8_diagonal_gfni512=3,
                per16_gfni512=34,per16_boolean=54,per16_record_shuffles=64,
                excludes='loads/stores/routing/diagonal expansion/spills/outer',
                native_implemented=False,whole_code_certified=False)


def fullstate_matrices(normalized=True):
    """Full16 Cauchy extension: width8/state8, same eta2 interface."""
    t=tuple(tuple(INV[i^(128+j)] for j in range(16)) for i in range(16))
    b=tuple(row[:8] for row in t[:8]);a0=tuple(row[8:] for row in t[:8])
    c0=tuple(row[:8] for row in t[8:]);e=tuple(row[8:] for row in t[8:])
    e_sum=xor_sum(e[0]);ei=scale(e,INV[MUL[e_sum][e_sum]])
    assert matmul(e,ei)==identity(8)
    c=matmul(ei,c0)
    if normalized:
        l=INV[xor_sum(a0[0])];r=INV[xor_sum(c[0])]
        b=scale(b,MUL[l][r]);a0=scale(a0,l);c=scale(c,r)
    return b,scale(a0,INV[ETA]),c,ETA


def fullstate_extended(normalized=True):
    f,a,c,eta=fullstate_matrices(normalized)
    a=scale(a,eta)
    return tuple(f[i]+a[i] for i in range(8))+tuple(c[i]+identity(8)[i] for i in range(8))


def fullstate_packet(values,state,dx,ds,adjoint=False):
    """Three eight-point circuits; effective eta is folded into A."""
    from gf8_fullstate_cost import convolution19
    f,a,c,eta=fullstate_matrices()
    a=scale(a,eta)
    f_sum=xor_sum(f[0]);f_kernel=tuple(MUL[INV[f_sum]][v] for v in f[0])
    if not adjoint:
        values=tuple(MUL[d][v] for d,v in zip(dx,values))
        state=tuple(MUL[d][v] for d,v in zip(ds,state))
    # All three8x8 matrices are symmetric XOR-circulants.
    mixed=tuple(MUL[f_sum][v] for v in convolution19(f_kernel,values))
    emitted=convolution19(c[0] if adjoint else a[0],state)
    injected=convolution19(a[0] if adjoint else c[0],values)
    output=tuple(v^e for v,e in zip(mixed,emitted))
    next_state=tuple(v^e for v,e in zip(state,injected))
    if adjoint:
        output=tuple(MUL[d][v] for d,v in zip(dx,output))
        next_state=tuple(MUL[d][v] for d,v in zip(ds,next_state))
    return output,next_state


def fullstate_four_packet_adjoint(packets,state,input_diagonals,state_diagonals):
    from gf8_fullstate_cost import four_convolutions
    f,a,c,eta=fullstate_matrices()
    a=scale(a,eta);f_sum=xor_sum(f[0])
    f_kernel=tuple(MUL[INV[f_sum]][v] for v in f[0])
    moments=four_convolutions(a[0],packets)
    # Destructively replace retained raw inputs by their feedthrough moments
    # before creating any state history. This avoids a third8-ZMM bank.
    mixed=tuple(tuple(MUL[f_sum][v] for v in values)
                for values in four_convolutions(f_kernel,packets))
    history=[]
    for moment,diagonal in zip(moments,state_diagonals):
        history.append(state)
        state=tuple(MUL[d][h^v] for d,h,v in zip(diagonal,state,moment))
    emissions=four_convolutions(c[0],history)
    outputs=tuple(tuple(MUL[d][v^e] for d,v,e in zip(diagonal,values,emitted))
                  for diagonal,values,emitted in zip(input_diagonals,mixed,emissions))
    return outputs,state


def fullstate_cost():
    return dict(width=8,state=8,batch_packets=4,batch_output_records=32,
                batch_input_mix_gfni512=27,batch_state_maps_gfni512=38,
                batch_diagonal_gfni512=16,batch_gfni512=81,batch_boolean=127,
                batch_record_shuffles=96,per16_gfni512=40.5,per16_boolean=63.5,
                per16_record_shuffles=48,proposed_peak_live_zmm_at_most=29,
                native_implemented=False,whole_code_certified=False)


def fullstate_factored_packet(values,state,dx,ds,adjoint=False):
    """Same8/8 map, using F=242*(effectiveA)*C to remove its third map."""
    from gf8_fullstate_cost import convolution19
    _,a,c,eta=fullstate_matrices();a=scale(a,eta)
    if not adjoint:
        values=tuple(MUL[d][v] for d,v in zip(dx,values))
        state=tuple(MUL[d][v] for d,v in zip(ds,state))
    moment=convolution19(a[0] if adjoint else c[0],values)
    emission_input=tuple(h^MUL[242][m] for h,m in zip(state,moment))
    output=convolution19(c[0] if adjoint else a[0],emission_input)
    next_state=tuple(h^m for h,m in zip(state,moment))
    if adjoint:
        output=tuple(MUL[d][v] for d,v in zip(dx,output))
        next_state=tuple(MUL[d][v] for d,v in zip(ds,next_state))
    return output,next_state


def fullstate_factored_four_packet_adjoint(packets,state,input_diagonals,state_diagonals):
    from gf8_fullstate_cost import four_convolutions
    _,a,c,eta=fullstate_matrices();a=scale(a,eta)
    moments=four_convolutions(a[0],packets)
    histories=[]
    for moment,diagonal in zip(moments,state_diagonals):
        histories.append(tuple(h^MUL[242][m] for h,m in zip(state,moment)))
        state=tuple(MUL[d][h^m] for d,h,m in zip(diagonal,state,moment))
    emitted=four_convolutions(c[0],histories)
    return (tuple(tuple(MUL[d][v] for d,v in zip(diagonal,values))
                  for diagonal,values in zip(input_diagonals,emitted)),state)


def fullstate_factored_cost():
    return dict(width=8,state=8,batch_packets=4,batch_output_records=32,
                batch_state_maps_gfni512=38,batch_history_scale_gfni512=8,
                batch_diagonal_gfni512=16,batch_gfni512=62,batch_boolean=90,
                batch_record_shuffles=64,per16_gfni512=31,per16_boolean=45,
                per16_record_shuffles=32,proposed_peak_live_zmm_at_most=20,
                native_implemented=False,whole_code_certified=False)


def largestate_matrices(normalized=True):
    """Full24 Cauchy extension, width8/state16, same eta2 interface."""
    t=tuple(tuple(INV[i^(128+j)] for j in range(24)) for i in range(24))
    b=tuple(row[:8] for row in t[:8]);a0=tuple(row[8:] for row in t[:8])
    c0=tuple(row[:8] for row in t[8:]);e=tuple(row[8:] for row in t[8:])
    e_sum=xor_sum(e[0]);ei=scale(e,INV[MUL[e_sum][e_sum]])
    c=matmul(ei,c0)
    if normalized:
        beta=(xor_sum(c[0]),xor_sum(c[8]))
        l=INV[MUL[xor_sum(a0[0][:8])][beta[0]]]
        b=scale(b,l)
        a0=tuple(tuple(MUL[l][MUL[v][beta[j//8]]] for j,v in enumerate(row)) for row in a0)
        c=tuple(tuple(MUL[INV[beta[i//8]]][v] for v in row) for i,row in enumerate(c))
    return b,scale(a0,INV[ETA]),c,ETA


def largestate_extended(normalized=True):
    f,a,c,eta=largestate_matrices(normalized);a=scale(a,eta)
    return tuple(f[i]+a[i] for i in range(8))+tuple(c[i]+identity(16)[i] for i in range(16))


def largestate_packet(values,state,dx,ds,adjoint=False):
    """Factored8/16 packet: F=177*Aeff0*C0, secondA block sums110."""
    from gf8_fullstate_cost import convolution19
    _,a,c,eta=largestate_matrices();a=scale(a,eta)
    a0=a[0][:8];a1=tuple(MUL[INV[110]][v] for v in a[0][8:])
    c0,c1=c[0],c[8]
    if not adjoint:
        values=tuple(MUL[d][v] for d,v in zip(dx,values))
        state=tuple(MUL[d][v] for d,v in zip(ds,state))
    if adjoint:
        m0=convolution19(a0,values)
        m1=tuple(MUL[110][v] for v in convolution19(a1,values))
        e0=convolution19(c0,tuple(h^MUL[177][m] for h,m in zip(state[:8],m0)))
        e1=convolution19(c1,state[8:])
    else:
        m0,m1=convolution19(c0,values),convolution19(c1,values)
        e0=convolution19(a0,tuple(h^MUL[177][m] for h,m in zip(state[:8],m0)))
        e1=tuple(MUL[110][v] for v in convolution19(a1,state[8:]))
    output=tuple(x^y for x,y in zip(e0,e1))
    next_state=tuple(h^m for h,m in zip(state,m0+m1))
    if adjoint:
        output=tuple(MUL[d][v] for d,v in zip(dx,output))
        next_state=tuple(MUL[d][v] for d,v in zip(ds,next_state))
    return output,next_state


def largestate_four_packet_adjoint(packets,state,input_diagonals,state_diagonals):
    from gf8_fullstate_cost import four_convolutions
    _,a,c,eta=largestate_matrices();a=scale(a,eta)
    a0=a[0][:8];a1=tuple(MUL[INV[110]][v] for v in a[0][8:])
    moments0=four_convolutions(a0,packets)
    # The raw-input bank is overwritten by these second-block moments.
    moments1=tuple(tuple(MUL[110][v] for v in row) for row in four_convolutions(a1,packets))
    histories0=[];histories1=[]
    h0,h1=state[:8],state[8:]
    for m0,m1,diagonal in zip(moments0,moments1,state_diagonals):
        histories0.append(tuple(h^MUL[177][m] for h,m in zip(h0,m0)))
        histories1.append(h1)
        h0=tuple(MUL[d][h^m] for d,h,m in zip(diagonal[:8],h0,m0))
        h1=tuple(MUL[d][h^m] for d,h,m in zip(diagonal[8:],h1,m1))
    emissions0=four_convolutions(c[0],histories0)
    emissions1=four_convolutions(c[8],histories1)
    outputs=tuple(tuple(MUL[d][x^y] for d,x,y in zip(diagonal,e0,e1))
                  for diagonal,e0,e1 in zip(input_diagonals,emissions0,emissions1))
    return outputs,tuple(h0)+tuple(h1)


def largestate_cost():
    return dict(width=8,state=16,batch_packets=4,batch_output_records=32,
                batch_unit_circulants_gfni512=76,batch_nonunit_scale_gfni512=8,
                batch_history_scale_gfni512=8,batch_diagonal_gfni512=24,
                batch_gfni512=116,batch_boolean=180,batch_record_shuffles=128,
                per16_gfni512=58,per16_boolean=90,per16_record_shuffles=64,
                proposed_peak_live_zmm_at_most=30,native_implemented=False,
                whole_code_certified=False)


def check_largestate():
    f,a,c,eta=largestate_matrices();ae=scale(a,eta)
    returned=add(f,matmul(ae,c))
    counts=dict(birth=all_minors(f),return_map=all_minors(returned))
    assert counts==dict(birth=12869,return_map=12869)
    assert rank(largestate_extended())==24
    t=tuple(tuple(INV[i^(128+j)] for j in range(24)) for i in range(24))
    e=tuple(row[8:] for row in t[8:]);se=xor_sum(e[0])
    assert se==16
    ei=scale(e,INV[MUL[se][se]])
    assert matmul(e,ei)==matmul(ei,e)==identity(16)
    assert all(e[i][j]==e[0][i^j] for i in range(16) for j in range(16))
    left=tuple(identity(24)[:8])+tuple((0,)*8+row for row in ei)
    assert matmul(left,t)==largestate_extended(False)
    row_scale=(234,)*8+(INV[184],)*8+(INV[24],)*8
    column_scale=(1,)*8+(184,)*8+(24,)*8
    scaled=tuple(tuple(MUL[row_scale[i]][MUL[v][column_scale[j]]] for j,v in enumerate(row))
                 for i,row in enumerate(largestate_extended(False)))
    assert scaled==largestate_extended()
    assert (xor_sum(f[0]),xor_sum(ae[0][:8]),xor_sum(ae[0][8:]),xor_sum(c[0]),xor_sum(c[8]))==(177,1,110,1,1)
    assert f==scale(matmul(tuple(row[:8] for row in ae),c[:8]),177)
    assert all(ae[i][j]==ae[0][i^(j%8)+(j//8)*8] for i in range(8) for j in range(16))
    assert all(c[i][j]==c[(i//8)*8][(i%8)^j] for i in range(16) for j in range(8))
    rng=Random(240816)
    matrix=largestate_extended()
    rank_cases=0
    for _ in range(1024):
        z=tuple(rng.sample(range(8),rng.randrange(9)))
        i=tuple(rng.sample(range(8),rng.randrange(9)))
        h=tuple(rng.sample(range(16),rng.randrange(17)))
        columns=i+tuple(8+j for j in h)
        for returning in (False,True):
            rows=z+(tuple(range(8,24)) if returning else ())
            assert rank(tuple(tuple(matrix[u][v] for v in columns) for u in rows))==min(len(z)+16*returning,len(i)+len(h))
            rank_cases+=1
    scalar_cases=0
    for trial in range(24*256+32):
        values=[0]*24
        if trial<24*256:
            coordinate,scalar=divmod(trial,256);values[coordinate]=scalar
        else:values=[rng.randrange(256) for _ in range(24)]
        dx=tuple(rng.randrange(1,256) for _ in range(8));ds=tuple(rng.randrange(1,256) for _ in range(16))
        dense=tuple(tuple(MUL[v][d] for v,d in zip(row,dx+ds)) for row in matrix)
        for adjoint in (False,True):
            got=largestate_packet(values[:8],values[8:],dx,ds,adjoint)
            flat=matvec(tuple(zip(*dense)) if adjoint else dense,values)
            assert got==(flat[:8],flat[8:]);scalar_cases+=1
    for _ in range(128):
        packets=tuple(tuple(rng.randrange(256) for _ in range(8)) for _ in range(4))
        state=tuple(rng.randrange(256) for _ in range(16))
        dx=tuple(tuple(rng.randrange(1,256) for _ in range(8)) for _ in range(4))
        ds=tuple(tuple(rng.randrange(1,256) for _ in range(16)) for _ in range(4))
        got,terminal=largestate_four_packet_adjoint(packets,state,dx,ds)
        expected=[]
        for i in range(4):
            output,state=largestate_packet(packets[i],state,dx[i],ds[i],True)
            expected.append(output)
        assert got==tuple(expected) and terminal==state
    return dict(minors=counts,selected_rank_cases=rank_cases,
                packed_basis_scalar_random=scalar_cases,four_packet_batches=128)


def check_fullstate():
    f,a,c,eta=fullstate_matrices()
    returned=add(f,scale(matmul(a,c),eta))
    counts=dict(birth=all_minors(f),return_map=all_minors(returned),
                emission=all_minors(a),injection=all_minors(c))
    assert all(value==12869 for value in counts.values())
    assert rank(fullstate_extended())==16
    f0,a0,c0,_=fullstate_matrices(False)
    assert (xor_sum(f0[0]),xor_sum(scale(a0,2)[0]),xor_sum(c0[0]))==(15,161,46)
    assert xor_sum(f[0])==242 and xor_sum(scale(a,2)[0])==xor_sum(c[0])==1
    assert f==scale(matmul(scale(a,2),c),242)
    t=tuple(tuple(INV[i^(128+j)] for j in range(16)) for i in range(16))
    ei=scale(f0,INV[MUL[15][15]])
    left=tuple(identity(16)[:8])+tuple((0,)*8+row for row in ei)
    assert matmul(left,t)==fullstate_extended(False)
    scaled=tuple(tuple(MUL[124 if i<8 else 1][MUL[v][162 if j<8 else 1]]
                       for j,v in enumerate(row)) for i,row in enumerate(fullstate_extended(False)))
    assert scaled==fullstate_extended()
    rng=Random(160840)
    base=fullstate_extended()
    cases=0
    for trial in range(16*256+32):
        values=[0]*16
        if trial<16*256:
            coordinate,scalar=divmod(trial,256);values[coordinate]=scalar
        else:values=[rng.randrange(256) for _ in range(16)]
        dx=tuple(rng.randrange(1,256) for _ in range(8));ds=tuple(rng.randrange(1,256) for _ in range(8))
        dense=tuple(tuple(MUL[v][d] for v,d in zip(row,dx+ds)) for row in base)
        for adjoint in (False,True):
            got=fullstate_packet(values[:8],values[8:],dx,ds,adjoint)
            assert got==fullstate_factored_packet(values[:8],values[8:],dx,ds,adjoint)
            flat=matvec(tuple(zip(*dense)) if adjoint else dense,values)
            assert got==(flat[:8],flat[8:]);cases+=1
    for _ in range(128):
        packets=tuple(tuple(rng.randrange(256) for _ in range(8)) for _ in range(4))
        state=tuple(rng.randrange(256) for _ in range(8))
        dx=tuple(tuple(rng.randrange(1,256) for _ in range(8)) for _ in range(4))
        ds=tuple(tuple(rng.randrange(1,256) for _ in range(8)) for _ in range(4))
        got,terminal=fullstate_four_packet_adjoint(packets,state,dx,ds)
        assert (got,terminal)==fullstate_factored_four_packet_adjoint(packets,state,dx,ds)
        expected=[]
        for i in range(4):
            output,state=fullstate_packet(packets[i],state,dx[i],ds[i],True)
            expected.append(output)
        assert got==tuple(expected) and terminal==state
    return dict(minors=counts,packed_basis_scalar_random=cases,four_packet_batches=128)


def check():
    assert matmul(E,E_INVERSE)==matmul(E_INVERSE,E)==identity(4)
    assert E_SUM==2 and OUTPUT_SCALES==(100,8) and INPUT_SCALES==(200,16)
    f,a,c,eta=matrices()
    counts=dict(birth=all_minors(f),return_map=all_minors(return_matrix()),
                emission=all_minors(a),injection=all_minors(c))
    assert counts==dict(birth=12869,return_map=12869,emission=494,injection=494)
    assert rank(extended())==rank(extended(False))==12
    # Check the exact factorization which proves every joint-support rank.
    bottom=tuple((0,)*8+row for row in E_INVERSE)
    left=tuple(identity(12)[:8])+bottom
    assert matmul(left,T)==extended(False)
    rows=tuple(OUTPUT_SCALES[i//4] for i in range(8))+(1,)*4
    columns=tuple(INPUT_SCALES[i//4] for i in range(8))+(1,)*4
    scaled=tuple(tuple(MUL[rows[i]][MUL[v][columns[j]]] for j,v in enumerate(row))
                 for i,row in enumerate(extended(False)))
    assert scaled==extended()
    for row in (0,4):assert xor_sum(block_kernel(EFFECTIVE_A_NORMAL,row,0))==1
    for column in (0,4):assert xor_sum(block_kernel(C_NORMAL,0,column))==1
    rng=Random(120834)
    rank_cases=0
    # Every support-cardinality triple, plus arbitrary patterns. The universal
    # rank statement follows from the checked factorization and Cauchy formula.
    patterns=[]
    for z in range(9):
        for w in range(9):
            for h in range(5):
                patterns.append((tuple(rng.sample(range(8),z)),tuple(rng.sample(range(8),w)),
                                 tuple(rng.sample(range(4),h))))
    for _ in range(512):
        patterns.append(tuple(tuple(i for i in range(n) if rng.randrange(2)) for n in (8,8,4)))
    for z,i,h in patterns:
        for returning in (False,True):
            assert selected_rank(z,i,h,returning)==min(len(z)+4*returning,len(i)+len(h))
            rank_cases+=1
    packed_cases=0
    base=extended()
    for trial in range(12*256+32):
        values=[0]*12
        if trial<12*256:
            coordinate,scalar=divmod(trial,256);values[coordinate]=scalar
        else:values=[rng.randrange(256) for _ in range(12)]
        dx=tuple(rng.randrange(1,256) for _ in range(8));ds=tuple(rng.randrange(1,256) for _ in range(4))
        dense=tuple(tuple(MUL[v][d] for v,d in zip(row,dx+ds)) for row in base)
        for adjoint in (False,True):
            got=packet(values[:8],values[8:],dx,ds,adjoint)
            flat=matvec(tuple(zip(*dense)) if adjoint else dense,values)
            assert got==(flat[:8],flat[8:])
            packed_cases+=1
    return dict(minors=counts,joint_rank_cases=rank_cases,packed_basis_scalar_random=packed_cases)


if __name__=='__main__':
    print(json.dumps(dict(checks=check(),cost=cost(),matrices=matrices())))
