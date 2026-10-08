"""Exact packet maps for GF256 inner-design screens.

The general scaled-Cauchy family has C=A^T and four phases. The preferred
width-eight map repairs the old embedded-GF4 cubic circuit and has eight
phases, C!=A^T, and stronger return weights. Both families have CA=0.
This file establishes local algebra, not code distance.
No native encoder, build, or benchmark is involved.
"""

from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
import json
from math import comb


def multiply(a,b):
    value=0
    for _ in range(8):
        if b&1:value^=a
        a=((a<<1)^(0x11b if a&128 else 0))&255
        b>>=1
    return value


MUL=tuple(bytes(multiply(a,b) for b in range(256)) for a in range(256))
INV=(0,)+tuple(next(b for b in range(1,256) if MUL[a][b]==1) for a in range(1,256))


def xor_sum(values):
    result=0
    for value in values:result^=value
    return result


def product(values):
    result=1
    for value in values:result=MUL[result][value]
    return result


def power(value,exponent):
    result=1
    while exponent:
        if exponent&1:result=MUL[result][value]
        value=MUL[value][value]
        exponent>>=1
    return result


def geometry(width,state,phase=0):
    if width not in (8,16,32) or state not in (4,8,16) or 2*state>width or not 0<=phase<4:
        raise ValueError('expected width8/16/32, power-of-two state<=width/2, and phase0..3')
    labels=tuple((phase*width)^i for i in range(width))
    residues=tuple(128^j for j in range(state))
    derivative=product(range(1,state))
    return labels,residues,derivative


@lru_cache(maxsize=None)
def matrices(width,state,phase=0):
    labels,residues,derivative=geometry(width,state,phase)
    a=[]
    for point in labels:
        row_scale=MUL[product(point^v for v in residues)][INV[derivative]]
        a.append(tuple(MUL[row_scale][INV[point^v]] for v in residues))
    a=tuple(a)
    return a,tuple(zip(*a))


def matmul(left,right):
    return tuple(tuple(xor_sum(MUL[a][b] for a,b in zip(row,column))
                       for column in zip(*right)) for row in left)


def matvec(matrix,vector):
    return tuple(xor_sum(MUL[a][b] for a,b in zip(row,vector)) for row in matrix)


def rank(matrix):
    a=[list(row) for row in matrix]
    if not a:return 0
    answer=0
    for column in range(len(a[0])):
        pivot=next((i for i in range(answer,len(a)) if a[i][column]),None)
        if pivot is None:continue
        a[answer],a[pivot]=a[pivot],a[answer]
        scale=INV[a[answer][column]]
        a[answer]=[MUL[scale][x] for x in a[answer]]
        for i in range(answer+1,len(a)):
            scale=a[i][column]
            if scale:a[i]=[x^MUL[scale][y] for x,y in zip(a[i],a[answer])]
        answer+=1
        if answer==len(a):break
    return answer


def determinant(matrix):
    size=len(matrix)
    if any(len(row)!=size for row in matrix):raise ValueError('square matrix required')
    a=[list(row) for row in matrix]
    result=1
    for column in range(size):
        pivot=next((i for i in range(column,size) if a[i][column]),None)
        if pivot is None:return 0
        a[column],a[pivot]=a[pivot],a[column]
        result=MUL[result][a[column][column]]
        inverse=INV[a[column][column]]
        for i in range(column+1,size):
            scale=MUL[a[i][column]][inverse]
            for j in range(column,size):a[i][j]^=MUL[scale][a[column][j]]
    return result


def cauchy_minor(width,state,phase,rows,columns):
    """Cauchy determinant formula, independent of Gaussian elimination."""
    labels,residues,derivative=geometry(width,state,phase)
    u=tuple(labels[i] for i in rows)
    v=tuple(residues[j] for j in columns)
    if len(u)!=len(v):raise ValueError('square minor required')
    numerator=product(a^b for a,b in combinations(u,2))
    numerator=MUL[numerator][product(a^b for a,b in combinations(v,2))]
    denominator=product(a^b for a in u for b in v)
    scales=product(MUL[product(a^b for b in residues)][INV[derivative]] for a in u)
    return MUL[MUL[numerator][scales]][INV[denominator]]


def returned_matrix(width,state,phase=0,eta=2):
    if eta not in range(2,256):raise ValueError('eta must be nonzero and not one')
    a,c=matrices(width,state,phase)
    gram=matmul(a,c)
    return tuple(tuple(int(i==j)^MUL[eta][gram[i][j]] for j in range(width)) for i in range(width))


def proved_return_minima(width,state):
    geometry(width,state)
    if width==2*state:
        return tuple(max(1,min(state+1,2*(state+1-w))) if w<=state else 1
                     for w in range(1,width+1))
    return (width-state+1,)+tuple(max(1,width-state+1-w) for w in range(2,width+1))


def packed_cost(width,state):
    geometry(width,state)
    maps=width*state//2
    diagonal=(width+state)//4
    boolean=width*state//4+(width+state)//4
    return dict(fixed_map_gfni512=maps,input_and_state_diagonal_gfni512=diagonal,
                packet_gfni512=maps+diagonal,packet_xor_or_xor3=boolean,
                reused_record_broadcasts=width+state,
                gfni512_per_output_record=str(F(maps+diagonal,width)),
                eta_folded_into_constants=True,
                scope='algebraic dense four-record packing; excludes loads/stores/materialization/spills')


# These are the actual AES-field subfields, not the low integer intervals.
GF4=(0,1,188,189)
GF16_REPRESENTATIVES=(0,92,176,236)
GF16_PHASES=(0,2,6,4)
CUBIC_NODES=(16,17,172,173)
CUBIC_RETURN_MINIMA=(8,6,5,4,3,2,2,1)
CUBIC_CONSTRUCTION='gf256-repaired-cubic-w8-s4-eta2-c3-regional-nonzero-diagonals'


def cubic_geometry(phase=0):
    if phase not in range(8):raise ValueError('cubic phase must be0..7')
    parent,half=divmod(phase,2)
    labels=tuple(GF16_PHASES[parent]^representative^x
                 for representative in GF16_REPRESENTATIVES[2*half:2*half+2]
                 for x in GF4)
    return labels,CUBIC_NODES


@lru_cache(maxsize=None)
def cubic_matrices(phase=0,injection_shift=3):
    """Preferred map uses A=cube+1 and C=cube+3; old C used shift0."""
    if injection_shift not in range(256):raise ValueError('GF256 shift required')
    labels,nodes=cubic_geometry(phase)
    a=tuple(tuple(power(u^v,3)^1 for v in nodes) for u in labels)
    c=tuple(tuple(power(u^v,3)^injection_shift for u in labels) for v in nodes)
    return a,c


def cubic_returned_matrix(phase=0):
    a,c=cubic_matrices(phase)
    product_ac=matmul(a,c)
    return tuple(tuple(int(i==j)^MUL[2][product_ac[i][j]] for j in range(8))
                 for i in range(8))


def zeta4(values):
    x0,x1,x2,x3=values
    lo,hi=x0^x1,x2^x3
    return lo^hi,x1^x3,hi,x3


def five_product_circulant(coefficients,values):
    """Four-record XOR convolution with coefficient sum one, five products."""
    one,g1,g2,g3=zeta4(coefficients)
    if one!=1:raise ValueError('five-product circuit needs coefficient sum one')
    z0,z1,z2,z3=zeta4(values)
    w0=z0
    w1=z1^MUL[g1][z0]
    w2=z2^MUL[g2][z0]
    w3=z3^MUL[g1][z2]^MUL[g2][z1]^MUL[g3][z0]
    return zeta4((w0,w1,w2,w3))


def transpose_four_records(rows):
    """Four ZMMs, four 128-bit records each; eight VSHUFI32X4 operations."""
    if len(rows)!=4 or any(len(row)!=4 for row in rows):raise ValueError('4x4 record array required')
    def shuffle(a,b,immediate):
        return (a[immediate&3],a[(immediate>>2)&3],b[(immediate>>4)&3],b[(immediate>>6)&3])
    a,b,c,d=rows
    ab0,ab1=shuffle(a,b,0x44),shuffle(a,b,0xee)
    cd0,cd1=shuffle(c,d,0x44),shuffle(c,d,0xee)
    return (shuffle(ab0,cd0,0x88),shuffle(ab0,cd0,0xdd),
            shuffle(ab1,cd1,0x88),shuffle(ab1,cd1,0xdd))


def cubic_packed_maps(values,state,phase=0,adjoint=False):
    """Model four independent 128-bit lanes in each of five GFNI512 sites.

    Forward lanes are C0*x0, C1*x1, A0*s, A1*s. Adjoint swaps A and C.
    Record packing and unpacking are exact transposes, not scalar broadcasts.
    The returned pair is (eight emitted records, four accumulated records).
    """
    if len(values)!=8 or len(state)!=4:raise ValueError('eight inputs, four state records required')
    a,c=cubic_matrices(phase)
    emission,injection=(tuple(zip(*c)),tuple(zip(*a))) if adjoint else (a,c)
    kernels=(tuple(injection[j][0] for j in range(4)),
             tuple(injection[j][4] for j in range(4)),emission[0],emission[4])
    lane_inputs=(tuple(values[:4]),tuple(values[4:]),tuple(state),tuple(state))
    packed_inputs=transpose_four_records(lane_inputs)
    transformed=[]
    for lane in range(4):
        transformed.append(five_product_circulant(kernels[lane],
                           tuple(packed_inputs[j][lane] for j in range(4))))
    # A native schedule transposes these four output ZMMs back to records.
    packed_outputs=tuple(zip(*transformed))
    unpacked=transpose_four_records(packed_outputs)
    return unpacked[2]+unpacked[3],tuple(x^y for x,y in zip(unpacked[0],unpacked[1]))


def cubic_step(values,state,input_diagonal,state_diagonal,phase=0,adjoint=False):
    if adjoint:
        emitted,injected=cubic_packed_maps(values,state,phase,adjoint=True)
        x=tuple(MUL[d][v^e] for d,v,e in zip(input_diagonal,values,emitted))
        s=tuple(MUL[d][h^MUL[2][v]] for d,h,v in zip(state_diagonal,state,injected))
        return x,s
    x=tuple(MUL[d][v] for d,v in zip(input_diagonal,values))
    r=tuple(MUL[d][v] for d,v in zip(state_diagonal,state))
    emitted,injected=cubic_packed_maps(x,tuple(MUL[2][v] for v in r),phase)
    return tuple(v^e for v,e in zip(x,emitted)),tuple(v^e for v,e in zip(r,injected))


def cubic_packed_cost():
    return dict(fixed_map_gfni512=5,input_diagonal_gfni512=2,
                state_diagonal_gfni512=1,eta_gfni512=1,packet_gfni512=9,
                map_xor_or_xor3=12,packet_xor_or_xor3=15,
                record_transpose_shuffles_at_most=16,
                fixed_coefficient_zmm_tables_per_phase=3,
                gfni512_per_output_record='9/8',
                scope='explicit arithmetic and record-transpose schedule; excludes memory, spills, and matrix expansion')


def has_full_support_kernel(matrix,columns):
    """Exact when columns<256: a union of columns proper hyperplanes cannot cover."""
    if not 0<columns<256:raise ValueError('need0<columns<field size')
    if any(len(row)!=columns for row in matrix):raise ValueError('column mismatch')
    r=rank(matrix)
    return r<columns and all(rank(tuple(tuple(v for j,v in enumerate(row) if j!=drop)
                                           for row in matrix))==r for drop in range(columns))


def exact_return_minima(matrix):
    """Exhaust every support, using exact-support kernel feasibility, not samples."""
    width=len(matrix)
    if not 0<width<=8 or any(len(row)!=width for row in matrix):
        raise ValueError('bounded exhaustive helper requires square width<=8')
    result=[]
    for weight in range(1,width+1):
        found=False
        for output_weight in range(1,width+1):
            for inputs in combinations(range(width),weight):
                for outputs in combinations(range(width),output_weight):
                    allowed=set(outputs)
                    restricted=tuple(tuple(matrix[i][j] for j in inputs)
                                     for i in range(width) if i not in allowed)
                    if has_full_support_kernel(restricted,weight):
                        found=True
                        break
                if found:break
            if found:break
        result.append(output_weight)
    return tuple(result)


def cubic_profile(include_matrices=False):
    result=dict(construction=CUBIC_CONSTRUCTION,field='AES GF256 polynomial0x11b',width=8,state=4,eta=2,
                phase_labels=[list(cubic_geometry(p)[0]) for p in range(8)],
                residues=list(CUBIC_NODES),
                normalization='A[u,j]=(u+v_j)^3+1; C[j,u]=(u+v_j)^3+3',injection_constant=3,
                return_weight_lower_bounds=list(CUBIC_RETURN_MINIMA),
                exact_return_minima=list(CUBIC_RETURN_MINIMA),
                square_minors_per_map_per_phase=494,packed_cost=cubic_packed_cost(),
                whole_code_certified=False)
    if include_matrices:result['phase_matrices']=[cubic_matrices(p) for p in range(8)]
    return result


def binary_rank(values):
    pivots={}
    for value in values:
        while value:
            bit=value.bit_length()-1
            if bit in pivots:value^=pivots[bit]
            else:
                pivots[bit]=value
                break
    return len(pivots)


def profile(width,state,include_matrices=False):
    labels,residues,derivative=geometry(width,state)
    kernels=[list(matrices(width,state,p)[0][first])
             for p in range(4) for first in range(0,width,state)]
    result=dict(field='AES GF256 polynomial0x11b',width=width,state=state,eta=2,
                phase_labels=[list(geometry(width,state,p)[0]) for p in range(4)],
                residues=list(residues),derivative=derivative,
                normalization='A[u,j]=Q(u)/(Q_prime*(u+v_j)); C=A^T',
                return_weight_lower_bounds=list(proved_return_minima(width,state)),
                return_event='nonzero refreshed state and zero successor; equivalently Cx!=0 and r=Cx',
                square_minors_per_phase=sum(comb(width,k)*comb(state,k) for k in range(1,state+1)),
                coefficient_blocks='state-by-state symmetric XOR-circulants, row sum one',
                kernels=kernels,kernel_GF2_ranks=[binary_rank(kernel) for kernel in kernels],
                packed_cost=packed_cost(width,state),whole_code_certified=False)
    if (width,state)==(8,4):result['structured_packed_cost']=cubic_packed_cost()
    if include_matrices:result['phase_matrices']=[matrices(width,state,p)[0] for p in range(4)]
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profile',type=int,nargs=2,action='append')
    parser.add_argument('--matrices',action='store_true')
    parser.add_argument('--cubic',action='store_true',help='preferred width8/state4 repair of the old GF4 circuit')
    args=parser.parse_args()
    if args.cubic:print(json.dumps(cubic_profile(args.matrices)),flush=True)
    else:
        for width,state in args.profile or [(8,4),(16,8),(32,8),(32,16)]:
            print(json.dumps(profile(width,state,args.matrices)),flush=True)
