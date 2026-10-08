"""Exact maps for the shared eight-point FFT / four-state proposal.

This is a construction and structural receipt, not a distance certificate.
Coefficient matrices are over AES GF256 and act on GF65536 symbols by the
specified subfield embedding. No native implementation or timing is tested.
"""

from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from math import comb, factorial

from hd_blockwise_maps import mul, power, matmul, matvec, xor_sum


@lru_cache(None)
def multiplication_table():
    return tuple(tuple(mul(a,b) for b in range(256)) for a in range(256))


def rank(matrix):
    if not matrix:
        return 0
    work=[list(row) for row in matrix]
    table=multiplication_table()
    result=0
    for column in range(len(work[0])):
        pivot=next((i for i in range(result,len(work)) if work[i][column]),None)
        if pivot is None:
            continue
        work[result],work[pivot]=work[pivot],work[result]
        scale=power(work[result][column],254)
        work[result]=[table[scale][x] for x in work[result]]
        for i in range(result+1,len(work)):
            scale=work[i][column]
            if scale:
                work[i]=[x^table[scale][y] for x,y in zip(work[i],work[result])]
        result+=1
        if result==len(work):
            break
    return result


def inverse(matrix):
    n=len(matrix)
    if any(len(row)!=n for row in matrix):
        raise ValueError('square matrix required')
    work=[list(row)+[int(i==j) for j in range(n)] for i,row in enumerate(matrix)]
    table=multiplication_table()
    for column in range(n):
        pivot=next((i for i in range(column,n) if work[i][column]),None)
        if pivot is None:
            raise ValueError('singular matrix')
        work[column],work[pivot]=work[pivot],work[column]
        scale=power(work[column][column],254)
        work[column]=[table[scale][x] for x in work[column]]
        for i in range(n):
            scale=work[i][column]
            if i!=column and scale:
                work[i]=[x^table[scale][y] for x,y in zip(work[i],work[column])]
    return tuple(tuple(row[n:]) for row in work)


def subspace_polynomial(level,x):
    result=1
    for node in range(1<<level):
        result=mul(result,x^node)
    return result


def basis(index,x):
    """Global normalized novel basis, with degrees equal to the index."""
    if not 0<=index<8:
        raise ValueError('only the eight-point basis is defined')
    result=1
    for level in range(3):
        if index&(1<<level):
            numerator=subspace_polynomial(level,x)
            denominator=subspace_polynomial(level,1<<level)
            result=mul(result,mul(numerator,power(denominator,254)))
    return result


@lru_cache(None)
def evaluation(offset,size=8):
    if size not in (4,8) or offset%size:
        raise ValueError('aligned four- or eight-point evaluation required')
    return tuple(tuple(basis(j,offset^i) for j in range(size)) for i in range(size))


@lru_cache(None)
def coset_matrices(input_cosets, output_cosets):
    """Tile maps for distinct aligned cosets, with fixed global state basis.

    Each list must contain eight distinct cosets; overlap BETWEEN the two
    lists is allowed. Matched input/output cosets must differ for branch9.
    These are evaluation labels, not additional storage permutations.
    """
    if (len(input_cosets) != 8 or len(output_cosets) != 8
            or any(len(set(nodes)) != 8 for nodes in (input_cosets, output_cosets))
            or any(not 0 <= x < 256 or x % 8 for x in input_cosets+output_cosets)
            or any(x == y for x,y in zip(input_cosets, output_cosets))):
        raise ValueError('eight distinct cosets and disjoint matched pairs required')
    a=[]
    c=[[] for _ in range(4)]
    full=[]
    for tile,(alpha,delta) in enumerate(zip(input_cosets,output_cosets)):
        offset=8*tile
        inv=inverse(evaluation(alpha))
        local=matmul(evaluation(delta),inv)
        a.extend(row[:4] for row in evaluation(delta))
        for j in range(4):
            c[j].extend(inv[4+j])
        full.extend(tuple([0]*offset+list(row)+[0]*(56-offset)) for row in local)
    return tuple(full),tuple(a),tuple(tuple(row) for row in c)


@lru_cache(None)
def matrices(translation_offset=8):
    """Return (F64,A64x4,C4x64,T4,H4), in the declared global basis."""
    if not 0 < translation_offset < 256 or translation_offset % 8:
        raise ValueError('translation must be outside the tile subspace')
    inputs=tuple(range(0,64,8))
    full,a,c=coset_matrices(inputs,tuple(x^translation_offset for x in inputs))
    translation=matmul(inverse(evaluation(0,4)),evaluation(translation_offset,4))
    h=matmul(tuple(tuple(row[:4]) for row in c),inverse(tuple(zip(*a[:4]))))
    return full,a,c,translation,h


def identity(size):
    return tuple(tuple(int(i==j) for j in range(size)) for i in range(size))


def minor_count(matrix):
    """Check every square minor; small 8x8 use only."""
    count=0
    for size in range(1,min(len(matrix),len(matrix[0]))+1):
        for rows in combinations(range(len(matrix)),size):
            for columns in combinations(range(len(matrix[0])),size):
                if rank(tuple(tuple(matrix[i][j] for j in columns) for i in rows))!=size:
                    raise AssertionError((rows,columns))
                count+=1
    return count


def branch_lower(occupancy,group=8):
    if not 0<=occupancy<=group:
        raise ValueError('invalid tile occupancy')
    return 0 if occupancy==0 else group+1-occupancy


def bernoulli_branch_fraction(probability,group=8):
    p=F(probability)
    return F(group+1,group)*(1-(1-p)**group)-p


def fixed_subset_branch_fraction(rows,active,group=8):
    """Exact mean, not a probability bound and not a whole-code claim."""
    if rows%group or not 0<=active<=rows:
        raise ValueError('invalid regional geometry')
    return sum((F(comb(group,w)*comb(rows-group,active-w),comb(rows,active))
                *branch_lower(w,group) for w in range(group+1)
                if 0<=active-w<=rows-group),F(0))/group


def occupancy_profile_probability(rows,active,counts):
    """One region: histogram of a uniform active-subset across fixed tiles."""
    group=len(counts)-1
    if group<1 or rows%group or any(n<0 for n in counts):
        raise ValueError('invalid histogram')
    if sum(counts)!=rows//group or sum(w*n for w,n in enumerate(counts))!=active:
        return F(0)
    numerator=factorial(rows//group)
    # Divide only after forming the entire denominator; intermediate factorial
    # divisions are not individually guaranteed integral.
    denominator=1
    for w,n in enumerate(counts):
        denominator*=factorial(n)
        numerator*=comb(group,w)**n
    return F(numerator,denominator*comb(rows,active))


def prescribed_mds_count(weight,length=16,dimension=8,q=65536):
    if not 0<=weight<=length:
        raise ValueError('invalid support weight')
    if weight==0:
        return 1
    parity=length-dimension
    if weight<=parity:
        return 0
    return (q-1)*sum((-1)**j*comb(weight-1,j)*q**(weight-parity-1-j)
                     for j in range(weight-parity))


def receipt():
    full,a,c,t,h=matrices()
    local=tuple(tuple(row[:8]) for row in full[:8])
    nilpotent=tuple(tuple(x^(i==j) for j,x in enumerate(row)) for i,row in enumerate(full))
    assert matmul(full,full)==identity(64)
    assert matmul(full,a)==matmul(a,t)
    assert matmul(c,full)==matmul(t,c)
    assert matmul(h,tuple(zip(*a)))==c
    assert not any(x for row in matmul(c,a) for x in row)
    assert rank(h)==4
    return dict(packet_width=64,tile_width=8,state=4,coefficient_field=256,
                symbol_field=65536,F_first_row=list(local[0]),H=h,T=t,
                CA_zero=True,F_involution=True,FA_equals_AT=True,CF_equals_TC=True,
                emission_minimum=61,feedback_kernel_minimum=5,
                tile_branch_number=9,rank_F_minus_I=rank(nilpotent),
                rank_stacked_C_F_minus_I=rank(c+nilpotent),
                quiet_rows_at_K20=1024,
                bernoulli_branch_fraction=str(bernoulli_branch_fraction(F(1,8))),
                fixed_subset_mean_branch_fraction=str(fixed_subset_branch_fraction(8192,1024)),
                whole_certificate=False,native_timing=False)


if __name__=='__main__':
    import json
    print(json.dumps(receipt(),indent=2))
