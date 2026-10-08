"""Exact full-Cauchy width4/state4 maps for the speed candidate.

Field: GF256 with AES polynomial0x11b. State-row normalization of the
8-by-8 Cauchy matrix preserves the selected-output/all-return rank law.
The scalar input/output normalizations preserve Hamming supports.
"""
from hashlib import sha256
from itertools import combinations

from gf8_packet_maps import MUL,INV,matmul,matvec,rank,xor_sum
from gf8_mds_packet import identity,scale,add,all_minors


WIDTH=STATE=4
ETA=2
ROWS=tuple(range(8))
COLUMNS=tuple(range(128,136))
CAUCHY=tuple(tuple(INV[u^v] for v in COLUMNS) for u in ROWS)
B=tuple(row[:4] for row in CAUCHY[:4])
A0=tuple(row[4:] for row in CAUCHY[:4])
C0=tuple(row[:4] for row in CAUCHY[4:])
E=tuple(row[4:] for row in CAUCHY[4:])
E_SUM=xor_sum(E[0])
E_INVERSE=scale(E,INV[MUL[E_SUM][E_SUM]])
C_RAW=matmul(E_INVERSE,C0)
OUTPUT_SCALE=INV[xor_sum(A0[0])]
INPUT_SCALE=INV[xor_sum(C_RAW[0])]
F=scale(B,MUL[OUTPUT_SCALE][INPUT_SCALE])
A_EFFECTIVE=scale(A0,OUTPUT_SCALE)
C=scale(C_RAW,INPUT_SCALE)
A=scale(A_EFFECTIVE,INV[ETA])
AC=matmul(A_EFFECTIVE,C)
FACTORED_LAMBDA=MUL[F[0][0]][INV[AC[0][0]]]


def matrices():
    return F,A,C,ETA


def extended():
    return tuple(F[i]+A_EFFECTIVE[i] for i in range(4))+tuple(C[i]+identity(4)[i] for i in range(4))


def matrix_digest():
    return sha256(bytes(v for matrix in (F,A,C) for row in matrix for v in row)+bytes([ETA])).hexdigest()


def packet(values,state,dx,ds,adjoint=False):
    if any(len(v)!=4 for v in (values,state,dx,ds)) or not all(dx) or not all(ds):
        raise ValueError('four live diagonal entries and four coordinates required')
    t=extended()
    if adjoint:
        t=tuple(zip(*t))
        image=matvec(t,tuple(values)+tuple(state))
        return (tuple(MUL[a][b] for a,b in zip(image[:4],dx)),
                tuple(MUL[a][b] for a,b in zip(image[4:],ds)))
    source=tuple(MUL[a][b] for a,b in zip(values,dx))+tuple(MUL[a][b] for a,b in zip(state,ds))
    image=matvec(t,source)
    return tuple(image[:4]),tuple(image[4:])


def check():
    if matmul(E,E_INVERSE)!=identity(4):
        raise ArithmeticError('state normalization inverse mismatch')
    if F!=scale(AC,FACTORED_LAMBDA):
        raise ArithmeticError('feedthrough factorization mismatch')
    if any(xor_sum(row)!=1 for matrix in (A_EFFECTIVE,C) for row in matrix):
        raise ArithmeticError('state maps not normalized to unit sum')
    for matrix in (F,A_EFFECTIVE,C,add(F,AC)):
        if any(matrix[i][j]!=matrix[0][i^j] for i in range(4) for j in range(4)):
            raise ArithmeticError('four-point XOR-circulant mismatch')
        all_minors(matrix)
    t=extended()
    tested=0
    for outputs in range(16):
        selected=tuple(j for j in range(4) if outputs>>j&1)
        for support in range(256):
            cols=tuple(j for j in range(8) if support>>j&1)
            for returned in (False,True):
                rows=selected+(tuple(range(4,8)) if returned else ())
                actual=rank(tuple(tuple(t[i][j] for j in cols) for i in rows))
                if actual!=min(len(rows),len(cols)):
                    raise ArithmeticError('selected-output/all-return rank mismatch')
                tested+=1
    return dict(selected_rank_cases=tested,matrix_sha256=matrix_digest(),
                f_row=list(F[0]),effective_a_row=list(A_EFFECTIVE[0]),c_row=list(C[0]),
                eta=ETA,factored_lambda=FACTORED_LAMBDA)


if __name__=='__main__':
    import json
    print(json.dumps(check()))
