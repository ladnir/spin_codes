"""Whole-code directed certificate for the fixed seven-of-eight GF256 code.

Replays the exact packet identity and the three setup-failure components.
Distance refers to retained output symbols. The certificate is for the
padded rate3/4 code, hence also its shortened requested-K20 subcode.
No performance result is asserted by this helper.
"""
from __future__ import annotations
from fractions import Fraction as F
import json

from gf8_mds_packet import fullstate_matrices,fullstate_extended,identity,scale
from gf8_packet_maps import INV,MUL,matmul,determinant
from gf8_punctured_profile import metadata
from gf8_punctured_sparse_certificate import certify as sparse_certificate
from gf8_punctured_middle_certificate import certify as middle_certificate
from gf8_punctured_dense import certify as dense_certificate


MATRIX_SHA256='6e3cf3e5d6b5876838afc5b2d7230f931974f363f9f67b402843f677521450a1'


def validate_maps():
    """Check the exact Cauchy factorization behind all selected-rank claims."""
    g=metadata()
    if g['matrix_sha256']!=MATRIX_SHA256:
        raise ArithmeticError('frozen packet matrix digest changed')
    nodes_rows=tuple(range(16))
    nodes_columns=tuple(range(128,144))
    if len(set(nodes_rows+nodes_columns))!=32:
        raise ArithmeticError('Cauchy row and column nodes overlap')
    cauchy=tuple(tuple(INV[i^j] for j in nodes_columns) for i in nodes_rows)
    f0,a0,c0,eta=fullstate_matrices(False)
    e=tuple(row[8:] for row in cauchy[8:])
    e_inverse=scale(f0,INV[MUL[15][15]])
    if matmul(e,e_inverse)!=identity(8) or matmul(e_inverse,e)!=identity(8):
        raise ArithmeticError('bottom-state normalization is not invertible')
    left=tuple(identity(16)[:8])+tuple((0,)*8+row for row in e_inverse)
    if matmul(left,cauchy)!=fullstate_extended(False):
        raise ArithmeticError('packet is not the asserted normalized Cauchy matrix')
    raw=fullstate_extended(False)
    normalized=tuple(tuple(MUL[124 if i<8 else 1][MUL[v][162 if j<8 else 1]]
                           for j,v in enumerate(row)) for i,row in enumerate(raw))
    if normalized!=fullstate_extended(True):
        raise ArithmeticError('additional nonzero row/column scalings changed')
    f,a,c,normalized_eta=fullstate_matrices(True)
    if eta!=2 or normalized_eta!=2 or not determinant(f):
        raise ArithmeticError('dense full causal inverse hypothesis failed')
    kept=tuple(g['kept_outputs'])
    if kept!=tuple(range(7)) or not determinant(tuple(tuple(f[i][j] for j in kept) for i in kept)):
        raise ArithmeticError('retained-output causal pivot is singular')
    return dict(matrix_sha256=g['matrix_sha256'],cauchy_row_nodes=list(nodes_rows),
                cauchy_column_nodes=list(nodes_columns),exact_state_normalization=True,
                exact_nonzero_row_column_scalings=True,full_feedthrough_invertible=True,
                retained_output_pivot_invertible=True,
                selected_minor_justification='Cauchy determinant and exact checked block factorization')


def component_upper(receipt,kind):
    bits=receipt.get('margin_bits_lower')
    if not isinstance(bits,int) or bits<40:
        raise ValueError('invalid component margin')
    if kind in ('sparse','dense'):
        try:
            p,q=(int(value,16) for value in receipt['failure_upper_dyadic_hex'])
            upper=F(p,q)
        except (KeyError,ValueError,ZeroDivisionError) as error:
            raise ValueError('invalid exact component failure bound') from error
    else:
        try:
            upper=F(receipt['failure_upper'])
        except (KeyError,ValueError,ZeroDivisionError) as error:
            raise ValueError('invalid decimal component failure bound') from error
    if not 0<upper<=F(1,1<<bits):
        raise ValueError('component failure bound does not support its stated margin')
    return upper


def combine(sparse,middle,dense):
    expected=metadata()
    for receipt in (sparse,middle,dense):
        if receipt.get('complete') is not True or receipt.get('rigorous') is not True:
            raise ValueError('every component must be complete and directed')
        for key,value in expected.items():
            if receipt.get(key)!=value:
                raise ValueError(f'component metadata mismatch: {key}')
    if expected['matrix_sha256']!=MATRIX_SHA256:
        raise ValueError('frozen packet matrix digest changed')
    first_active,last_active=sparse['active_range']
    first_middle,last_middle=middle['weights']
    first_dense,last_dense=dense['weights']
    outer_distance=expected['outer'][0]-expected['outer'][1]+1
    if first_active!=1 or not 1<=last_active<=expected['rows']:
        raise ValueError('invalid sparse active-row cover')
    if not 1<=first_middle<=last_middle<first_dense<=last_dense==expected['input_n']:
        raise ValueError('invalid middle/dense inverse-input ordering')
    if last_middle+1!=first_dense:
        raise ValueError('middle/dense weight gap')
    if outer_distance*(last_active+1)<first_middle:
        raise ValueError('sparse/middle weight gap')
    lost=expected['input_n']-expected['output_n']
    if (dense.get('omitted_symbols')!=lost or dense.get('completion_factor_log2')!=8*lost
            or dense.get('dense_inner_hypothesis')!=
            'invertible full causal inner; independent nonzero input diagonals'):
        raise ValueError('dense completion factor or inverse hypothesis mismatch')
    if not (10*expected['output_threshold']<expected['output_n']
            <=10*(expected['output_threshold']+1)):
        raise ValueError('distance threshold is not the strict one-tenth cutoff')
    upper=sum((component_upper(receipt,kind) for receipt,kind in
               ((sparse,'sparse'),(middle,'middle'),(dense,'dense'))),F(0))
    bits=min(r['margin_bits_lower'] for r in (sparse,middle,dense))
    while upper>F(1,1<<bits):
        bits-=1
    if bits<40:
        raise ArithmeticError('whole-code failure margin below forty bits')
    return dict(**expected,complete=True,rigorous=True,whole_code=True,
                minimum_distance=expected['output_threshold']+1,
                margin_bits_lower=bits,failure_upper_dyadic_exponent=bits,
                failure_upper_dyadic_hex=['0x1',hex(1<<bits)],
                component_margins_bits=[r['margin_bits_lower'] for r in (sparse,middle,dense)],
                sparse_active_range=sparse['active_range'],middle_weights=middle['weights'],
                dense_weights=dense['weights'],minimum_outer_distance=outer_distance,
                first_remaining_outer_weight=outer_distance*(last_active+1),
                completion_factor_log2=8*lost,
                probability_space='once-sampled independent uniform regional permutations and nonzero GF256 diagonals',
                scalar_extension='same scalar code on all sixteen byte lanes; no sixteen-fold setup-failure union',
                certified_domain='padded code; shortened requested-dimension subcode inherits the bound',
                native_timing_claim=False)


def certify():
    maps=validate_maps()
    sparse=sparse_certificate()
    middle=middle_certificate()
    dense=dense_certificate()
    result=combine(sparse,middle,dense)
    result['local_map_checks']=maps
    result['component_partition_sha256']=[sparse['partition_sha256'],middle['partition_sha256'],dense['witness_sha256']]
    return result


if __name__=='__main__':
    print(json.dumps(certify()),flush=True)
