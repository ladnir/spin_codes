"""Join the independently directed components of the smaller speed code."""
from fractions import Fraction as F

from gf8_packet_maps import MUL,matmul,determinant
from gf8_mds_packet import identity
from gf8_punctured_certificate import component_upper
from gf8_speed_profile import metadata,EXPECTED_MATRIX_SHA256
from gf8_speed_maps import (check,ROWS,COLUMNS,CAUCHY,E,E_INVERSE,OUTPUT_SCALE,
                           INPUT_SCALE,extended,matrices)


def validate_maps():
    if len(set(ROWS+COLUMNS))!=16 or not OUTPUT_SCALE or not INPUT_SCALE:
        raise ArithmeticError('invalid Cauchy nodes or normalization')
    if matmul(E,E_INVERSE)!=identity(4) or matmul(E_INVERSE,E)!=identity(4):
        raise ArithmeticError('state normalization is not invertible')
    left=identity(8)[:4]+tuple((0,)*4+row for row in E_INVERSE)
    raw=matmul(left,CAUCHY)
    scaled=tuple(tuple(MUL[OUTPUT_SCALE if i<4 else 1][MUL[INPUT_SCALE if j<4 else 1][v]]
                       for j,v in enumerate(row)) for i,row in enumerate(raw))
    if scaled!=extended():
        raise ArithmeticError('Cauchy block factorization mismatch')
    f,_,_,eta=matrices()
    if eta!=2 or not determinant(f) or not determinant(tuple(row[:3] for row in f[:3])):
        raise ArithmeticError('full or retained causal pivot singular')
    result=check()
    if result['matrix_sha256']!=EXPECTED_MATRIX_SHA256:
        raise ArithmeticError('matrix digest changed')
    result.update(exact_cauchy_factorization=True,full_feedthrough_invertible=True,
                  retained_output_pivot_invertible=True,cauchy_row_nodes=list(ROWS),
                  cauchy_column_nodes=list(COLUMNS))
    return result


def combine(sparse,middle,dense):
    dimension=sparse['outer'][1]
    g=metadata(dimension,sparse['k_log'])
    checks=validate_maps()
    for receipt in (sparse,middle,dense):
        if receipt.get('complete') is not True or receipt.get('rigorous') is not True:
            raise ValueError('complete directed components required')
        for key,value in g.items():
            if receipt.get(key)!=value:
                raise ValueError(f'component metadata mismatch:{key}')
    a,b=sparse['active_range'];lo,hi=middle['weights'];tail,end=dense['weights']
    row_distance=65-dimension
    if not(a==1 and 1<=b<=g['rows'] and 1<=lo<=hi<tail<=end==g['input_n']
           and hi+1==tail and row_distance*(b+1)>=lo):
        raise ValueError('invalid component coverage or a gap')
    lost=g['input_n']-g['output_n']
    if (dense.get('omitted_symbols')!=lost or dense.get('completion_factor_log2')!=8*lost
            or dense.get('dense_inner_hypothesis')!=
            'invertible full causal inner; independent nonzero input diagonals'):
        raise ValueError('dense completion factor or inverse hypothesis mismatch')
    if not 10*g['output_threshold']<g['output_n']<=10*(g['output_threshold']+1):
        raise ValueError('strict symbol-distance threshold mismatch')
    upper=sum((component_upper(r,kind) for r,kind in
               ((sparse,'sparse'),(middle,'middle'),(dense,'dense'))),F(0))
    bits=min(r['margin_bits_lower'] for r in (sparse,middle,dense))
    while upper>F(1,1<<bits):
        bits-=1
    if bits<40:
        raise ArithmeticError('whole-code setup-failure margin below40bits')
    return dict(**g,complete=True,rigorous=True,whole_code=True,
        minimum_distance=g['output_threshold']+1,margin_bits_lower=bits,
        failure_upper_dyadic_hex=['0x1',hex(1<<bits)],
        component_margins_bits=[r['margin_bits_lower'] for r in (sparse,middle,dense)],
        sparse_active_range=[a,b],middle_weights=[lo,hi],dense_weights=[tail,end],
        minimum_outer_distance=row_distance,first_remaining_outer_weight=row_distance*(b+1),
        completion_factor_log2=8*lost,local_map_checks=checks,
        scalar_extension='same scalar code on all sixteen byte lanes; no sixteen-fold setup-failure union',
        probability_space='once-sampled independent uniform regional permutations and nonzero GF256 diagonals',
        certified_domain='padded code; shortened requested-dimension subcode inherits the bound',
        native_timing_claim=False)


def certify(dimension=40,k_log=20,progress=None):
    from gf8_speed_sparse_certificate import certify as sparse
    from gf8_speed_middle_certificate import certify as middle
    from gf8_speed_dense_certificate import certify as dense
    g=metadata(dimension,k_log)
    first=(3*g['input_n']+99)//100
    last=(first-1)//(65-dimension)
    a=sparse(last=last,dimension=dimension,k_log=k_log,progress=progress)
    b=middle(dimension,k_log)
    c=dense(dimension,k_log)
    result=combine(a,b,c)
    result['component_partition_sha256']=[a['partition_sha256'],b['partition_sha256'],c['witness_sha256']]
    return result


if __name__=='__main__':
    import json
    print(json.dumps(certify()),flush=True)
