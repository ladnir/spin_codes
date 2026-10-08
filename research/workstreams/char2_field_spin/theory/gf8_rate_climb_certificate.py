"""Join independently replayed components for a new higher-rate code.

No frozen certificate is retagged. Every receipt must carry the exact new
geometry and probability-space metadata, and all probability bounds are
added as exact rationals before making a whole-code claim.
"""
from fractions import Fraction as F

from gf8_mds_packet import fullstate_matrices
from gf8_packet_maps import determinant
from gf8_punctured_certificate import validate_maps as validate_frozen_maps,component_upper
from gf8_rate_climb_profile import metadata as fixed_metadata
from gf8_periodic_climb_profile import metadata as periodic_metadata


def combine(sparse,middle,dense):
    dimension=sparse['outer'][1]
    if 'output_pattern' in sparse:
        pattern=tuple(sparse['output_pattern'])
        g=periodic_metadata(dimension,pattern,sparse['k_log'])
    else:
        pattern=(sparse['outputs'],)
        g=fixed_metadata(dimension,pattern[0],sparse['k_log'])
    checks=validate_frozen_maps()
    f,_,_,_=fullstate_matrices(True)
    for outputs in set(pattern):
        if not determinant(tuple(tuple(f[i][j] for j in range(outputs)) for i in range(outputs))):
            raise ArithmeticError('retained-output causal pivot failed')
    for receipt in (sparse,middle,dense):
        if receipt.get('complete') is not True or receipt.get('rigorous') is not True:
            raise ValueError('complete directed component required')
        for key,value in g.items():
            if receipt.get(key)!=value:raise ValueError(f'component metadata mismatch: {key}')
    a,b=sparse['active_range'];lo,hi=middle['weights'];tail,end=dense['weights']
    row_distance=129-dimension
    if not (a==1 and 1<=b<=g['rows'] and 1<=lo<=hi<tail<=end==g['input_n']
            and hi+1==tail and row_distance*(b+1)>=lo):
        raise ValueError('component coverage gap or invalid endpoint')
    lost=g['input_n']-g['output_n']
    if dense.get('omitted_symbols')!=lost or dense.get('completion_factor_log2')!=8*lost:
        raise ValueError('missing omitted-output completion factor')
    if dense.get('dense_inner_hypothesis')!='invertible full causal inner; independent nonzero input diagonals':
        raise ValueError('dense inner hypothesis mismatch')
    if not 10*g['output_threshold']<g['output_n']<=10*(g['output_threshold']+1):
        raise ValueError('distance threshold mismatch')
    upper=sum((component_upper(r,kind) for r,kind in
               ((sparse,'sparse'),(middle,'middle'),(dense,'dense'))),F(0))
    bits=min(r['margin_bits_lower'] for r in (sparse,middle,dense))
    while upper>F(1,1<<bits):bits-=1
    if bits<40:raise ArithmeticError('whole-code margin below40bits')
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


def certify_fixed(dimension=80,outputs=6):
    from gf8_rate_climb_sparse import certify as sparse
    from gf8_rate_climb_middle_certificate import certify as middle
    from gf8_rate_climb_dense import certify as dense
    joins={(80,6):F(4,5),(67,5):F(9,10)}
    if (dimension,outputs) not in joins:
        raise ValueError('no frozen middle witnesses for this fixed profile')
    g=fixed_metadata(dimension,outputs)
    first=(3*g['input_n']+99)//100
    active=(first-1)//(129-dimension)
    join=joins[dimension,outputs]
    a=sparse(last=active,dimension=dimension,outputs=outputs)
    b=middle(dimension,outputs,F(3,100),join)
    c=dense(dimension,outputs,join)
    result=combine(a,b,c)
    result['component_partition_sha256']=[a['partition_sha256'],b['partition_sha256'],c['witness_sha256']]
    return result


def certify_periodic(dimension=74,pattern=(5,6)):
    from gf8_periodic_climb_sparse import certify as sparse
    from gf8_periodic_climb_middle_certificate import certify as middle
    from gf8_periodic_climb_dense import certify as coarse_dense
    from gf8_exact_spectrum_dense import certify_profile as exact_dense
    pattern=tuple(pattern)
    joins={(84,(6,6,6,7)):F(17,20),(74,(5,6)):F(9,10)}
    if (dimension,pattern) not in joins:
        raise ValueError('no frozen middle witnesses for this periodic profile')
    g=periodic_metadata(dimension,pattern)
    first=(3*g['input_n']+99)//100
    active=(first-1)//(129-dimension)
    join=joins[dimension,pattern]
    a=sparse(last=active,dimension=dimension,pattern=pattern)
    b=middle(dimension,pattern,F(3,100),join)
    c=(exact_dense(g,join) if dimension==74
       else coarse_dense(dimension,pattern,join))
    result=combine(a,b,c)
    result['component_partition_sha256']=[a['partition_sha256'],b['partition_sha256'],c['witness_sha256']]
    if 'dense_method' in c:
        result['dense_method']=c['dense_method']
    return result


if __name__=='__main__':
    import argparse,json
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profile',choices=('keep6','keep5','periodic84','periodic74'),default='keep6')
    args=parser.parse_args()
    if args.profile in ('keep6','keep5'):
        result=certify_fixed(*((80,6) if args.profile=='keep6' else (67,5)))
    else:
        result=certify_periodic(*((84,(6,6,6,7)) if args.profile=='periodic84' else (74,(5,6))))
    print(json.dumps(result),flush=True)
