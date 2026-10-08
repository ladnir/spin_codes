"""Replay and join the complete higher-rate K20 GRS[128,87] certificate.

The CLI recomputes each component. Geometry and integer coverage must match,
then exact Fractions sum their outward upper bounds. Fixed source-coordinate
permutations and arbitrary source shortening preserve the resulting distance.
"""
from decimal import Context, Decimal, ROUND_CEILING
from fractions import Fraction as F
import json

from field128_sparse_certificate import certify as sparse_certificate
from field128_density_certificate import certify as middle_certificate
from field128_dense_certificate import certify as dense_certificate
from field256_whole_certificate import margin_bits


def combine(sparse,middle,dense):
    for r in (sparse,middle,dense):
        if (r['outer']!=[128,87] or r['q']!=256 or r['width']!=16 or
            r['state']!=8 or r['k_log']!=20 or r['rows']!=12096):
            raise ValueError('component construction mismatch')
    for r in (sparse,middle):
        if not r['rigorous'] or r['k']!=1052352 or r['n']!=1548288 or r['output_threshold']!=154828:
            raise ValueError('inconsistent geometry or nonrigorous component')
    if not sparse['complete_sparse_range'] or sparse['active_range']!=[1,1105]:
        raise ValueError('incomplete active-row cover')
    if (middle['minimum_input_weight']!=46449 or middle['maximum_input_weight']!=371589 or
        middle['covered_weights']!=325141):
        raise ValueError('incomplete middle-weight cover')
    if (not dense['dense_rigorous'] or not dense['cutoff_cover_complete'] or
        dense['padded_source']!=1052352 or dense['output_length']!=1548288 or
        dense['forbidden_output_threshold']!=154828 or dense['dense_weights']!=[371590,1548288] or
        dense['dense_margin_bits_lower']!=13690 or dense['failure_upper_power_of_two']!=-13690):
        raise ValueError('incomplete dense-weight cover')
    if 42*(sparse['active_range'][1]+1)<middle['minimum_input_weight']:
        raise ValueError('gap between active-row and weight coverage')
    if middle['maximum_input_weight']+1!=dense['dense_weights'][0]:
        raise ValueError('gap between middle and dense weights')
    a,b=map(int,sparse['failure_upper_dyadic'])
    bounds=(F(a,b),F(middle['failure_upper']),F(1,2**13690))
    if any(not 0<x<1 for x in bounds):
        raise ArithmeticError('invalid failure upper bound')
    upper=sum(bounds,F(0));bits=margin_bits(upper)
    if bits<40:
        raise ArithmeticError('combined union fails forty bits')
    display=Context(prec=80,rounding=ROUND_CEILING).divide(Decimal(upper.numerator),Decimal(upper.denominator))
    return dict(outer=[128,87],q=256,width=16,state=8,k_log=20,
                requested_k=2**20,padded_k=1052352,n=1548288,rows=12096,
                design_rate='87/128',effective_requested_rate='128/189',
                relative_distance='1/10',minimum_distance=154829,
                whole_code_margin_bits_lower=bits,failure_upper=str(display),
                component_margin_bits_lower=[margin_bits(x) for x in bounds],
                sparse_active_range=[1,1105],sparse_intervals=sparse['interval_count'],
                sparse_partition_sha256=sparse['partition_sha256'],
                middle_weights=[46449,371589],middle_cells=middle['checked_cells'],
                middle_partition_sha256=middle['partition_sha256'],dense_weights=[371590,1548288],
                rigorous=True,whole_code=True,binary_field_degrees=[8,32,64,128],
                extension_degrees_over_gf256=[1,4,8,16],
                scope='ideal independent regional permutations and GF256 nonzero diagonals',
                source_order_invariance='fixed source-coordinate permutations and arbitrary shortening')


def certify(progress=None):
    sparse=sparse_certificate()
    if progress:progress('sparse',sparse['certified_integer_bits'])
    middle=middle_certificate()
    if progress:progress('middle',middle['class_margin_bits_lower'])
    dense=dense_certificate()
    if progress:progress('dense',dense['dense_margin_bits_lower'])
    return combine(sparse,middle,dense)


if __name__=='__main__':
    print(json.dumps(certify(lambda name,bits:print(json.dumps(dict(component=name,
                    margin_bits_lower=bits)),flush=True))),flush=True)
