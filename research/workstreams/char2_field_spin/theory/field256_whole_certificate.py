"""Replay and join the three K20 GRS256/176 distance bounds.

All numerical acceptance is directed or exact. This certifies the ideal
width16/state8 setup ensemble, not a pseudorandom implementation or timing.
Each component is replayed; no saved result is accepted by the CLI.
"""

from decimal import Context, Decimal, ROUND_CEILING
from fractions import Fraction as F
import json

from field256_coverage_audit import audit as dense_certificate
from state8_density_certificate import certify as middle_certificate
from state8_sparse_certificate import certify as sparse_certificate


def margin_bits(upper):
    if not 0 < upper < 1:
        raise ArithmeticError('expected a positive probability bound below one')
    bits = upper.denominator.bit_length()-upper.numerator.bit_length()
    while upper > F(2)**(-bits):
        bits -= 1
    while upper <= F(2)**(-bits-1):
        bits += 1
    return bits


def combine(sparse, middle, dense):
    """Check geometry/coverage and sum successful component upper bounds."""
    for receipt in (sparse,middle,dense):
        if (receipt['outer'] != [256,176] or receipt['q'] != 256
            or receipt['width'] != 16 or receipt['state'] != 8
            or receipt['k_log'] != 20 or receipt['rows'] != 6016):
            raise ValueError('component construction mismatch')
    for receipt in (sparse,middle):
        if (not receipt['rigorous'] or receipt['k'] != 1058816
            or receipt['n'] != 1540096 or receipt['output_threshold'] != 154009):
            raise ValueError('component is not a matching rigorous replay')
    if (not sparse['complete_sparse_range'] or sparse['active_range'] != [1,600]
        or sparse['interval_count'] != 80):
        raise ValueError('incomplete sparse active-row cover')
    if (middle['minimum_input_weight'] != 46203
        or middle['maximum_input_weight'] != 369623
        or middle['covered_weights'] != 323421):
        raise ValueError('incomplete intermediate weight cover')
    if (not dense['dense_rigorous'] or not dense['cutoff_cover_complete']
        or dense['padded_source'] != 1058816 or dense['output_length'] != 1540096
        or dense['forbidden_output_threshold'] != 154009
        or dense['dense_weights'] != [369624,1540096]
        or dense['dense_margin_bits_lower'] != 12072):
        raise ValueError('incomplete dense weight cover')
    # For A>=601 the MDS constituent gives W>=81A. Otherwise sparse applies.
    if 81*(sparse['active_range'][1]+1) < middle['minimum_input_weight']:
        raise ValueError('gap between sparse and weight-based covers')
    if middle['maximum_input_weight']+1 != dense['dense_weights'][0]:
        raise ValueError('gap between intermediate and dense covers')
    numerator,denominator = map(int,sparse['failure_upper_dyadic'])
    bounds = (F(numerator,denominator),F(middle['failure_upper']),F(1,2**12072))
    if any(not 0 < bound < 1 for bound in bounds):
        raise ArithmeticError('invalid component bound')
    upper = sum(bounds,F(0))
    bits = margin_bits(upper)
    if bits < 40:
        raise ArithmeticError('whole-code union fails forty-bit target')
    upward = Context(prec=80,rounding=ROUND_CEILING)
    display = upward.divide(Decimal(upper.numerator),Decimal(upper.denominator))
    return dict(outer=[256,176],q=256,width=16,state=8,k_log=20,
                requested_k=2**20,padded_k=1058816,n=1540096,rows=6016,
                design_rate='11/16',effective_requested_rate=str(F(2**20,1540096)),
                relative_distance='1/10',minimum_distance=154010,
                whole_code_margin_bits_lower=bits,failure_upper=str(display),
                component_margin_bits_lower=[margin_bits(bound) for bound in bounds],
                sparse_active_range=[1,600],middle_weights=[46203,369623],
                dense_weights=[369624,1540096],middle_cells=middle['checked_cells'],
                middle_partition_sha256=middle['partition_sha256'],
                rigorous=True,whole_code=True,scalar_extension_degrees=[8,32,64,128],
                scope='ideal independent regional permutations and GF256 nonzero diagonals',
                native_performance_measured=False)


def certify(progress=None):
    sparse = sparse_certificate()
    if progress:
        progress('sparse',sparse['certified_integer_bits'])
    middle = middle_certificate()
    if progress:
        progress('middle',middle['class_margin_bits_lower'])
    dense = dense_certificate()
    if progress:
        progress('dense',dense['dense_margin_bits_lower'])
    return combine(sparse,middle,dense)


if __name__ == '__main__':
    print(json.dumps(certify(lambda name,bits: print(
        json.dumps(dict(component=name,margin_bits_lower=bits)),flush=True))),flush=True)
