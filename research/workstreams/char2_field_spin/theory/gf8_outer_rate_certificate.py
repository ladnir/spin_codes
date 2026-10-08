"""Join independently verified GF256 rate-sweep class certificates.

Only geometry, exact probability arithmetic and coverage are handled here.
Calling combine does not verify an arbitrary external receipt's mathematics;
the intended inputs are the matching local directed/exact replay functions.
Frozen dimension-87 certificates are not modified.
"""
from decimal import Context, Decimal, ROUND_CEILING
from fractions import Fraction as F

from field256_whole_certificate import margin_bits


CONSTRUCTION = 'symmetric-cauchy-w16-s8-eta2-regional-nonzero-diagonals'


def metadata(k_log, dimension):
    if (not isinstance(k_log, int) or not 1 <= k_log <= 30
            or not isinstance(dimension, int) or not 1 <= dimension < 128):
        raise ValueError('invalid size or outer dimension')
    requested = 1 << k_log
    rows = ((requested + 64*dimension-1)//(64*dimension))*64
    n = 128*rows
    return dict(outer=[128, dimension], q=256, width=16, state=8,
                k_log=k_log, requested_k=requested, padded_k=dimension*rows,
                rows=rows, n=n, output_threshold=(n+9)//10-1,
                construction=CONSTRUCTION)


def probability(receipt):
    if 'failure_upper_dyadic' in receipt:
        numerator, denominator = map(int, receipt['failure_upper_dyadic'])
        value = F(numerator, denominator)
    elif 'failure_upper_dyadic_exponent' in receipt:
        exponent = receipt['failure_upper_dyadic_exponent']
        if not isinstance(exponent, int) or exponent <= 0:
            raise ValueError('invalid dyadic exponent')
        value = F(1, 2**exponent)
    else:
        # A directed upper decimal is an exact rational, not a binary float.
        endpoint = receipt['failure_upper']
        if not isinstance(endpoint, str):
            raise ValueError('directed endpoint must be exact text')
        value = F(endpoint)
    if not 0 < value < 1:
        raise ValueError('component probability must be strictly between zero and one')
    return value


def combine(k_log, dimension, sparse, middle, dense, target_bits=40):
    if not isinstance(target_bits, int) or target_bits < 40:
        raise ValueError('rate sweep requires at least forty bits')
    expected = metadata(k_log, dimension)
    for receipt in (sparse, middle, dense):
        if (any(receipt.get(key) != value for key, value in expected.items())
                or receipt.get('complete') is not True
                or receipt.get('rigorous') is not True):
            raise ValueError('component ensemble, geometry or rigor mismatch')
    first_active, last_active = sparse['active_range']
    middle_first, middle_last = middle['weights']
    dense_first, dense_last = dense['weights']
    if (first_active != 1 or not 1 <= last_active <= expected['rows']
            or not 1 <= middle_first <= middle_last < expected['n']
            or not 1 <= dense_first <= dense_last == expected['n']
            or (129-dimension)*(last_active+1) < middle_first
            or middle_last+1 != dense_first):
        raise ValueError('sparse/middle/dense cover has a gap or invalid endpoints')
    bounds = tuple(probability(r) for r in (sparse, middle, dense))
    total = sum(bounds, F(0))
    if total > F(1, 2**target_bits):
        raise ArithmeticError('whole-code probability sum misses target')
    display = Context(prec=80, rounding=ROUND_CEILING).divide(
        Decimal(total.numerator), Decimal(total.denominator))
    return dict(**expected, design_rate=str(F(dimension, 128)),
                effective_requested_rate=str(F(expected['requested_k'], expected['n'])),
                relative_distance='1/10', minimum_distance=expected['output_threshold']+1,
                whole_code_margin_bits_lower=margin_bits(total),
                failure_upper=str(display),
                component_margin_bits_lower=[margin_bits(b) for b in bounds],
                sparse=sparse, middle=middle, dense=dense, rigorous=True, whole_code=True,
                scope='ideal independent regional permutations and GF256 nonzero diagonals')


def certify_k18_dimension88(progress=None):
    """Replay every component of the newly closed, naturally padded profile."""
    from gf8_outer_rate_sparse import certify_cover
    from gf8_outer_rate_middle import certify as certify_middle
    from gf8_outer_rate_dense import certify as certify_dense
    sparse = certify_cover(375, 18, 88, precision=4096, cell_bits=64)
    if progress:
        progress('sparse', sparse)
    middle = certify_middle(18, 88, first=F(1, 25))
    if progress:
        progress('middle', middle)
    dense = certify_dense(88, 18)
    if progress:
        progress('dense', dense)
    return combine(18, 88, sparse, middle, dense)


if __name__ == '__main__':
    import json
    def show_progress(name, receipt):
        print(json.dumps(dict(component=name,
            margin_bits_lower=margin_bits(probability(receipt)),
            covered_range=receipt.get('active_range', receipt.get('weights')))), flush=True)
    result = certify_k18_dimension88(show_progress)
    # Keep routine terminal output compact; the return value retains components.
    summary = {key:value for key,value in result.items() if key not in ('sparse','middle','dense')}
    print(json.dumps(summary), flush=True)
