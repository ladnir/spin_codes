"""Six local quiet-fiber counting comparisons, not a distance screen.

The new full-rank fiber bound retains exact nonzero RS row enumerators.
The old comparison keeps a fixed labelled A, excludes wholly zero rows,
and retains total pre-H weight >= a*(n-k+1); it relaxes the individual
row floors.  Both use the two-moment-zero branch only.  This helper builds
no regional operator, changes no witness and cannot certify distance.
"""
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from math import comb, log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_two_sum_closed_kernel import zero_two_moments


DIRECTORY = Path(__file__).resolve().parent
BANK = 'hd_two_sum_gf16_source_strata_anchors_float.json'
BANK_SHA = '5050033694bb6a1ef7852aa309ee3abca252c9740fc713a83aa8945a380fb1cf'


def row_weight_coefficients(q=256, length=16, dimension=8):
    """Exact nonnegative weight enumerator of NONZERO RS[n,k] words."""
    if not (type(q) is int and type(length) is int and type(dimension) is int
            and q >= 2 and not q & (q-1) and 1 <= dimension < length <= q):
        raise ValueError('Proper RS parameters over GF(2^b) required')
    distance = length-dimension+1
    coefficients = [0]*(length+1)
    for w in range(distance, length+1):
        prescribed = (q-1)*sum((-1)**j*comb(w-1, j)*q**(w-distance-j)
                               for j in range(w-distance+1))
        coefficients[w] = comb(length, w)*prescribed
    if min(coefficients) < 0 or sum(coefficients) != q**dimension-1:
        raise ArithmeticError('Exact MDS enumerator failed normalization')
    return tuple(coefficients)


def weighted_row_enumerator(markers, *, q=256, dimension=8):
    """Exact marked nonzero-row sum via elementary symmetric polynomials."""
    markers = tuple(Fraction(x) for x in markers)
    if any(x < 0 for x in markers):
        raise ValueError('Nonnegative markers required')
    length = len(markers)
    coefficients = row_weight_coefficients(q, length, dimension)
    elementary = [Fraction(1)]+[Fraction(0)]*length
    for count, marker in enumerate(markers, 1):
        for w in range(count, 0, -1):
            elementary[w] += marker*elementary[w-1]
    return sum((Fraction(coefficients[w], comb(length, w))*elementary[w]
                for w in range(1, length+1)), Fraction(0))


def log_polynomial(coefficients, log_marker):
    terms = [log2(c)+w*log_marker for w, c in enumerate(coefficients) if c]
    return float(np.logaddexp2.reduce(terms)) if terms else -np.inf


def fiber_log_upper(arity, log_marker, *, q=256, length=16, dimension=8,
                    completion_rows=5):
    """Uniform bound on good-setup zero fibers with every RS row nonzero."""
    k = arity-completion_rows
    if k <= 0 or completion_rows*dimension < 2*length:
        raise ValueError('Positive row-marginal order and full-rank completion required')
    coefficients = row_weight_coefficients(q, length, dimension)
    row_sum = log_polynomial(coefficients, log_marker*arity/k)
    return ((arity*dimension-2*length)*log2(q)
            +k*(row_sum-dimension*log2(q)))


def quiet_word_logs(arity, log_marker, alpha=1., *, q=256, length=16, dimension=8):
    """Positive column comparison with labelled row activity kept exactly.

    c previously seen rows, h new rows and j old rows give multiplicity
    C(a-c,h)C(c,j).  Each exact column support has Z2(h+j) assignments.
    There is no additional C(a,h+j).  All arithmetic after the exact
    integer counts is FLOAT and positive; no signed FLOAT inclusion sum.
    """
    row_weight_coefficients(q, length, dimension)
    if not (1 <= arity <= min(q, 16) and 0 < alpha <= 1):
        raise ValueError('Valid labelled row count and fractional exponent required')
    maximum = arity*length
    values = np.full((arity+1, maximum+1), -np.inf)
    values[0, 0] = 0.
    counts = [zero_two_moments(w, q) for w in range(arity+1)]
    transitions = []
    for c in range(arity+1):
        for h in range(arity-c+1):
            for j in range(c+1):
                w = h+j
                if counts[w]:
                    multiplicity = comb(arity-c, h)*comb(c, j)
                    transitions.append((c, c+h, w, log2(multiplicity)
                                        +alpha*log2(counts[w])+w*log_marker))
    for column in range(length):
        updated = np.full_like(values, -np.inf)
        end = arity*column+1
        for old, new, weight, value in transitions:
            target = updated[new, weight:weight+end]
            np.logaddexp2(target, values[old, :end]+value, out=target)
        values = updated
    floor = arity*(length-dimension+1)
    normalizer = alpha*(length-dimension)*arity*log2(q-1)
    return dict(all_rows_seen=True, individual_row_floors=False,
                aggregate_floor=floor,
                log2_upper=float(np.logaddexp2.reduce(values[arity, floor:]))-normalizer,
                log2_upper_without_aggregate_floor=float(np.logaddexp2.reduce(values[arity]))-normalizer)


def compare(arity, activity_log, alpha):
    """Compare one fixed activity marker under both fractional orderings.

    To keep B^weight AFTER alpha, the unpowered fiber uses t=B^(1/alpha).
    This defines an activity-count diagnostic only: it does not assert
    that a fractional regional envelope can be rooted to recover R1.
    """
    alpha = float(alpha)
    unpowered_marker = activity_log/alpha
    fiber = alpha*fiber_log_upper(arity, unpowered_marker)
    mean = quiet_word_logs(arity, unpowered_marker, 1.)
    column_fractional = quiet_word_logs(arity, activity_log, alpha)
    return dict(arity=arity, log2_activity_marker=float(activity_log), alpha=alpha,
                unpowered_log2_marker=unpowered_marker,
                new_fiber_then_alpha_log2=fiber,
                old_mean_then_alpha_log2=alpha*mean['log2_upper'],
                old_column_alpha_log2=column_fractional['log2_upper'],
                new_minus_old_mean_bits=fiber-alpha*mean['log2_upper'],
                new_minus_old_column_alpha_bits=fiber-column_fractional['log2_upper'],
                new_all_rows_nonzero_rs=True, old_all_rows_seen=True,
                old_aggregate_floor=mean['aggregate_floor'], old_individual_row_floors=False,
                compared_scope='fixed A, all sixteen column moments zero',
                regional_envelope_recovered=False, distance_certificate=False)


def run():
    data = (DIRECTORY/BANK).read_bytes()
    if sha256(data).hexdigest() != BANK_SHA:
        raise ValueError('Frozen activity source changed')
    banks = json.loads(data)['banks']
    expected = ((18, '691/700'), (24, '172/175'), (32, '171/175'))
    if tuple((b['active_groups'], b['mark']) for b in banks) != expected or any(
            b['alpha'] != '1/5' or len(b['log2_slopes']) != 4 for b in banks):
        raise ValueError('Only the three frozen activity slopes are allowed')
    started = perf_counter()
    rows = [dict(active_groups=b['active_groups'], mark=b['mark'],
                 **compare(a, b['log2_slopes'][0], Fraction(b['alpha'])))
            for b in banks for a in (8, 9)]
    return dict(source= BANK, source_sha256=BANK_SHA, requested_controls=6,
                regional_banks_built=0, fitted_parameters=0, directed=False,
                whole_code_certificate=False, rows=rows, elapsed_seconds=perf_counter()-started)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a frozen local control')
    result = run()
    result['implementation_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as output:
        json.dump(result, output, indent=2, allow_nan=False)
        output.write('\n')
    print(json.dumps(result, allow_nan=False))
