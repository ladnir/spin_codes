"""Exact character-sum data for a small-field batched inner.

The code fixes C=A^T and uniform nonzero refreshed state R. Its input
measure assigns weight v to each nonzero symbol and weight 1 to zero.
Output weight is marked by z. Fourier orthogonality enforces C X=R,
giving the joint return moment without independent upper/lower envelopes.

Histogram construction uses exact integer counts. Fraction evaluation is
exact; floating evaluation is a research screen, NOT a certificate.
No native-encoder benchmark is performed here.
"""

import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import json
import math
from numbers import Integral

import numpy as np

from batched_transfer import multiply


def validate_geometry(q, state, width):
    if (not all(isinstance(value, Integral) for value in (q, state, width))
            or q not in (4, 16)
            or not 1 <= state <= width <= q
            or int(q)**(2*int(state)) > 2**25):
        raise ValueError('supported geometry requires GF4/GF16, 1<=state<=width<=q, and q^(2*state)<=2^25')


def evaluation_arguments(v, z):
    """Validate marks before any signed sum or floating-point power."""
    for value in (v, z):
        if isinstance(value, (Integral, Fraction)):
            continue
        try:
            finite = math.isfinite(value)
        except (TypeError, ValueError, OverflowError) as error:
            raise ValueError('marks must be finite real numbers') from error
        if not finite:
            raise ValueError('marks must be finite real numbers')
    exact = isinstance(v, Fraction) or isinstance(z, Fraction)
    if exact:
        v, z = Fraction(v), Fraction(z)
    if v < 0 or not 0 < z <= 1:
        raise ValueError('require finite v>=0 and 0<z<=1')
    return v, z, exact


def require_finite_transfer(result):
    if any(not math.isfinite(entry) for row in result for entry in row):
        raise ArithmeticError('nonfinite floating transfer; use Fraction marks or rescale the screen')


def exact_sum(values):
    # Empty transition polynomials must remain rational zero.  Starting with
    # integer zero would produce float zero at the subsequent division.
    return sum(values, Fraction(0))


def field_tables(q):
    if q not in (4, 16):
        raise ValueError('supported exact audit fields: GF4 and GF16')
    modulus = 7 if q == 4 else 0x13
    product = np.array([[multiply(a, b, modulus) for b in range(q)]
                        for a in range(q)], dtype=np.uint8)
    trace = []
    for a in range(q):
        value, term = 0, a
        for _ in range(q.bit_length()-1):
            value ^= term
            term = int(product[term, term])
        if value not in (0, 1):
            raise ArithmeticError('invalid absolute trace')
        trace.append(value)
    return product, np.array(trace, dtype=np.uint8)


def vectors_and_codewords(q, state, width):
    validate_geometry(q, state, width)
    product, trace = field_tables(q)
    vectors = np.zeros((q**state, state), dtype=np.uint8)
    numbers = np.arange(q**state, dtype=np.uint32)
    for j in range(state):
        vectors[:, j] = (numbers // q**j) % q
    columns = np.ones((width, state), dtype=np.uint8)
    for j in range(1, state):
        columns[:, j] = product[columns[:, j-1], np.arange(width)]
    codewords = np.zeros((q**state, width), dtype=np.uint8)
    for j in range(state):
        codewords ^= product[vectors[:, j, None], columns[None, :, j]]
    return product, trace, vectors, codewords


@lru_cache(maxsize=None)
def build_histogram(q=16, state=3, width=16):
    product, trace, vectors, codewords = vectors_and_codewords(q, state, width)
    signed = Counter()
    all_zero = codewords == 0
    zero_count = all_zero.sum(axis=1)
    # Bounded temporary arrays; each batch handles 64 refreshed states.
    for begin in range(1, q**state, 64):
        end = min(begin+64, q**state)
        n00 = np.zeros((end-begin, q**state), dtype=np.uint32)
        minus = np.zeros_like(n00)
        for j in range(width):
            n00 += all_zero[begin:end, j, None] & all_zero[None, :, j]
            minus += trace[product[codewords[begin:end, j, None], codewords[None, :, j]]]
        n0n = zero_count[begin:end, None]-n00
        nn0 = zero_count[None, :]-n00
        dot = np.zeros_like(n00, dtype=np.uint8)
        for j in range(state):
            dot ^= product[vectors[begin:end, j, None], vectors[None, :, j]]
        negative = trace[dot].astype(bool)
        # Each count fits five bits, including the all-zero count of sixteen.
        keys = n00 | (n0n.astype(np.uint32) << 5) | (nn0.astype(np.uint32) << 10) | (minus << 15)
        for is_negative, sign in ((False, 1), (True, -1)):
            values, counts = np.unique(keys[negative == is_negative], return_counts=True)
            for key, count in zip(values, counts):
                signed[int(key)] += sign*int(count)
    entries = []
    for key, count in sorted(signed.items()):
        if count:
            n00, n0n, nn0, minus = ((key >> shift) & 31 for shift in (0, 5, 10, 15))
            plus = width-n00-n0n-nn0-minus
            if plus < 0:
                raise ArithmeticError('inconsistent character categories')
            entries.append([n00, n0n, nn0, plus, minus, count])
    return {'q': q, 'state': state, 'width': width,
            'pairs': (q**state-1)*q**state,
            'all_state_zero_counts': dict(sorted(Counter(map(int, zero_count)).items())),
            'nonzero_state_zero_counts': dict(sorted(Counter(map(int, zero_count[1:])).items())),
            'signed_categories': entries}


def moment(v, z, histogram):
    """Return the exact weighted 2x2 transfer, in the input numeric type.

    Use Fraction v,z for exact evaluation. Floats are diagnostic only.
    The entries sum over all input values; they are not normalized by the
    total input mass (1+(q-1)*v)^width.
    """
    q, state, width = (histogram[k] for k in ('q', 'state', 'width'))
    validate_geometry(q, state, width)
    v, z, exact = evaluation_arguments(v, z)
    add = exact_sum if exact else math.fsum
    zero = (1+(q-1)*v*z, 1-v*z)
    birth = add(int(count)*zero[0]**int(h)*zero[1]**(width-int(h))
                for h, count in histogram['all_state_zero_counts'].items()) / q**state
    row0 = (1+(q-1)*v*z)**width
    emitted_nonzero = z+v*(1+(q-2)*z)
    row1 = add(int(count)*zero[0]**int(h)*emitted_nonzero**(width-int(h))
               for h, count in histogram['nonzero_state_zero_counts'].items()) / (q**state-1)
    factors = (zero[0], zero[1], emitted_nonzero, z+v*(1-2*z), z-v)
    terms = []
    for *counts, signed_count in histogram['signed_categories']:
        value = signed_count
        for factor, count in zip(factors, counts):
            value *= factor**count
        terms.append(value)
    joint = add(terms) / (q**state*(q**state-1))
    result = [[birth, row0-birth], [joint, row1-joint]]
    if exact:
        if any(entry < 0 for row in result for entry in row):
            raise ArithmeticError('negative exact marked probability')
    else:
        require_finite_transfer(result)
        scale = max(row0, row1)
        if any(entry < -1e-10*scale for row in result for entry in row):
            raise ArithmeticError('unstable signed floating evaluation')
        result = [[max(0.0, entry) for entry in row] for row in result]
    return result


def aggregate(v, z, state=3):
    """GF16 screen interface; all floating sums use positive coefficients."""
    try:
        result = np.asarray(positive_aggregate(v, z, state), dtype=float)
    except OverflowError as error:
        raise ArithmeticError('floating transfer overflow; rescale the screen') from error
    require_finite_transfer(result)
    return result


def polynomial_product(left, right):
    result = Counter()
    for (w, h), a in left.items():
        for (v, j), b in right.items():
            result[w+v, h+j] += a*b
    return {key: value for key, value in result.items() if value}


@lru_cache(maxsize=None)
def integer_enumerator(q=16, state=3, width=16):
    """Positive integer IO counts, eliminating signed numerical cancellation.

    Row0 counts X for R=0. Row1 counts (R,X) for every R!=0, so normalize
    that row by q^state-1. Coefficients index input/output nonzero weights.
    The Fourier numerator must divide exactly by q^state coefficientwise.
    """
    hist = build_histogram(q, state, width)
    factors = ({(0, 0): 1, (1, 1): q-1},
               {(0, 0): 1, (1, 1): -1},
               {(0, 1): 1, (1, 0): 1, (1, 1): q-2},
               {(0, 1): 1, (1, 0): 1, (1, 1): -2},
               {(0, 1): 1, (1, 0): -1})
    powers = []
    for factor in factors:
        table = [{(0, 0): 1}]
        for _ in range(width):
            table.append(polynomial_product(table[-1], factor))
        powers.append(table)
    returned = Counter()
    for *counts, signed_count in hist['signed_categories']:
        poly = {(0, 0): signed_count}
        for table, count in zip(powers, counts):
            if count:
                poly = polynomial_product(poly, table[count])
        returned.update(poly)
    for key, value in list(returned.items()):
        quotient, remainder = divmod(value, q**state)
        if remainder or quotient < 0:
            raise ArithmeticError('character numerator has invalid integer count')
        returned[key] = quotient
    birth, all_nonzero = Counter(), Counter()
    for h, count in hist['all_state_zero_counts'].items():
        poly = polynomial_product(powers[0][h], powers[1][width-h])
        for key, value in poly.items():
            birth[key] += count*value
    for key, value in list(birth.items()):
        quotient, remainder = divmod(value, q**state)
        if remainder or quotient < 0:
            raise ArithmeticError('birth numerator has invalid integer count')
        birth[key] = quotient
    for h, count in hist['nonzero_state_zero_counts'].items():
        poly = polynomial_product(powers[0][h], powers[2][width-h])
        for key, value in poly.items():
            all_nonzero[key] += count*value
    all_zero = powers[0][width]
    stay0 = {key: all_zero.get(key, 0)-birth.get(key, 0) for key in set(all_zero) | set(birth)}
    stay1 = {key: all_nonzero.get(key, 0)-returned.get(key, 0) for key in set(all_nonzero) | set(returned)}
    result = [[birth, stay0], [returned, stay1]]
    if any(value < 0 for row in result for poly in row for value in poly.values()):
        raise ArithmeticError('negative IO enumerator coefficient')
    result = [[{key: value for key, value in poly.items() if value} for poly in row] for row in result]
    # This identity checks every input-weight class, including both state
    # destinations.  It is stronger than checking just the total mass.
    for row, polynomials in enumerate(result):
        totals = [0] * (width+1)
        for polynomial in polynomials:
            for (weight, output_weight), count in polynomial.items():
                if not 0 <= weight <= width or not 0 <= output_weight <= width:
                    raise ArithmeticError('IO enumerator weight outside the lane range')
                totals[weight] += count
        for weight, count in enumerate(totals):
            expected = math.comb(width, weight)*(q-1)**weight*(q**state-1 if row else 1)
            if count != expected:
                raise ArithmeticError('IO enumerator input-weight mass mismatch')
    return result


def positive_aggregate(v, z, state=3, q=16, width=16):
    """Evaluate the exact positive integer enumerator, with no signed sums."""
    v, z, exact = evaluation_arguments(v, z)
    polynomials = integer_enumerator(q, state, width)
    try:
        vp, zp = [v**w for w in range(width+1)], [z**h for h in range(width+1)]
        add = exact_sum if exact else math.fsum
        result = [[add(count*vp[w]*zp[h] for (w, h), count in poly.items()) /
                   (q**state-1 if row else 1) for poly in polynomials[row]] for row in range(2)]
    except OverflowError as error:
        raise ArithmeticError('floating transfer overflow; use Fraction marks or rescale the screen') from error
    if not exact:
        require_finite_transfer(result)
    return result


@lru_cache(maxsize=None)
def enumerator_by_weight(q=16, state=3, width=16):
    """Group positive counts once for repeated sparse-weight evaluations."""
    polynomials = integer_enumerator(q, state, width)
    return tuple(tuple(tuple(tuple((h, count) for (w, h), count in polynomial.items() if w == weight)
                             for polynomial in row) for row in polynomials)
                 for weight in range(width+1))


def conditional_weight_transfer(weight, z, state=3, q=16, width=16):
    """Return T_weight(z) for a uniform mask and uniform nonzero values.

    Row zero conditions on R=0. Row one averages uniform R!=0. Therefore
    each coefficient count is divided by binom(width,weight)*(q-1)^weight,
    and row one is additionally divided by q^state-1. No input mark v
    remains. Fraction z gives exact rational entries.
    """
    validate_geometry(q, state, width)
    if not isinstance(weight, Integral) or not 0 <= weight <= width:
        raise ValueError('input weight must be an integer between zero and width')
    _, z, exact = evaluation_arguments(0, z)
    return evaluate_conditional_weight(weight, z, exact, state, q, width)


def evaluate_conditional_weight(weight, z, exact, state, q, width, powers=None):
    powers = powers if powers is not None else [z**h for h in range(width+1)]
    add = exact_sum if exact else math.fsum
    input_mass = math.comb(width, weight)*(q-1)**weight
    grouped = enumerator_by_weight(q, state, width)[weight]
    result = [[add(count*powers[h] for h, count in polynomial) /
               (input_mass*(q**state-1 if row else 1))
               for polynomial in grouped[row]] for row in range(2)]
    if not exact:
        require_finite_transfer(result)
    return result


def conditional_transfers(z, state=3, q=16, width=16):
    """Return every T_w(z), w=0..width, reusing powers for sparse screens."""
    validate_geometry(q, state, width)
    _, z, exact = evaluation_arguments(0, z)
    powers = [z**h for h in range(width+1)]
    return [evaluate_conditional_weight(weight, z, exact, state, q, width, powers)
            for weight in range(width+1)]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--q', type=int, choices=(4, 16), default=16)
    parser.add_argument('--state', type=int, default=3)
    parser.add_argument('--width', type=int, default=16)
    parser.add_argument('--histogram', action='store_true')
    parser.add_argument('--v', default='1/10')
    parser.add_argument('--z', default='1/2')
    args = parser.parse_args()
    hist = build_histogram(args.q, args.state, args.width)
    if args.histogram:
        print(json.dumps(hist, indent=2))
    else:
        transfer = positive_aggregate(Fraction(args.v), Fraction(args.z), args.state, args.q, args.width)
        print(json.dumps({'q': args.q, 'state': args.state, 'width': args.width,
                          'pairs': hist['pairs'], 'categories': len(hist['signed_categories']),
                          'transfer': [[str(v) for v in row] for row in transfer]}, indent=2))
