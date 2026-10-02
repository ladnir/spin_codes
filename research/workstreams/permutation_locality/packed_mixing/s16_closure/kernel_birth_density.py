"""Map-specific capped-density representations of newborn state measures.

Let m_s=E[z^wt(X) 1{CX=s}] for a fixed input law.  For any d>=0, the
nonzero measure splits into min(m_s,d) plus (m_s-d)_+.  The first part is
dominated by uniform-density coordinate (2**bits-1)*d.  Summing the second
part within each expansion-weight class gives valid class coordinates.
Every candidate below is a complete representation of this same measure.
Choosing a whole candidate row is valid; taking entrywise minima is not.

Outward Walsh inversion uses exact signed int64 arithmetic on dyadic
interval enclosures.  Proposal-only floating inversion is kept separate.
"""
from math import comb
import numpy as np
from flint import arb, arb_poly
import birth_classes

aq, up = birth_classes.aq, birth_classes.up


def attach(data):
    S = 1 << data['bits']; W = data['windows']; radix = W+1
    chars = np.arange(S, dtype=np.uint32)
    codes = np.zeros(S, dtype=np.uint64)
    powers = np.array([radix**j for j in range(5)], dtype=np.uint64)
    for offset in range(0, 4*W, 4):
        weights = sum(np.bitwise_count(chars & c)&1 for c in data['columns'][offset:offset+4])
        codes += powers[weights]
    records = np.array([sum(int(n)*radix**j for j, n in enumerate(row))
                        for row in data['records']], dtype=np.uint64)
    order = np.argsort(records)
    positions = np.searchsorted(records[order], codes)
    if np.any(positions >= len(records)) or not np.array_equal(records[order[positions]], codes):
        raise ArithmeticError('character-to-profile lookup failed')
    indices = order[positions]
    if not np.array_equal(np.bincount(indices, minlength=len(records)), data['multiplicities']):
        raise ArithmeticError('character profile multiplicities disagree')
    image_weights = np.array([int(x).bit_count() for x in data['map_images']], dtype=np.int16)
    return dict(data, birth_character_indices=indices, birth_state_weights=image_weights)


def walsh(values):
    result = np.array(values, copy=True)
    n = len(result); half = 1
    if not n or n & (n-1):
        raise ValueError('power-of-two first dimension required')
    while half < n:
        blocks = result.reshape((-1, 2*half, *result.shape[1:]))
        left, right = blocks[:, :half].copy(), blocks[:, half:].copy()
        blocks[:, :half], blocks[:, half:] = left+right, left-right
        half *= 2
    return result


def float_law(data, values):
    return np.maximum(0., walsh(np.asarray(values)[data['birth_character_indices']])/(1 << data['bits']))


def upper_laws(data, profile_values):
    """Return upper numerators and a separate dyadic divisor per column."""
    S = 1 << data['bits']; precision = min(44, 61-data['bits']); budget = 1 << precision
    rows = list(profile_values)
    if len(rows) != len(data['records']) or not rows or not rows[0]:
        raise ValueError('one nonempty value row per exact character profile required')
    width = len(rows[0]); lower = np.zeros((len(rows), width), dtype=np.int64)
    zero_character = int(data['birth_character_indices'][0])
    scales = []
    for value in rows[zero_character]:
        # The zero-character Fourier coefficient is the total tilted mass.
        # Every other coefficient has absolute value at most that mass.
        # Scaling each occupancy separately avoids an absolute noise floor.
        upper = abs(value).upper()
        if not upper.is_finite() or not upper > 0:
            raise ArithmeticError('positive finite total Fourier mass required')
        mantissa, exponent = upper.man_exp()
        mass_exponent = int(exponent)+abs(int(mantissa)).bit_length()
        fractional_bits = precision-mass_exponent
        if not 0 <= fractional_bits <= 4096:
            raise ArithmeticError('unsupported Fourier scaling exponent')
        scales.append(1 << fractional_bits)
    errors = [0]*width
    for index, (row, multiplicity) in enumerate(zip(rows, data['multiplicities'])):
        if len(row) != width:
            raise ValueError('matching occupancy ranges required')
        for j, value in enumerate(row):
            scale = scales[j]
            if not value.is_finite():
                raise ArithmeticError('finite Fourier coefficient required')
            lo = int((value*scale).lower().floor().unique_fmpz())
            hi = int((value*scale).upper().ceil().unique_fmpz())
            # Each normalized Fourier coefficient has absolute value <=1.
            # An interval wider than this generous dyadic rounding budget
            # must be retried at greater Arb precision, not silently accepted.
            if max(abs(lo), abs(hi)) > budget+4 or not 0 <= hi-lo <= 4:
                raise ArithmeticError('insufficient precision for bounded dyadic Fourier inversion')
            lower[index, j] = lo
            errors[j] += int(multiplicity)*(hi-lo)
    if max(errors) >= 1 << 61 or S*(budget+4) >= 1 << 62:
        raise OverflowError('signed Walsh accumulation bound exceeded')
    transformed = walsh(lower[data['birth_character_indices']])
    result = np.maximum(0, transformed+np.array(errors, dtype=np.int64)[None, :])
    if any(sum(map(int, result[:, j])) >= 1 << 62 for j in range(width)):
        raise OverflowError('class residual accumulation bound exceeded')
    return result, [S*scale for scale in scales]


def conditional_laws(data, z):
    W = data['windows']
    values = [((1+z)**(4-r)*(1-z)**r-1)/15 for r in range(5)]
    powers = [[arb_poly([1, value])**j for j in range(W+1)] for value in values]
    rows = []
    divisors = [comb(W, j) for j in range(W+1)]
    for profile in data['records']:
        polynomial = arb_poly([1])
        for table, count in zip(powers, profile):
            polynomial *= table[int(count)]
        rows.append([polynomial[j]/divisors[j] for j in range(W+1)])
    return upper_laws(data, rows)


def thresholds(values):
    positive = values[values > 0]
    if not len(positive):
        return []
    # These points only choose valid caps. No quantile claim enters a bound.
    sorted_values = np.sort(positive)
    return sorted(set(sorted_values[min(len(sorted_values)-1, len(sorted_values)*j//8)].item()
                      for j in (0, 2, 4, 6, 7, 8)))


def float_rows(data, baseline, law, *, class_labels=None, state_labels=None):
    yield baseline
    L = (1 << data['bits'])-1
    weights = data['birth_state_weights'][1:] if state_labels is None else np.asarray(state_labels)[1:]
    classes = data['birth_class_levels'] if class_labels is None else class_labels
    for cap in thresholds(law[1:]):
        row = baseline.copy(); row[2] = L*cap
        residual = np.maximum(0., law[1:]-cap)
        for target, level in enumerate(classes, 3):
            row[target] = min(baseline[target], float(residual[weights == level].sum()))
        yield row


def outward_rows(data, baseline, numerators, denominator, *, class_labels=None, state_labels=None):
    yield baseline
    L = (1 << data['bits'])-1
    weights = data['birth_state_weights'][1:] if state_labels is None else np.asarray(state_labels)[1:]
    classes = data['birth_class_levels'] if class_labels is None else class_labels
    for cap in thresholds(numerators[1:]):
        cap = int(cap)
        row = list(baseline); row[2] = up(arb(L*cap)/denominator)
        residual = np.maximum(0, numerators[1:]-cap)
        for target, level in enumerate(classes, 3):
            value = int(residual[weights == level].sum())
            row[target] = min(baseline[target], up(arb(value)/denominator))
        yield row


def continuation(matrix):
    """A floating positive potential chooses rows, never proves inequalities."""
    values, vectors = np.linalg.eig(matrix)
    index = int(np.argmax(values.real))
    potential = np.abs(vectors[:, index].real)
    if not np.all(np.isfinite(potential)) or not potential.max() > 0:
        return np.ones(matrix.shape[0])
    return np.maximum(potential/potential.max(), 1e-200)


def refine_local(data, local, z, activity):
    """Select one complete birth row per occupancy using a common potential."""
    W = data['windows']; n = local[0].nrows()
    arrays = np.array([[[float(matrix[i, j]) for j in range(n)] for i in range(n)] for matrix in local])
    p = float(activity)
    if not 0 <= p <= 1:
        raise ValueError('valid proposal activity required')
    weights = np.array([comb(W, j)*p**j*(1-p)**(W-j) for j in range(W+1)])
    potential = continuation(np.tensordot(weights, arrays, axes=1))
    numerators, denominators = conditional_laws(data, z)
    result = []
    for j, matrix in enumerate(local):
        baseline = [matrix[0, k] for k in range(n)]
        options = list(outward_rows(data, baseline, numerators[:, j], denominators[j]))
        selected = min(options, key=lambda row: sum(float(x)*v for x, v in zip(row, potential)))
        value = matrix*arb(1)
        for target, mass in enumerate(selected):
            value[0, target] = mass
        result.append(value)
    return result
