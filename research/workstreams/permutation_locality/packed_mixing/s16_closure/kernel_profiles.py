"""Uniform-GL envelopes partitioned by exact expansion packet profiles.

The first three coordinates retain zero mass, arbitrary nonzero mass, and
uniform-density mass. Each further coordinate bounds arbitrary mass on one
exact histogram of the 32 expansion-packet weights. This partition is
separate from the retained total-weight birth classes; their identifiers
are never reinterpreted or passed to legacy birth-class operators.

Input laws are exchangeable across packet positions. Thus the emitted
weight moment is identical for every state in an expansion-profile class.
The update remains fresh uniform GL, followed by raw-input feedback.
"""
from fractions import Fraction as Q
from math import comb
import numpy as np
from flint import arb, arb_mat
import kernel_maps
import kernel_birth_density as density
import occupancy_kernel
import birth_classes
import scalar_cover as sc

aq, up = kernel_maps.aq, kernel_maps.up


def attach(data):
    images = kernel_maps.authenticate(data)
    if data['distribution'] != 'uniform_gl':
        raise ValueError('profile kernel requires exact fresh uniform GL')
    if 'birth_character_indices' not in data:
        data = density.attach(data)
    W = data['windows']; radix = W+1
    powers = np.array([radix**j for j in range(5)], dtype=np.uint64)
    low = np.array([x & ((1 << 64)-1) for x in images], dtype=np.uint64)
    high = np.array([x >> 64 for x in images], dtype=np.uint64)
    codes = np.zeros(len(images), dtype=np.uint64)
    for packet in range(W):
        half = low if packet < 16 else high
        weights = np.bitwise_count((half >> (4*(packet % 16))) & 15)
        codes += powers[weights]
    profiles = np.array([sum(int(n)*radix**j for j, n in enumerate(row))
                         for row in data['histograms']], dtype=np.uint64)
    order = np.argsort(profiles)
    positions = np.searchsorted(profiles[order], codes[1:])
    if np.any(positions >= len(profiles)) or not np.array_equal(profiles[order[positions]], codes[1:]):
        raise ArithmeticError('state-to-expansion-profile lookup failed')
    indices = np.empty(len(images), dtype=np.int32); indices[0] = -1
    indices[1:] = order[positions]
    if not np.array_equal(np.bincount(indices[1:], minlength=len(profiles)), data['histogram_multiplicities']):
        raise ArithmeticError('expansion-profile partition has incorrect multiplicities')
    return dict(data, profile_partition=True, profile_state_indices=indices,
                profile_count=len(profiles), profile_kernel='exact-packet-profile-uniform-gl-1')


def prepare(feedback='weight5', seed=0, refresh='uniform', updates=8, birth_density='capped'):
    if refresh != 'uniform':
        raise ValueError('profile kernel implements only exact uniform GL')
    data, record = kernel_maps.prepare(feedback, seed, refresh, updates, birth_density)
    return attach(data), record


def validate(data):
    kernel_maps.authenticate(data)
    if (data['distribution'] != 'uniform_gl' or data.get('profile_partition') is not True
            or data.get('profile_count') != len(data['histograms'])
            or len(data.get('profile_state_indices', [])) != 1 << data['bits']):
        raise ValueError('fresh explicit expansion-profile partition required')


def _floating_law(probabilities):
    law = np.asarray(probabilities, dtype=float)
    if law.shape != (5,) or np.any(law < 0) or not np.isclose(law.sum(), 1., rtol=0., atol=1e-14):
        raise ValueError('five normalized packet-weight probabilities required')
    p = 1-float(law[0])
    if not np.allclose(law[1:], [comb(4, w)*p/15 for w in range(1, 5)], rtol=1e-12, atol=1e-15):
        raise ValueError('active packets must be uniform nonzero GF16 labels')
    return p


def floating(data, probabilities, tilt):
    validate(data)
    p = _floating_law(probabilities); z = np.exp(-float(tilt))
    S = 1 << data['bits']; L = S-1; count = data['profile_count']; n = 3+count
    values = np.array([((1+z)**(4-r)*(1-z)**r-1)/15 for r in range(5)])
    fourier = np.prod((1-p+p*values)[None, :]**data['records'], axis=1)
    law = density.float_law(data, fourier)
    normal = np.prod((1-p+p*np.array([1., -1/15, -1/15, -1/15, -1/15]))[None, :]**data['records'], axis=1)
    nonzero = max(0., min(1., 1-float(data['multiplicities']@normal)/S))
    def moments(value):
        factors = (1-16*p/15)*value**np.arange(5)+(p/15)*(1+value)**4
        return np.prod(factors[None, :]**data['histograms'], axis=1)
    first, second = moments(z), moments(z*z)
    totals = [float(first.max()), float(data['histogram_multiplicities']@first)/L, *first]
    seconds = [float(second.max()), float(data['histogram_multiplicities']@second)/L, *second]
    matrix = np.zeros((n, n)); matrix[0, 0] = law[0]
    matrix[0, 3:] = np.bincount(data['profile_state_indices'][1:], weights=law[1:], minlength=count)
    for source, (total, total2) in enumerate(zip(totals, seconds), 1):
        matrix[source, 0] = min(total, np.sqrt(max(0., total2*nonzero)))/L
        matrix[source, 2] = total
    if data.get('birth_density', 'classes') == 'classes':
        return matrix
    options = []
    for row in density.float_rows(data, matrix[0], law, class_labels=range(count), state_labels=data['profile_state_indices']):
        option = matrix.copy(); option[0] = row; options.append(option)
    return min(options, key=sc.log_power)


def _rows(data, first, second, nonzero, zero, numerators, denominator):
    """Build one complete outward operator before optional cap selection."""
    count = data['profile_count']; n = 3+count; L = (1 << data['bits'])-1
    mean = up(sum((int(c)*v for c, v in zip(data['histogram_multiplicities'], first)), arb(0))/L)
    mean2 = up(sum((int(c)*v for c, v in zip(data['histogram_multiplicities'], second)), arb(0))/L)
    totals = [max(map(up, first)), mean, *map(up, first)]
    seconds = [max(map(up, second)), mean2, *map(up, second)]
    rows = [[arb(0)]*n for _ in range(n)]; rows[0][0] = zero
    for profile in range(count):
        mass = int(numerators[data['profile_state_indices'] == profile].sum())
        rows[0][3+profile] = up(arb(mass)/denominator)
    for source, (total, total2) in enumerate(zip(totals, seconds), 1):
        rows[source][0] = up(min(total, up((total2*nonzero).sqrt()))/L)
        rows[source][2] = total
    return rows


def _options(data, rows, numerators, denominator):
    if data.get('birth_density', 'classes') == 'classes':
        return [rows[0]]
    return list(density.outward_rows(data, rows[0], numerators, denominator,
        class_labels=range(data['profile_count']), state_labels=data['profile_state_indices']))


def outward_at_z(data, probabilities, z):
    validate(data); p = kernel_maps._activity(probabilities)
    if not 0 < z <= 1:
        raise ValueError('output weight in (0,1] required')
    S = 1 << data['bits']
    fourier, _ = birth_classes.base.fourier_parts(data, p, z, 0)
    normal, _ = birth_classes.base.fourier_parts(data, p, arb(1), 0)
    numerators, denominators = density.upper_laws(data, [[v] for v in fourier])
    zero = max(arb(0), up(sum((int(c)*v for c, v in zip(data['multiplicities'], fourier)), arb(0))/S))
    nonzero = min(arb(1), max(arb(0), up(1-sum((int(c)*v for c, v in zip(data['multiplicities'], normal)), arb(0))/S)))
    def moments(value):
        factors = [aq(1-16*p/15)*value**w+aq(p/15)*(1+value)**4 for w in range(5)]
        result = []
        for profile in data['histograms']:
            moment = arb(1)
            for factor, multiplicity in zip(factors, profile):
                moment *= factor**int(multiplicity)
            result.append(moment)
        return result
    rows = _rows(data, moments(z), moments(z*z), nonzero, zero, numerators[:, 0], denominators[0])
    options = []
    for row in _options(data, rows, numerators[:, 0], denominators[0]):
        value = arb_mat([row, *rows[1:]])
        array = np.array([[float(value[i, j]) for j in range(value.ncols())] for i in range(value.nrows())])
        options.append((sc.log_power(array), value))
    return min(options, key=lambda item: item[0])[1]


def outward(data, probabilities, tilt):
    if Q(tilt) <= 0:
        raise ValueError('positive output tilt required')
    return outward_at_z(data, probabilities, (-aq(Q(tilt))).exp())


def outward_local_at_z(data, z, activity=Q(1, 2)):
    validate(data)
    if not 0 < z <= 1 or not 0 <= Q(activity) <= 1:
        raise ValueError('valid output weight and row-selection activity required')
    W = data['windows']; count = data['profile_count']; n = 3+count
    numerators, denominators = density.conditional_laws(data, z)
    first = [occupancy_kernel.polynomial(h, z) for h in data['histograms']]
    second = [occupancy_kernel.polynomial(h, z*z) for h in data['histograms']]
    all_rows = []
    for j in range(W+1):
        zero = min(up(arb(int(numerators[0, j]))/denominators[j]),
                   up(z**j*aq(Q(data['zero_probabilities'][j]))))
        rows = _rows(data, [v[j] for v in first], [v[j] for v in second],
            aq(1-Q(data['zero_probabilities'][j])), zero, numerators[:, j], denominators[j])
        all_rows.append(rows)
    if data.get('birth_density', 'classes') == 'classes':
        return list(map(arb_mat, all_rows))
    p = float(activity)
    masses = np.array([comb(W, j)*p**j*(1-p)**(W-j) for j in range(W+1)])
    arrays = np.array([[[float(v) for v in row] for row in rows] for rows in all_rows])
    potential = density.continuation(np.tensordot(masses, arrays, axes=1))
    result = []
    for j, rows in enumerate(all_rows):
        selected = min(_options(data, rows, numerators[:, j], denominators[j]),
            key=lambda row: sum(float(x)*v for x, v in zip(row, potential)))
        result.append(arb_mat([selected, *rows[1:]]))
    return result
