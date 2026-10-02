"""Authenticated fixed-map kernels for S16 closure work.

For each epoch, emit A*a+X and then set a'=T*a+C*X.  T is fresh and
independent of X and the entering state.  ``uniform_gl`` samples T uniformly
from GL(bits,2); ``transvections`` uses the retained r-update distribution.
The uniform branch has exactly zero lazy mass, not a large-r approximation.

The coordinates are zero mass, arbitrary nonzero mass, uniform-density
mass, and the retained expansion-weight classes.  Numeric proposals do not
certify an outer domain.  Outward calls check only their stated local input.
"""
from fractions import Fraction as Q
import hashlib
import json
from math import comb
from pathlib import Path
import sys

import numpy as np
from flint import arb, arb_mat

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE.parents[1]/'gf16_packets'))
import s16_maps
import birth_classes
import birth_classes_probe
import refresh_kernel

aq, up = birth_classes.aq, birth_classes.up


def _identity(images, columns, bits):
    rows = [int(images[1 << j]) for j in range(bits)]
    payload = dict(bits=bits, width=len(columns), expansion_rows=list(map(hex, rows)),
                   feedback_columns=list(map(int, columns)))
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def authenticate(data):
    """Bind fresh exact preparation to the declared immutable linear maps."""
    if data.get('distribution') not in ('uniform_gl', 'transvections'):
        raise ValueError('explicit state-update distribution required')
    images = data.get('map_images')
    if not isinstance(images, tuple) or len(images) != 1 << data['bits']:
        raise ValueError('prepared immutable expansion images required')
    if data.get('map_sha256') != _identity(images, data['columns'], data['bits']):
        raise ArithmeticError('prepared map identity changed')
    if data['windows']*4 != len(data['columns']):
        raise ArithmeticError('prepared packet geometry changed')
    if data['distribution'] == 'uniform_gl' and 'updates' in data:
        raise ValueError('uniform GL has no transvection update parameter')
    return images


def prepare_maps(images, columns, bits=16, updates=8, distribution='uniform_gl', birth_density='capped'):
    """Regenerate all exact profiles from explicit maps; accept no saved caps."""
    if distribution not in ('uniform_gl', 'transvections'):
        raise ValueError('unknown state-update distribution')
    if birth_density not in ('classes', 'capped'):
        raise ValueError('birth density must be classes or capped')
    images, columns = tuple(map(int, images)), list(map(int, columns))
    if s16_maps.binary_rank(columns) != bits:
        raise ValueError('full-rank feedback required')
    data = birth_classes.prepare(images, columns, bits, updates)
    data.update(map_images=images, map_sha256=_identity(images, columns, bits),
                distribution=distribution, birth_density=birth_density)
    if distribution == 'uniform_gl':
        del data['updates']
        if birth_density == 'capped':
            import kernel_birth_density
            data = kernel_birth_density.attach(data)
    authenticate(data)
    return data


def prepare(feedback='weight5', seed=0, refresh='uniform', updates=8, birth_density='capped'):
    """Return freshly prepared S16 data and a fully specified map record."""
    if refresh not in ('uniform', 'transvections'):
        raise ValueError('refresh must be uniform or transvections')
    images, columns, record = s16_maps.candidate(feedback, seed)
    distribution = 'uniform_gl' if refresh == 'uniform' else 'transvections'
    data = prepare_maps(images, columns, 16, updates, distribution, birth_density)
    record = dict(record, map_sha256=data['map_sha256'], distribution=distribution,
                  updates=None if refresh == 'uniform' else updates)
    return data, record


def _activity(probabilities):
    probabilities = list(map(Q, probabilities))
    if len(probabilities) != 5 or min(probabilities) < 0 or sum(probabilities) != 1:
        raise ValueError('five normalized packet-weight probabilities required')
    p = 1-probabilities[0]
    if probabilities[1:] != [Q(comb(4, w), 15)*p for w in range(1, 5)]:
        raise ValueError('active packets must be uniform nonzero GF16 labels')
    return p


def floating(data, probabilities, tilt):
    """Floating search proposal; the uniform branch is evaluated directly."""
    authenticate(data)
    if data['distribution'] == 'transvections':
        return birth_classes.floating(data, probabilities, tilt)
    law = np.asarray(probabilities, dtype=float)
    if law.shape != (5,) or np.any(law < 0) or not np.isclose(law.sum(), 1., rtol=0., atol=1e-14):
        raise ValueError('five normalized packet-weight probabilities required')
    p = 1-float(law[0]); z = np.exp(-float(tilt))
    if not np.allclose(law[1:], [comb(4, w)*p/15 for w in range(1, 5)], rtol=1e-12, atol=1e-15):
        raise ValueError('active packets must be uniform nonzero GF16 labels')
    W = data['windows']; S = 1 << data['bits']; L = S-1
    values = np.array([((1+z)**(4-r)*(1-z)**r-1)/15 for r in range(5)])
    fourier = np.prod((1-p+p*values)[None, :]**data['records'], axis=1)
    zero = max(0., float(data['multiplicities']@fourier)/S)
    normal = np.prod((1-p+p*np.array([1., -1/15, -1/15, -1/15, -1/15]))[None, :]**data['records'], axis=1)
    nonzero = max(0., min(1., 1-float(data['multiplicities']@normal)/S))
    def moments(value):
        factors = (1-16*p/15)*value**np.arange(5)+(p/15)*(1+value)**4
        return np.prod(factors[None, :]**data['histograms'], axis=1)
    first, second = moments(z), moments(z*z)
    levels = data['birth_class_levels']; n = 3+len(levels)
    matrix = np.zeros((n, n)); matrix[0, 0] = zero
    matrix[0, 3:] = np.maximum(0, data['birth_class_float']@fourier/S)
    selections = [np.arange(len(first)), None,
                  *[np.flatnonzero(data['image_histogram_weights'] == level) for level in levels]]
    for source, selected in enumerate(selections, 1):
        if selected is None:
            total = float(data['histogram_multiplicities']@first)/L
            total2 = float(data['histogram_multiplicities']@second)/L
        else:
            total, total2 = float(first[selected].max()), float(second[selected].max())
        matrix[source, 0] = min(total, np.sqrt(max(0., total2*nonzero)))/L
        matrix[source, 2] = total
    if data.get('birth_density', 'classes') == 'capped':
        import kernel_birth_density as density
        law = density.float_law(data, fourier)
        options = []
        for row in density.float_rows(data, matrix[0], law):
            option = matrix.copy(); option[0] = row; options.append(option)
        return min(options, key=birth_classes.sc.log_power)
    return matrix


def outward_at_z(data, probabilities, z):
    authenticate(data)
    if data['distribution'] == 'transvections':
        return birth_classes.outward_at_z(data, probabilities, z)
    p = _activity(probabilities)
    if not 0 < z <= 1:
        raise ValueError('output weight in (0,1] required')
    W = data['windows']; S = 1 << data['bits']; L = S-1
    values, _ = birth_classes.base.fourier_parts(data, p, z, 0)
    normal, _ = birth_classes.base.fourier_parts(data, p, arb(1), 0)
    zero = max(arb(0), up(sum((int(n)*v for n, v in zip(data['multiplicities'], values)), arb(0))/S))
    nonzero = min(arb(1), max(arb(0), up(1-sum((int(n)*v for n, v in zip(data['multiplicities'], normal)), arb(0))/S)))
    def moments(value):
        factors = [aq(1-16*p/15)*value**w+aq(p/15)*(1+value)**4 for w in range(5)]
        result = []
        for profile in data['histograms']:
            mass = arb(1)
            for factor, count in zip(factors, profile):
                mass *= factor**int(count)
            result.append(mass)
        return result
    first, second = moments(z), moments(z*z)
    levels = list(map(int, data['birth_class_levels'])); n = 3+len(levels)
    rows = [[arb(0)]*n for _ in range(n)]
    rows[0][0] = zero
    rows[0][3:] = birth_classes._class_masses_from_values(data, values)
    selections = [range(len(first)), None,
                  *[[i for i, w in enumerate(data['image_histogram_weights']) if w == level] for level in levels]]
    for source, selected in enumerate(selections, 1):
        if selected is None:
            total = up(sum((int(n)*v for n, v in zip(data['histogram_multiplicities'], first)), arb(0))/L)
            total2 = up(sum((int(n)*v for n, v in zip(data['histogram_multiplicities'], second)), arb(0))/L)
        else:
            total, total2 = max(up(first[i]) for i in selected), max(up(second[i]) for i in selected)
        rows[source][0] = up(min(total, up((total2*nonzero).sqrt()))/L)
        rows[source][2] = total
    if data.get('birth_density', 'classes') == 'capped':
        import kernel_birth_density as density
        numerators, denominators = density.upper_laws(data, [[v] for v in values])
        options = []
        for row in density.outward_rows(data, rows[0], numerators[:, 0], denominators[0]):
            option = arb_mat([row, *rows[1:]])
            proposed = np.array([[float(option[i, j]) for j in range(n)] for i in range(n)])
            options.append((birth_classes.sc.log_power(proposed), option))
        return min(options, key=lambda item: item[0])[1]
    return arb_mat(rows)


def outward(data, probabilities, tilt):
    if Q(tilt) <= 0:
        raise ValueError('positive output tilt required')
    return outward_at_z(data, probabilities, (-aq(Q(tilt))).exp())
