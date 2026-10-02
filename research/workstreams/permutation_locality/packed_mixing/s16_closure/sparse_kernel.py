"""Fixed-occupancy envelopes for a fresh uniform invertible state map.

For input x and entering state q, the epoch emits x+Aq and leaves Tq+Cx.
T is sampled independently and uniformly from GL(s,2) for each epoch.
The input has a uniform j-subset of packet positions and independent
uniform nonzero four-bit labels. These are bounds, not a code certificate.

Coordinates are zero mass, arbitrary nonzero mass, a measure dominated
pointwise by a multiple of uniform nonzero mass, and arbitrary measures
restricted to the retained expansion-weight classes. For q != 0, Tq is
uniform nonzero independently of the emitted weight. Translation by Cx
gives each nonzero destination probability at most 1/(2**s-1). A return
to zero additionally requires Cx != 0. No return is possible at j=0.
"""
from fractions import Fraction as Q
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
for path in (ROOT, ROOT/'gf16_packets'):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from flint import arb, arb_mat
import occupancy_birth_classes as births
import occupancy_kernel as occupancy

aq, up = occupancy.aq, occupancy.up


def validate(data, z):
    if (data.get('distribution') != 'uniform_gl' or 'updates' in data
            or type(data.get('bits')) is not int or not 1 <= data['bits'] <= 24
            or type(data.get('windows')) is not int or not 1 <= data['windows'] <= 32
            or not z.is_finite() or not 0 < z <= 1):
        raise ValueError('explicit fresh uniform-GL distribution and valid geometry/weight required')
    levels = list(map(int, data['birth_class_levels']))
    weights = list(map(int, data['image_histogram_weights']))
    if (not levels or levels != sorted(set(levels)) or min(levels) <= 0
            or set(weights) != set(levels)
            or len(weights) != len(data['histograms'])
            or len(data['zero_probabilities']) != data['windows']+1
            or any(not 0 <= Q(p) <= 1 for p in data['zero_probabilities'])
            or Q(data['zero_probabilities'][0]) != 1
            or sum(map(int, data['histogram_multiplicities'])) != (1 << data['bits'])-1):
        raise ValueError('complete image classes and exact zero-feedback probabilities required')
    return levels, weights


def outward_at_z(data, z):
    """Return one outward matrix for each exact packet occupancy 0..W."""
    z = arb(z)
    levels, weights = validate(data, z)
    W, L = data['windows'], (1 << data['bits'])-1
    masses = births.class_masses(data, z, include_zero=True)
    profiles = [occupancy.polynomial(h, z) for h in data['histograms']]
    squares = [occupancy.polynomial(h, z*z) for h in data['histograms']]
    selected = [[k for k, weight in enumerate(weights) if weight == level] for level in levels]
    n, result = 3+len(levels), []
    for j in range(W+1):
        total = [up(row[j]) for row in profiles]
        second = [up(row[j]) for row in squares]
        mean = up(sum((int(c)*v for c, v in zip(data['histogram_multiplicities'], total)), arb(0))/L)
        mean2 = up(sum((int(c)*v for c, v in zip(data['histogram_multiplicities'], second)), arb(0))/L)
        nonzero = aq(1-Q(data['zero_probabilities'][j]))
        rows = [[arb(0)]*n for _ in range(n)]
        rows[0][0], rows[0][3:] = masses[j][0], masses[j][1:]

        def refresh_row(index, bound, bound2):
            # Cauchy controls the feedback-nonzero restriction without claiming
            # independence between feedback and the emitted output weight.
            fresh = min(bound, up((bound2*nonzero).sqrt()))
            rows[index][0] = arb(0) if j == 0 else up(fresh/L)
            rows[index][2] = bound

        refresh_row(1, max(total), max(second))
        refresh_row(2, mean, mean2)
        for i, indices in enumerate(selected, 3):
            refresh_row(i, max(total[k] for k in indices), max(second[k] for k in indices))
        matrix = arb_mat(rows)
        if any(not matrix[i,k].is_finite() or not matrix[i,k] >= 0
               for i in range(n) for k in range(n)):
            raise ArithmeticError('finite nonnegative uniform-GL envelope required')
        result.append(matrix)
    return result


def outward(data, tilt):
    if Q(tilt) <= 0:
        raise ValueError('positive exact tilt required')
    return outward_at_z(data, (-aq(Q(tilt))).exp())
