"""Floating exact-occupancy composition of the symmetric state-rank cap cover.

The packet contract dominates only nonnegative continuations nondecreasing
in the nonzero rank cap. Positive products and sums preserve that cone.
No exact-rank distribution or whole-code certificate is claimed.
"""

import argparse
import json
import math

import numpy as np

from cauchy_firstmoment_screen import geometry, log_choose
from mask_regional_screen import normalized_power
from state_rank_cap_transfer import local_matrix
from superregular_return_screen import bound as two_state_bound


WITNESSES = {
    8: (6.0198647674, 5.2348778115),
    16: (5.7892288799, 2.9328870797),
    32: (5.7272037908, 2.2433313929),
    64: (5.6671559436, 1.9649272430),
    128: (5.618724367961515, 1.8572075805232853),
}


def region(active, y, tau, rows):
    if rows % 16 or not 0 <= active <= rows:
        raise ValueError('invalid regional geometry')
    local = np.asarray([local_matrix(w, y, tau) for w in range(min(active,16)+1)], dtype=float)
    coefficients = normalized_power(local, 16, rows//16, active)
    marks = np.asarray([log_choose(active,w)+w*math.log(tau) for w in range(active+1)])
    scale = float(marks.max())
    matrix = np.einsum('i,ijk->jk', np.exp(marks-scale), coefficients)
    peak = float(matrix.max())
    if peak <= 0 or not math.isfinite(peak):
        raise ArithmeticError('invalid regional scale')
    return matrix/peak, scale+math.log(peak)


def bound(active, point, k_log=20, length=64, dimension=44):
    if active <= 0:
        raise ValueError('positive active-row count required')
    g = geometry(k_log, length, dimension)
    logtilt, logtau = map(float, point)
    logy = -math.exp(logtilt)*active/g['n']
    if -g['rows']*logy > 600:
        raise ValueError('regional attenuation exceeds this floating screen guard')
    matrix, scale = region(active, math.exp(logy), math.exp(logtau), g['rows'])
    values = np.ones(9)
    accumulated = 0.
    for _ in range(length):
        # Compose backwards so the continuation cone can be checked directly.
        values = matrix @ values
        peak = float(values.max())
        if not peak > 0 or not math.isfinite(peak):
            raise ArithmeticError('invalid continuation scale')
        values /= peak
        if any(values[r] > values[r+1]+1e-12 for r in range(1,8)):
            raise ArithmeticError('nonzero-cap continuation lost monotonicity')
        accumulated += scale+math.log(peak)
    if not values[0] > 0:
        raise ArithmeticError('zero initial-state continuation underflow')
    return (log_choose(g['rows'],active)-((length-dimension)*active+1)*logtau
            -g['threshold']*logy+accumulated+math.log(float(values[0])))/math.log(2)


def compare(active, point, k_log=20, length=64, dimension=44):
    return dict(active=active, k_log=k_log, outer=[length, dimension], state=8,
                logtilt=point[0], logtau=point[1],
                rank_cap_margin_bits=-bound(active,point,k_log,length,dimension),
                two_state_margin_bits=-two_state_bound(active,point,k_log,length,dimension,8),
                rigorous=False, whole_code=False, method='nine-state monotone rank caps')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--active', type=int, nargs='+', default=[16,32,64,128])
    parser.add_argument('--point', type=float, nargs=2)
    parser.add_argument('--k-log', type=int, default=20)
    parser.add_argument('--outer', type=int, nargs=2, default=[64,44])
    args = parser.parse_args()
    for active in args.active:
        print(json.dumps(compare(active, args.point or WITNESSES[active], args.k_log,
                                 *args.outer)), flush=True)
