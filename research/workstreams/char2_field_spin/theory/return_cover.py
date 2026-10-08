"""Whole-cover screen adding independent state-return equations.

At a zero-state checkpoint choose up to four remaining diagonal columns
as MDS pivots. Each consumes its randomness once and contracts the alive
rank potential by 1/tau. This is an exploratory floating cover, not a
directed certificate; the exact local inequality has separate tests.
"""

import argparse
from fractions import Fraction
from itertools import product
import json
import math
from random import Random

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

from regional_sparse_cover import (LOG_BINOMIAL, RETURN_TYPES, gf4_mul,
                                   gf4_rank)


def return_screen_objective(variables, alpha, delta=.1, length=32, parity=10):
    log_x, log_y, log_tau, log_tilt, log_eta = variables
    log_union = logsumexp([-log_tau+log_x, log_y, log_x+log_y])
    zero_zero = logsumexp([0]+[LOG_BINOMIAL[w]+w*(log_x+log_y)-4*log_tau
                              for w in range(5, 17)])
    zero_nonzero = logsumexp([LOG_BINOMIAL[w]+w*(log_x+log_y) for w in range(1, 17)])
    nonzero_nonzero = logsumexp([LOG_BINOMIAL[w]+w*log_union for w in range(13, 17)])
    # Return pivots use only diagonals at I intersect J. I minus J was
    # already consumed by the output-zero phase.
    both = RETURN_TYPES[:, 2]-RETURN_TYPES[:, 1]
    valid = (both > 0) | (RETURN_TYPES[:, 3] != 1)
    nonzero_zero = logsumexp((RETURN_TYPES[:, 0]
        -(RETURN_TYPES[:, 1]+np.minimum(both, 4))*log_tau
        +RETURN_TYPES[:, 2]*log_x+RETURN_TYPES[:, 3]*log_y)[valid])
    shift = max(zero_zero, zero_nonzero, nonzero_zero, nonzero_nonzero)
    a, b, c, d = [math.exp(value-shift) for value in
                  (zero_zero, zero_nonzero, nonzero_zero, nonzero_nonzero)]
    log_radius = shift+math.log((a+d+math.sqrt((a-d)**2+4*b*c))/2)
    log_row = logsumexp([length*math.log1p(-alpha),
        -parity*log_tau-(parity+1)*log_eta
        +length*np.logaddexp(math.log1p(-alpha),
                             log_tau+log_tilt+log_eta+math.log(alpha))])
    return (log_radius/16-alpha*log_x-delta*log_y
            -alpha*log_tilt+log_row/length)/math.log(2)


def screen(alpha):
    choices = [minimize(return_screen_objective, [-1, -2, lt, 0, 0], args=(alpha,),
                method='L-BFGS-B', bounds=[(-25, 10), (-20, 0),
                        (0, math.log(255)), (-15, 15), (0, 15)])
               for lt in (0, 1, 2, 4, math.log(255))]
    best = min(choices, key=lambda result: result.fun)
    return dict(alpha=alpha, log2_envelope_per_N=float(best.fun),
                logarithmic_parameters=best.x.tolist(), rigorous=False,
                scope='complete support cover with return pivots; finite factors omitted')


def exact_return_checks():
    """Check arbitrary input-form dependencies, using exact probabilities."""
    rng = Random(13003222)
    count = 0
    for width in range(1, 4):
        columns = [(1, value) for value in range(width)]
        for _ in range(80):
            forms = []
            for _ in range(width):
                value = (0, 0)
                while not any(value):
                    value = tuple(rng.randrange(4) for _ in range(2))
                forms.append(value)
            offset = [[rng.randrange(4) for _ in range(2)] for _ in range(2)]
            for tau in (Fraction(1), Fraction(2), Fraction(3)):
                total = Fraction(0)
                for diagonals in product((1, 2, 3), repeat=width):
                    equations = [row[:] for row in offset]
                    for column, diagonal, form in zip(columns, diagonals, forms):
                        for i in range(2):
                            for j in range(2):
                                equations[i][j] ^= gf4_mul(column[i], gf4_mul(diagonal, form[j]))
                    rank = gf4_rank(equations)
                    if any(gf4_rank(equations+[list(form)]) == rank for form in forms):
                        continue
                    total += tau**(-rank)
                assert total/3**width <= tau**(-min(width, 2))
                count += 1
    return count


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--screen', action='store_true')
    parser.add_argument('--alpha', type=float, nargs='+', default=(.001, .01, .04, .1, .15, .2, .2399))
    args = parser.parse_args()
    if args.check:
        print(json.dumps(dict(exact_return_potential_cases=exact_return_checks())), flush=True)
    if args.screen:
        for alpha in args.alpha:
            print(json.dumps(screen(alpha)), flush=True)
