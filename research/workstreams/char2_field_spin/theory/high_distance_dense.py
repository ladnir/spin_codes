"""Floating dense-class screens for the high-distance field exploration.

Uses exact MDS support counts and the established positive-excess bound.
The input and state diagonals are uniform over the larger field's nonzero
elements. This is NOT valid for small coefficient subsets merely by changing
field_bits. A point-density screen is not interval coverage or a certificate.
"""
from fractions import Fraction as F
from functools import lru_cache
import argparse
import json
import math

from scipy.optimize import minimize_scalar
from scipy.special import logsumexp

from gf8_exact_spectrum_dense import positive_excess
from gf8_punctured_screen import profile


def logqminus1(bits):
    if bits not in (8, 16, 32, 64, 128):
        raise ValueError('unsupported extension of AES GF256')
    return bits * math.log(2) + math.log1p(-2.**(-bits))


def ball_bits_per_symbol(distance, bits):
    d = float(F(distance))
    if not 0 < d <= .5:
        raise ValueError('distance must be in (0,1/2]')
    return (d*logqminus1(bits)-d*math.log(d)-(1-d)*math.log1p(-d))/math.log(2)


@lru_cache(maxsize=64)
def excess_logs(length, dimension, bits):
    return tuple((w, math.log(math.comb(length, w))+
                  math.log(e.numerator)-math.log(e.denominator))
                 for w, e in enumerate(positive_excess(length, dimension, 1 << bits)) if e)


def log_membership_row(alpha, logtilt, length, dimension, bits):
    if not 0 < alpha <= 1 or logtilt < 0:
        raise ValueError('positive density and tilt at least one required')
    baseline = (-(length-dimension)*bits*math.log(2)
                +length*math.log1p(math.expm1(logtilt)*alpha))
    tail = []
    for w, value in excess_logs(length, dimension, bits):
        if alpha == 1 and w < length:
            continue
        tail.append(value+w*logtilt + ((length-w)*math.log1p(-alpha) if w < length else 0))
    return float(logsumexp([baseline, *tail]))


def margin(alpha, logtilt, bits=16, distance=F(2,5), length=64,
           dimension=16, width=4, outputs=2, k_log=20):
    g = profile(width, outputs, length, dimension, k_log)
    membership = (length*math.log(g['rows']+1)
                  +g['rows']*log_membership_row(alpha, logtilt, length, dimension, bits)
                  -alpha*g['input_n']*logtilt)/math.log(2)
    completion = bits*(g['input_n']-g['output_n'])
    ball = g['output_n']*ball_bits_per_symbol(distance, bits)
    return -membership-completion-ball


def search(alpha, bits=16, distance=F(2,5), length=64, dimension=16,
           width=4, outputs=2, k_log=20):
    objective = lambda x: -margin(alpha, x, bits, distance, length, dimension, width, outputs, k_log)
    optimum = minimize_scalar(objective, bounds=(0, max(20.,bits*math.log(2))), method='bounded',
                              options={'xatol': 1e-11})
    witness = min((0., float(optimum.x)), key=objective)
    g = profile(width, outputs, length, dimension, k_log)
    return dict(field_bits=bits, distance=str(F(distance)), outer=[length, dimension],
                width=width, outputs=outputs, input_density=alpha,
                requested_k=2**k_log, input_n=g['input_n'], output_n=g['output_n'],
                q=1 << bits, coefficient_count=(1 << bits)-1,
                diagonal_distribution='independent uniform nonzero elements of the full coefficient field',
                completion_factor_log2=bits*(g['input_n']-g['output_n']),
                ball_bits_per_symbol=ball_bits_per_symbol(distance,bits),
                rate=str(F(dimension*width, length*outputs)), logtilt=witness,
                margin_bits=-objective(witness), rigorous=False, whole_code=False,
                scope='one dense inverse-input density; output completion union included')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bits', type=int, default=16)
    p.add_argument('--distance', type=F, default=F(2,5))
    p.add_argument('--outer', type=int, nargs=2, default=(64,16))
    p.add_argument('--width', type=int, default=4)
    p.add_argument('--outputs', type=int, default=2)
    p.add_argument('--density', type=float, nargs='+', default=(.8,.9,.95,1.))
    a=p.parse_args()
    for density in a.density:
        print(json.dumps(search(density,a.bits,a.distance,*a.outer,a.width,a.outputs)),flush=True)
