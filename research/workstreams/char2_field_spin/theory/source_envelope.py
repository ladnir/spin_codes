"""Nonempty-source envelope for fresh-packet field-SPIN (floating screen).

For normalized regional moments R_h, dominate every entry by v**h M.
The resulting row enumerator retains the nonzero-systematic condition:
((1+15*v)**k-1)**a. No coverage or directed certificate is claimed.
"""

import argparse
import json
import math

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import logsumexp

import random_packet_screen as packet
from field_rate_tradeoff import log_choose


def log_nonzero_binomial(dimension, logv):
    base = dimension*np.logaddexp(0., math.log(15)+logv)
    return float(base+math.log(-math.expm1(-base)))


def dominate(regions, logv):
    return np.max(regions-np.arange(len(regions))[:, None, None]*logv, axis=0)


def source_envelope(regions, active, length, dimension, logv, parity_logv=None):
    """Log upper bound on the exact nonempty-systematic source moment."""
    if parity_logv is None:
        parity_logv = logv
    systematic = dominate(regions, logv)
    parity = dominate(regions, parity_logv)
    vector = np.array([0., -math.inf])
    for region in range(length):
        matrix = systematic if region < dimension else parity
        vector = logsumexp(vector[:, None]+matrix, axis=0)
    scale = active*(log_nonzero_binomial(dimension, logv)
                    +(length-dimension)*(np.logaddexp(0., math.log(15)+parity_logv)-math.log(16)))
    return float(logsumexp(vector)+scale)


def envelope_bound(active, z, k_log=20, length=64, dimension=40, width=16, state=5, identity=True):
    geometry = packet.geometry(k_log, length, dimension, width)
    regions = packet.log_regions(z, active, geometry['rows'], width, state, identity)
    result = minimize_scalar(lambda lv: source_envelope(regions, active, length, dimension, lv),
                             bounds=(-25, 2), method='bounded')
    value = (log_choose(geometry['rows'], active)-math.log(15)
             -geometry['threshold']*math.log(z)+result.fun)
    return dict(margin_bits=-value/math.log(2), v=math.exp(result.x), z=z,
                rigorous=False, method='nonempty-source-regional-envelope')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outer', nargs=2, type=int, default=(64, 40))
    parser.add_argument('--active', nargs='+', type=int, default=(16, 64, 128))
    parser.add_argument('--state', type=int, default=5)
    parser.add_argument('--z', type=float, default=.99)
    args = parser.parse_args()
    for active in args.active:
        print(json.dumps(dict(active=active, **envelope_bound(active, args.z,
                         length=args.outer[0], dimension=args.outer[1], state=args.state))), flush=True)
