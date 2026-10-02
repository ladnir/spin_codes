"""Bounded, exact refinements of the shared-group support counts.

The first refinement averages exponential basis weights; floating arithmetic
selects a rational tilt only. The second continues positive primal/dual
shortening identities to a stated iteration limit. Neither changes the code.
"""
import argparse
from fractions import Fraction as Q
from math import gcd, log2, prod
from functools import reduce
from pathlib import Path
from time import monotonic
import json

import numpy as np
from scipy.optimize import brentq
from scipy.special import logsumexp
import shared_support
from basis_lattice import improve_caps
from bch_joint_support import rank_total, support_caps
from shortening_moments import improve
from dual_moments import dual_shell_caps, improve_from_dual
from dual_shortening import complement_caps


def basis_laplace_cap(spectrum, g, h, u, tilt):
    """Jensen bound on the number of rank-h ordered g-tuples at support <=u.

    Every h-dimensional subspace has B ordered bases. Their mean total
    basis weight is h*2^(h-1)*support/(2^h-1). Convexity lower-bounds the
    sum of tilt^(basis weight). Nonzero words in C^h upper-bound that sum
    over all independent bases. Multiply by the number of spanning g-tuples
    per subspace. Interpolation exploits the lattice of possible weights.
    """
    tilt = Q(tilt)
    if (not 1 <= h <= g or not 0 <= u < len(spectrum) or not 0 < tilt <= 1
            or any(type(v) is not int or v < 0 for v in spectrum)):
        raise ValueError('valid ranks, support, spectrum caps and tilt in (0,1] required')
    weights = [w for w, c in enumerate(spectrum) if w and c]
    if not weights:
        return 0
    lattice = reduce(gcd, weights)
    mean = Q(h*(1 << (h-1))*u, (1 << h)-1)
    if mean < h*weights[0]:
        return 0
    lo = (mean.numerator // (mean.denominator*lattice))*lattice
    fraction = (mean-lo)/lattice
    lower = tilt**lo*((1-fraction)+fraction*tilt**lattice)
    polynomial = Q(0)
    for c in reversed([0]+spectrum[1:]):
        polynomial = polynomial*tilt+c
    bases = prod((1 << h)-(1 << j) for j in range(h))
    embeddings = prod((1 << g)-(1 << j) for j in range(h))
    upper = embeddings*polynomial**h/(bases*lower)
    return upper.numerator // upper.denominator


def improve_basis(caps, spectrum):
    weights = np.array([w for w, c in enumerate(spectrum) if w and c])
    logs = np.array([log2(spectrum[w])*np.log(2) for w in weights])
    revised = [row[:] for row in caps]
    for h, row in enumerate(revised, 1):
        for u, previous in enumerate(row):
            if not previous:
                continue
            mean = h*(1 << (h-1))*u/((1 << h)-1)
            def derivative(lam):
                scores = logs-lam*weights
                return mean-h*float(np.exp(scores-logsumexp(scores)) @ weights)
            if derivative(0) >= 0:
                lam = 0.
            elif derivative(4) <= 0:
                lam = 4.
            else:
                lam = brentq(derivative, 0., 4.)
            tilt = Q(max(1, round(np.exp(-lam)*(1 << 20))), 1 << 20)
            row[u] = min(previous, basis_laplace_cap(spectrum, len(caps), h, u, tilt))
        for u in range(len(row)-2, -1, -1):
            row[u] = min(row[u], row[u+1])
    return revised


def refine_pair(ranks, spectrum, iterations=4, record=None, lp_through=104):
    """Regenerate all bounds; recorded numerical caps are never inputs."""
    if (type(iterations) is not int or not 1 <= iterations <= 30
            or type(lp_through) is not int or not 104 <= lp_through <= 192):
        raise ValueError('one through thirty iterations required')
    ranks = improve_basis(ranks, spectrum)
    if record: record('basis-laplace', ranks)
    from shortening_polynomial import improve_dimensions
    from dual_shortening import improve_dimensions as dual_dimensions
    dims = dual_dimensions(improve_dimensions(shared_support.dimension_caps(last_lp=lp_through)), last_lp=lp_through)
    ddims = complement_caps(shared_support.dimension_caps(256, 128, 30, last_lp=lp_through), dims, 128)
    if lp_through > 104:
        caps = improve_caps(support_caps(spectrum, g=4, dimensions=dims), spectrum)
        ranks = [[min(a, b) for a, b in zip(old, new)] for old, new in zip(ranks, caps)]
        ranks, _ = improve(ranks, dims)
        if record: record(f'longer-shortening-lp-{lp_through}', ranks)
    dspectrum = dual_shell_caps()
    dual = improve_caps(support_caps(dspectrum, g=4, dimensions=ddims), dspectrum)
    dual, _ = improve(dual, ddims)
    for iteration in range(1, iterations+1):
        before = ranks, dual
        ranks, _ = improve_from_dual(ranks, dual, 128)
        dual, _ = improve_from_dual(dual, ranks, 128)
        ranks, _ = improve(ranks, dims)
        dual, _ = improve(dual, ddims)
        if record: record(f'dual-containment-{iteration}', ranks)
        if (ranks, dual) == before:
            print('EXACT FIXED POINT at iteration', iteration, flush=True)
            break
    return ranks, dual, dims, ddims


def refine(ranks, spectrum, iterations=4, record=None, lp_through=104):
    return refine_pair(ranks, spectrum, iterations, record, lp_through)[0]


def run(iterations, output, lp_through=104):
    start = monotonic()
    baseline, ranks = shared_support.shared_counts(True)
    spectrum = shared_support.authenticated_caps()
    stages = []
    def record(name, revised):
        values = [sum(row[u] for row in revised) for u in range(257)]
        gain = lambda u: log2(baseline[u])-log2(values[u]) if values[u] else None
        points = {str(u): gain(u) for u in (96, 104, 112, 120, 128, 136, 144, 160, 192)}
        stages.append(dict(method=name, seconds=monotonic()-start, gains_bits=points,
                           counts=list(map(str, values)), rank_cdfs=[[str(v) for v in row] for row in revised]))
        print('EXACT COUNT STAGE', name, 'gains versus baseline', points, flush=True)
        if output:
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_text(json.dumps(dict(schema='shared-gf16-count-refinements-1',
                baseline=list(map(str, baseline)), lp_through=lp_through, stages=stages), indent=2)+'\n')
        return values
    refine(ranks, spectrum, iterations, record, lp_through)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--iterations', type=int, default=8)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--lp-through', type=int, default=104)
    args = parser.parse_args()
    run(args.iterations, args.output, args.lp_through)
