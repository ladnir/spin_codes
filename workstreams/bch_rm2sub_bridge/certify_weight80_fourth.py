"""Bound the mixed fourth-order dependence of two shared-permutation rows.

The optimization proposes only two rational dual coefficients. An exact
hinge majorant certifies the bound independently of optimizer correctness.
This is not a bound on the full zero-state second moment.
"""
import math
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from flint import arb, ctx
import bridge as base
import christoffel_caps as christoffel
import screen_exponential_modes as cap_source


def run():
    caps, sources = cap_source.latest_caps()
    family_path = base.HERE/'generated/weight80_family.json'
    family = base.read(family_path)
    for name, digest in family['local_sha256'].items():
        assert base.sha(base.HERE/name) == digest
    a = family['shell_size_lower']
    ctx.prec = 512
    christoffel.build(verify=True)
    from audit_bch_q1_full_arb import rational
    weights = [0] + list(range(38, 161, 2))
    upper = [F(1, a)]
    intersections = [None]
    for c in weights[1:]:
        i, j = c//2, 80-c//2
        kernel = sum((arb(christoffel.kraw(c, r, i)**2 *
                          christoffel.kraw(256-c, s, j)**2) /
                      (math.comb(c, r)*math.comb(256-c, s))
                      for r in range(15) for s in range(15-r)), arb(0)).lower()
        assert kernel > 0
        intersection = min(caps[80], math.comb(c, i)*math.comb(256-c, j),
                           math.floor(rational((arb(2)**128/kernel).upper())))
        intersections.append(intersection)
        upper.append(min(F(1), F(caps[c]*intersection, a*a)))
    objective = [(55-c//2)**2 for c in weights]
    solution = linprog(-np.array(objective, dtype=float),
                      A_eq=np.array([[1]*len(weights), weights], dtype=float),
                      b_eq=[1, 110], bounds=[(0, float(b)) for b in upper],
                      method='highs')
    assert solution.success
    intercept, slope = [F.from_float(float(-v)) for v in solution.eqlin.marginals]
    hinges = [max(F(0), f-intercept-slope*c) for c, f in zip(weights, objective)]
    assert all(intercept+slope*c+h >= f and h >= 0
               for c, f, h in zip(weights, objective, hinges))
    second = intercept + 110*slope + sum(b*h for b, h in zip(upper, hinges))
    n, p = 256, F(5, 16)
    v = p*(1-p)
    pair_mean = F(80*79, 256*255)
    elementary = pair_mean*(p-pair_mean)
    baseline = F(n*n, n-1)*v*v
    tau = min(elementary, (second-baseline)/(n*(n-1)))
    assert tau >= 0
    normalized = tau/(2620*v*v)
    assert normalized < F(1, 5000)
    rows = [dict(difference_weight=c, intersection_upper=i,
                 probability_upper=base.encode(b), dual_hinge=base.encode(h))
            for c, i, b, h in zip(weights, intersections, upper, hinges)]
    result = dict(status='EXACT_WEIGHT80_MIXED_FOURTH_BOUND',
                  occupation=2620, row_weight=80, rows=rows,
                  dual_intercept=base.encode(intercept), dual_slope=base.encode(slope),
                  centered_inner_product_second_moment_upper=base.encode(second),
                  off_diagonal_covariance_variance_upper=base.encode(tau),
                  elementary_variance_upper=base.encode(elementary),
                  normalized_shared_edge_mixed_fourth_upper=base.encode(normalized),
                  full_second_moment_certified=False,
                  local_sha256={str(path.relative_to(base.HERE)):base.sha(path)
                                for path in [Path(__file__), Path(christoffel.__file__),
                                             family_path]+sources})
    path = base.HERE/'generated/weight80_fourth.json'
    if path.exists():
        assert result == base.read(path)
    else:
        base.write_new(path, result)
    print('Exact fourth bound: E(inner product squared) <=', float(second),
          'tau <=', float(tau), 'elementary tau <=', float(elementary),
          'normalized mixed fourth <=', float(normalized), flush=True)


if __name__ == '__main__':
    run()
