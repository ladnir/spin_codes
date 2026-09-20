"""Per-A fixed-total multipliers and convex envelopes; numerical screening."""
import math
from pathlib import Path
import numpy as np
from scipy.optimize import linprog, minimize_scalar
from scipy.special import gammaln, logsumexp
from scipy.sparse import diags, vstack, csr_matrix
import bridge as base


def run():
    output = base.HERE/'generated/weighted_adaptive_moment_screen.json'
    assert not output.exists()
    curve_path = base.HERE/'generated/weighted_overlap_curves_screen.json'
    saved = base.read(curve_path)
    for name, digest in saved['local_sha256'].items():
        assert base.sha(base.HERE/name) == digest
    curves = [row for row in saved['curves'] if row['multiplier'] in (0, -.0625, -.25, -.5, -1, -2)]
    length = 898
    curvature = diags([-np.ones(length-2), 2*np.ones(length-2), -np.ones(length-2)],
                      [0, 1, 2], shape=(length-2, length), format='csr')
    slopes = csr_matrix(([1., -1., -1., 1.], ([0, 0, 1, 1], [0, 1, length-2, length-1])), shape=(2, length))
    constraints = vstack([curvature, slopes])
    k = np.arange(208801)
    factorial = gammaln(k+1)
    law_path = base.HERE/'generated/johnson_convex_law.json'
    law = base.read(law_path)
    atoms = np.array([r['overlap'] for r in law['rows']])
    logmass = np.log([float(base.decode(r['probability'])) for r in law['rows']])
    r = np.arange(2611)
    def logcomb(n, j):
        return gammaln(n+1)-gammaln(j+1)-gammaln(n-j+1)
    support_logmass = logcomb(2610, r)+logcomb(5579, 2610-r)-logcomb(8189, 2610)
    rows = []
    for a in [10000, 20800, 30000, 40000, 60000, 80000, 100000, 150000, 189440, 200000]:
        def outer(theta):
            return logsumexp(support_logmass+r*logsumexp(logmass+atoms*theta))-a*theta
        mass_fit = minimize_scalar(outer, bounds=(-6, 6), method='bounded')
        mass_bound = min(0, mass_fit.fun)
        trials, best = [], None
        mean = a/256
        objective = np.exp(-mean+k[:length]*math.log(mean)-factorial[:length])+1e-10
        for curve in curves:
            g = np.array(curve['core_log_bounds'])
            for right in [0., .5, 1., 1.5, 2., 3., 4., 6.]:
                rhs = np.r_[np.zeros(length-2), .7, right]
                fit = linprog(objective, A_ub=constraints, b_ub=rhs,
                              bounds=list(zip(g, [None]*length)), method='highs')
                assert fit.success
                psi = fit.x
                assert max(float(np.max(g-psi)), float(np.max(constraints@psi-rhs))) < 1e-5
                extension = np.r_[psi, psi[-1]+right*np.arange(1, a+2-length)]
                coefficient = extension-factorial[:a+1]
                def inner(theta):
                    return gammaln(a+1)-a*(math.log(256)+theta)+256*logsumexp(coefficient+k[:a+1]*theta)
                tilt = minimize_scalar(inner, bounds=(-12, 12), method='bounded')
                trial = dict(multiplier=curve['multiplier'], extension_slope=right,
                             inner_upper_log_screen=float(tilt.fun), inner_log_tilt=float(tilt.x))
                trials.append(trial)
                if best is None or tilt.fun < best['inner_upper_log_screen']:
                    best = dict(trial, proposed_log_envelope=psi.tolist())
        row = dict(total_core_overlap=a, best=best, trials=trials,
                   actual_mass_bound_log_screen=float(mass_bound), mass_log_tilt=float(mass_fit.x),
                   combined_summand_bound_bits_screen=float((best['inner_upper_log_screen']+mass_bound)/math.log(2)))
        rows.append(row)
        print('Weighted adaptive A', a, 'bound bits', row['combined_summand_bound_bits_screen'],
              'b', best['multiplier'], 'slope', best['extension_slope'], flush=True)
    base.write_new(output, dict(status='FIXED_TOTAL_WEIGHTED_ADAPTIVE_MOMENT_SCREEN_ONLY', rows=rows,
        full_second_moment_certified=False, actual_second_moment_estimated=False,
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in [Path(__file__), curve_path, law_path]}))


if __name__ == '__main__':
    run()
