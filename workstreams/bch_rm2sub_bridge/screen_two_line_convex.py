"""Three-variable optimization of a global convex envelope; screen only."""
import math
from pathlib import Path
import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.special import gammaln, logsumexp
import bridge as base


def run():
    output = base.HERE/'generated/two_line_convex_screen.json'
    assert not output.exists()
    source = base.HERE/'generated/global_overlap_optimized_screen.json'
    g = np.array(base.read(source)['correction_robust_log_bound_screen'])
    k = np.arange(2610*80+1)
    q = np.minimum(k, 897)
    lo, hi = 740**2/8192, 900**2/8192
    central = .00012*np.maximum((q-hi)**2, (q+3-lo)**2)
    central += np.maximum(k-897, 0)*(.00024*(900-lo))
    factorial = gammaln(k+1)
    def envelope(left, right):
        intercept_l = max(g[:40]+left*np.arange(40))
        intercept_r = max(g[158:]-right*np.arange(158, 898))
        return np.maximum(central, np.maximum(intercept_l-left*k, intercept_r+right*k))
    target = 20800
    def objective(parameters):
        left, right, theta = parameters
        psi = envelope(left, right)
        return gammaln(target+1)-target*(math.log(256)+theta)+256*logsumexp(
            psi[:target+1]-factorial[:target+1]+k[:target+1]*theta)
    fit = minimize(objective, [.35, 1.8, math.log(target/256)], method='Powell',
                   bounds=[(.1, 3), (1, 3), (2, 6)], options={'xtol':1e-9, 'ftol':1e-10})
    psi = envelope(*fit.x[:2])
    assert np.max(g-psi[:898]) < 1e-9
    assert np.min(np.diff(psi[:898], n=2)) > -1e-9
    law_path = base.HERE/'generated/overlap_convex_law.json'
    law = base.read(law_path)
    atoms = np.array([r['overlap'] for r in law['rows']])
    mass = np.array([float(base.decode(r['probability'])) for r in law['rows']])*2610/8189
    mass[0] += 1-2610/8189
    logmass = np.log(mass)
    rows = []
    for a in [10000, 15000, 19000, 20000, 20800, 22000, 25000, 30000, 40000, 60000, 80000, 100000, 150000, 200000]:
        def inner(theta):
            return gammaln(a+1)-a*(math.log(256)+theta)+256*logsumexp(psi[:a+1]-factorial[:a+1]+k[:a+1]*theta)
        def outer(theta):
            return 2610*logsumexp(logmass+atoms*theta)-a*theta
        first = minimize_scalar(inner, bounds=(-12, 12), method='bounded')
        second = minimize_scalar(outer, bounds=(-6, 6), method='bounded')
        rows.append(dict(total_core_overlap=a, multinomial_bound_log_screen=float(first.fun),
                         outer_coefficient_bound_log_screen=float(min(0, second.fun)),
                         combined_summand_bound_bits_screen=float((first.fun+min(0, second.fun))/math.log(2))))
    base.write_new(output, dict(status='TWO_LINE_CONVEX_SECOND_MOMENT_SCREEN_ONLY',
        left_slope_magnitude=float(fit.x[0]), right_slope=float(fit.x[1]),
        proposed_log_envelope=psi[:898].tolist(), rows=rows,
        full_second_moment_certified=False, actual_second_moment_estimated=False,
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in [Path(__file__), source, law_path]}))
    print('Two-line envelope slopes', fit.x[:2].tolist(), 'selected summand bound bits',
          [(r['total_core_overlap'], round(r['combined_summand_bound_bits_screen'], 3)) for r in rows], flush=True)


if __name__ == '__main__':
    run()
