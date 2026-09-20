"""Per-total-overlap convex envelopes plus actual-MGF bounds; discovery only."""
import math
from pathlib import Path
import numpy as np
from scipy.optimize import linprog, minimize_scalar
from scipy.special import gammaln, logsumexp
from scipy.sparse import diags, vstack, csr_matrix
import bridge as base


def run():
    output = base.HERE/'generated/adaptive_low_overlap_screen.json'
    assert not output.exists()
    source = base.HERE/'generated/global_overlap_optimized_screen.json'
    coverage_path = base.HERE/'generated/extended_overlap_low_coverage.json'
    coverage = base.read(coverage_path)
    assert coverage['actual_overlaps'] == list(range(161))
    actual = np.array([r['log_ratio_upper_screen'] for r in base.read(source)['rows']])
    for k in range(161):
        actual[k] = .00012*max((k-740**2/8192)**2, (k-900**2/8192)**2)
    g = np.array([max(actual[k:k+4]) for k in range(898)])
    length = len(g)
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
    for a in [10000, 15000, 19000, 20000, 20800, 22000, 25000, 30000, 40000, 60000, 80000, 100000, 150000, 200000]:
        def outer(theta):
            return logsumexp(support_logmass+r*logsumexp(logmass+atoms*theta))-a*theta
        mass_bound = min(0, minimize_scalar(outer, bounds=(-6, 6), method='bounded').fun)
        trials, best_psi, best_bound = [], None, math.inf
        for right in [0., .5, 1., 1.5, 2., 2.5, 3.]:
            rhs = np.r_[np.zeros(length-2), .7, right]
            def oracle(weight):
                fit = linprog(weight+1e-13, A_ub=constraints, b_ub=rhs,
                              bounds=list(zip(g, [None]*length)), method='highs')
                assert fit.success
                return fit.x
            def extended(psi):
                return np.r_[psi, psi[-1]+right*np.arange(1, a+2-length)]
            def inner(psi, theta):
                return gammaln(a+1)-a*(math.log(256)+theta)+256*logsumexp(extended(psi)-factorial[:a+1]+k[:a+1]*theta)
            def best_tilt(psi):
                return minimize_scalar(lambda theta:inner(psi, theta), bounds=(-12, 12), method='bounded')
            mean = a/256
            psi = oracle(np.exp(-mean+k[:length]*math.log(mean)-factorial[:length]))
            for _ in range(12):
                tilt = best_tilt(psi)
                logweights = extended(psi)-factorial[:a+1]+k[:a+1]*tilt.x
                weights = np.exp(logweights-logsumexp(logweights))
                gradient = weights[:length].copy()
                gradient[-1] += sum(weights[length:])
                vertex = oracle(gradient)
                if 256*float(gradient@(psi-vertex)) < 1e-6:
                    break
                direction = vertex-psi
                step = minimize_scalar(lambda eta:inner(psi+eta*direction, tilt.x), bounds=(0, 1), method='bounded')
                psi += step.x*direction
            tilt = best_tilt(psi)
            assert max(float(np.max(g-psi)), float(np.max(constraints@psi-rhs))) < 1e-5
            trials.append(dict(extension_slope=right, inner_upper_log_screen=float(tilt.fun)))
            if tilt.fun < best_bound:
                best_bound, best_psi, chosen_slope = float(tilt.fun), psi, right
        row = dict(total_core_overlap=a, combined_summand_bound_bits_screen=float((best_bound+mass_bound)/math.log(2)),
                   actual_mass_bound_log_screen=float(mass_bound), inner_bound_log_screen=best_bound,
                   selected_extension_slope=chosen_slope, proposed_log_envelope=best_psi.tolist(), trials=trials)
        rows.append(row)
        print('Adaptive overlap', a, 'combined bound bits', row['combined_summand_bound_bits_screen'],
              'slope', chosen_slope, flush=True)
    base.write_new(output, dict(status='ADAPTIVE_ENVELOPE_AND_ACTUAL_MGF_SCREEN_ONLY', rows=rows,
        full_second_moment_certified=False, actual_second_moment_estimated=False,
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in [Path(__file__), source, coverage_path, law_path]}))


if __name__ == '__main__':
    run()
