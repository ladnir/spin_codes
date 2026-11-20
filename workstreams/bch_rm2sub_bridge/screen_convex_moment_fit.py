"""Optimize the exponential moment of a global convex envelope; screen only.

Frank-Wolfe updates preserve linear feasibility. Each update uses an LP
and scalar line search. All outputs remain numerical proposals.
"""
import math
from pathlib import Path
import numpy as np
from scipy.optimize import linprog, minimize_scalar
from scipy.special import gammaln, logsumexp
from scipy.sparse import diags, vstack, csr_matrix
import bridge as base


def run():
    output = base.HERE/'generated/convex_moment_fit_screen.json'
    assert not output.exists()
    source = base.HERE/'generated/global_overlap_optimized_screen.json'
    g = np.array(base.read(source)['correction_robust_log_bound_screen'])
    length, target, extension_slope = len(g), 20800, 1.8
    curvature = diags([-np.ones(length-2), 2*np.ones(length-2), -np.ones(length-2)],
                      [0, 1, 2], shape=(length-2, length), format='csr')
    slopes = csr_matrix(([1., -1., -1., 1.], ([0, 0, 1, 1], [0, 1, length-2, length-1])),
                        shape=(2, length))
    constraints = vstack([curvature, slopes])
    rhs = np.r_[np.zeros(length-2), .7, extension_slope]
    k = np.arange(2610*80+1)
    factorial = gammaln(k+1)
    def oracle(weights):
        fit = linprog(weights+1e-14, A_ub=constraints, b_ub=rhs,
                      bounds=list(zip(g, [None]*length)), method='highs')
        assert fit.success
        return fit.x
    mean = target/256
    psi = oracle(np.exp(-mean+k[:length]*np.log(mean)-factorial[:length]))
    def extension(values, a):
        return np.r_[values, values[-1]+extension_slope*np.arange(1, max(1, a+2-length))][:a+1]
    def inner(values, a, theta):
        return gammaln(a+1)-a*(math.log(256)+theta)+256*logsumexp(
            extension(values, a)-factorial[:a+1]+k[:a+1]*theta)
    def best_tilt(values, a):
        return minimize_scalar(lambda theta:inner(values, a, theta), bounds=(-12, 12), method='bounded',
                               options={'xatol':1e-10})
    history = []
    for iteration in range(101):
        tilt = best_tilt(psi, target)
        logweights = extension(psi, target)-factorial[:target+1]+k[:target+1]*tilt.x
        weights = np.exp(logweights-logsumexp(logweights))
        gradient = weights[:length].copy()
        gradient[-1] += sum(weights[length:])
        vertex = oracle(gradient)
        gap = 256*float(gradient@(psi-vertex))
        if iteration % 10 == 0 or iteration == 100:
            history.append(dict(iteration=iteration, target_bound_bits=float(tilt.fun/math.log(2)), gap=gap))
            print('Convex moment fit', history[-1], flush=True)
        if gap < 1e-7 or iteration == 100:
            break
        direction = vertex-psi
        step = minimize_scalar(lambda eta:inner(psi+eta*direction, target, tilt.x),
                               bounds=(0, 1), method='bounded', options={'xatol':1e-10})
        psi += step.x*direction
    assert max(float(np.max(g-psi)), float(np.max(constraints@psi-rhs))) < 1e-6
    law_path = base.HERE/'generated/overlap_convex_law.json'
    law = base.read(law_path)
    atoms = np.array([r['overlap'] for r in law['rows']])
    mass = np.array([float(base.decode(r['probability'])) for r in law['rows']])*2610/8189
    mass[0] += 1-2610/8189
    rows = []
    for a in [10000, 15000, 19000, 20000, 20800, 22000, 25000, 30000, 40000, 60000, 80000, 100000, 150000, 200000]:
        first = best_tilt(psi, a)
        second = minimize_scalar(lambda theta:2610*logsumexp(np.log(mass)+atoms*theta)-a*theta,
                                 bounds=(-6, 6), method='bounded')
        rows.append(dict(total_core_overlap=a, multinomial_bound_log_screen=float(first.fun),
                         outer_coefficient_bound_log_screen=float(min(0, second.fun)),
                         combined_summand_bound_bits_screen=float((first.fun+min(0, second.fun))/math.log(2))))
    base.write_new(output, dict(status='EXPONENTIAL_MOMENT_OPTIMIZED_CONVEX_ENVELOPE_SCREEN_ONLY',
        proposed_log_envelope=psi.tolist(), affine_extension_slope=extension_slope,
        optimization_history=history, rows=rows, full_second_moment_certified=False,
        actual_second_moment_estimated=False,
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in [Path(__file__), source, law_path]}))
    print('Selected combined summand bound bits:',
          [(r['total_core_overlap'], round(r['combined_summand_bound_bits_screen'], 3)) for r in rows], flush=True)


if __name__ == '__main__':
    run()
