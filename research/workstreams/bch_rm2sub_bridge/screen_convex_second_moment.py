"""Discovery-only convex envelope and two scalar Chernoff layers.

This evaluates proposed upper bounds, not actual second moments. A large
value demonstrates slack in this proposal, not actual setup failure.
"""
import math
from pathlib import Path
import numpy as np
from scipy.optimize import linprog, minimize_scalar
from scipy.special import gammaln, logsumexp
from scipy.sparse import diags, vstack, csr_matrix
import bridge as base


def run():
    path = base.HERE/'generated/convex_second_moment_screen.json'
    assert not path.exists()
    source = base.HERE/'generated/global_overlap_optimized_screen.json'
    proposal = base.read(source)
    g = np.array(proposal['correction_robust_log_bound_screen'])
    length = len(g)
    curvature = diags([-np.ones(length-2), 2*np.ones(length-2), -np.ones(length-2)],
                      [0, 1, 2], shape=(length-2, length), format='csr')
    slope = csr_matrix(([-1., 1.], ([0, 0], [length-2, length-1])), shape=(1, length))
    k = np.arange(length)
    mean = 2610**2/8189*25/256
    objective = np.exp(-mean+k*np.log(mean)-gammaln(k+1))+1e-8
    fit = linprog(objective, A_ub=vstack([curvature, slope]),
                  b_ub=np.r_[np.zeros(length-2), 3.], bounds=list(zip(g, [None]*length)), method='highs')
    assert fit.success
    psi = fit.x
    maximum_violation = max(float(np.max(g-psi)), float(np.max(curvature@psi)))
    assert maximum_violation < 1e-6
    final_slope = psi[-1]-psi[-2]
    core_max = 2610*80
    all_k = np.arange(core_max+1)
    extension = np.r_[psi, psi[-1]+final_slope*np.arange(1, core_max+2-length)]
    log_coefficients = extension-gammaln(all_k+1)
    law_source = base.HERE/'generated/overlap_convex_law.json'
    law = base.read(law_source)
    atoms = np.array([r['overlap'] for r in law['rows']])
    masses = np.array([float(base.decode(r['probability'])) for r in law['rows']])
    p = 2610/8189
    masses *= p
    masses[0] += 1-p
    log_mass = np.log(masses)
    rows = []
    for a in [10000, 15000, 19000, 20000, 20800, 22000, 25000, 30000, 40000, 60000, 80000, 100000, 150000, 200000]:
        def inner(theta):
            return gammaln(a+1)-a*(math.log(256)+theta)+256*logsumexp(log_coefficients[:a+1]+all_k[:a+1]*theta)
        def outer(theta):
            return 2610*logsumexp(log_mass+atoms*theta)-a*theta
        first = minimize_scalar(inner, bounds=(-12, 12), method='bounded', options={'xatol':1e-10})
        second = minimize_scalar(outer, bounds=(-6, 6), method='bounded', options={'xatol':1e-10})
        total = first.fun+min(0, second.fun)
        rows.append(dict(total_core_overlap=a, multinomial_bound_log_screen=float(first.fun),
                         outer_coefficient_bound_log_screen=float(min(0, second.fun)),
                         combined_summand_bound_bits_screen=float(total/math.log(2)),
                         inner_log_tilt=float(first.x), outer_log_tilt=float(second.x)))
    base.write_new(path, dict(status='CONVEX_SECOND_MOMENT_PROPOSAL_SCREEN_ONLY',
        proposed_log_envelope=psi.tolist(), affine_extension_slope=float(final_slope),
        maximum_lp_constraint_violation=maximum_violation, rows=rows,
        full_second_moment_certified=False, actual_second_moment_estimated=False,
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in [Path(__file__), source, law_source]}))
    print('Convex proposal selected summand upper-bound bits:',
          [(r['total_core_overlap'], round(r['combined_summand_bound_bits_screen'], 3)) for r in rows], flush=True)


if __name__ == '__main__':
    run()
