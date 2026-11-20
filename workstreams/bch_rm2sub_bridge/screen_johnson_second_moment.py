"""Screen the improved convex row law with the exact support-intersection law."""
import math
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import gammaln, logsumexp
import bridge as base


def run():
    output = base.HERE/'generated/johnson_second_moment_screen.json'
    assert not output.exists()
    envelope_path = base.HERE/'generated/convex_moment_fit_screen.json'
    envelope = base.read(envelope_path)
    psi = np.array(envelope['proposed_log_envelope'])
    k = np.arange(2610*80+1)
    extension = np.r_[psi, psi[-1]+envelope['affine_extension_slope']*np.arange(1, len(k)+1-len(psi))]
    coefficient = extension-gammaln(k+1)
    law_path = base.HERE/'generated/johnson_convex_law.json'
    law = base.read(law_path)
    atoms = np.array([r['overlap'] for r in law['rows']])
    logmass = np.log([float(base.decode(r['probability'])) for r in law['rows']])
    r = np.arange(2611)
    def logcomb(n, j):
        return gammaln(n+1)-gammaln(j+1)-gammaln(n-j+1)
    support_logmass = logcomb(2610, r)+logcomb(5579, 2610-r)-logcomb(8189, 2610)
    assert abs(logsumexp(support_logmass)) < 1e-9
    rows = []
    for a in [10000, 15000, 19000, 20000, 20800, 22000, 25000, 30000, 40000, 60000, 80000, 100000, 150000, 200000]:
        def inner(theta):
            return gammaln(a+1)-a*(math.log(256)+theta)+256*logsumexp(coefficient[:a+1]+k[:a+1]*theta)
        def outer(theta):
            log_mgf = logsumexp(logmass+atoms*theta)
            return logsumexp(support_logmass+r*log_mgf)-a*theta
        first = minimize_scalar(inner, bounds=(-12, 12), method='bounded')
        second = minimize_scalar(outer, bounds=(-6, 6), method='bounded')
        rows.append(dict(total_core_overlap=a, multinomial_bound_log_screen=float(first.fun),
                         hypergeometric_coefficient_bound_log_screen=float(min(0, second.fun)),
                         combined_summand_bound_bits_screen=float((first.fun+min(0, second.fun))/math.log(2))))
    base.write_new(output, dict(status='JOHNSON_CONVEX_LAW_SECOND_MOMENT_SCREEN_ONLY', rows=rows,
        full_second_moment_certified=False, actual_second_moment_estimated=False,
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in [Path(__file__), envelope_path, law_path]}))
    print('Improved-law combined summand bound bits:',
          [(r['total_core_overlap'], round(r['combined_summand_bound_bits_screen'], 3)) for r in rows], flush=True)


if __name__ == '__main__':
    run()
