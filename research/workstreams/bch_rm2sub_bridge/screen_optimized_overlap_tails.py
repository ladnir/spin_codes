"""Probe positive-tilt improvement at uncovered overlap types; no certificate."""
import math
from pathlib import Path
import numpy as np
from scipy.special import gammaln
from flint import arb_poly, ctx
import bridge as base
import tightened_occupancy as tight
import pair_positive_tilt as tilt


def run():
    path = base.HERE/'generated/optimized_overlap_tails_screen.json'
    assert not path.exists()
    fitter = tilt.PairTilt()
    m = np.array([6638., 736., 736., 82.])
    theta = np.log(m[1:]/m[0])
    value, gradient = fitter.value_gradient(theta, m)
    for j in range(3):
        direction = np.zeros(3)
        direction[j] = 1e-5
        numerical = (fitter.value_gradient(theta+direction, m)[0]-
                     fitter.value_gradient(theta-direction, m)[0])/2e-5
        assert abs(numerical-gradient[j]) < 1e-5
    ctx.prec = 192
    spectrum = tight.kernel_spectrum()
    kernel = arb_poly([spectrum.get(j, 0) for j in range(129)])**64
    singles = {x:float((kernel[x]/math.comb(8192, x)).log()) for x in (740, 820, 900)}
    rows = []
    for x, y in [(740, 740), (820, 820), (900, 900), (740, 900)]:
        for k in [0, 20, 39, 161, 200, 400, 600, 740]:
            m = np.array([8192-x-y+k, y-k, x-k, k], dtype=float)
            p, objective, _ = fitter.fit(m)
            bound = objective-gammaln(8193)+sum(gammaln(m+1))-singles[x]-singles[y]
            rows.append(dict(first_weight=x, second_weight=y, overlap=k,
                             optimized_chernoff_log_ratio_screen=float(bound),
                             tilt_probabilities=p.tolist()))
    base.write_new(path, dict(status='OPTIMIZED_OVERLAP_TAIL_SCREEN_ONLY', rows=rows,
        analytic_gradient_finite_difference_test_passed=True, outward_certified=False,
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in
                      [Path(__file__), Path(tilt.__file__), Path(tight.__file__),
                       base.HERE/'generated/t128_s15_pair_spectrum.json']}))
    print('Optimized positive bounds, x=y740:',
          [(r['overlap'], round(r['optimized_chernoff_log_ratio_screen'], 6))
           for r in rows if r['first_weight'] == r['second_weight'] == 740], flush=True)


if __name__ == '__main__':
    run()
