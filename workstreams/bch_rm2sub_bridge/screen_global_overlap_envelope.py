"""Discovery-only global region-ratio bounds, including both overlap tails.

Positive Chernoff tilts cover every feasible x,y,k in the central marginal
window. The retained atlas improves k40..160. Nothing here is outward.
"""
import math
from pathlib import Path
import numpy as np
from scipy.special import gammaln
from flint import arb_poly, ctx
import bridge as base
import tightened_occupancy as tight
from screen_pair_type_bound import load_pairs, evaluate


def run():
    path = base.HERE/'generated/global_overlap_envelope_screen.json'
    assert not path.exists()
    powers, counts = load_pairs()
    ctx.prec = 192
    spectrum = tight.kernel_spectrum()
    kernel = arb_poly([spectrum.get(j, 0) for j in range(129)])**64
    weights = np.arange(740, 901, 2)
    single = np.array([float((kernel[int(j)]/math.comb(8192, int(j))).log()) for j in weights])
    # Nine narrow marginal blocks; all unordered blocks cover the transpose.
    blocks = [np.arange(i, min(i+10, 81)) for i in range(0, 81, 10)]
    lf = gammaln(np.arange(8193)+1)
    maxima = np.full(901, -np.inf)
    witnesses = [None]*901
    for bi, first in enumerate(blocks):
        for second in blocks[bi:]:
            maximum_k = int(min(weights[first[-1]], weights[second[-1]]))
            ks = np.arange(maximum_k+1)
            centers_x = np.maximum(float(np.mean(weights[first])), ks)
            centers_y = np.maximum(float(np.mean(weights[second])), ks)
            cells = np.array([8192-centers_x-centers_y+ks, centers_y-ks,
                              centers_x-ks, ks], dtype=float)
            cells = np.maximum(cells, .25)
            p = cells/cells.sum(axis=0)
            log_kernel = 64*np.log(evaluate(p, powers, counts))
            logp = np.log(p)
            x, y = np.meshgrid(weights[first], weights[second], indexing='ij')
            singles = single[first, None]+single[None, second]
            trivial = np.minimum(single[first, None], single[None, second])-singles
            for k in ks:
                feasible = (x >= k) & (y >= k)
                m = [8192-x-y+k, np.maximum(y-k, 0), np.maximum(x-k, 0),
                     np.full(x.shape, k)]
                bound = log_kernel[k]-lf[8192]-singles
                for cell, lp in zip(m, logp[:, k]):
                    bound = bound+lf[cell]-cell*lp
                bound = np.minimum(bound, trivial)
                if 40 <= k <= 160:
                    bound = np.minimum(bound, .00012*(k-x*y/8192)**2)
                bound = np.where(feasible, bound, -np.inf)
                index = np.unravel_index(np.argmax(bound), bound.shape)
                maximum = float(bound[index])
                if maximum > maxima[k]:
                    maxima[k] = maximum
                    witnesses[k] = [int(x[index]), int(y[index]), int(k)]
        print('Global overlap screen marginal block', bi, 'complete', flush=True)
    assert np.all(np.isfinite(maxima))
    corrected = [float(max(maxima[k:min(k+4, 901)])) for k in range(898)]
    base.write_new(path, dict(status='GLOBAL_REGION_OVERLAP_SCREEN_ONLY',
        all_actual_overlaps=list(range(901)), marginal_weights=weights.tolist(),
        rows=[dict(overlap=k, log_ratio_upper_screen=float(maxima[k]), worst_type=witnesses[k])
              for k in range(901)],
        correction_robust_log_bound_screen=corrected, outward_certified=False,
        full_second_moment_certified=False,
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in
                      [Path(__file__), Path(tight.__file__), base.HERE/'screen_pair_type_bound.py',
                       base.HERE/'generated/t128_s15_pair_spectrum.json',
                       base.HERE/'generated/overlap_atlas_coverage.json']}))
    print('All overlaps screened:', [(k, round(float(maxima[k]), 6))
                                      for k in (0, 20, 39, 40, 80, 120, 160, 161, 200, 300, 400, 600, 740, 900)], flush=True)


if __name__ == '__main__':
    run()
