"""Optimized positive-tilt global overlap screen, retaining replay proposals."""
import math
from pathlib import Path
import numpy as np
from scipy.special import gammaln
from flint import arb_poly, ctx
import bridge as base
import tightened_occupancy as tight
import pair_positive_tilt as tilt


def run():
    path = base.HERE/'generated/global_overlap_retiled_screen.json'
    assert not path.exists()
    fitter = tilt.PairTilt()
    ctx.prec = 192
    spectrum = tight.kernel_spectrum()
    kernel = arb_poly([spectrum.get(j, 0) for j in range(129)])**64
    weights = np.arange(740, 901, 2)
    single = np.array([float((kernel[int(j)]/math.comb(8192, int(j))).log()) for j in weights])
    blocks = [np.arange(i, min(i+17, 81)) for i in range(0, 81, 17)]
    lf = gammaln(np.arange(8193)+1)
    maxima = np.full(901, -np.inf)
    witnesses, proposals = [None]*901, []
    for bi, first in enumerate(blocks):
        for bj in range(bi, len(blocks)):
            second = blocks[bj]
            maximum_k = int(min(weights[first[-1]], weights[second[-1]]))
            x, y = np.meshgrid(weights[first], weights[second], indexing='ij')
            singles = single[first, None]+single[None, second]
            trivial = np.minimum(single[first, None], single[None, second])-singles
            saved_tilts = []
            initial = None
            for k in range(maximum_k+1):
                cx = float(np.mean(weights[first][weights[first] >= k]))
                cy = float(np.mean(weights[second][weights[second] >= k]))
                center = np.array([8192-cx-cy+k, cy-k, cx-k, k])
                p, objective, initial = fitter.fit(center, initial)
                logp = np.log(p)
                log_kernel = objective+center@logp
                saved_tilts.append(p[1:].tolist())
                feasible = (x >= k) & (y >= k)
                m = [8192-x-y+k, np.maximum(y-k, 0), np.maximum(x-k, 0),
                     np.full(x.shape, k)]
                bound = log_kernel-lf[8192]-singles
                for cell, lp in zip(m, logp):
                    bound = bound+lf[cell]-cell*lp
                bound = np.minimum(bound, trivial)
                if 40 <= k <= 160:
                    bound = np.minimum(bound, .00012*(k-x*y/8192)**2)
                bound = np.where(feasible, bound, -np.inf)
                index = np.unravel_index(np.argmax(bound), bound.shape)
                maximum = float(bound[index])
                if maximum > maxima[k]:
                    maxima[k] = maximum
                    witnesses[k] = [int(x[index]), int(y[index]), k]
            proposals.append(dict(first_weights=weights[first].tolist(),
                                  second_weights=weights[second].tolist(),
                                  nonzero_tilt_probabilities_by_overlap=saved_tilts))
            print('Optimized global overlap block', bi, bj, 'complete', flush=True)
    assert np.all(np.isfinite(maxima))
    base.write_new(path, dict(status='OPTIMIZED_GLOBAL_REGION_OVERLAP_SCREEN_ONLY',
        rows=[dict(overlap=k, log_ratio_upper_screen=float(maxima[k]), worst_type=witnesses[k])
              for k in range(901)],
        correction_robust_log_bound_screen=[float(max(maxima[k:min(k+4, 901)])) for k in range(898)],
        tilt_proposals=proposals, outward_certified=False, full_second_moment_certified=False,
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in
                      [Path(__file__), Path(tilt.__file__), Path(tight.__file__),
                       base.HERE/'generated/t128_s15_pair_spectrum.json',
                       base.HERE/'generated/overlap_atlas_replay.json']}))
    print('Optimized global screen:', [(k, round(float(maxima[k]), 6))
                                       for k in (0, 20, 39, 40, 80, 120, 160, 161, 200, 300, 400, 600, 740, 900)], flush=True)


if __name__ == '__main__':
    run()
