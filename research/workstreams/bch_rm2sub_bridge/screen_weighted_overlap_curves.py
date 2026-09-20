"""Screen all region types with multipliers preserving core marginal totals."""
import math
from pathlib import Path
import numpy as np
from flint import arb_poly, ctx
from scipy.special import gammaln
import bridge as base
import tightened_occupancy as tight
from screen_pair_type_bound import load_pairs, evaluate


def run():
    output = base.HERE/'generated/weighted_overlap_curves_screen.json'
    assert not output.exists()
    proposals_path = base.HERE/'generated/global_overlap_optimized_screen.json'
    proposals = base.read(proposals_path)
    for name, digest in proposals['local_sha256'].items():
        assert base.sha(base.HERE/name) == digest
    atlas_path = base.HERE/'generated/overlap_envelope_201_coverage.json'
    atlas = base.read(atlas_path)
    for name, digest in atlas['local_sha256'].items():
        assert base.sha(base.HERE/name) == digest
    intercepts = list(map(lambda v:float(base.decode(v)), atlas['intercept_by_overlap']))
    ctx.prec = 192
    spectrum = tight.kernel_spectrum()
    kernel = arb_poly([spectrum.get(j, 0) for j in range(129)])**64
    singles = {x:float((kernel[x]/math.comb(8192, x)).log()) for x in range(740, 901, 2)}
    multipliers = np.array([0., -1/1024, -1/256, -1/64, -1/16, -1/8, -1/4, -1/2, -1., -2.])
    maxima = np.full((len(multipliers), 901), -np.inf)
    worst = np.zeros((len(multipliers), 901, 2), dtype=np.int64)
    factorial = gammaln(np.arange(8193)+1)
    powers, counts = load_pairs()
    for block, proposal in enumerate(proposals['tilt_proposals']):
        xs, ys = proposal['first_weights'], proposal['second_weights']
        x, y = np.meshgrid(xs, ys, indexing='ij')
        single = np.array([[singles[a]+singles[b] for b in ys] for a in xs])
        trivial = np.minimum(np.array([singles[a] for a in xs])[:, None],
                             np.array([singles[b] for b in ys])[None, :])-single
        weighted_offset = -multipliers[:, None, None]*(x+y-1631.25)[None, :, :]
        p3 = np.array(proposal['nonzero_tilt_probabilities_by_overlap']).T
        p = np.vstack([1-p3.sum(axis=0), p3])
        assert np.all(p > 0)
        logp = np.log(p)
        logkernel = 64*np.log(evaluate(p, powers, counts))
        for k in range(p.shape[1]):
            cells = [8192-x-y+k, np.maximum(y-k, 0), np.maximum(x-k, 0), np.full(x.shape, k)]
            bound = logkernel[k]-factorial[8192]-single
            for cell, lp in zip(cells, logp[:, k]):
                bound = bound+factorial[cell]-cell*lp
            bound = np.minimum(bound, trivial)
            if k <= 200:
                bound = np.minimum(bound, intercepts[k]+.00012*(k-x*y/8192)**2)
            bounds = np.where(((x >= k) & (y >= k))[None, :, :], bound[None, :, :]+weighted_offset, -np.inf)
            flat = bounds.reshape(len(multipliers), -1)
            indices = np.argmax(flat, axis=1)
            values = flat[np.arange(len(multipliers)), indices]
            changed = values > maxima[:, k]
            maxima[changed, k] = values[changed]
            worst[changed, k, 0] = x.ravel()[indices[changed]]
            worst[changed, k, 1] = y.ravel()[indices[changed]]
        print('Weighted type-curve block', block, 'complete', flush=True)
    assert np.all(np.isfinite(maxima))
    curves = [dict(multiplier=float(b), actual_log_bounds=values.tolist(), worst_marginals=points.tolist(),
                   core_log_bounds=[float(max(values[k:k+4])) for k in range(898)])
              for b, values, points in zip(multipliers, maxima, worst)]
    base.write_new(output, dict(status='FIXED_MARGINAL_TOTAL_WEIGHTED_REGION_SCREEN_ONLY',
        core_total_weight=208800, mean_core_region_weight=815.625,
        bound_form='log R <= b*(x+y-1631.25)+log h(core_overlap), b<=0', curves=curves,
        outward_certified=False, full_second_moment_certified=False,
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in
                      [Path(__file__), proposals_path, atlas_path, Path(tight.__file__),
                       base.HERE/'screen_pair_type_bound.py', base.HERE/'generated/t128_s15_pair_spectrum.json']}))
    print('Weighted core log bound at740:', [(r['multiplier'], r['core_log_bounds'][740]) for r in curves], flush=True)


if __name__ == '__main__':
    run()
