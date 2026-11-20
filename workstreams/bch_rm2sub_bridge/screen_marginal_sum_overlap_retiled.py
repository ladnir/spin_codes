"""Retain the marginal sum in the full region-ratio screen; no certificate."""
import math
from pathlib import Path
import numpy as np
from flint import arb_poly, ctx
from scipy.special import gammaln
import bridge as base
import tightened_occupancy as tight
from screen_pair_type_bound import load_pairs, evaluate


def run():
    output = base.HERE/'generated/marginal_sum_overlap_retiled_screen.json'
    assert not output.exists()
    proposal_path = base.HERE/'generated/global_overlap_retiled_screen.json'
    proposals = base.read(proposal_path)
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
    table = np.full((901, 161), -np.inf)
    factorial = gammaln(np.arange(8193)+1)
    powers, counts = load_pairs()
    for proposal in proposals['tilt_proposals']:
        xs, ys = proposal['first_weights'], proposal['second_weights']
        x, y = np.meshgrid(xs, ys, indexing='ij')
        index = ((x+y-1480)//2).ravel()
        single = np.array([[singles[a]+singles[b] for b in ys] for a in xs])
        trivial = np.minimum(np.array([singles[a] for a in xs])[:, None],
                             np.array([singles[b] for b in ys])[None, :])-single
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
            feasible = (x >= k) & (y >= k)
            np.maximum.at(table[k], index, np.where(feasible, bound, -np.inf).ravel())
    core = np.maximum.reduce([table[d:d+898] for d in range(4)])
    assert np.all(np.isfinite(np.max(core, axis=1)))
    rows = [[float(v) if np.isfinite(v) else None for v in row] for row in core]
    base.write_new(output, dict(status='MARGINAL_SUM_AND_CORE_OVERLAP_REGION_SCREEN_ONLY',
        marginal_sums=list(range(1480, 1801, 2)), core_overlaps=list(range(898)), log_bounds=rows,
        core_mean_marginal_sum=1631.25, correction_robust=True,
        outward_certified=False, full_second_moment_certified=False,
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in
                      [Path(__file__), proposal_path, atlas_path, Path(tight.__file__),
                       base.HERE/'screen_pair_type_bound.py', base.HERE/'generated/t128_s15_pair_spectrum.json']}))
    print('Retained marginal sums for all898 core overlaps:', int(np.isfinite(core).sum()), 'cells', flush=True)


if __name__ == '__main__':
    run()
