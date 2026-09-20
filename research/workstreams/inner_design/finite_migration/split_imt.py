"""Directed seven-state IMT bounds with an explicit all-one-row count.

For Q active rows, h all-one rows put exactly h ones in every region.
Only d=Q-h ordinary rows use the nonconstant-band envelope. The frozen
positive fold is dimension-independent; the terminal kernel here is seven
states rather than the historical four-state RM2Sub kernel.
"""
import argparse
from fractions import Fraction as F
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import gammaln
from flint import arb, ctx

import ladder
import constant_split_grid as fold
import scalable_split_grid as selected
import sparse_ranges
import sparse_bch

model = ladder.model
up = fold.up
UNSET = selected.UNSET
POLICY = 'imt-seven-state;all-one-count-exact;ordinary-nonconstant-bands-v1'


def terminal_256(matrices, exponents):
    """Directed fixed-width powers; never use a BLAS reduction for a bound."""
    assert matrices.ndim == 3 and matrices.shape[1:] == (7, 7)
    assert exponents.shape == (len(matrices),)
    assert np.isfinite(matrices).all() and (matrices >= 0).all()
    assert (matrices <= 2).all() and (np.abs(exponents) < 10**9).all()
    a, powers = matrices.copy(), exponents.copy()
    for _ in range(8):
        out = up(a[:, :, 0, None] * a[:, None, 0, :])
        out = up(out + up(a[:, :, 1, None] * a[:, None, 1, :]))
        out = up(out + up(a[:, :, 2, None] * a[:, None, 2, :]))
        out = up(out + up(a[:, :, 3, None] * a[:, None, 3, :]))
        out = up(out + up(a[:, :, 4, None] * a[:, None, 4, :]))
        out = up(out + up(a[:, :, 5, None] * a[:, None, 5, :]))
        out = up(out + up(a[:, :, 6, None] * a[:, None, 6, :]))
        maximum = out.max(axis=(1, 2))
        assert (maximum > 0).all() and np.isfinite(maximum).all()
        _, shift = np.frexp(maximum)
        a = up(np.ldexp(out, -shift[:, None, None]))
        powers = 2 * powers + shift
    total = up(a[:, 0, 0] + a[:, 0, 1])
    total = up(total + a[:, 0, 2])
    total = up(total + a[:, 0, 3])
    total = up(total + a[:, 0, 4])
    total = up(total + a[:, 0, 5])
    total = up(total + a[:, 0, 6])
    return total, powers


def choose(engine, region, ordinary):
    """Binary64 witness discovery only; stable even above 1024 rows."""
    assert ordinary > 0 and len(region) == ordinary + 1
    logs = np.array([[float(v.log()) if v > 0 else -math.inf for v in row]
                     for row in region]).reshape(-1, 7, 7)
    floating = np.exp(logs - float(logs.max()))
    counts = model.outer.inputs.caps_module.caps()
    j = np.arange(ordinary + 1)
    log_counts = gammaln(ordinary + 1) - gammaln(j + 1) - gammaln(ordinary - j + 1)
    result = []
    for band in model.BANDS[:-1]:
        weights = np.array(band)
        costs = np.array([math.log(counts[w]) - math.log(math.comb(256, w)) for w in band])

        def objective(theta):
            p = 1 / (1 + math.exp(-theta))
            lp, lq = math.log(p), math.log1p(-p)
            density = float(np.max(costs - weights * lp - (256 - weights) * lq))
            probabilities = np.exp(log_counts + j * lp + (ordinary - j) * lq)
            matrix = np.einsum('i,ijk->jk', probabilities, floating)
            return sparse_bch.terminal_log(matrix) + ordinary * density

        fit = minimize_scalar(objective, bounds=(-8, 12), method='bounded')
        assert fit.success and math.isfinite(fit.fun)
        result.append(F.from_float(1 / (1 + math.exp(-float(fit.x)))))
    return result


def evaluate(engine, region, ps, lo, hi, tilt, best, owners, witness, replay=False):
    assert engine.n == 7 and len(region) == hi + 1
    assert best.shape == owners.shape == (hi - lo + 1, hi + 1)
    assert 1 <= lo <= hi <= engine.length
    assert len(ps) == 12 and all(0 < p < 1 for p in ps)
    assert model.BANDS[-1] == (256,) and len(model.BANDS) == 13
    costs = engine.costs(list(ps) + [F(1)])[:-1]
    roots = [model.up((c.log() / 256).exp()) for c in costs]
    probabilities = list(map(model.number, ps))
    left = up(np.array([float(model.up(r * (1 - p))) for r, p in zip(roots, probabilities)]))
    right = up(np.array([float(model.up(r * p)) for r, p in zip(roots, probabilities)]))
    keep = fold.sparse.hull.indices(left, right)
    correction = model.up((engine.cutoff * model.number(tilt).exp()).exp())
    man, exp = correction.man_exp()
    correction_exp = int(exp) + int(man).bit_length()
    correction_man = float(up(float(correction * arb(2)**(-correction_exp))))
    initial = sparse_ranges.initial(region, 7)
    locations = [math.comb(engine.length, q) for q in range(hi + 1)]
    ordinary_count, visited = 1, 0
    for d, matrices, powers in fold.folds(*initial, left[keep], right[keep]):
        qs = np.arange(max(lo, d), hi + 1, dtype=np.int64)
        hs = qs - d
        mask = owners[qs - lo, hs] == witness if replay else best[qs - lo, hs] > -80
        qs, hs = qs[mask], hs[mask]
        if len(qs):
            terminal, terminal_exp = terminal_256(matrices[hs], powers[hs])
            terminal = up(terminal * correction_man)
            previous_q, binomial = -1, 0
            for q0, h0, value, exponent in zip(qs, hs, terminal, terminal_exp):
                q, h = int(q0), int(h0)
                binomial = binomial * q // (q - d) if q == previous_q + 1 and q > d else math.comb(q, d)
                count = locations[q] * binomial * ordinary_count
                power = fold.dyadic_ceiling(value, int(exponent) + correction_exp, count)
                if replay:
                    assert power <= int(best[q - lo, h]), (q, h, power, int(best[q - lo, h]))
                elif power < best[q - lo, h]:
                    best[q - lo, h], owners[q - lo, h] = power, witness
                previous_q = q
            visited += len(qs)
        ordinary_count *= 12
    return visited


def margins(best, lo, hi):
    powers = selected.pack(best, lo, hi)
    values = fold.bounds(powers, lo, hi)
    passing = [q for q, value in zip(range(lo, hi + 1), values) if value < F(2)**-50]
    return values, passing


def run(output, exponent, lo, hi, jobs, verify=False):
    saved = model.base.read(output) if verify else None
    if saved:
        model.authenticate(saved)
        assert saved['status'] == 'IMT_SPLIT_GRID' and saved['policy'] == POLICY
        exponent, (lo, hi) = saved['exponent'], saved['interval']
    elif output.exists():
        raise FileExistsError('Use a fresh output path')
    ctx.prec = 512 if verify else 256
    engine = ladder.Engine(exponent)
    assert 1 <= lo <= hi <= engine.length and (hi - lo + 1) * (hi + 1) <= 10**7
    if saved:
        assert saved['instance'] == engine.identity()
        best = selected.unpack(saved['powers'], lo, hi)
        owners = np.array(saved['owners'], dtype=np.int64)
        assert owners.shape == best.shape
        jobs = saved['witnesses']
        assert all(np.all((owners[q - lo, :q + 1] >= 0) & (owners[q - lo, :q + 1] < len(jobs)))
                   for q in range(lo, hi + 1))
    else:
        best, owners = selected.table(lo, hi), np.full((hi - lo + 1, hi + 1), -1, dtype=np.int64)
    witnesses, visited = [], 0
    started = time.monotonic()
    for index, job in enumerate(jobs):
        tilt = F(job['tilt'])
        q, h = job['anchor']
        assert -12 <= tilt <= 3 and 0 <= h < q <= hi
        region = engine.region(tilt, hi)
        ps = list(map(model.base.decode, job['probabilities'])) if saved else choose(engine, region[h:q + 1], q - h)
        visited += evaluate(engine, region, ps, lo, hi, tilt, best, owners, index, verify)
        witnesses.append(dict(tilt=str(tilt), anchor=[q, h], probabilities=list(map(model.base.encode, ps))))
        values, passing = margins(best, lo, hi)
        print('split', ctx.prec, exponent, index + 1, str(tilt), q, h,
              'passing', len(passing), '/', hi - lo + 1,
              'worst-case', int(best[0, :lo + 1].max()) if lo == hi else '',
              'seconds', round(time.monotonic() - started, 2), flush=True)
        engine.region.cache_clear()
        engine.epoch.cache_clear()
    values, passing = margins(best, lo, hi)
    if saved:
        assert visited == sum(q + 1 for q in range(lo, hi + 1))
        assert passing == saved['passing_50_bits']
        model.base.write_new(output.with_name(output.stem + '_replay.json'), dict(
            status='IMT_SPLIT_GRID_512_BIT_REPLAY_PASSED', producer_sha256=model.base.sha(output),
            cases_checked=visited, instance=engine.identity(), interval=[lo, hi], full_distance_proved=False))
    else:
        model.base.write_new(output, dict(status='IMT_SPLIT_GRID', policy=POLICY,
            exponent=exponent, instance=engine.identity(), interval=[lo, hi],
            witnesses=witnesses, powers=selected.pack(best, lo, hi), owners=owners.tolist(),
            passing_50_bits=passing, full_distance_proved=False, source_sha256=ladder.candidate.sources()))
    if lo == hi:
        value = values[0]
        print('occupancy margin', math.log2(value.denominator) - math.log2(value.numerator), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--m', type=int, choices=(16, 18), default=16)
    p.add_argument('--first', type=int, default=218)
    p.add_argument('--last', type=int, default=218)
    p.add_argument('--witness', action='append', help='log-lambda-index:Q:h; index is divided by 10')
    p.add_argument('--verify', action='store_true')
    a = p.parse_args()
    if not a.verify and not a.witness:
        p.error('Provide witness proposals')
    jobs = []
    for text in a.witness or []:
        tilt, q, h = text.split(':')
        jobs.append(dict(tilt=str(F(tilt) / 10), anchor=[int(q), int(h)]))
    run(a.output.resolve(), a.m, a.first, a.last, jobs, a.verify)
