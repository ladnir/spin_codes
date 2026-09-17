"""Matched-map Q1 study of transvection rounds and message length.

Screening is binary64 only. Certification reconstructs selected witnesses in
Arb; replay uses 512 bits and linear region products. No higher-occupancy
certificate transfers between mixer laws. 'refresh' means the ideal marginal
limit, not zero transvections and not a certified full-distance instance.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import math
from pathlib import Path

import numpy as np
from flint import ctx
import q1_slack as audit

ladder, model = audit.ladder, audit.model
HERE = Path(__file__).resolve().parent


def epsilon(rounds):
    if rounds == 'refresh':
        return F(0)
    r = int(rounds)
    if r < 1:
        raise ValueError('Rounds must be positive; refresh is a separate limit')
    return F(1, 2**r)


def transfers(spectrum, injection_weights, cancellation_weights, z, eps):
    """Twelve-state positive transfer, also usable with exact rational inputs."""
    levels = sorted(spectrum)
    t, m = len(injection_weights), sum(spectrum.values())
    old_n, n = len(levels)+2, 2*len(levels)+2
    eta = 1-eps
    assert 0 <= eps <= 1 and len(cancellation_weights) == t
    counts = Counter(injection_weights)
    zero, one = ([z*0 for _ in range(n*n)] for _ in range(2))
    zero[0] = z*0+1
    f0 = {v: z**v for v in levels}
    f1 = {v: (v*z**(v-1)+(t-v)*z**(v+1))/t for v in levels}
    d0, d1 = max(f0.values()), max(f1.values())
    zero[n+1] = eps*d0
    one[n+1] = eps*d1
    one[n] = eps*max(z**w for w in cancellation_weights)/t + eta*d1/m
    for i, v in enumerate(levels):
        c, h = i+2, old_n+i
        one[h] = z*counts[v]/t
        zero[c*n+c] = eps*f0[v]
        zero[h*n+h] = eps*f0[v]
        one[c*n+1] = eps*f1[v]
        one[h*n+1] = eps*f1[v]
        cancels = [z**w for iw, w in zip(injection_weights, cancellation_weights) if iw == v]
        one[c*n] = eps*sum(cancels, z*0)/(t*spectrum[v]) + eta*f1[v]/m
        one[h*n] = eps*max(cancels, default=z*0)/t + eta*f1[v]/m
        for j, w in enumerate(levels):
            dest = j+2
            zero[n+dest] = eta*d0*spectrum[w]/m
            one[n+dest] = eta*d1*spectrum[w]/m
            for source in (c, h):
                zero[source*n+dest] += eta*f0[v]*spectrum[w]/m
                one[source*n+dest] += eta*f1[v]*spectrum[w]/m
    return tuple(zero), tuple(one), n


class Engine(ladder.Engine):
    def __init__(self, exponent, rounds):
        super().__init__(exponent)
        self.rounds = str(rounds)
        self.eps = epsilon(self.rounds)
        self.n = 2*len(self.levels)+2
        assert len(set(self.columns)) == len(self.columns) and all(self.columns)
        image = lambda q: sum(((q & a).bit_count() & 1) << i for i, a in enumerate(self.a_columns))
        images = [image(q) for q in self.columns]
        self.injection_weights = [v.bit_count() for v in images]
        self.cancellation_weights = [(v ^ (1 << i)).bit_count() for i, v in enumerate(images)]

    def identity(self):
        result = super().identity()
        result['inner']['transvection_rounds'] = None if self.rounds == 'refresh' else int(self.rounds)
        result['inner']['mixing_law'] = 'ideal_uniform_nonzero_marginal' if self.rounds == 'refresh' else 'independent_transvection_product'
        result['setup'] = ('independent row/region permutations; '+
                           ('ideal uniform nonzero refresh marginal' if self.rounds == 'refresh' else self.rounds+' independent transvections per epoch')+
                           '; zero start; output before update; persistent state; no flush')
        return result

    @lru_cache(maxsize=128)
    def epoch(self, log_lam):
        z = (-model.number(log_lam).exp()).exp()
        zero, one, n = transfers(self.spectrum, self.injection_weights,
                                 self.cancellation_weights, z, model.number(self.eps))
        assert n == self.n
        return tuple(map(model.up, zero)), tuple(map(model.up, one))

    def screen_coefficients(self, log_lam):
        lam = math.exp(log_lam)
        zero, one, n = transfers(self.spectrum, self.injection_weights,
                                 self.cancellation_weights, math.exp(-lam), float(self.eps))
        zero, one = np.array(zero).reshape(n,n), np.array(one).reshape(n,n)
        rz, ra = np.eye(n), np.zeros((n,n))
        remaining = self.length//128
        while remaining:
            if remaining & 1:
                ra, rz = ra@zero+rz@one, rz@zero
            remaining >>= 1
            if remaining:
                one, zero = one@zero+zero@one, zero@zero
        ra /= self.length//128
        current = np.zeros((257,n))
        current[0,0] = 1
        # Normalize binomial coefficients at each step; all summands positive.
        for count in range(1,257):
            nxt = np.zeros_like(current)
            weights = np.arange(count)
            nxt[:count] = (current[:count]@rz)*((count-weights)/count)[:,None]
            nxt[1:count+1] += (current[:count]@ra)*((weights+1)/count)[:,None]
            current = nxt
        logs = np.log(np.maximum(current.sum(axis=1), 1e-300)) + math.log(self.length)+self.cutoff*lam
        return {w: float(logs[w]) for w in model.base.WEIGHTS}


def screen(output):
    if output.exists():
        raise FileExistsError('Use a fresh output path')
    cells = []
    for m in (16,18,20):
        for rounds in ('1','2','3','4','8','refresh'):
            engine = Engine(m, rounds)
            best = {w: math.inf for w in model.base.WEIGHTS}
            witnesses = {}
            for index in range(-20,21):
                tilt = math.log(3/engine.length)+index/16
                values = engine.screen_coefficients(tilt)
                for w in best:
                    if values[w] < best[w]:
                        best[w], witnesses[w] = values[w], str(F.from_float(tilt))
            values = {w: F.from_float(math.exp(max(-690., min(0., best[w])))) for w in best}
            upper, _, rest = model.base.bch_bound(values)
            # Low shells carry the objective; additional original-grid tilts
            # cover the small remaining-shell tail during outward evaluation.
            selected = sorted({witnesses[w] for w in (38,40,42,44,46,48)}, key=lambda x: float(F(x)))
            row = dict(exponent=m, rounds=rounds, instance=engine.identity(),
                       diagnostic_q1_bits=audit.bits(upper), selected_tilts=selected,
                       remaining_shell_fraction=float(rest/upper))
            cells.append(row)
            print(m, rounds, row['diagnostic_q1_bits'], 'witnesses',len(selected),flush=True)
    model.base.write_new(output, dict(status='BINARY64_IMT_MIXING_ROUNDS_Q1_SCREEN', cells=cells,
        full_distance_proved=False, source_sha256=ladder.candidate.sources()))


def certify(output, source, cases, verify=False):
    ctx.prec = 512 if verify else 256
    if not verify and output.exists():
        raise FileExistsError('Use a fresh output path')
    saved = model.base.read(output) if verify else None
    if saved:
        model.authenticate(saved)
        rows = saved['cells']
    else:
        proposal = model.base.read(source)
        model.authenticate(proposal)
        rows = [r for r in proposal['cells'] if not cases or f"{r['exponent']}:{r['rounds']}" in cases]
        if cases:
            assert len(rows) == len(cases)
    results = []
    for row in rows:
        engine = Engine(row['exponent'],row['rounds'])
        assert row['instance'] == engine.identity()
        tilts = row['tilts'] if saved else row['selected_tilts']
        best = None
        for tilt in tilts:
            current = audit.coefficients(engine,F(tilt),linear=verify)
            best = current if best is None else audit.minimum(best,current)
        if saved:
            accepted = {int(w): model.base.decode(v) for w,v in row['coefficients'].items()}
            assert set(accepted) == set(best) and all(best[w] <= accepted[w] for w in best)
            best = accepted
        upper, _, rest = model.base.bch_bound(best)
        result = dict(exponent=engine.exponent,rounds=engine.rounds,instance=engine.identity(),
                      tilts=tilts,coefficients={str(w):model.base.encode(v) for w,v in best.items()},
                      q1_upper=model.base.encode(upper),q1_margin_bits=audit.bits(upper),
                      remaining_shell_fraction=float(rest/upper))
        if saved:
            assert result == row
        results.append(result)
        print(ctx.prec,engine.exponent,engine.rounds,result['q1_margin_bits'],flush=True)
    if saved:
        model.base.write_new(output.with_name(output.stem+'_replay.json'),dict(
            status='IMT_MIXING_ROUNDS_Q1_512_BIT_LINEAR_REPLAY_PASSED',
            producer_sha256=model.base.sha(output),cells_checked=len(results),full_distance_proved=False))
    else:
        model.base.write_new(output,dict(status='OUTWARD_IMT_MIXING_ROUNDS_Q1',cells=results,
            full_distance_proved=False,source_sha256={**ladder.candidate.sources(),
                source.relative_to(model.ROOT).as_posix():model.base.sha(source)}))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=('screen','certify'))
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--source',type=Path)
    p.add_argument('--cases',nargs='*')
    p.add_argument('--verify',action='store_true')
    args = p.parse_args()
    if args.mode == 'screen':
        screen(args.output.resolve())
    else:
        if not args.verify and not args.source:
            p.error('certify requires --source unless replaying')
        certify(args.output.resolve(),args.source.resolve() if args.source else None,args.cases,args.verify)
