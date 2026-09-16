"""Outward point checks of the fixed high-band cover backend."""
import argparse
from fractions import Fraction as F
from pathlib import Path

from flint import arb,ctx
import high_band_dense as dense

model = dense.model


def run(seed,output,verify=False):
    saved = model.base.read(output) if verify else None
    if saved:
        model.authenticate(saved)
        assert saved['status'] == 'IMT_FIXED_HIGH_BAND_POINTS'
    else:
        assert not output.exists()
    ctx.prec = 512 if verify else 256
    prior = model.base.read(seed)
    model.authenticate(prior)
    checker = dense.Checker(16)
    assert prior['instance'] == checker.engine.identity()
    rows = []
    for i,point in enumerate(prior['points']):
        q,v,w = point['q'],F(point['coordinate']),point['witness']
        upper = checker.split_bound(q,q,v,v,w)
        if saved:
            old = saved['points'][i]
            assert (old['q'],old['coordinate'],old['witness']) == (q,str(v),w)
            assert upper <= model.number(model.base.decode(old['log_upper']))
        rows.append(dict(q=q,coordinate=str(v),witness=w,margin_bits=float(-upper/arb(2).log()),
                         log_upper=model.base.encode(model.exact(upper))))
        print('fixed high-band',ctx.prec,q,str(v),rows[-1]['margin_bits'],flush=True)
    if saved:
        assert len(saved['points']) == len(rows)
        model.base.write_new(output.with_name(output.stem+'_replay.json'),dict(
            status='IMT_FIXED_HIGH_BAND_POINTS_512_REPLAY_PASSED',
            producer_sha256=model.base.sha(output),full_distance_proved=False))
    else:
        sources = dense.parent.parent.short_dense.mixed_dense.ladder.candidate.sources()
        sources[seed.relative_to(model.ROOT).as_posix()] = model.base.sha(seed)
        model.base.write_new(output,dict(status='IMT_FIXED_HIGH_BAND_POINTS',instance=checker.engine.identity(),
            points=rows,cutoff=str(dense.CUTOFF),xi=str(dense.XI),full_distance_proved=False,source_sha256=sources))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--seed',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--verify',action='store_true')
    a = p.parse_args()
    run(a.seed.resolve(),a.output.resolve(),a.verify)
