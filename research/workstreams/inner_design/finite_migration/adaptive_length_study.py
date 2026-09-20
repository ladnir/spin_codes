"""Reproducible Figure 1 study: smaller steps, one transvection, fixed maps.

All numerical outputs are binary64 Q1 diagnostics. The BCH-256 baseline uses
the implemented maps; its alternatives are explicitly separate configurations.
"""
import argparse
from contextlib import contextmanager
import json
import math
from pathlib import Path

import numpy as np
import adaptive_step as core

HERE = Path(__file__).resolve().parent


def cases():
    return [(b,t,m) for b in (64,128,256)
            for m in (range(12,27,2) if b != 256 else range(14,25,2))
            for t in ((8,16,32,64) if b != 256 else (16,32,64,128))]


def record_for(b,t):
    if (b,t) == (256,128):
        engine = core.mixing.Engine(16,'1')
        return dict(t=128,s=19,spectrum=engine.spectrum,levels=engine.levels,
                    cancellation=list(map(list,zip(engine.injection_weights,engine.cancellation_weights))),
                    expansion_columns=engine.a_columns,feedback_columns=engine.columns,
                    transvection_rounds=1)
    return core.inner(t,min(19 if b == 256 else 20,core.capacity(t)))


def baseline_t(b):
    return 128 if b == 256 else 64


def sources():
    paths = {Path(core.__file__).resolve(),Path(__file__).resolve(),
             Path(core.nc.__file__).resolve(),Path(core.grid.__file__).resolve(),
             Path(core.mixing.__file__).resolve(),Path(core.grid.wm.__file__).resolve()}
    paths.update(core.grid.outer(b)[1] for b in (64,128))
    result = core.grid.ladder.candidate.sources()
    result.update({p.relative_to(core.grid.ROOT).as_posix():core.grid.maps.sha(p) for p in paths})
    return result


def run(output):
    if output.exists(): raise FileExistsError('Use a fresh output')
    cells, maps = [], {}
    for b,t,m in cases():
        record = record_for(b,t)
        maps[f'b{b}_t{t}'] = record
        for refresh in ([False,True] if t == baseline_t(b) else [False]):
            row = core.evaluate(record,b,m,refresh,sharp=True)
            cells.append(row)
            print(b,t,record['s'],m,'refresh' if refresh else 'IMT',round(row['q1_margin_bits'],6),flush=True)
    core.grid.ladder.model.base.write_new(output,dict(status='BINARY64_IMT_ADAPTIVE_LENGTH_STUDY',
        cells=cells,maps=maps,full_distance_proved=False,source_sha256=sources(),
        selection='largest Q1 margin among tested steps; ties choose larger t; not a runtime optimum'))


@contextmanager
def remove_lazy_returns():
    """Ablate terms of the positive envelope; result is NOT a valid Q1 bound."""
    original = core.epoch_logs
    def modified(record,tilts,refresh=False,sharp=False):
        assert sharp and not refresh
        zero,one = original(record,tilts,refresh,sharp)
        levels = sorted(record['spectrum'])
        lam = np.exp(tilts)
        t = record['t']
        # Only the uniform-refresh part of return-to-zero is retained.
        f1 = np.array([np.logaddexp(math.log(v/t)-lam*(v-1),
                                   math.log1p(-v/t)-lam*(v+1)) for v in levels]).T
        fresh = math.log(.5/((1 << record['s'])-1))
        one[:,1,0] = fresh+f1.max(axis=1)
        for i in range(len(levels)):
            one[:,i+2,0] = fresh+f1[:,i]
            one[:,len(levels)+2+i,0] = fresh+f1[:,i]
        return zero,one
    core.epoch_logs = modified
    try: yield
    finally: core.epoch_logs = original


def diagnose(output):
    if output.exists(): raise FileExistsError('Use a fresh output')
    cells = []
    for b,t,m in ((128,64,12),(128,64,16),(256,128,16),(256,128,20)):
        record = record_for(b,t)
        row = dict(b=b,t=t,s=record['s'],exponent=m)
        row['coarse'] = core.evaluate(record,b,m)['q1_margin_bits']
        row['sharp'] = core.evaluate(record,b,m,sharp=True)['q1_margin_bits']
        row['refresh'] = core.evaluate(record,b,m,refresh=True,sharp=True)['q1_margin_bits']
        with remove_lazy_returns():
            row['ablated_score_not_a_bound'] = core.evaluate(record,b,m,sharp=True)['q1_margin_bits']
        cells.append(row)
        print(row,flush=True)
    core.grid.ladder.model.base.write_new(output,dict(status='IMT_LENGTH_ROOT_CAUSE_DIAGNOSTIC',
        cells=cells,ablation_is_valid_bound=False,full_distance_proved=False,source_sha256=sources()))


def verify(source,output):
    if output.exists(): raise FileExistsError('Use a fresh output')
    saved = core.grid.ladder.model.base.read(source)
    core.grid.ladder.model.authenticate(saved)
    assert saved['status'] == 'BINARY64_IMT_ADAPTIVE_LENGTH_STUDY'
    assert not saved['full_distance_proved']
    expected = {(b,t,m,r) for b,t,m in cases()
                for r in ([False,True] if t==baseline_t(b) else [False])}
    assert {(r['b'],r['t'],r['exponent'],r['refresh']) for r in saved['cells']} == expected
    assert len(saved['cells']) == len(expected)
    seen, worst = set(),0.
    for row in saved['cells']:
        b,t,m = row['b'],row['t'],row['exponent']
        key = f'b{b}_t{t}'
        record = record_for(b,t)
        canonical = json.loads(json.dumps(record))
        assert saved['maps'][key] == canonical
        seen.add(key)
        # Recompute every cell using log-domain vector products. This differs
        # from the producer's normalized positive coefficient multiplication.
        actual = core.evaluate(record,b,m,row['refresh'],sharp=True,log_replay=True)
        for field,value in actual.items():
            if isinstance(value,float):
                difference = abs(value-row[field])
                assert math.isfinite(value) and difference < 1e-7,(b,t,m,field,difference)
                worst = max(worst,difference)
            else: assert value == row[field],(field,value,row[field])
        print('replayed',b,t,m,row['refresh'],flush=True)
    assert seen == set(saved['maps'])
    core.grid.ladder.model.base.write_new(output,dict(status='VERIFIED_BINARY64_IMT_ADAPTIVE_LENGTH_STUDY',
        producer_sha256=core.grid.maps.sha(source),cells_checked=len(expected),maps_checked=len(seen),
        log_domain_cells_checked=len(expected),largest_difference=worst,full_distance_proved=False,
        source_sha256=sources()))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=('run','diagnose','verify'))
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--source',type=Path)
    a = p.parse_args()
    if a.mode == 'verify': verify(a.source.resolve(),a.output.resolve())
    elif a.mode == 'diagnose': diagnose(a.output.resolve())
    else: run(a.output.resolve())
