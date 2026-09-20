"""Check the complete IMT Q1 grid and repeat selected cells in log arithmetic."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
import parameter_q1 as grid


def verify(path, output):
    model = grid.ladder.model
    saved = model.base.read(path)
    model.authenticate(saved)
    assert saved['status'] == 'BINARY64_IMT_PARAMETER_Q1_ONLY'
    assert saved['geometry_complete'] and saved['full_distance_proved'] is False
    assert saved['expansion_chain_seed'] == grid.SEED_A and saved['feedback_chain_seed'] == grid.SEED_B
    expected = grid.geometries()
    assert [(r['b'],r['t'],r['s'],r['exponent']) for r in saved['rows']] == expected
    tilts = np.array(saved['tilts'])
    np.testing.assert_array_equal(tilts,np.arange(-180,1,dtype=float)/10)
    # Independent log-domain products cover both outers, all epoch widths,
    # shortest and largest message lengths, and small and large state sizes.
    checks = {(64,64,7,20),(128,64,20,12),(128,64,20,26),
              (128,128,8,20),(64,128,20,20),(128,256,9,20),(128,256,20,20)}
    seen_maps = set()
    worst_difference = 0.
    for row in saved['rows']:
        b,t,s,m = key = row['b'],row['t'],row['s'],row['exponent']
        name = f't{t}_s{s}'
        record = grid.inner(t,s)
        canonical = json.loads(json.dumps(record))
        assert saved['maps'][name] == canonical
        seen_maps.add(name)
        actual = grid.screen(record,b,m,tilts)
        assert row['full_margin_bits'] is None and row['full_distance_proved'] is False
        for field,value in actual.items():
            if type(value) is float:
                assert math.isfinite(value) and abs(value-row[field]) < 1e-8
            else:
                assert row[field] == value
        if key in checks:
            lam = np.exp(tilts)
            rz,ra = grid.wm.regions(*grid.wm.transfers(record,lam,1),row['outer_rows']//t)
            moments = grid.wm.coefficients(rz,ra-math.log(row['outer_rows']//t),b)
            best = np.minimum(0.,moments+row['bad_weight']*lam[:,None]).min(axis=0)
            counts,_ = grid.outer(b)
            margin = -float(np.logaddexp.reduce([math.log(row['outer_rows']*n)+best[w]
                                                for w,n in counts.items()]))/math.log(2)
            difference = abs(margin-row['q1_margin_bits'])
            assert difference < 1e-7
            worst_difference = max(worst_difference,difference)
            print('log-domain check',key,difference,flush=True)
    assert seen_maps == set(saved['maps'])
    model.base.write_new(output,dict(status='VERIFIED_BINARY64_IMT_Q1_GRID',
        producer_sha256=model.base.sha(path), geometries_checked=len(expected),
        exact_maps_reconstructed=len(seen_maps),log_domain_cells_checked=len(checks),
        largest_log_domain_difference_bits=worst_difference,full_distance_proved=False,
        source_sha256=grid.ladder.candidate.sources()))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a = p.parse_args()
    verify(a.input.resolve(),a.output.resolve())
