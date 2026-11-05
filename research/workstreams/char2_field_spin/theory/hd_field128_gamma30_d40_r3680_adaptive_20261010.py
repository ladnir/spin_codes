"""Freeze positive r3680 marker proposals, then recompute every FLOAT cap."""

from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

import hd_field128_gamma_midpoint_adaptive as adaptive


def main():
    here = Path(__file__).resolve().parent
    source_path = here/'hd_field128_gamma30_d40_r3680_fresh_midpoint_float_20261010.json'
    raw = source_path.read_bytes()
    digest = 'dd9a570987c2a4a0ce80174cc6cc5ef26b44d2d976585efdc3e02656625f5da5'
    if sha256(raw).hexdigest() != digest:
        raise RuntimeError('Frozen geometric calibration changed')
    source = json.loads(raw)
    if (source['active_groups'] != 3680 or source['state_dimension'] != 30
            or source['target_distance'] != '2/5'
            or source['construction'] != adaptive.provider.gamma.construction(30)
            or source['profile_pairs_evaluated'] != 1841**2
            or not source['all_profiles_retained']):
        raise ValueError('Unexpected complete original-weight source')
    witness = source['witness']
    markers = [[str(F(2.**float(F(value)))) for value in row]
               for row in witness['marker_log_odds_exact']]
    if len(markers) != 76 or any(len(row) != 2 or any(F(v) <= 0 for v in row) for row in markers):
        raise ArithmeticError('Invalid positive marker conversion')
    pins = adaptive.source_pins()
    for name, value in source['source_sha256'].items():
        if name in pins and pins[name] != value:
            raise RuntimeError('Source closure conflict: '+name)
        pins[name] = value
    pins[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()

    def check():
        if source_path.read_bytes() != raw:
            raise RuntimeError('Calibration changed')
        for name, value in pins.items():
            if sha256((here/name).read_bytes()).hexdigest() != value:
                raise RuntimeError('Pinned source changed: '+name)
        for parent in source['source']['sources']:
            if sha256(Path(parent['path']).read_bytes()).hexdigest() != parent['sha256']:
                raise RuntimeError('Geometric-mark source changed')

    check()
    out = here/'hd_field128_gamma30_d40_r3680_adaptive_20261010'
    out.mkdir(exist_ok=False)
    proposal = dict(schema='gamma-original-weight-positive-literal-proposal-v1',
        construction=source['construction'], active_groups=3680, state_dimension=30,
        field_order=1 << 128, target_distance='2/5',
        source=dict(path=str(source_path), sha256=digest), source_sha256=dict(sorted(pins.items())),
        witness=dict(z=witness['z'], alpha=witness['alpha'], markers=markers, grid=128, state_dimension=30),
        marker_count=76, source_marker_logs_exact=witness['marker_log_odds_exact'],
        conversion='exp2 of each exact binary64 log proposal, frozen as an exact positive binary64 rational',
        source_role='Positive marker proposals only; every numerical cap is rebuilt',
        old_numerical_bounds_reused=False, no_weight_interpolation=True,
        directed=False, selected_weight_certificate=False, whole_code_certificate=False)
    proposal_pin = adaptive.provider.search.write_new(out/'proposal.json', proposal)
    result = adaptive.run(out/'proposal.json', 3680, max_seconds=90., max_fits=32,
                          max_calls=128, target_bits=4096)
    check()
    result_pin = adaptive.provider.search.write_new(out/'adaptive.json', result)
    last = result['last_complete_step']
    print(json.dumps(dict(proposal=proposal_pin, result=result_pin, status=result['status'],
        markers=result['marker_count'], winners=result['possible_winner_only_count'],
        score=last['log2_upper_screen'] if last else None, pending=result['pending_step'])), flush=True)


if __name__ == '__main__':
    main()
