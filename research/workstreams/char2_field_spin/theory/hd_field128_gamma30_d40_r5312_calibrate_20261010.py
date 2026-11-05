"""Bounded original-weight r5312 FLOAT fitting from pinned positive literals."""

from hashlib import sha256
import json
from pathlib import Path

import hd_field128_gamma_midpoint_dense_adaptive as adaptive


HERE = Path(__file__).resolve().parent
SOURCE = HERE/'hd_field128_gamma30_d40_r5888_winner_replay_20261010/proposal.json'
SOURCE_SHA256 = '43cea88db107cfc364ab5989978269fa97b17d0068e32ca5d24aa930887e4754'
OUTPUT = HERE/'hd_field128_gamma30_d40_r5312_fresh_adaptive_20261010'


def main():
    raw = SOURCE.read_bytes()
    if sha256(raw).hexdigest() != SOURCE_SHA256:
        raise RuntimeError('Frozen literal parent changed')
    if sha256(Path(adaptive.__file__).read_bytes()).hexdigest() != 'f96e9cf781b40828a6aa83880879f904909a6411fe02003cfc0a81d3c9331e61':
        raise RuntimeError('Reviewed adaptive helper changed')
    source = json.loads(raw)
    if (source['active_groups'] != 5888 or source['state_dimension'] != 30
            or source['target_distance'] != '2/5'
            or source['construction'] != adaptive.provider.gamma.construction(30)
            or len(source['witness']['markers']) != 17):
        raise ValueError('Unexpected exact literal source scope')
    pins = adaptive.source_pins()
    for name, digest in source['source_sha256'].items():
        if name in pins and pins[name] != digest:
            raise RuntimeError('Source closure conflict: '+name)
        pins[name] = digest
    pins[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()

    def check():
        if SOURCE.read_bytes() != raw:
            raise RuntimeError('Literal parent changed during calibration')
        for name, digest in pins.items():
            if sha256((HERE/name).read_bytes()).hexdigest() != digest:
                raise RuntimeError('Pinned helper changed: '+name)

    check()
    OUTPUT.mkdir(exist_ok=False)
    proposal = dict(schema='gamma-original-weight-positive-literal-proposal-v1',
        construction=source['construction'], field_order=1 << 128,
        state_dimension=30, target_distance='2/5', active_groups=5312,
        source=dict(path=str(SOURCE), sha256=SOURCE_SHA256),
        source_sha256=dict(sorted(pins.items())), source_calibration_active_groups=5888,
        new_requested_active_groups=5312, witness=source['witness'], marker_count=17,
        source_role='Exact positive literals only; no source numerical bounds transferred',
        old_numerical_bounds_reused=False, no_weight_interpolation=True,
        directed=False, selected_weight_certificate=False, whole_code_certificate=False)
    proposal_pin = adaptive.provider.search.write_new(OUTPUT/'initial_proposal.json', proposal)
    result = adaptive.run(OUTPUT/'initial_proposal.json', 5312, max_seconds=90.,
                          max_fits=32, max_calls=128, target_bits=80)
    completed = result['last_complete_step']
    if completed and completed['profile_pairs_evaluated'] != 1441**2:
        raise ArithmeticError('Incomplete original r5312 profile domain')
    check()
    result_pin = adaptive.provider.search.write_new(OUTPUT/'adaptive.json', result)
    print(json.dumps(dict(proposal=proposal_pin, result=result_pin,
        status=result['status'], complete_score=completed['log2_upper_screen'] if completed else None,
        completed_markers=result['marker_count'], pending=result['pending_step'],
        winner_marker_count=result['possible_winner_only_count'],
        elapsed_seconds=result['elapsed_seconds'])), flush=True)
    return result, result_pin


if __name__ == '__main__':
    main()
