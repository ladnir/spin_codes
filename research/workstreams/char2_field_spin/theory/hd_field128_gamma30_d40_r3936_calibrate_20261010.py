"""Fresh original-weight FLOAT calibration from exact positive r3680 literals."""

from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path

import hd_field128_gamma_midpoint_adaptive as adaptive


def main():
    here = Path(__file__).resolve().parent
    parent = here/'hd_field128_gamma30_d40_r3680_winners_replay_20261010/proposal.json'
    digest = '23dd34aecb5b6b8e4562be03c89fd22904eeb582398b9c2e62283939c12e0141'
    raw = parent.read_bytes()
    if sha256(raw).hexdigest() != digest:
        raise RuntimeError('Frozen positive-literal parent changed')
    source = json.loads(raw)
    if (source['active_groups'] != 3680 or source['state_dimension'] != 30
            or source['target_distance'] != '2/5'
            or source['construction'] != adaptive.provider.gamma.construction(30)
            or len(source['witness']['markers']) != 60):
        raise ValueError('Unexpected positive-literal parent')
    pins = adaptive.source_pins()
    for name, value in source['source_sha256'].items():
        if name in pins and pins[name] != value:
            raise RuntimeError('Source closure conflict: '+name)
        pins[name] = value
    pins[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()

    def check():
        if parent.read_bytes() != raw:
            raise RuntimeError('Parent changed')
        for name, value in pins.items():
            if sha256((here/name).read_bytes()).hexdigest() != value:
                raise RuntimeError('Pinned source changed: '+name)

    check()
    out = here/'hd_field128_gamma30_d40_r3936_calibration_20261010'
    out.mkdir(exist_ok=False)
    proposal = dict(schema='gamma-original-weight-positive-literal-proposal-v1',
        construction=source['construction'], active_groups=3936, state_dimension=30,
        field_order=1 << 128, target_distance='2/5',
        source=dict(path=str(parent), sha256=digest), source_sha256=dict(sorted(pins.items())),
        witness=deepcopy(source['witness']), marker_count=60,
        source_role='Exact positive literals only; every numerical cap is rebuilt at weight3936',
        old_numerical_bounds_reused=False, no_weight_interpolation=True,
        directed=False, selected_weight_certificate=False, whole_code_certificate=False)
    proposal_pin = adaptive.provider.search.write_new(out/'proposal.json', proposal)
    result = adaptive.run(out/'proposal.json', 3936, max_seconds=90., max_fits=17,
                          max_calls=128, target_bits=4096)
    check()
    result_pin = adaptive.provider.search.write_new(out/'adaptive.json', result)
    last = result['last_complete_step']
    print(json.dumps(dict(proposal=proposal_pin, result=result_pin, status=result['status'],
        markers=result['marker_count'], winners=result['possible_winner_only_count'],
        score=last['log2_upper_screen'] if last else None, pending=result['pending_step'])), flush=True)


if __name__ == '__main__':
    main()
