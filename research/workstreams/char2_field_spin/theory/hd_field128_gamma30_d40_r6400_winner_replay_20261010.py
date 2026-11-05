"""Exact literal selection and fresh complete original-r6400 FLOAT replay."""

import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'

from hashlib import sha256
import json
from pathlib import Path

import hd_field128_gamma_midpoint_dense_adaptive as adaptive


def main():
    here, runner = Path(__file__).resolve().parent, Path(__file__).resolve()
    source_path = here/'hd_field128_gamma30_d40_r6400_strengthened_float_20261010.json'
    raw = source_path.read_bytes()
    digest = sha256(raw).hexdigest()
    if digest != '07150979e76724b4d7e77530fc60f65d04fdd37c75a0b5c688440722c0cf1f10':
        raise RuntimeError('Frozen complete adaptive source changed')
    source = json.loads(raw)
    completed = source['last_complete_step']
    indices = source['winning_marker_indices_at_last_complete']
    markers = source['witness']['markers']
    if (source['active_groups'] != 6400 or source['state_dimension'] != 30
            or source['target_distance'] != '2/5' or source['pending_step'] is not None
            or source['construction'] != adaptive.provider.gamma.construction(30)
            or completed['profile_pairs_evaluated'] != 897**2
            or not completed['complete_selected_weight_screen']
            or completed['marker_count'] != len(markers)
            or not indices or indices != sorted(set(indices))
            or any(type(i) is not int or not 0 <= i < len(markers) for i in indices)):
        raise ValueError('Unexpected completed original-r6400 source scope')
    pins = adaptive.source_pins()
    for name, value in source['source_sha256'].items():
        if name in pins and pins[name] != value:
            raise RuntimeError('Source closure conflict: '+name)
        pins[name] = value
    pins[runner.name] = sha256(runner.read_bytes()).hexdigest()
    def check():
        if source_path.read_bytes() != raw:
            raise RuntimeError('Completed source changed')
        for name, value in pins.items():
            if sha256((here/name).read_bytes()).hexdigest() != value:
                raise RuntimeError('Source closure changed: '+name)
    check()
    output = here/'hd_field128_gamma30_d40_r6400_winner_replay_20261010'
    output.mkdir(exist_ok=False)
    witness = source['witness']
    proposal = dict(schema='gamma-dense-selected-positive-literal-proposal-v1',
        construction=source['construction'], state_dimension=30, target_distance='2/5', active_groups=6400,
        source=dict(path=str(source_path), sha256=digest), source_sha256=dict(sorted(pins.items())),
        selected_source_marker_indices=indices, source_marker_count=len(markers), marker_count=len(indices),
        witness=dict(z=witness['z'], alpha=witness['alpha'], markers=[markers[i] for i in indices],
                     grid=witness['grid'], state_dimension=30),
        source_role='Exact positive literal selection only; no numerical cap or score inherited',
        selection_role='Fresh complete first-winner selection, pending independent subset replay',
        old_numerical_bounds_reused=False, directed=False, selected_weight_certificate=False,
        whole_code_certificate=False, no_weight_interpolation=True)
    proposal_pin = adaptive.provider.search.write_new(output/'proposal.json', proposal)
    replay = adaptive.run(output/'proposal.json', 6400, max_seconds=90., max_fits=0,
                          max_calls=128, target_bits=4096)
    fresh = replay['last_complete_step']
    matches = (replay['pending_step'] is None and fresh is not None
               and fresh['profile_pairs_evaluated'] == 897**2
               and fresh['regional_logs_sha256'] == completed['regional_logs_sha256']
               and fresh['log2_upper_screen'] == completed['log2_upper_screen'])
    replay['independent_replay_comparison'] = dict(source=proposal['source'],
        expected_regional_logs_sha256=completed['regional_logs_sha256'],
        expected_log2_upper_screen=completed['log2_upper_screen'],
        complete_score_and_regional_hash_match=matches, comparison_is_not_a_directed_certificate=True)
    check()
    replay_pin = adaptive.provider.search.write_new(output/'replay.json', replay)
    summary = dict(schema='gamma-original-weight-winner-replay-links-FLOAT-v1',
        active_groups=6400, proposal=proposal_pin, replay=replay_pin,
        source=proposal['source'], source_sha256=dict(sorted(pins.items())),
        marker_count=len(indices), complete_score=completed['log2_upper_screen'],
        exact_replay_match=matches, directed=False, selected_weight_certificate=False,
        whole_code_certificate=False, no_weight_interpolation=True)
    check()
    pin = adaptive.provider.search.write_new(output/'summary.json', summary)
    print(json.dumps(dict(summary=pin, marker_count=len(indices), exact_replay_match=matches,
                          score=summary['complete_score'])), flush=True)
    if not matches:
        raise ArithmeticError('Fresh selected literal replay differs from completed source')


if __name__ == '__main__':
    main()
