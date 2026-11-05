"""Bounded FLOAT recalibration from exact positive literals, then winner replay."""

from hashlib import sha256
import json
from pathlib import Path

import hd_field128_gamma_midpoint_dense_adaptive as adaptive


def main():
    here = Path(__file__).resolve().parent
    runner = Path(__file__).resolve()
    runner_hash = sha256(runner.read_bytes()).hexdigest()
    source_path = here/'hd_field128_gamma30_d40_r5120_winner18_replay_20261010/proposal.json'
    source_raw = source_path.read_bytes()
    source_hash = sha256(source_raw).hexdigest()
    if source_hash != 'b37a28b0a992ca313536260a25547e639d300f542b0814345c7a3789cc448c8d':
        raise RuntimeError('Frozen exact18 source changed')
    if sha256(Path(adaptive.__file__).read_bytes()).hexdigest() != 'f96e9cf781b40828a6aa83880879f904909a6411fe02003cfc0a81d3c9331e61':
        raise RuntimeError('Reviewed adaptive helper changed')
    source = json.loads(source_raw)
    if (source['active_groups'] != 5120 or source['state_dimension'] != 30
            or source['target_distance'] != '2/5'
            or source['construction'] != adaptive.provider.gamma.construction(30)
            or len(source['witness']['markers']) != 18):
        raise ValueError('Unexpected exact literal source scope')
    pins = adaptive.source_pins()
    for name, digest in source['source_sha256'].items():
        if name in pins and pins[name] != digest:
            raise RuntimeError('Source closure conflict: '+name)
        pins[name] = digest
    pins[runner.name] = runner_hash
    def check():
        if source_path.read_bytes() != source_raw:
            raise RuntimeError('Literal parent changed')
        for name, digest in pins.items():
            if sha256((here/name).read_bytes()).hexdigest() != digest:
                raise RuntimeError('Source closure changed: '+name)
    check()
    output = here/'hd_field128_gamma30_d40_r5888_dense_calibration_20261010'
    output.mkdir(exist_ok=False)
    proposal = dict(schema='gamma-original-weight-positive-literal-proposal-v1',
        construction=source['construction'], state_dimension=30, target_distance='2/5', active_groups=5888,
        source=dict(path=str(source_path), sha256=source_hash), source_sha256=dict(sorted(pins.items())),
        source_calibration_active_groups=5120, new_requested_active_groups=5888,
        witness=dict(source['witness']), marker_count=18,
        source_role='Exact positive literals only; source calibration is not transferred',
        old_numerical_bounds_reused=False, directed=False, selected_weight_certificate=False,
        whole_code_certificate=False, no_weight_interpolation=True)
    proposal_pin = adaptive.provider.search.write_new(output/'initial_proposal.json', proposal)
    result = adaptive.run(output/'initial_proposal.json', 5888, max_seconds=90., max_fits=17,
                          max_calls=128, target_bits=80)
    check()
    result_pin = adaptive.provider.search.write_new(output/'adaptive.json', result)
    completed = result['last_complete_step']
    report = dict(initial_proposal=proposal_pin, adaptive=result_pin, status=result['status'],
        complete_score=completed['log2_upper_screen'] if completed else None,
        complete_marker_count=result['marker_count'], elapsed_seconds=result['elapsed_seconds'],
        source_sha256=dict(sorted(pins.items())), directed=False, whole_code_certificate=False)
    if completed is not None and completed['log2_upper_screen'] <= -80:
        indices = result['winning_marker_indices_at_last_complete']
        markers = result['witness']['markers']
        if (indices != sorted(set(indices)) or not indices
                or any(type(i) is not int or not 0 <= i < len(markers) for i in indices)
                or completed['profile_pairs_evaluated'] != 1153**2):
            raise ArithmeticError('Invalid completed winner selection')
        witness = result['witness']
        selected = dict(schema='gamma-dense-selected-positive-literal-proposal-v1',
            construction=source['construction'], state_dimension=30, target_distance='2/5', active_groups=5888,
            source=result_pin, source_sha256=dict(sorted(pins.items())),
            selected_source_marker_indices=indices, source_marker_count=len(markers), marker_count=len(indices),
            witness=dict(z=witness['z'], alpha=witness['alpha'], markers=[markers[i] for i in indices],
                         grid=witness['grid'], state_dimension=30),
            source_role='Exact positive literal selection only; no numerical bound inherited',
            selection_role='Saved first-winner indices are a heuristic pending fresh replay',
            old_numerical_bounds_reused=False, directed=False, selected_weight_certificate=False,
            whole_code_certificate=False, no_weight_interpolation=True)
        check()
        selected_pin = adaptive.provider.search.write_new(output/'winner_proposal.json', selected)
        replay = adaptive.run(output/'winner_proposal.json', 5888, max_seconds=90., max_fits=0,
                              max_calls=128, target_bits=80)
        fresh = replay['last_complete_step']
        matched = (fresh is not None and replay['pending_step'] is None
            and fresh['profile_pairs_evaluated'] == 1153**2
            and fresh['regional_logs_sha256'] == completed['regional_logs_sha256']
            and fresh['log2_upper_screen'] == completed['log2_upper_screen'])
        replay['independent_replay_comparison'] = dict(source=result_pin,
            expected_regional_logs_sha256=completed['regional_logs_sha256'],
            expected_log2_upper_screen=completed['log2_upper_screen'],
            complete_score_and_regional_hash_match=matched, comparison_is_not_a_directed_certificate=True)
        check()
        replay_pin = adaptive.provider.search.write_new(output/'winner_replay.json', replay)
        report.update(winner_proposal=selected_pin, winner_replay=replay_pin, exact_replay_match=matched,
                      winner_marker_count=len(indices))
        if not matched:
            adaptive.provider.search.write_new(output/'failed_replay_summary.json', report)
            raise ArithmeticError('Winner replay changed complete score or regional cap hash')
    check()
    summary_pin = adaptive.provider.search.write_new(output/'summary.json', report)
    print(json.dumps(dict(summary=summary_pin, **report)), flush=True)


if __name__ == '__main__':
    main()
