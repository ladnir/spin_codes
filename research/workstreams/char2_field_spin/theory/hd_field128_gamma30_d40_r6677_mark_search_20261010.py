"""One bounded fixed19-marker original-r6677 search and fresh winner replay."""

import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'

from hashlib import sha256
import json
from pathlib import Path

import hd_field128_gamma_midpoint_dense_adaptive as adaptive
import hd_field128_gamma_midpoint_dense_mark_search as search


def main():
    here, runner = Path(__file__).resolve().parent, Path(__file__).resolve()
    runner_hash = sha256(runner.read_bytes()).hexdigest()
    source_path = here/'hd_field128_gamma30_d40_r6677_endpoint_calibration_20261010/proposal.json'
    raw = source_path.read_bytes()
    source_hash = sha256(raw).hexdigest()
    if source_hash != '73b90c0a6d0fcfe064466fca718bab98083f9f818709b4e44ef020473638936a':
        raise RuntimeError('Reviewed exact19 parent changed')
    source = json.loads(raw)
    if (source['active_groups'] != 6677 or source['state_dimension'] != 30
            or source['target_distance'] != '2/5'
            or source['construction'] != adaptive.provider.gamma.construction(30)
            or len(source['witness']['markers']) != 19):
        raise ValueError('Unexpected complete literal parent scope')
    pins = adaptive.source_pins()
    for added in (search.source_pins(), source['source_sha256']):
        for name, digest in added.items():
            if name in pins and pins[name] != digest:
                raise RuntimeError('Source closure conflict: '+name)
            pins[name] = digest
    pins[runner.name] = runner_hash
    def check():
        if source_path.read_bytes() != raw:
            raise RuntimeError('Exact literal parent changed')
        for name, digest in pins.items():
            if sha256((here/name).read_bytes()).hexdigest() != digest:
                raise RuntimeError('Source closure changed: '+name)
    def matches(row, previous):
        return (row is not None and row['profile_pairs_evaluated'] == 758**2
                and row['regional_logs_sha256'] == previous['regional_logs_sha256']
                and row['log2_upper_screen'] == previous['log2_upper_screen'])
    check()
    output = here/'hd_field128_gamma30_d40_r6677_mark_search_20261010'
    output.mkdir(exist_ok=False)
    trials = {}
    def checkpoint(row):
        check()
        row['source_sha256'][runner.name] = runner_hash
        pin = adaptive.provider.search.write_new(output/f"evaluation_{row['evaluation_index']:03d}.json", row)
        trials[row['evaluation_index']] = pin
        print(json.dumps(dict(stage='search-trial', evaluation=row['evaluation_index'],
            score=row['log2_upper_screen'], elapsed_seconds=row['elapsed_seconds'])), flush=True)
        return pin
    result = search.run_search(source_path, 6677, max_evaluations=40,
                               max_seconds=120., checkpoint=checkpoint)
    result['source_sha256'][runner.name] = runner_hash
    check()
    result_pin = adaptive.provider.search.write_new(output/'search_summary.json', result)
    best = result['best']
    report = dict(schema='gamma-original-weight-search-winner-links-FLOAT-v1', active_groups=6677,
        construction=source['construction'], state_dimension=30, target_distance='2/5', search=result_pin,
        source=dict(path=str(source_path), sha256=source_hash), source_sha256=dict(sorted(pins.items())),
        status='no-complete-search-trial', directed=False, selected_weight_certificate=False,
        whole_code_certificate=False, no_weight_interpolation=True)
    if best is not None:
        best_pin = trials[best['evaluation_index']]
        report.update(best_trial=best_pin, complete_score=best['log2_upper_screen'])
        baseline = adaptive.run(best_pin['path'], 6677, max_seconds=90., max_fits=0,
                                max_calls=128, target_bits=4096)
        check()
        baseline_pin = adaptive.provider.search.write_new(output/'baseline_replay.json', baseline)
        baseline_matches = baseline['pending_step'] is None and matches(baseline['last_complete_step'], best)
        report.update(baseline=baseline_pin, fresh_baseline_matches_search=baseline_matches)
        if not baseline_matches:
            adaptive.provider.search.write_new(output/'failed_baseline_summary.json', report)
            raise ArithmeticError('Fresh complete baseline differs from exact best search trial')
        indices = baseline['winning_marker_indices_at_last_complete']
        markers, witness = baseline['witness']['markers'], baseline['witness']
        if (not indices or indices != sorted(set(indices))
                or any(type(i) is not int or not 0 <= i < len(markers) for i in indices)):
            raise ArithmeticError('Invalid complete first-winner selection')
        proposal = dict(schema='gamma-dense-selected-positive-literal-proposal-v1',
            construction=source['construction'], state_dimension=30, target_distance='2/5', active_groups=6677,
            source=baseline_pin, source_sha256=dict(sorted(pins.items())),
            selected_source_marker_indices=indices, source_marker_count=len(markers), marker_count=len(indices),
            witness=dict(z=witness['z'], alpha=witness['alpha'], markers=[markers[i] for i in indices],
                         grid=witness['grid'], state_dimension=30),
            source_role='Exact positive literal selection only; no numerical cap or score inherited',
            selection_role='Fresh complete first-winner selection, pending independent subset replay',
            old_numerical_bounds_reused=False, directed=False, selected_weight_certificate=False,
            whole_code_certificate=False, no_weight_interpolation=True)
        proposal_pin = adaptive.provider.search.write_new(output/'proposal.json', proposal)
        replay = adaptive.run(output/'proposal.json', 6677, max_seconds=90., max_fits=0,
                              max_calls=128, target_bits=4096)
        matched = replay['pending_step'] is None and matches(replay['last_complete_step'], best)
        replay['independent_replay_comparison'] = dict(source=best_pin,
            expected_regional_logs_sha256=best['regional_logs_sha256'],
            expected_log2_upper_screen=best['log2_upper_screen'],
            complete_score_and_regional_hash_match=matched, comparison_is_not_a_directed_certificate=True)
        check()
        replay_pin = adaptive.provider.search.write_new(output/'replay.json', replay)
        report.update(proposal=proposal_pin, replay=replay_pin, marker_count=len(indices),
                      exact_replay_match=matched, status='complete-replay' if matched else 'replay-mismatch')
        if not matched:
            adaptive.provider.search.write_new(output/'failed_replay_summary.json', report)
            raise ArithmeticError('Fresh subset replay differs from complete best trial')
    check()
    pin = adaptive.provider.search.write_new(output/'summary.json', report)
    print(json.dumps(dict(summary=pin, score=report.get('complete_score'), markers=report.get('marker_count'),
                          exact_replay_match=report.get('exact_replay_match'))), flush=True)


if __name__ == '__main__':
    main()
