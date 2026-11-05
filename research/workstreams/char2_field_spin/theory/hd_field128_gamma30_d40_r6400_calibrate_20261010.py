"""Bounded original-r6400 FLOAT search, adaptive markers, and fresh winner replay."""

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
    source_path = here/'hd_field128_gamma30_d40_r5888_winner_replay_20261010/proposal.json'
    raw = source_path.read_bytes()
    source_hash = sha256(raw).hexdigest()
    if source_hash != '43cea88db107cfc364ab5989978269fa97b17d0068e32ca5d24aa930887e4754':
        raise RuntimeError('Frozen exact17 source changed')
    expected = {
        adaptive: 'f96e9cf781b40828a6aa83880879f904909a6411fe02003cfc0a81d3c9331e61',
        search: '99ca1cde5d29e20eef67f2aa975bec2fc3f8d5e27ac3a621a5685fbe5072ab05'}
    for module, digest in expected.items():
        if sha256(Path(module.__file__).read_bytes()).hexdigest() != digest:
            raise RuntimeError('Reviewed helper changed: '+module.__name__)
    source = json.loads(raw)
    if (source['active_groups'] != 5888 or source['state_dimension'] != 30
            or source['target_distance'] != '2/5'
            or source['construction'] != adaptive.provider.gamma.construction(30)
            or len(source['witness']['markers']) != 17):
        raise ValueError('Unexpected exact literal source scope')
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
        return (row is not None and row['profile_pairs_evaluated'] == 897**2
                and row['regional_logs_sha256'] == previous['regional_logs_sha256']
                and row['log2_upper_screen'] == previous['log2_upper_screen'])
    check()
    output = here/'hd_field128_gamma30_d40_r6400_calibration_20261010'
    output.mkdir(exist_ok=False)
    seed = dict(schema='gamma-original-weight-positive-literal-proposal-v1',
        construction=source['construction'], state_dimension=30, target_distance='2/5', active_groups=6400,
        source=dict(path=str(source_path), sha256=source_hash), source_sha256=dict(sorted(pins.items())),
        source_calibration_active_groups=5888, new_requested_active_groups=6400,
        witness=dict(source['witness']), marker_count=17,
        source_role='Exact positive literals only; no source calibration or numerical bound transferred',
        old_numerical_bounds_reused=False, directed=False, selected_weight_certificate=False,
        whole_code_certificate=False, no_weight_interpolation=True)
    seed_pin = adaptive.provider.search.write_new(output/'seed.json', seed)
    completed_pins = {}
    def checkpoint(row):
        check()
        row['source_sha256'][runner.name] = runner_hash
        pin = adaptive.provider.search.write_new(output/f"evaluation_{row['evaluation_index']:03d}.json", row)
        completed_pins[row['evaluation_index']] = pin
        print(json.dumps(dict(stage='search-trial', evaluation=row['evaluation_index'],
            score=row['log2_upper_screen'], largest=row['largest_profile'],
            elapsed_seconds=row['elapsed_seconds'])), flush=True)
        return pin
    searched = search.run_search(output/'seed.json', 6400, max_evaluations=40,
                                 max_seconds=120., checkpoint=checkpoint)
    searched['source_sha256'][runner.name] = runner_hash
    check()
    search_pin = adaptive.provider.search.write_new(output/'search_summary.json', searched)
    best = searched['best']
    report = dict(schema='gamma-original-weight-calibration-links-FLOAT-v1', active_groups=6400,
        construction=source['construction'], state_dimension=30, target_distance='2/5', seed=seed_pin,
        search=search_pin, source_sha256=dict(sorted(pins.items())), openblas_num_threads='1',
        directed=False, selected_weight_certificate=False, whole_code_certificate=False,
        no_weight_interpolation=True, status='no-complete-search-trial')
    if best is not None:
        best_pin = completed_pins[best['evaluation_index']]
        report.update(best_trial=best_pin, search_best_score=best['log2_upper_screen'])
        result = adaptive.run(best_pin['path'], 6400, max_seconds=90., max_fits=17,
                              max_calls=128, target_bits=80)
        check()
        result_pin = adaptive.provider.search.write_new(output/'adaptive.json', result)
        baseline_matches = bool(result['history']) and matches(result['history'][0], best)
        report.update(adaptive=result_pin, status=result['status'],
                      fresh_adaptive_baseline_matches_search=baseline_matches)
        completed = result['last_complete_step']
        if not baseline_matches:
            adaptive.provider.search.write_new(output/'failed_baseline_summary.json', report)
            raise ArithmeticError('Fresh adaptive baseline differs from exact selected search trial')
        report['complete_score'] = completed['log2_upper_screen']
        if completed['log2_upper_screen'] <= -80:
            indices = result['winning_marker_indices_at_last_complete']
            markers, witness = result['witness']['markers'], result['witness']
            if (not indices or indices != sorted(set(indices))
                    or any(type(i) is not int or not 0 <= i < len(markers) for i in indices)):
                raise ArithmeticError('Invalid complete first-winner selection')
            proposal = dict(schema='gamma-dense-selected-positive-literal-proposal-v1',
                construction=source['construction'], state_dimension=30, target_distance='2/5', active_groups=6400,
                source=result_pin, source_sha256=dict(sorted(pins.items())),
                selected_source_marker_indices=indices, source_marker_count=len(markers), marker_count=len(indices),
                witness=dict(z=witness['z'], alpha=witness['alpha'], markers=[markers[i] for i in indices],
                             grid=witness['grid'], state_dimension=30),
                source_role='Exact positive literal selection only; no numerical cap or score inherited',
                selection_role='Fresh complete first-winner selection, pending independent subset replay',
                old_numerical_bounds_reused=False, directed=False, selected_weight_certificate=False,
                whole_code_certificate=False, no_weight_interpolation=True)
            proposal_pin = adaptive.provider.search.write_new(output/'proposal.json', proposal)
            replay = adaptive.run(output/'proposal.json', 6400, max_seconds=90., max_fits=0,
                                  max_calls=128, target_bits=80)
            matched = replay['pending_step'] is None and matches(replay['last_complete_step'], completed)
            replay['independent_replay_comparison'] = dict(source=result_pin,
                expected_regional_logs_sha256=completed['regional_logs_sha256'],
                expected_log2_upper_screen=completed['log2_upper_screen'],
                complete_score_and_regional_hash_match=matched, comparison_is_not_a_directed_certificate=True)
            check()
            replay_pin = adaptive.provider.search.write_new(output/'replay.json', replay)
            report.update(proposal=proposal_pin, replay=replay_pin, marker_count=len(indices),
                          exact_replay_match=matched)
            if not matched:
                adaptive.provider.search.write_new(output/'failed_replay_summary.json', report)
                raise ArithmeticError('Fresh winner subset changed complete score or regional cap hash')
    check()
    final_pin = adaptive.provider.search.write_new(output/'summary.json', report)
    print(json.dumps(dict(summary=final_pin, status=report['status'], score=report.get('complete_score'),
                          marker_count=report.get('marker_count'), exact_replay_match=report.get('exact_replay_match'))), flush=True)


if __name__ == '__main__':
    main()
