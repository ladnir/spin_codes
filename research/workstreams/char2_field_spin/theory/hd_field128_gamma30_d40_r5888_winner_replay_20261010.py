"""Fresh FLOAT first-winner selection and independent exact-literal subset replay."""

from hashlib import sha256
import json
from pathlib import Path

import hd_field128_gamma_midpoint_dense_adaptive as adaptive


def main():
    here = Path(__file__).resolve().parent
    runner = Path(__file__).resolve()
    runner_hash = sha256(runner.read_bytes()).hexdigest()
    folder = here/'hd_field128_gamma30_d40_r5888_dense_mark_search_20261010'
    source_path, summary_path = folder/'evaluation_038.json', folder/'summary.json'
    source_raw, summary_raw = source_path.read_bytes(), summary_path.read_bytes()
    source_hash, summary_hash = sha256(source_raw).hexdigest(), sha256(summary_raw).hexdigest()
    if source_hash != 'f20fe8bcc909007924eaa519845d8e0c0118dcf2f661ad57fca48a2162624cae':
        raise RuntimeError('Frozen best trial changed')
    if summary_hash != '265f74954d774ff6c9f1275e4ecf55ced15137d7f90520d39d7c2e3a9d01e7ea':
        raise RuntimeError('Frozen search summary changed')
    if sha256(Path(adaptive.__file__).read_bytes()).hexdigest() != 'f96e9cf781b40828a6aa83880879f904909a6411fe02003cfc0a81d3c9331e61':
        raise RuntimeError('Reviewed adaptive source changed')
    source, summary = json.loads(source_raw), json.loads(summary_raw)
    if (source != summary['best'] or summary['completed_trials'] != 40
            or source['evaluation_index'] != 38 or source['active_groups'] != 5888
            or source['state_dimension'] != 30 or source['target_distance'] != '2/5'
            or source['construction'] != adaptive.provider.gamma.construction(30)
            or source['profile_pairs_evaluated'] != 1153**2
            or not source['complete_selected_weight_screen']
            or source['log2_upper_screen'] != min(row['log2_upper_screen'] for row in summary['history'])
            or len(source['witness']['markers']) != 35):
        raise ValueError('Unexpected complete best-trial selection')
    pins = adaptive.source_pins()
    for name, digest in source['source_sha256'].items():
        if name in pins and pins[name] != digest:
            raise RuntimeError('Source conflict: '+name)
        pins[name] = digest
    pins[runner.name] = runner_hash
    def check():
        if source_path.read_bytes() != source_raw or summary_path.read_bytes() != summary_raw:
            raise RuntimeError('Source artifact changed')
        for name, digest in pins.items():
            if sha256((here/name).read_bytes()).hexdigest() != digest:
                raise RuntimeError('Source closure changed: '+name)
    def same_score(row):
        return (row is not None and row['profile_pairs_evaluated'] == 1153**2
                and row['regional_logs_sha256'] == source['regional_logs_sha256']
                and row['log2_upper_screen'] == source['log2_upper_screen'])
    check()
    output = here/'hd_field128_gamma30_d40_r5888_winner_replay_20261010'
    output.mkdir(exist_ok=False)
    baseline = adaptive.run(source_path, 5888, max_seconds=90., max_fits=0, max_calls=128, target_bits=80)
    baseline['source_sha256'][runner.name] = runner_hash
    check()
    baseline_pin = adaptive.provider.search.write_new(output/'winner_reconstruction.json', baseline)
    if baseline['pending_step'] is not None or not same_score(baseline['last_complete_step']):
        raise ArithmeticError('Fresh baseline does not reproduce selected trial score and cap hash')
    indices = baseline['winning_marker_indices_at_last_complete']
    markers = baseline['witness']['markers']
    if (not indices or indices != sorted(set(indices))
            or any(type(i) is not int or not 0 <= i < len(markers) for i in indices)):
        raise ArithmeticError('Invalid saved first-winner indices')
    witness = source['witness']
    proposal = dict(schema='gamma-dense-selected-positive-literal-proposal-v1',
        construction=source['construction'], state_dimension=30, target_distance='2/5', active_groups=5888,
        source=baseline_pin, original_trial=dict(path=str(source_path), sha256=source_hash),
        search_summary=dict(path=str(summary_path), sha256=summary_hash), source_sha256=dict(sorted(pins.items())),
        selected_source_marker_indices=indices, source_marker_count=35, marker_count=len(indices),
        witness=dict(z=witness['z'], alpha=witness['alpha'], markers=[markers[i] for i in indices],
                     grid=witness['grid'], state_dimension=30),
        source_role='Exact positive literal selection only; no numerical cap or score inherited',
        selection_role='Fresh first-winner indices are a heuristic pending independent subset replay',
        old_numerical_bounds_reused=False, directed=False, selected_weight_certificate=False,
        whole_code_certificate=False, no_weight_interpolation=True)
    proposal_pin = adaptive.provider.search.write_new(output/'proposal.json', proposal)
    replay = adaptive.run(output/'proposal.json', 5888, max_seconds=90., max_fits=0,
                          max_calls=128, target_bits=80)
    matched = replay['pending_step'] is None and same_score(replay['last_complete_step'])
    replay['independent_replay_comparison'] = dict(source=baseline_pin,
        original_trial=dict(path=str(source_path), sha256=source_hash),
        expected_regional_logs_sha256=source['regional_logs_sha256'],
        expected_log2_upper_screen=source['log2_upper_screen'], complete_score_and_regional_hash_match=matched,
        comparison_is_not_a_directed_certificate=True)
    check()
    replay_pin = adaptive.provider.search.write_new(output/'replay.json', replay)
    report = dict(schema='gamma-dense-winner-subset-replay-summary-FLOAT-v1', active_groups=5888,
        source_sha256=dict(sorted(pins.items())), original_trial=dict(path=str(source_path), sha256=source_hash),
        winner_reconstruction=baseline_pin, proposal=proposal_pin, replay=replay_pin,
        selected_source_marker_indices=indices, marker_count=len(indices), exact_replay_match=matched,
        log2_upper_screen=replay['last_complete_step']['log2_upper_screen'] if replay['last_complete_step'] else None,
        directed=False, selected_weight_certificate=False, whole_code_certificate=False, no_weight_interpolation=True)
    check()
    report_pin = adaptive.provider.search.write_new(output/'summary.json', report)
    print(json.dumps(dict(summary=report_pin, proposal=proposal_pin, replay=replay_pin,
        exact_match=matched, marker_count=len(indices), score=report['log2_upper_screen'])), flush=True)
    if not matched:
        raise ArithmeticError('Fresh subset replay did not preserve complete score and regional cap hash')


if __name__ == '__main__':
    main()
