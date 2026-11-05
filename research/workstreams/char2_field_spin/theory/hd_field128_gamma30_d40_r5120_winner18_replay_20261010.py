"""Select exact marker literals, then independently rebuild the complete FLOAT sum."""

from hashlib import sha256
import json
from pathlib import Path

import hd_field128_gamma_midpoint_dense_adaptive as adaptive


def main():
    here = Path(__file__).resolve().parent
    runner = Path(__file__).resolve()
    runner_hash = sha256(runner.read_bytes()).hexdigest()
    source_path = here/'hd_field128_gamma30_d40_r5120_dense_adaptive_float_20261010.json'
    source_raw = source_path.read_bytes()
    source_hash = sha256(source_raw).hexdigest()
    if source_hash != '7f2dc681c668c30e4462e9a286f29f03da26fb4dea276a69a1ab5135241b2ac8':
        raise RuntimeError('Frozen117-marker source changed')
    if sha256(Path(adaptive.__file__).read_bytes()).hexdigest() != 'f96e9cf781b40828a6aa83880879f904909a6411fe02003cfc0a81d3c9331e61':
        raise RuntimeError('Reviewed adaptive source changed')
    source = json.loads(source_raw)
    indices = source['winning_marker_indices_at_last_complete']
    markers = source['witness']['markers']
    if (source['active_groups'] != 5120 or source['state_dimension'] != 30 or source['target_distance'] != '2/5'
            or source['construction'] != adaptive.provider.gamma.construction(30)
            or source['pending_step'] is not None or source['unscored_applied_markers']
            or len(markers) != 117 or source['last_complete_step']['marker_count'] != 117
            or len(indices) != 18 or indices != sorted(set(indices))
            or any(type(i) is not int or not 0 <= i < 117 for i in indices)):
        raise ValueError('Unexpected frozen complete marker selection')
    pins = adaptive.source_pins()
    for name, digest in source['source_sha256'].items():
        if name in pins and pins[name] != digest:
            raise RuntimeError('Source conflict: '+name)
        pins[name] = digest
    pins[runner.name] = runner_hash
    def check():
        if source_path.read_bytes() != source_raw:
            raise RuntimeError('Source artifact changed')
        for name, digest in pins.items():
            if sha256((here/name).read_bytes()).hexdigest() != digest:
                raise RuntimeError('Source closure changed: '+name)
    check()
    output = here/'hd_field128_gamma30_d40_r5120_winner18_replay_20261010'
    output.mkdir(exist_ok=False)
    witness = source['witness']
    proposal = dict(schema='gamma-dense-selected-positive-literal-proposal-v1',
        construction=source['construction'], state_dimension=30, target_distance='2/5', active_groups=5120,
        source=dict(path=str(source_path), sha256=source_hash), source_sha256=dict(sorted(pins.items())),
        selected_source_marker_indices=indices, source_marker_count=117, marker_count=18,
        witness=dict(z=witness['z'], alpha=witness['alpha'], markers=[markers[i] for i in indices],
                     grid=witness['grid'], state_dimension=30),
        source_role='Exact positive literal selection only; numerical caps and scores are not inherited',
        selection_role='Saved first-winner indices are a proposal heuristic until fresh replay',
        old_numerical_bounds_reused=False, directed=False, selected_weight_certificate=False,
        whole_code_certificate=False, scope='Untested18-marker exact proposal at original weight5120')
    proposal_pin = adaptive.provider.search.write_new(output/'proposal.json', proposal)
    result = adaptive.run(output/'proposal.json', 5120, max_seconds=90., max_fits=0, max_calls=128, target_bits=80)
    completed = result['last_complete_step']
    old = source['last_complete_step']
    matched = (completed is not None and result['pending_step'] is None
               and completed['profile_pairs_evaluated'] == 1537**2
               and completed['regional_logs_sha256'] == old['regional_logs_sha256']
               and completed['log2_upper_screen'] == old['log2_upper_screen'])
    result.update(exact_literal_proposal=proposal_pin,
        independent_replay_comparison=dict(source=dict(path=str(source_path), sha256=source_hash),
            expected_regional_logs_sha256=old['regional_logs_sha256'], expected_log2_upper_screen=old['log2_upper_screen'],
            complete_score_and_regional_hash_match=matched, comparison_is_not_a_directed_certificate=True))
    check()
    result_pin = adaptive.provider.search.write_new(output/'replay.json', result)
    print(json.dumps(dict(proposal=proposal_pin, replay=result_pin, exact_match=matched,
                          score=completed['log2_upper_screen'] if completed else None,
                          status=result['status'], elapsed_seconds=result['elapsed_seconds'])), flush=True)
    if not matched:
        raise RuntimeError('Fresh18-marker complete replay did not match the frozen117-marker score and hash')


if __name__ == '__main__':
    main()
