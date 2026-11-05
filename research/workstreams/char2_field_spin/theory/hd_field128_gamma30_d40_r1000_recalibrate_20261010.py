"""Two authorized FLOAT prune/refit cycles, then at most one local mark search."""

from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path

import hd_field128_gamma_midpoint_prune as prune
import hd_field128_gamma_midpoint_mark_search as search


HERE = Path(__file__).resolve().parent
INITIAL = HERE/'hd_field128_gamma30_d40_r1000_trial39_pruned_enriched_float_20261010.json'
OUTPUT = HERE/'hd_field128_gamma30_d40_r1000_recalibration_cycles_20261010'
write_new = search.provider.search.write_new


def source_pins():
    pins = search.source_pins()
    for path in (Path(prune.__file__), Path(__file__)):
        pins[path.name] = sha256(path.read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def projection(source_path, pins):
    """Bind a selected completed round and its exact witness prefix explicitly."""
    raw = Path(source_path).read_bytes()
    value = json.loads(raw)
    if value['requested_weights'] != [1000] or value['state_dimension'] != 30 or value['target_distance'] != '2/5':
        raise ValueError('Expected explicitly calibrated gamma30/d40 r1000 enrichment')
    row = value['last_complete_round']
    if not row or not row['complete_all_requested_weights'] or row['pending_weights']:
        raise ValueError('No complete selected-weight round to project')
    if len(row['completed']) != 1 or row['completed'][0]['active_groups'] != 1000:
        raise ValueError('Ambiguous selected-weight calibration')
    count, witness = row['marker_count'], deepcopy(value['witness'])
    witness['markers'] = witness['markers'][:count]
    witness['marker_log_odds_exact'] = witness['marker_log_odds_exact'][:count]
    if len(witness['markers']) != count:
        raise ValueError('Missing completed-round marker prefix')
    closure = dict(pins)
    for name, digest in value['source_sha256'].items():
        if name in closure and closure[name] != digest:
            raise ValueError('Projection source conflict: '+name)
        closure[name] = digest
    return dict(schema='gamma-r1000-completed-round-exact-proposal-projection-v1', active_groups=1000,
                requested_weights=[1000], state_dimension=30, target_distance='2/5', construction=value['construction'],
                witness=witness, marker_count=count, source=dict(path=str(Path(source_path).resolve()), sha256=sha256(raw).hexdigest()),
                source_sha256=dict(sorted(closure.items())), source_completed_round=row['round'],
                source_had_pending_round=value['pending_round'] is not None,
                projection_rule='Explicit active_groups1000 from requested_weights and complete round; exact witness prefix only',
                no_numerical_caps_or_bounds_imported=True, directed=False, proposal_only=True, whole_code_certificate=False)


def main():
    OUTPUT.mkdir(exist_ok=False)
    pins, source_path, cycles, last_score = source_pins(), INITIAL, [], None
    for number in (1, 2):
        projected = projection(source_path, pins)
        proposal_path = OUTPUT/f'cycle_{number:02d}_source_proposal.json'
        proposal_pin = write_new(proposal_path, projected)
        pruned = prune.run(proposal_path, 1000)
        pruned_path = OUTPUT/f'cycle_{number:02d}_pruned.json'
        pruned_pin = write_new(pruned_path, pruned)
        step = dict(cycle=number, proposal=proposal_pin, pruned=pruned_pin,
                    retained_markers=pruned['marker_count'], score_before=pruned['rebuilt_pruned']['log2_upper_screen'])
        cycles.append(step)
        print(json.dumps(dict(stage='pruned', **step)), flush=True)
        if pruned['marker_count'] > 124:
            step['status'] = 'insufficient-four-marker-capacity'
            last_score = step['score_before']
            break
        enriched = search.enrich.run(pruned_path, [1000], target_bits=80, max_seconds=90., max_rounds=12, max_markers=128)
        source_path = OUTPUT/f'cycle_{number:02d}_enriched.json'
        step['enriched'] = write_new(source_path, enriched)
        complete = enriched['last_complete_round']
        last_score = complete['completed'][0]['log2_upper_screen'] if complete else None
        step.update(status=enriched['status'], score_after=last_score, marker_count=enriched['marker_count'])
        if last_score is None or last_score <= -80.:
            break
    mark_search = None
    if last_score is not None and last_score > -80.:
        projected = projection(source_path, pins)
        mark_path = OUTPUT/'local_mark_source_proposal.json'
        mark_pin = write_new(mark_path, projected)
        directory = OUTPUT/'local_mark_search'
        directory.mkdir(exist_ok=False)
        def checkpoint(row):
            pin = write_new(directory/f"evaluation_{row['evaluation_index']:03d}.json", row)
            print(json.dumps(dict(stage='mark-trial', evaluation=row['evaluation_index'],
                                  score=row['log2_upper_screen'], largest=row['largest_profile'])), flush=True)
            return pin
        result = search.run_search(mark_path, 1000, max_evaluations=40, max_seconds=120., checkpoint=checkpoint)
        mark_search = dict(source=mark_pin, result=write_new(directory/'summary.json', result),
                           best_score=None if result['best'] is None else result['best']['log2_upper_screen'],
                           optimizer=result['optimizer'])
    if source_pins() != pins:
        raise RuntimeError('Orchestration source closure changed')
    result = dict(schema='gamma-r1000-bounded-recalibration-cycles-FLOAT-v1', active_groups=1000,
                  cycles=cycles, local_mark_search=mark_search, source_sha256=pins,
                  target_log2_score=-80, maximum_prune_refit_cycles=2, maximum_enrichment_rounds_per_cycle=12,
                  maximum_enrichment_seconds_per_cycle=90, maximum_search_evaluations=40, maximum_search_seconds=120,
                  directed=False, selected_weight_certificate=False, whole_code_certificate=False,
                  scope='Separate explicit r1000 FLOAT calibration steps; no coefficient certificate or interpolated coverage')
    print(json.dumps(dict(stage='saved', artifact=write_new(OUTPUT/'summary.json', result))), flush=True)


if __name__ == '__main__':
    main()
