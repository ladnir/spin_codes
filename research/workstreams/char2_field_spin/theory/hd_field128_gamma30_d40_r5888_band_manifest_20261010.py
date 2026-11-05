"""File-backed links to individually audited bands, not a new aggregate bound."""

from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

import hd_field128_gamma_profile_certificate as profiles


def main():
    here, source = Path(__file__).resolve().parent, Path(__file__).resolve()
    source_hash = sha256(source.read_bytes()).hexdigest()
    folder = here/'hd_field128_gamma30_d40_r5888_winner17_directed256_20261010'
    ranges = ((5345, 5407), (5408, 5471), (5472, 5535), (5536, 5599), (5600, 5663),
              (5664, 5727), (5728, 5791), (5792, 5855), (5856, 5919), (5920, 5983))
    inputs = {}
    def load(path):
        path = path.resolve()
        raw = path.read_bytes()
        pin = dict(path=str(path), sha256=sha256(raw).hexdigest())
        inputs[path] = raw
        return json.loads(raw), pin
    proposal, proposal_pin = load(here/'hd_field128_gamma30_d40_r5888_winner_replay_20261010/proposal.json')
    if proposal_pin['sha256'] != '43cea88db107cfc364ab5989978269fa97b17d0068e32ca5d24aa930887e4754':
        raise ValueError('Exact17 proposal changed')
    bank, bank_pin = load(folder/'bank.json')
    if any(bank['witness'][key] != proposal['witness'][key] for key in ('z', 'alpha', 'markers')):
        raise ValueError('Bank literal mismatch')
    events, covered = [], set()
    for first, last in ranges:
        directory = folder/f'band_{first:04d}_{last:04d}'
        record, record_pin = load(directory/'summary.json')
        audit, audit_pin = load(directory/'independent_audit.json')
        weights = list(range(first, last+1))
        if (record['requested_weight_range'] != [first, last] or record['requested_weights'] != weights
                or record['status'] != 'complete' or record['complete_requested_band'] is not True
                or record['directed'] is not True or record['covers_all_matching_profiles'] is not True
                or record['meets_requested_band_target'] is not True or record['whole_code_certificate'] is not False
                or record['zero'] is not False or F(record['log2_upper_rational']) > -54
                or covered.intersection(weights)
                or record['context']['construction'] != proposal['construction']
                or record['context']['state_dimension'] != 30 or record['context']['target_distance'] != '2/5'
                or any(record['context']['witness'][key] != proposal['witness'][key] for key in ('z', 'alpha', 'markers'))):
            raise ValueError('Invalid or overlapping completed band')
        if (audit['schema'] != 'gamma-midpoint-weight-band-independent-audit-v1'
                or audit['audit_result'] != 'PASS' or audit['input_summary'] != record_pin
                or audit['requested_weight_range'] != [first, last]
                or audit['complete_requested_band'] is not True
                or audit['log2_upper_rational'] != record['log2_upper_rational']
                or audit['coefficient_replay_in_this_audit'] is not False
                or audit['source_and_input_bytes_checked_at_exit'] is not True):
            raise ValueError('Missing matching completed-band audit')
        covered.update(weights)
        events.append(dict(weight_range=[first, last], summary=record_pin, independent_audit=audit_pin,
            log2_upper_rational=record['log2_upper_rational'], tree_nodes=record['tree_nodes'],
            accepted_leaves=len(record['accepted_leaves']), required_event_failure_bits=54))
    excluded, excluded_pin = load(folder/'band_5344_5407/summary.json')
    if (excluded['status'] != 'partial' or excluded['complete_requested_band'] is not False
            or excluded['log2_upper_rational'] is not None):
        raise ValueError('Expected preserved partial original band')
    if covered != set(range(5345, 5984)):
        raise ValueError('Unexpected manifest weight union')
    result = dict(schema='gamma-audited-band-link-manifest-v1', construction=proposal['construction'],
        state_dimension=30, target_distance='2/5', bank=bank_pin, proposal=proposal_pin,
        individually_audited_events=events, event_count=len(events), explicit_covered_range=[5345, 5983],
        covered_weight_count=len(covered), nonoverlapping_weight_ranges=True, no_weight_interpolation=True,
        excluded_partial_events=[dict(summary=excluded_pin, requested_weight_range=[5344, 5407],
            reason='One r5344 profile has insufficient bound; original whole-band receipt remains partial')],
        not_claimed_by_this_manifest=[[5280, 5344]],
        source_sha256={source.name: source_hash}, all_linked_input_bytes_rechecked=True,
        joint_probability_claim=False, index_update_performed=False, whole_code_certificate=False,
        scope='File-backed list of individually audited events; no new coefficient or arithmetic replay')
    if sha256(source.read_bytes()).hexdigest() != source_hash or any(path.read_bytes() != raw for path, raw in inputs.items()):
        raise RuntimeError('Manifest input/source changed')
    pin = profiles.write_new(folder/'audited_band_manifest.json', result)
    print(json.dumps(dict(manifest=pin, events=len(events), covered_weights=len(covered), range=[5345, 5983])), flush=True)


if __name__ == '__main__':
    main()
