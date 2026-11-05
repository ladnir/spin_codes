"""File-backed individually audited upper-gap events; no new aggregate bound."""

from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from shutil import disk_usage

import hd_field128_gamma_profile_certificate as profiles


def main():
    here, source = Path(__file__).resolve().parent, Path(__file__).resolve()
    source_hash = sha256(source.read_bytes()).hexdigest()
    groups = (
        ('hd_field128_gamma30_d40_r5984_winner13_directed256_20261010',
         'hd_field128_gamma30_d40_r5984_endpoint_calibration_20261010/proposal.json',
         '74406acabc8a132172642bab535c58bdfcd3c711b591f781e9998e7d9caee4a0',
         ((5984, 6047), (6048, 6111), (6112, 6175))),
        ('hd_field128_gamma30_d40_r6400_winner12_directed256_20261010',
         'hd_field128_gamma30_d40_r6400_winner_replay_20261010/proposal.json',
         'd7719fe4a2fafbc3a926bbb7189524094cb8fc45821e308cb2d676d941d7826f',
         ((6176, 6239), (6240, 6303), (6304, 6367), (6368, 6431), (6432, 6495))),
        ('hd_field128_gamma30_d40_r6677_winner19_directed256_20261010',
         'hd_field128_gamma30_d40_r6677_mark_search_20261010/proposal.json',
         '01d34522cd16e96c325f6347ea67e724d18bed3b94e1f6971a500609c0e186dd',
         ((6496, 6559), (6560, 6623), (6624, 6677))))
    inputs = {}
    def load(path):
        path = path.resolve()
        raw = path.read_bytes()
        inputs[path] = raw
        return json.loads(raw), dict(path=str(path), sha256=sha256(raw).hexdigest())
    events, banks, covered, construction = [], [], set(), None
    for folder_name, proposal_name, expected_hash, ranges in groups:
        folder = here/folder_name
        proposal, proposal_pin = load(here/proposal_name)
        if proposal_pin['sha256'] != expected_hash:
            raise ValueError('Exact reviewed literal proposal changed')
        if construction is None:
            construction = proposal['construction']
        if (proposal['construction'] != construction or proposal['state_dimension'] != 30
                or proposal['target_distance'] != '2/5'):
            raise ValueError('Different construction in upper-gap manifest')
        bank, bank_pin = load(folder/'bank.json')
        if any(bank['witness'][key] != proposal['witness'][key] for key in ('z', 'alpha', 'markers')):
            raise ValueError('Actual bank and literal proposal differ')
        banks.append(dict(bank=bank_pin, proposal=proposal_pin))
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
                    or record['context']['construction'] != construction
                    or record['context']['state_dimension'] != 30 or record['context']['target_distance'] != '2/5'
                    or any(record['context']['witness'][key] != proposal['witness'][key]
                           for key in ('z', 'alpha', 'markers'))):
                raise ValueError('Invalid, partial, or overlapping completed band')
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
                bank=bank_pin, proposal=proposal_pin, log2_upper_rational=record['log2_upper_rational'],
                tree_nodes=record['tree_nodes'], accepted_leaves=len(record['accepted_leaves']),
                required_event_failure_bits=54))
    if covered != set(range(5984, 6678)) or len(events) != 11:
        raise ValueError('Unexpected upper-gap weight union')
    result = dict(schema='gamma-audited-band-link-manifest-v1', construction=construction,
        state_dimension=30, target_distance='2/5', banks=banks,
        individually_audited_events=events, event_count=len(events), explicit_covered_range=[5984, 6677],
        covered_weight_count=len(covered), nonoverlapping_weight_ranges=True, no_weight_interpolation=True,
        excluded_partial_events=[], source_sha256={source.name: source_hash},
        all_linked_input_bytes_rechecked=True, joint_probability_claim=False,
        index_update_performed=False, whole_code_certificate=False,
        scope='File-backed list of individually audited events; no new coefficient or arithmetic replay')
    if sha256(source.read_bytes()).hexdigest() != source_hash or any(path.read_bytes() != raw for path, raw in inputs.items()):
        raise RuntimeError('Manifest input/source changed')
    if disk_usage(here).free < 750_000_000:
        raise OSError('Publication guard: less than750MB disk space remains')
    pin = profiles.write_new(here/'hd_field128_gamma30_d40_upper_gap_5984_6677_manifest_20261010.json', result)
    print(json.dumps(dict(manifest=pin, events=len(events), covered_weights=len(covered), range=[5984, 6677])), flush=True)


if __name__ == '__main__':
    main()
