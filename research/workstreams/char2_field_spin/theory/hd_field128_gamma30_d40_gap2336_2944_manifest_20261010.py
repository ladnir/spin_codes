"""File-backed individually audited assigned-gap events; no new aggregate bound."""

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
        ('hd_field128_gamma30_d40_r2336_gap_retry02_directed256_20261010',
         'hd_field128_gamma30_d40_r2336_gap_calibration_20261010/proposal.json',
         '60dbab49585c3a5b4988c2b783fcf0c2c98b23961b923be2490d30c81c194a47',
         ((2257, 2303), (2304, 2367), (2368, 2431))),
        ('hd_field128_gamma30_d40_r2944_gap_directed256_20261010',
         'hd_field128_gamma30_d40_r2944_gap_calibration_20261010/proposal.json',
         '8a0f5f235bcc87cd83e18d62f428ad67fbf41ab3d3c829bb5fbc2993aa0a0120',
         ((2848, 2911), (2912, 2975), (2976, 3007))))
    inputs = {}
    def load(path):
        path = path.resolve()
        raw = path.read_bytes()
        inputs[path] = raw
        return json.loads(raw), dict(path=str(path), sha256=sha256(raw).hexdigest())
    failure, failure_pin = load(here/'hd_field128_gamma30_d40_r2336_gap_directed256_20261010/band_2304_2367/failure.json')
    if (failure_pin['sha256'] != '5d7ba28272946fdc8664c96b6714c709348bde8b2d8b1767dfa672de1e238dee'
            or failure['complete_requested_band'] is not False or failure['log2_upper_rational'] is not None
            or failure['requested_weight_range'] != [2304, 2367]
            or failure['failure']['message'] != 'Checkpoint guard: less than750MB disk space remains'):
        raise ValueError('Unexpected excluded resource-stop receipt')
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
            raise ValueError('Different construction in assigned-gap manifest')
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
    if covered != set(range(2257, 2432)) | set(range(2848, 3008)) or len(events) != 6:
        raise ValueError('Unexpected assigned-gap weight union')
    result = dict(schema='gamma-audited-band-link-manifest-v1', construction=construction,
        state_dimension=30, target_distance='2/5', banks=banks,
        individually_audited_events=events, event_count=len(events),
        explicit_covered_ranges=[[2257, 2431], [2848, 3007]],
        covered_weight_count=len(covered), nonoverlapping_weight_ranges=True, no_weight_interpolation=True,
        excluded_partial_events=[dict(receipt=failure_pin, reason='Resource guard; no directed band bound')],
        source_sha256={source.name: source_hash},
        all_linked_input_bytes_rechecked=True, joint_probability_claim=False,
        index_update_performed=False, whole_code_certificate=False,
        scope='File-backed list of individually audited events; no new coefficient or arithmetic replay')
    if sha256(source.read_bytes()).hexdigest() != source_hash or any(path.read_bytes() != raw for path, raw in inputs.items()):
        raise RuntimeError('Manifest input/source changed')
    if disk_usage(here).free < 750_000_000:
        raise OSError('Publication guard: less than750MB disk space remains')
    pin = profiles.write_new(here/'hd_field128_gamma30_d40_gaps2257_2431_2848_3007_manifest_20261010.json', result)
    print(json.dumps(dict(manifest=pin, events=len(events), covered_weights=len(covered),
        ranges=[[2257, 2431], [2848, 3007]])), flush=True)


if __name__ == '__main__':
    main()
