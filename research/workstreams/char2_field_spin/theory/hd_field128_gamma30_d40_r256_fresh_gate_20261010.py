"""One authorized bounded directed replay; no optimization or whole-code claim."""

from fractions import Fraction as F
from hashlib import sha256
import json
from math import isfinite
from pathlib import Path
from time import monotonic

import hd_field128_gamma_all_pair_certificate as provider
import hd_field128_lift_alias_cover as restored


DIRECTORY = Path(__file__).resolve().parent
SOURCE = DIRECTORY / 'hd_field128_gamma30_d40_r256_initial_20261010.json'
OUTPUT = DIRECTORY / 'hd_field128_gamma30_d40_r256_fresh_bank_20261010'


def write_new(name, record):
    with (OUTPUT / name).open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(record, stream, indent=2, allow_nan=False)
        stream.write('\n')


def positive_proposal(log_literal):
    if isinstance(log_literal, float):
        raise ValueError('Expected an exact saved log literal')
    value = 2. ** float(F(log_literal))
    if not isfinite(value) or value <= 0:
        raise ValueError('Conversion did not produce a finite positive proposal')
    return F(value)


def main():
    started = monotonic()
    source_bytes = SOURCE.read_bytes()
    source = json.loads(source_bytes)
    witness = source['witness']
    if (source.get('active_groups') != 256 or source.get('state_dimension') != 30
            or source.get('construction') != provider.local.construction(30)
            or witness.get('target_distance') != '2/5'
            or witness.get('field_order') != provider.FIELD_ORDER
            or witness.get('pair_swaps') is not True
            or witness.get('outer_measure') != 'invertible-toeplitz'):
        raise ValueError('Unexpected source construction or calibration')
    z, alpha = F(witness['z']), F(witness['alpha'])
    markers = tuple(tuple(positive_proposal(x) for x in row)
                    for row in witness['marker_log_odds_exact'])
    if (not markers or any(len(row) != 2 for row in markers)
            or len(set(markers)) != len(markers)):
        raise ValueError('Converted markers must be distinct positive pairs')
    ratio = positive_proposal(witness['log2_potential_ratio_exact'])
    ratios = (ratio / 2, ratio, ratio * 2)
    runner = Path(__file__).resolve()
    pins = dict(provider=provider.source_pins(), band=restored.source_pins(),
                runner={runner.name: sha256(runner.read_bytes()).hexdigest()})

    def check_inputs():
        if SOURCE.read_bytes() != source_bytes:
            raise RuntimeError('Pinned source proposal changed')
        if (provider.source_pins() != pins['provider']
                or restored.source_pins() != pins['band']
                or sha256(runner.read_bytes()).hexdigest() != pins['runner'][runner.name]):
            raise RuntimeError('Pinned replay closure changed')

    provenance = dict(path=str(SOURCE), sha256=sha256(source_bytes).hexdigest(),
                      calibration_active_groups=256,
                      role='exact proposal literals only; source numerical bounds not reused')
    proposal = dict(schema='field128-gamma-exact-proposal-v1',
                    construction=provider.local.construction(30), state_dimension=30,
                    target_distance='2/5', field_order=provider.FIELD_ORDER,
                    calibration_active_groups=256, source=provenance,
                    conversion='F(2. ** float(F(saved exact log literal))); actual binary64 result frozen as a rational',
                    conversion_is_exact_exponential_claim=False,
                    source_marker_log_odds_exact=witness['marker_log_odds_exact'],
                    source_log2_potential_ratio_exact=witness['log2_potential_ratio_exact'],
                    state_potential_ratio=str(ratio), potential_candidates=list(map(str, ratios)),
                    marker_count=len(markers), old_numerical_ceilings_reused=False,
                    directed=False, whole_code_certificate=False,
                    witness=dict(z=str(z), alpha=str(alpha), markers=[list(map(str, row)) for row in markers],
                                 grid=128, precision=256, packet_mantissa_bits=80),
                    replay_source_sha256=pins,
                    requested_band=[241, 272], build_seconds=600., band_seconds=120.,
                    max_weights=32, failure_bits=64, support_method='runs')
    OUTPUT.mkdir(exist_ok=False)
    check_inputs()
    write_new('proposal.json', proposal)

    def progress(row):
        if row.get('stage') != 'directed-interval-all-pair-alias' or row['marker'] % 8 == 0 or row['marker'] == len(markers)-1:
            print(json.dumps(row, allow_nan=False), flush=True)

    try:
        bank = provider.AliasBank(z, alpha, markers, state=30, distance=F(2, 5),
                                  precision=256, grid=128, bits=80, max_seconds=600., progress=progress)
        check_inputs()
        bank_record = dict(bank.metadata(), witness_proposal=proposal,
                           bank_build_elapsed_seconds=monotonic()-started,
                           scope='fresh directed bank metadata and exact witnesses; no selected-weight or whole-code claim')
        write_new('bank.json', bank_record)
        print(json.dumps(dict(stage='bank-saved', marker_count=len(markers), elapsed_seconds=monotonic()-started)), flush=True)
    except Exception as error:
        write_new('summary.json', dict(status='bank-build-failed', directed=False,
                    whole_code_certificate=False, complete_requested_range=False,
                    pending_weight_ranges=[[241, 272]], log2_upper_rational=None,
                    witness_proposal=proposal, replay_source_sha256=pins,
                    failure=dict(error_type=type(error).__name__, message=str(error))))
        raise

    try:
        result = restored.band_receipt(bank, 241, 272, ratios, method='runs', max_weights=32,
                                       max_seconds=120., failure_bits=64, progress=progress)
        check_inputs()
        result.update(witness_proposal=provenance,
                      proposal_sha256=sha256((OUTPUT/'proposal.json').read_bytes()).hexdigest(),
                      bank_metadata_sha256=sha256((OUTPUT/'bank.json').read_bytes()).hexdigest(),
                      replay_source_sha256=pins, total_elapsed_seconds=monotonic()-started)
        write_new('band_0241_0272.json', result)
        summary = {key: result[key] for key in ('construction', 'target_distance', 'state_dimension',
                    'requested_weight_range', 'weights_evaluated', 'every_integer_weight_evaluated',
                    'completed_weight_range', 'pending_weight_ranges', 'complete_requested_range',
                    'completed_union_log2_upper_rational', 'log2_upper_rational',
                    'requested_range_meets_target', 'weights_above_target', 'failure_bits',
                    'elapsed_seconds', 'total_elapsed_seconds', 'failure')}
        summary.update(status='complete' if result['complete_requested_range'] else 'partial',
                       directed=True, whole_code_certificate=False, old_numerical_ceilings_reused=False,
                       witness_proposal=provenance, replay_source_sha256=pins,
                       artifact_sha256={name: sha256((OUTPUT/name).read_bytes()).hexdigest()
                           for name in ('proposal.json', 'bank.json', 'band_0241_0272.json')})
        write_new('summary.json', summary)
        print(json.dumps({key: summary[key] for key in ('status', 'weights_evaluated',
              'completed_weight_range', 'pending_weight_ranges', 'log2_upper_rational',
              'requested_range_meets_target', 'weights_above_target', 'failure',
              'total_elapsed_seconds', 'artifact_sha256')}, allow_nan=False), flush=True)
    except Exception as error:
        write_new('summary.json', dict(status='band-replay-failed', directed=False,
                    whole_code_certificate=False, complete_requested_range=False,
                    log2_upper_rational=None, witness_proposal=provenance,
                    replay_source_sha256=pins,
                    failure=dict(error_type=type(error).__name__, message=str(error))))
        raise


if __name__ == '__main__':
    main()
