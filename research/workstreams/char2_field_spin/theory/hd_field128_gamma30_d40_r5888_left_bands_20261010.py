"""Five explicit bounded directed left bands with the retained exact17 bank."""

from hashlib import sha256
import gc
import json
from pathlib import Path
from time import monotonic

import hd_field128_gamma_midpoint_weight_bands as wb


HERE = Path(__file__).resolve().parent
SOURCE = Path(__file__).resolve()
SOURCE_HASH = sha256(SOURCE.read_bytes()).hexdigest()
PROPOSAL = HERE/'hd_field128_gamma30_d40_r5888_winner_replay_20261010/proposal.json'
PROPOSAL_HASH = '43cea88db107cfc364ab5989978269fa97b17d0068e32ca5d24aa930887e4754'
OUT = HERE/'hd_field128_gamma30_d40_r5888_winner17_directed256_20261010'
RANGES = ((5344, 5407), (5408, 5471), (5472, 5535), (5536, 5599), (5600, 5663))


def run(bank, parent_check):
    pins = wb.source_pins()
    pins[SOURCE.name] = SOURCE_HASH
    proposal, actual = json.loads(PROPOSAL.read_bytes()), bank.metadata()['witness']
    if any(actual[key] != proposal['witness'][key] for key in ('z', 'alpha', 'markers')):
        raise ValueError('Retained bank differs from exact17 proposal')
    def check():
        parent_check()
        wb.check_source_pins(pins)
        if sha256(PROPOSAL.read_bytes()).hexdigest() != PROPOSAL_HASH:
            raise RuntimeError('Literal proposal changed')
    check()
    reports = []
    for first, last in RANGES:
        check()
        output = OUT/f'band_{first:04d}_{last:04d}'
        output.mkdir(exist_ok=False)
        evaluator, started = None, monotonic()
        try:
            print(json.dumps(dict(stage='band-start', first=first, last=last)), flush=True)
            evaluator = wb.Evaluator(bank, first, last, max_seconds=180.,
                witness_source=dict(path=str(PROPOSAL), sha256=PROPOSAL_HASH))
            result = evaluator.run_cover(max_nodes=512, max_seconds=60., cell_failure_bits=80, aggregate_bits=54,
                context_checkpoint=lambda row: wb.profiles.write_new(output/'context.json', row),
                checkpoint=lambda row: wb.profiles.write_new(output/f"node_{row['sequence']:06d}.json", row))
            check()
            result['orchestration_source_sha256'] = {SOURCE.name: SOURCE_HASH}
            result['total_preparation_and_cover_seconds'] = monotonic()-started
            pin = wb.profiles.write_new(output/'summary.json', result)
            report = dict(requested_weight_range=[first, last], summary=pin, status=result['status'],
                complete=result['complete_requested_band'], passing=result['meets_requested_band_target'],
                log2_upper_rational=result['log2_upper_rational'], nodes=result['tree_nodes'],
                elapsed_seconds=result['total_preparation_and_cover_seconds'])
        except (TimeoutError, MemoryError, OSError) as error:
            check()
            report = dict(requested_weight_range=[first, last], status='partial', complete=False,
                passing=False, log2_upper_rational=None,
                failure=dict(error_type=type(error).__name__, message=str(error)),
                directed=False, whole_code_certificate=False, elapsed_seconds=monotonic()-started,
                orchestration_source_sha256={SOURCE.name: SOURCE_HASH})
            try:
                report['checkpoint'] = wb.profiles.write_new(output/'failure.json', report)
            except OSError:
                print(json.dumps(report), flush=True)
                raise
        finally:
            evaluator = None
            gc.collect()
            check()
        reports.append(report)
        print(json.dumps(report), flush=True)
    check()
    summary = dict(schema='explicit-directed-band-run-links-v1', requested_ranges=[list(row) for row in RANGES],
        bands=reports, source_sha256=pins, witness_source=dict(path=str(PROPOSAL), sha256=PROPOSAL_HASH),
        joint_probability_claim=False, no_weight_interpolation=True, whole_code_certificate=False)
    pin = wb.profiles.write_new(OUT/'left_extension_links.json', summary)
    print(json.dumps(dict(stage='left-extensions-saved', summary=pin)), flush=True)
    return summary
