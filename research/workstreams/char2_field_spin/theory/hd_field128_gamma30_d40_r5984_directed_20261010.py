"""One exact13-marker directed bank and bounded explicit5984..6047 cover."""

import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'

from fractions import Fraction as F
from hashlib import sha256
import gc
import json
from shutil import disk_usage
from pathlib import Path
from time import monotonic

import hd_field128_gamma_all_pair_certificate as p
import hd_field128_gamma_midpoint_weight_bands as wb


HERE = Path(__file__).resolve().parent
RUNNER = Path(__file__).resolve()
RUNNER_HASH = sha256(RUNNER.read_bytes()).hexdigest()
PROPOSAL = HERE/'hd_field128_gamma30_d40_r5984_endpoint_calibration_20261010/proposal.json'
PROPOSAL_HASH = '74406acabc8a132172642bab535c58bdfcd3c711b591f781e9998e7d9caee4a0'
OUT = HERE/'hd_field128_gamma30_d40_r5984_winner13_directed256_20261010'
def write_new(path, value):
    if disk_usage(HERE).free < 750_000_000:
        raise OSError('Checkpoint guard: less than750MB disk space remains')
    return wb.profiles.write_new(path, value)


def initialize():
    raw = PROPOSAL.read_bytes()
    if sha256(raw).hexdigest() != PROPOSAL_HASH:
        raise RuntimeError('Frozen exact13 proposal changed')
    proposal = json.loads(raw)
    if (proposal['construction'] != p.local.construction(30) or proposal['state_dimension'] != 30
            or proposal['active_groups'] != 5984 or proposal['target_distance'] != '2/5'
            or proposal['marker_count'] != 13 or len(proposal['witness']['markers']) != 13):
        raise ValueError('Unexpected literal proposal scope')
    pins = wb.source_pins()
    for name, digest in proposal['source_sha256'].items():
        if name in pins and pins[name] != digest:
            raise RuntimeError('Source pin conflict: '+name)
        pins[name] = digest
    pins[RUNNER.name] = RUNNER_HASH
    def check():
        if PROPOSAL.read_bytes() != raw:
            raise RuntimeError('Frozen literal proposal changed')
        wb.check_source_pins(pins)
    check()
    OUT.mkdir(exist_ok=False)
    witness = proposal['witness']
    started = monotonic()
    bank = p.AliasBank(F(witness['z']), F(witness['alpha']),
        [tuple(map(F, row)) for row in witness['markers']], state=30, distance=F(2, 5),
        grid=128, precision=256, bits=80, max_seconds=180., max_table_bytes=64 << 20,
        progress=lambda row: print(json.dumps(row), flush=True))
    check()
    metadata = dict(bank.metadata(), witness_source=dict(path=str(PROPOSAL), sha256=PROPOSAL_HASH),
        calibration_active_groups=5984, orchestration_source_sha256={RUNNER.name: RUNNER_HASH},
        proposal_source_sha256=proposal['source_sha256'], old_numerical_bounds_reused=False,
        bank_build_seconds=monotonic()-started)
    bank_pin = write_new(OUT/'bank.json', metadata)
    print(json.dumps(dict(stage='bank-complete', bank=bank_pin,
                          elapsed_seconds=metadata['bank_build_seconds'])), flush=True)
    return bank, check


def first_band(bank, check):
    output = OUT/'band_5984_6047'
    output.mkdir(exist_ok=False)
    started, evaluator = monotonic(), None
    check()
    try:
        evaluator = wb.Evaluator(bank, 5984, 6047, max_seconds=180.,
            witness_source=dict(path=str(PROPOSAL), sha256=PROPOSAL_HASH),
            progress=lambda row: print(json.dumps(row), flush=True))
        result = evaluator.run_cover(max_nodes=512, max_seconds=60., cell_failure_bits=80, aggregate_bits=54,
            context_checkpoint=lambda row: write_new(output/'context.json', row),
            checkpoint=lambda row: write_new(output/f"node_{row['sequence']:06d}.json", row))
        check()
        result['orchestration_source_sha256'] = {RUNNER.name: RUNNER_HASH}
        result['total_preparation_and_cover_seconds'] = monotonic()-started
        pin = write_new(output/'summary.json', result)
        print(json.dumps(dict(stage='band-saved', summary=pin, status=result['status'],
            complete=result['complete_requested_band'], passing=result['meets_requested_band_target'],
            log2_upper_rational=result['log2_upper_rational'], nodes=result['tree_nodes'],
            elapsed_seconds=result['total_preparation_and_cover_seconds'])), flush=True)
        return result
    except (TimeoutError, MemoryError, OSError) as error:
        check()
        failure = dict(schema='gamma-directed-band-orchestration-failure-v1',
            requested_weight_range=[5984, 6047], complete_requested_band=False,
            log2_upper_rational=None, directed=False, whole_code_certificate=False,
            failure=dict(error_type=type(error).__name__, message=str(error)),
            witness_source=dict(path=str(PROPOSAL), sha256=PROPOSAL_HASH),
            orchestration_source_sha256={RUNNER.name: RUNNER_HASH}, elapsed_seconds=monotonic()-started)
        try:
            print(json.dumps(dict(failure=write_new(output/'failure.json', failure))), flush=True)
        except OSError:
            print(json.dumps(failure), flush=True)
        return failure
    finally:
        evaluator = None
        gc.collect()
        check()


if __name__ == '__main__':
    bank, check = initialize()
    last_result = first_band(bank, check)
    print('Retained actual13-marker bank as global bank; all Evaluator caches released.', flush=True)
